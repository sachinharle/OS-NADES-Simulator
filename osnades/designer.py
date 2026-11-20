"""Multicomponent NADES design: find the component mole fractions that dissolve the most solute.

Workflow (used by the web app's "NADES workflow" page and ``osnades design``)::

    solute + components (+ bounds) ──► check σ-surfaces / fusion data
        ──► screen a composition lattice (batched openCOSMO-RS: γ∞ of the solute, activity of each
            component → is the solvent liquid at every temperature?)
        ──► rank by an objective over the temperature range
        ──► recompute the best candidates with the exact (iterative) SLE solubility.

A *component* is a formula unit that is weighed in: a neutral molecule (``glycerol``) or a salt that
COSMO-RS treats as ions (``choline chloride`` → choline⁺ + Cl⁻). Mole fractions in this module are
component fractions z (salt = one unit); openCOSMO-RS sees the species fractions.

Liquidity of the solvent: component i crystallises when its activity reaches the ideal-solid limit,
ln a_i ≥ −ΔH_fus,i/R (1/T − 1/T_m,i). With a_i referenced to pure liquid i (salts: pure molten salt,
as on the eutectic page). Components without fusion data are not checked. COSMO-RS over-stabilises
choline chloride mixtures, so for salts this check is optimistic.
"""
from __future__ import annotations

import itertools
import math
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from .library import Library
from .thermo import Engine, SolventSystem, ideal_rhs, presets

OBJECTIVES = {
    "mean_x": "Highest solubility over the whole range (geometric mean of x)",
    "min_T": "Highest solubility at the lowest temperature",
    "g_per_kg": "Highest mass solubility (mean g solute / kg solvent)",
    "cooling": "Best for cooling crystallisation (largest g/kg gain from T_min to T_max)",
}


# --------------------------------------------------------------------------- components
@dataclass
class DesignComponent:
    name: str
    species: dict[str, float]  # COSMO species -> stoichiometric coefficient
    lo: float = 0.0  # bounds on the component mole fraction
    hi: float = 1.0
    dh_fus: float | None = None  # J/mol per formula unit (for the liquidity check)
    t_fus: float | None = None
    dcp_fus: float = 0.0
    fusion_source: str = ""

    @property
    def molar_mass(self) -> float:
        return getattr(self, "_mm", float("nan"))

    @property
    def has_fusion(self) -> bool:
        return bool(self.dh_fus and self.t_fus)


def component_choices(library: Library) -> list[str]:
    """Neutral library species, preset salts (if their ions are present) and preset species not yet computed."""
    names = [s.name for s in library.all() if s.charge == 0]
    for salt, d in presets().get("salts", {}).items():
        names.append(salt)
    for n in presets()["species"]:
        if n not in names and n not in ("choline", "chloride"):
            names.append(n)
    return list(dict.fromkeys(names))


def species_of(name: str) -> dict[str, float]:
    salts = presets().get("salts", {})
    if name in salts:
        return {k: float(v) for k, v in salts[name]["ions"].items()}
    if "+" in name:
        return {s.strip(): 1.0 for s in name.split("+")}
    return {name: 1.0}


def fusion_of(name: str, library: Library) -> tuple[float | None, float | None, float, str]:
    """(ΔH_fus, T_m, ΔCp_fus, source): library first, then presets (salts, well-known species)."""
    salts = presets().get("salts", {})
    if name in salts:
        d = salts[name]
        return d.get("dh_fus"), d.get("t_fus"), d.get("dcp_fus", 0.0), d.get("fusion_source", "preset")
    if library.has(name):
        sp = library.get(name)
        if sp.dh_fus and sp.t_fus:
            return sp.dh_fus, sp.t_fus, sp.dcp_fus or 0.0, sp.fusion_source or sp.notes or "library"
    d = presets()["species"].get(name, {})
    return d.get("dh_fus"), d.get("t_fus"), d.get("dcp_fus", 0.0), d.get("fusion_source", "")


def make_component(name: str, library: Library, lo: float = 0.0, hi: float = 1.0) -> DesignComponent:
    dh, tm, dcp, src = fusion_of(name, library)
    return DesignComponent(name, species_of(name), lo, hi, dh, tm, dcp or 0.0, src)


def missing_species(solute: str, components: list[DesignComponent], library: Library) -> list[str]:
    """Species whose σ-surface still has to be computed (in order)."""
    out = [] if library.has(solute) else [solute]
    for c in components:
        out += [s for s in c.species if not library.has(s) and s not in out]
    return out


# --------------------------------------------------------------------------- lattice
def lattice(components: list[DesignComponent], step: float) -> np.ndarray:
    """All component-fraction vectors on a simplex lattice (spacing ``step``) inside the bounds."""
    n = len(components)
    m = max(1, int(round(1.0 / step)))
    pts = []
    for c in itertools.combinations_with_replacement(range(n), m):
        z = np.bincount(c, minlength=n) / m
        pts.append(z)
    Z = np.unique(np.array(pts), axis=0) if n > 1 else np.ones((1, 1))
    lo = np.array([c.lo for c in components]) - 1e-9
    hi = np.array([c.hi for c in components]) + 1e-9
    return Z[np.all((Z >= lo) & (Z <= hi), axis=1)]


def lattice_size(n: int, step: float) -> int:
    m = max(1, int(round(1.0 / step)))
    return math.comb(n + m - 1, m)


def ratio_text(z: np.ndarray) -> str:
    """Compact molar ratio, smallest non-zero component = 1 (e.g. "1 : 2.1 : 0.4")."""
    nz = z[z > 1e-9]
    base = nz.min() if len(nz) else 1.0
    return " : ".join(f"{v / base:.3g}" if v > 1e-9 else "0" for v in z)


# --------------------------------------------------------------------------- designer
@dataclass
class DesignResult:
    solute: str
    components: list[DesignComponent]
    temps: np.ndarray
    objective: str
    screen: pd.DataFrame  # one row per lattice point: z_*, x_inf at each T, liquid flags, score
    best: pd.DataFrame  # top candidates with exact solubilities
    curves: pd.DataFrame  # exact solubility vs T for the top candidates
    notes: list[str] = field(default_factory=list)

    def best_system(self, rank: int = 0) -> SolventSystem:
        row = self.best.iloc[rank]
        z = np.array([row[f"z_{c.name}"] for c in self.components])
        return composition_system(self.components, z, name=row["name"])


def composition_system(components: list[DesignComponent], z: np.ndarray, name: str | None = None) -> SolventSystem:
    comps: dict[str, float] = {}
    for c, zi in zip(components, z):
        if zi <= 1e-12:
            continue
        for s, nu in c.species.items():
            comps[s] = comps.get(s, 0.0) + zi * nu
    label = name or " : ".join(c.name for c, zi in zip(components, z) if zi > 1e-12) + f" = {ratio_text(z[z > 1e-12])}"
    return SolventSystem(label, comps)


class Designer:
    def __init__(self, engine: Engine, solute: str, components: list[DesignComponent],
                 dh_fus: float, t_fus: float, dcp_fus: float = 0.0):
        if not components:
            raise ValueError("Add at least one solvent component")
        lib = engine.library
        self.engine, self.solute, self.components = engine, solute, components
        self.fus = (dh_fus, t_fus, dcp_fus or 0.0)
        names = [c.name for c in components]
        if len(set(names)) != len(names):
            raise ValueError("Each component may appear only once")
        self.solvent_species: list[str] = []
        for c in components:
            for s in c.species:
                if s == solute:
                    raise ValueError(f"The solute {solute!r} cannot also be a solvent component")
                if s in self.solvent_species:
                    raise ValueError(f"Species {s!r} appears in two components")
                self.solvent_species.append(s)
                if not lib.has(s):
                    raise KeyError(f"σ-surface of {s!r} is missing – compute it first")
            q = sum(lib.get(s).charge * nu for s, nu in c.species.items())
            if abs(q) > 1e-9:
                raise ValueError(f"Component {c.name!r} is not electroneutral (charge {q:+g})")
            c._mm = sum(lib.get(s).molar_mass * nu for s, nu in c.species.items())
        if sum(c.lo for c in components) > 1 + 1e-9:
            raise ValueError("The lower bounds add up to more than 1")
        self.species = [solute] + self.solvent_species
        # ν matrix: component × solvent species
        self.nu = np.array([[c.species.get(s, 0.0) for s in self.solvent_species] for c in components])
        self.m_solute = lib.get(solute).molar_mass

    # ------------------------------------------------------------------ helpers
    def species_x(self, Z: np.ndarray) -> np.ndarray:
        """Solvent species mole fractions for component fractions Z (rows)."""
        n = np.atleast_2d(Z) @ self.nu
        return n / n.sum(axis=1, keepdims=True)

    def solvent_molar_mass(self, Z: np.ndarray) -> np.ndarray:
        """Mass of solvent per mole of solvent *species* (matches the mole fraction basis of x)."""
        Z = np.atleast_2d(Z)
        mass = Z @ np.array([c._mm for c in self.components])
        return mass / (Z @ self.nu).sum(axis=1)

    def _lng(self, X: np.ndarray, T: float) -> np.ndarray:
        """ln γ (cosmo reference) for rows of X; rows where COSMOspace fails become NaN."""
        try:
            return self.engine.ln_gamma(self.species, X, T, refst="cosmo")
        except Exception:  # noqa: BLE001
            if len(X) == 1:
                return np.full((1, len(self.species)), np.nan)
            h = len(X) // 2
            return np.vstack([self._lng(X[:h], T), self._lng(X[h:], T)])

    # ------------------------------------------------------------------ screening
    def screen_T(self, Z: np.ndarray, T: float) -> dict[str, np.ndarray]:
        """Batched evaluation at one temperature: solute solubility (successive substitution from γ∞)
        and the liquidity margin of every component with fusion data."""
        nS, nC = len(self.solvent_species), len(self.components)
        Xs = self.species_x(Z)
        rows = [np.hstack([[0.0], x]) for x in Xs]  # solute at infinite dilution
        pure_solute = np.zeros(nS + 1)
        pure_solute[0] = 1.0
        pure_comp = []
        for i in range(nC):
            v = np.hstack([[0.0], self.nu[i] / self.nu[i].sum()])
            pure_comp.append(v)
        lg = self._lng(np.vstack(rows + [pure_solute] + pure_comp), T)
        nZ = len(Z)
        lg_mix, lg_ps, lg_pc = lg[:nZ], lg[nZ, 0], lg[nZ + 1:]
        rhs = ideal_rhs(T, *self.fus)
        lng_inf = lg_mix[:, 0] - lg_ps
        # SLE root ln x + ln γ(x) = rhs, solved for all compositions at once: g(ln x) on a coarse grid
        # of solute fractions (one batched openCOSMO-RS call), first sign change, linear interpolation.
        grid = np.concatenate([np.linspace(-20.0, -4.0, 5), np.linspace(-3.0, -0.02, 9)])
        xg = np.exp(grid)
        Xg = np.vstack([np.hstack([np.full((nZ, 1), v), (1 - v) * Xs]) for v in xg])
        lgg = (self._lng(Xg, T)[:, 0] - lg_ps).reshape(len(grid), nZ).T  # (nZ, n_grid)
        gv = grid[None, :] + lgg - rhs
        x = np.full(nZ, np.nan)
        for k in range(nZ):
            g = gv[k]
            if not np.all(np.isfinite(g)):
                continue
            if g[0] > 0:
                x[k] = np.exp(rhs - lng_inf[k])  # very low solubility: infinite-dilution limit
                continue
            sc = np.where((g[:-1] <= 0) & (g[1:] > 0))[0]
            if len(sc) == 0:
                x[k] = 1.0  # no saturation below x ≈ 1: (nearly) miscible
                continue
            a = sc[0]
            t = g[a] / (g[a] - g[a + 1])
            x[k] = float(np.exp(grid[a] + t * (grid[a + 1] - grid[a])))
        if rhs >= 0:
            x[:] = 1.0
        out = {"x": x, "ln_gamma_inf": lng_inf}
        # liquidity of the solute-free solvent
        margins = np.full((nZ, nC), np.inf)
        for i, c in enumerate(self.components):
            if not c.has_fusion:
                continue
            nu = self.nu[i]
            m = nu > 0
            x0 = nu / nu.sum()
            with np.errstate(divide="ignore", invalid="ignore"):
                lna = np.sum(nu[m] * (np.log(Xs[:, m]) + lg_mix[:, 1:][:, m]
                                      - np.log(x0[m]) - lg_pc[i, 1:][m]), axis=1)
            lim = ideal_rhs(T, c.dh_fus, c.t_fus, c.dcp_fus)
            margins[:, i] = np.where(Z[:, i] > 1e-12, lim - lna, np.inf)
        out["margin"] = margins.min(axis=1)
        out["limiting"] = np.array([self.components[j].name if np.isfinite(margins[k, j]) else ""
                                    for k, j in enumerate(margins.argmin(axis=1))])
        return out

    def score(self, xs: np.ndarray, gkg: np.ndarray, objective: str) -> np.ndarray:
        """xs, gkg: shape (n_points, n_T) with temperatures ascending."""
        with np.errstate(divide="ignore", invalid="ignore"):
            if objective == "min_T":
                return np.log(xs[:, 0])
            if objective == "g_per_kg":
                return np.nanmean(gkg, axis=1)
            if objective == "cooling":
                return gkg[:, -1] - gkg[:, 0]
            return np.nanmean(np.log(xs), axis=1)

    def _evaluate(self, Z: np.ndarray, temps: np.ndarray, progress=None, p0: float = 0.0, p1: float = 1.0,
                  label: str = "Screening"):
        nT = len(temps)
        xs, marg = np.zeros((len(Z), nT)), np.zeros((len(Z), nT))
        lim = np.empty((len(Z), nT), dtype=object)
        for j, T in enumerate(temps):
            if progress:
                progress(p0 + (p1 - p0) * j / nT, f"{label} {len(Z)} compositions at {T - 273.15:.1f} °C")
            r = self.screen_T(Z, float(T))
            xs[:, j], marg[:, j], lim[:, j] = r["x"], r["margin"], r["limiting"]
        return xs, marg, lim

    def _neighbours(self, z: np.ndarray, h: float) -> list[np.ndarray]:
        """Move a fraction h from one component to another (all pairs), inside the bounds."""
        lo = np.array([c.lo for c in self.components])
        hi = np.array([c.hi for c in self.components])
        out = []
        for i in range(len(z)):
            for j in range(len(z)):
                if i != j:
                    y = z.copy()
                    y[i] += h
                    y[j] -= h
                    if np.all(y >= lo - 1e-9) and np.all(y <= hi + 1e-9) and np.all(y >= -1e-12):
                        out.append(np.round(np.clip(y, 0, 1), 6))
        return out

    def run(self, temps, step: float = 0.1, objective: str = "mean_x", require_liquid: bool = True,
            top_n: int = 8, refine: bool = True, local_rounds: int = 2, progress=None) -> DesignResult:
        """Screen the lattice, then ``local_rounds`` of local search around the best points (step halved
        each round), then exact solubilities for the ``top_n`` candidates (if ``refine``)."""
        temps = np.sort(np.atleast_1d(np.asarray(temps, float)))
        Z = lattice(self.components, step)
        if len(Z) == 0:
            raise ValueError("No composition satisfies the bounds – widen them or use a finer step")
        notes: list[str] = []
        span = 0.6 if local_rounds else 0.7
        xs, marg, lim = self._evaluate(Z, temps, progress, 0.0, span)
        mm_s = np.array([self.m_solute])

        def scores(xs_, Z_, marg_):
            mm = self.solvent_molar_mass(Z_)[:, None]
            with np.errstate(divide="ignore", invalid="ignore"):
                g = 1000.0 * xs_ * mm_s[0] / ((1 - xs_) * mm)
            sc_ = self.score(xs_, g, objective)
            if require_liquid and (marg_ > 0).all(axis=1).any():
                sc_ = np.where((marg_ > 0).all(axis=1), sc_, -np.inf)
            return sc_

        h = step
        for r in range(local_rounds if len(self.components) > 1 else 0):
            h /= 2
            sc_now = scores(xs, Z, marg)
            seen = {tuple(np.round(z, 6)) for z in Z}
            new = []
            for k in np.argsort(-np.nan_to_num(sc_now, nan=-np.inf))[:3]:
                for y in self._neighbours(Z[k], h):
                    if tuple(y) not in seen:
                        seen.add(tuple(y))
                        new.append(y)
            if not new:
                break
            Zn = np.array(new)
            p0 = span + (0.7 - span) * r / local_rounds
            xn, mn, ln = self._evaluate(Zn, temps, progress, p0, p0 + (0.7 - span) / local_rounds,
                                        label=f"Refining (step {h:.3g}):")
            Z, xs, marg, lim = np.vstack([Z, Zn]), np.vstack([xs, xn]), np.vstack([marg, mn]), np.vstack([lim, ln])
        nT = len(temps)
        mm = self.solvent_molar_mass(Z)[:, None]
        with np.errstate(divide="ignore", invalid="ignore"):
            gkg = 1000.0 * xs * self.m_solute / ((1 - xs) * mm)
        liquid = marg > 0
        all_liquid = liquid.all(axis=1)
        sc = self.score(xs, gkg, objective)
        df = pd.DataFrame({f"z_{c.name}": Z[:, i] for i, c in enumerate(self.components)})
        df.insert(0, "ratio", [ratio_text(z) for z in Z])
        for j, T in enumerate(temps):
            df[f"x@{T:.2f}"] = xs[:, j]
            df[f"liquid@{T:.2f}"] = liquid[:, j]
        df["liquid_all_T"] = all_liquid
        df["limiting_component"] = [next((l for l, ok in zip(row_l, row_ok) if not ok and l), "")
                                    for row_l, row_ok in zip(lim, liquid)]
        df["score"] = sc
        unchecked = [c.name for c in self.components if not c.has_fusion]
        if unchecked:
            notes.append("No fusion data for " + ", ".join(unchecked) + ": they are assumed never to crystallise.")
        if any(len(c.species) > 1 for c in self.components):
            notes.append("COSMO-RS over-stabilises salt (ionic) components, so the liquidity check is optimistic "
                         "for them – verify the melting point of the chosen mixture experimentally.")
        cand = df[np.isfinite(df["score"])]
        if require_liquid:
            if not all_liquid.any():
                notes.append("No composition is liquid at every temperature of the range – showing the best "
                             "compositions regardless; raise T_min or add/allow other components (e.g. water).")
            else:
                cand = cand[cand["liquid_all_T"]]
        cand = cand.sort_values("score", ascending=False)
        top = cand.head(top_n)
        best_rows, curves = [], []
        for k, (idx, row) in enumerate(top.iterrows()):
            z = Z[idx]
            name = " : ".join(c.name for c, zi in zip(self.components, z) if zi > 1e-12) + \
                   f" = {ratio_text(z[z > 1e-12])}"
            sysm = composition_system(self.components, z, name=name)
            if progress:
                progress(0.7 + 0.3 * k / max(1, len(top)), f"Exact solubility for {name}")
            rec = {"rank": k + 1, "name": name, "ratio": row["ratio"], "liquid_all_T": bool(row["liquid_all_T"])}
            rec.update({f"z_{c.name}": float(z[i]) for i, c in enumerate(self.components)})
            xe = []
            for T in temps:
                if refine:
                    r = self.engine.solubility(self.solute, sysm, float(T), *self.fus)
                    xT, gT, note = r.x, r.g_per_kg, r.note
                else:
                    j = int(np.where(temps == T)[0][0])
                    xT, gT, note = xs[idx, j], gkg[idx, j], "screening value"
                xe.append((xT, gT))
                curves.append({"rank": k + 1, "name": name, "T": float(T), "x": xT, "g_per_kg": gT,
                               "liquid": bool(liquid[idx, int(np.where(temps == T)[0][0])]), "note": note})
                rec[f"x@{T - 273.15:.1f}°C"] = xT
                rec[f"g/kg@{T - 273.15:.1f}°C"] = gT
            xe_a = np.array(xe)
            rec["score"] = float(self.score(xe_a[None, :, 0], xe_a[None, :, 1], objective)[0])
            best_rows.append(rec)
        best = pd.DataFrame(best_rows)
        if len(best) and refine:
            best = best.sort_values("score", ascending=False).reset_index(drop=True)
            order = {n: i + 1 for i, n in enumerate(best["name"])}
            best["rank"] = range(1, len(best) + 1)
            curves = [dict(c, rank=order[c["name"]]) for c in curves]
        if progress:
            progress(1.0, "done")
        return DesignResult(self.solute, self.components, temps, objective, df, best, pd.DataFrame(curves), notes)


def ideal_curve(temps, dh_fus: float, t_fus: float, dcp_fus: float = 0.0) -> np.ndarray:
    return np.minimum(1.0, np.exp([ideal_rhs(float(T), dh_fus, t_fus, dcp_fus) for T in temps]))


# --------------------------------------------------------------------------- phase diagram of the solvent
@dataclass
class SolventPhase:
    """Solid–liquid phase diagram of a (solute-free) solvent with one composition of interest."""
    names: list[str]
    z: np.ndarray  # composition of interest (component mole fractions, same order as names)
    diagram: object  # eutectic.PhaseDiagram (2 components) or eutectic.TernaryDiagram (3)
    T_melt: float = float("nan")  # liquidus temperature of z (K)
    T_melt_below: float = float("nan")  # set instead when the liquidus lies below the checked range

    @property
    def ternary(self) -> bool:
        return len(self.names) == 3

    def liquid_at(self, T: float) -> bool | None:
        """Is z fully liquid at T? None if the liquidus could not be computed."""
        if np.isfinite(self.T_melt):
            return bool(self.T_melt < T)
        if np.isfinite(self.T_melt_below):
            return bool(self.T_melt_below <= T)
        return None

    def melt_text(self) -> str:
        if np.isfinite(self.T_melt):
            return f"{self.T_melt - 273.15:.1f} °C"
        if np.isfinite(self.T_melt_below):
            return f"< {self.T_melt_below - 273.15:.0f} °C"
        return "—"

    def eutectic_text(self) -> str:
        from .eutectic import ratio_label

        d = self.diagram
        T = (f"{d.T_e - 273.15:.1f} °C" if np.isfinite(d.T_e) else
             f"< {d.T_e_upper - 273.15:.0f} °C" if np.isfinite(d.T_e_upper) else "—")
        if self.ternary:
            comp = ratio_label(d.z_e) if np.all(np.isfinite(d.z_e)) else ""
            return f"{T} at {':'.join(self.names)} = {comp}" if comp else T
        return f"{T} at x_{self.names[0]} = {d.x_e:.3f}" if np.isfinite(d.x_e) else T


def phase_components(components: list[DesignComponent], z: np.ndarray) -> tuple[list[DesignComponent], np.ndarray]:
    """Components for the phase diagram: all of them when there are ≤ 3, else only those present in z."""
    z = np.asarray(z, float)
    if len(components) <= 3:
        return list(components), z
    keep = z > 1e-9
    return [c for c, k in zip(components, keep) if k], z[keep]


def solvent_phase_diagram(engine: Engine, components: list[DesignComponent], z, step: float = 0.1,
                          progress=None) -> SolventPhase:
    """Binary or ternary solid–liquid diagram of the solvent, plus the liquidus of composition z."""
    from .eutectic import Component, EutecticSolver, TernarySolver

    comps, zc = phase_components(components, z)
    if len(comps) not in (2, 3):
        raise ValueError(f"a phase diagram needs 2 or 3 components (the composition has {len(comps)})")
    missing = [c.name for c in comps if not c.has_fusion]
    if missing:
        raise ValueError("fusion data (ΔH_fus, T_m) missing for " + ", ".join(missing))
    ecs = [Component(c.name, c.species, float(c.dh_fus), float(c.t_fus), float(c.dcp_fus or 0.0),
                     c.fusion_source) for c in comps]
    out = SolventPhase([c.name for c in comps], zc, None)
    if len(comps) == 3:
        solver = TernarySolver(engine, ecs)
        T_lo = max(173.15, 0.5 * min(c.t_fus for c in ecs))
        out.diagram = solver.solve(step=step, T_lo=T_lo, progress=progress)
        pt = solver.evaluate(zc[None, :], T_lo).iloc[0]
        if np.isfinite(pt["T_liquidus"]):
            out.T_melt = float(pt["T_liquidus"])
        elif pt["below_range"]:
            out.T_melt_below = float(pt["T_upper"])
    else:
        solver = EutecticSolver(engine, *ecs)
        out.diagram = solver.solve(progress=progress)
        T_lo = max(50.0, 0.3 * min(c.t_fus for c in ecs))
        (tA, nA, bA), (tB, nB, bB) = (solver.liquidus(float(zc[0]), i, T_lo) for i in (0, 1))
        if np.isfinite([tA, tB]).any():
            out.T_melt = float(np.nanmax([tA, tB]))
        elif "below" in (nA, nB):
            out.T_melt_below = float(np.nanmax([bA, bB]))
    return out
