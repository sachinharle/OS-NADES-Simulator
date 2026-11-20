"""Solid–liquid phase diagram (eutectic) of a binary HBA + HBD mixture with openCOSMO-RS.

For each component i (immiscible pure solids, ideal solid phases)::

    ln a_i(x, T) = -ΔH_i/R (1/T - 1/T_m,i) + ΔCp_i/R (T_m,i/T - 1 - ln(T_m,i/T))

The liquidus of i is the temperature solving this equation at composition x; the mixture
liquidus is max(T_A, T_B) and the eutectic is where both branches meet.

A component may be a salt made of several COSMO species (choline chloride = choline⁺ + Cl⁻).
Its activity is the stoichiometric product of the ion activities, normalised to its own pure
(sub-cooled) molten liquid::

    ln a_i = Σ_k ν_ik [ln x_k + ln γ_k(mix) − ln x_k⁰ − ln γ_k(pure i)]

so ln a_i = 0 for pure liquid i, and for a neutral molecule it reduces to the usual
pure-liquid reference. "Ideal" curves use γ = 1 in the same expression, so the gap
between predicted and ideal eutectic shows how "deep" the eutectic is.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy.optimize import brentq

from .library import Library
from .thermo import R, Engine, ideal_rhs, presets


@dataclass
class Component:
    name: str
    species: dict[str, float]  # COSMO species -> stoichiometric coefficient
    dh_fus: float  # J/mol (per formula unit)
    t_fus: float  # K
    dcp_fus: float = 0.0
    source: str = ""

    @property
    def x0(self) -> np.ndarray:
        nu = np.array(list(self.species.values()), float)
        return nu / nu.sum()


def component_names(library: Library) -> list[str]:
    """Library neutral species plus preset salts whose ions are all in the library."""
    names = [s.name for s in library.all() if s.charge == 0]
    for salt, d in presets().get("salts", {}).items():
        if all(library.has(i) for i in d["ions"]):
            names.append(salt)
    return names


def fusion_defaults(name: str, library: Library) -> tuple[float | None, float | None, float, str]:
    """(ΔH_fus, T_m, ΔCp_fus, source) from the library, else the presets."""
    salts = presets().get("salts", {})
    if name in salts:
        d = salts[name]
        return d.get("dh_fus"), d.get("t_fus"), d.get("dcp_fus", 0.0), d.get("fusion_source", "preset")
    if library.has(name):
        sp = library.get(name)
        if sp.dh_fus and sp.t_fus:
            return sp.dh_fus, sp.t_fus, sp.dcp_fus or 0.0, "library"
    d = presets()["species"].get(name, {})
    return d.get("dh_fus"), d.get("t_fus"), d.get("dcp_fus", 0.0), d.get("fusion_source", d.get("notes", "preset"))


def resolve_component(name: str, library: Library, dh_fus: float | None = None, t_fus: float | None = None,
                      dcp_fus: float | None = None) -> Component:
    """Build a component from a library species, a preset salt, or ``"cation+anion"``."""
    salts = presets().get("salts", {})
    if name in salts:
        species = {k: float(v) for k, v in salts[name]["ions"].items()}
    elif "+" in name:
        species = {s.strip(): 1.0 for s in name.split("+")}
    else:
        species = {name: 1.0}
    for s in species:
        if not library.has(s):
            raise KeyError(f"Species {s!r} is not in the library (osnades add {s})")
    q = sum(library.get(s).charge * nu for s, nu in species.items())
    if abs(q) > 1e-9:
        raise ValueError(f"Component {name!r} is not electroneutral (charge {q:+g})")
    dh0, tm0, dcp0, src = fusion_defaults(name, library)
    dh = dh_fus if dh_fus is not None else dh0
    tm = t_fus if t_fus is not None else tm0
    dcp = dcp_fus if dcp_fus is not None else (dcp0 or 0.0)
    if not dh or not tm:
        raise ValueError(f"ΔH_fus and T_m are required for {name!r} (osnades set {name} --dh-fus .. --t-fus ..)")
    user = dh_fus is not None or t_fus is not None
    return Component(name, species, float(dh), float(tm), float(dcp), "user" if user else src)


def _t_from_lna(lna: float, c: Component) -> float:
    """Ideal-solid liquidus temperature for a given component activity."""
    if not np.isfinite(lna):
        return float("nan")
    if lna >= 0:
        return c.t_fus
    if not c.dcp_fus:
        return 1.0 / (1.0 / c.t_fus - R * lna / c.dh_fus)
    f = lambda T: ideal_rhs(T, c.dh_fus, c.t_fus, c.dcp_fus) - lna  # noqa: E731
    lo = 1.0
    return brentq(f, lo, c.t_fus, xtol=1e-6) if f(lo) < 0 else float("nan")


@dataclass
class PhaseDiagram:
    a: Component
    b: Component
    data: pd.DataFrame  # x_A, T_A, T_B, T_liquidus, T_ideal_A, T_ideal_B, T_ideal
    x_e: float
    T_e: float
    x_e_ideal: float
    T_e_ideal: float
    notes: list[str] = field(default_factory=list)
    T_e_upper: float = float("nan")  # set when the eutectic lies below the converged range

    @property
    def depression(self) -> float:
        """How far the predicted eutectic lies below the ideal one (K); > 0 = 'deep'."""
        return self.T_e_ideal - self.T_e

    def liquid_window(self, T: float) -> tuple[float, float] | None:
        """Composition range (x_A) that is fully liquid at temperature T (None if none)."""
        d = self.data
        ok = d["x_A"][(d["T_liquidus"] < T) | (d["below_range"] & (d["T_upper"] < T))]
        return (float(ok.min()), float(ok.max())) if len(ok) else None

    def summary(self) -> dict:
        return {"A": self.a.name, "B": self.b.name, "x_A_eutectic": self.x_e, "T_eutectic_K": self.T_e,
                "x_A_eutectic_ideal": self.x_e_ideal, "T_eutectic_ideal_K": self.T_e_ideal,
                "depression_K": self.depression, "T_eutectic_upper_bound_K": self.T_e_upper}


class EutecticSolver:
    def __init__(self, engine: Engine, a: Component, b: Component):
        if set(a.species) & set(b.species):
            raise ValueError("The two components must not share COSMO species")
        self.engine, self.a, self.b = engine, a, b
        self.species = list(a.species) + list(b.species)
        nA = np.array([a.species.get(s, 0.0) for s in self.species])
        nB = np.array([b.species.get(s, 0.0) for s in self.species])
        self.nu = (nA, nB)
        self.pure = (nA / nA.sum(), nB / nB.sum())

    def species_x(self, xA: float) -> np.ndarray:
        n = xA * self.nu[0] + (1 - xA) * self.nu[1]
        return n / n.sum()

    def ln_a(self, xA: float, T: float, which: int) -> float:
        """ln activity of component ``which`` (0 = A, 1 = B) at component fraction x_A."""
        x = self.species_x(xA)
        nu = self.nu[which]
        m = nu > 0
        lg = self.engine.ln_gamma(self.species, np.vstack([x, self.pure[which]]), T, refst="cosmo")
        return float(np.sum(nu[m] * (np.log(x[m]) + lg[0, m] - np.log(self.pure[which][m]) - lg[1, m])))

    def ln_a_ideal(self, xA: float, which: int) -> float:
        x, nu = self.species_x(xA), self.nu[which]
        m = nu > 0
        with np.errstate(divide="ignore"):
            return float(np.sum(nu[m] * (np.log(x[m]) - np.log(self.pure[which][m]))))

    def liquidus(self, xA: float, which: int, T_lo: float) -> tuple[float, str, float]:
        """(T, note, bound): bound is the lowest temperature checked when the liquidus lies below it."""
        c = (self.a, self.b)[which]
        xi = xA if which == 0 else 1 - xA
        if xi >= 1:
            return c.t_fus, "", float("nan")
        if xi <= 0:
            return float("nan"), "", float("nan")
        def f(T: float) -> float:
            try:
                return self.ln_a(xA, T, which) - ideal_rhs(T, c.dh_fus, c.t_fus, c.dcp_fus)
            except ValueError:  # COSMOspace not converged (typically ions at low T)
                return float("nan")

        f_hi = f(c.t_fus)
        if not np.isfinite(f_hi):
            return float("nan"), "openCOSMO-RS did not converge", float("nan")
        if f_hi >= 0:
            return float("nan"), "activity ≥ 1 at T_m (liquid–liquid demixing likely)", float("nan")
        # lowest evaluable temperature (ionic systems may not converge when very cold)
        lo, f_lo = T_lo, f(T_lo)
        while not np.isfinite(f_lo) and lo < c.t_fus - 10:
            lo = lo + 0.1 * (c.t_fus - lo) + 5
            f_lo = f(lo)
        if not np.isfinite(f_lo) or f_lo <= 0:
            return float("nan"), "below", lo
        try:
            return brentq(f, lo, c.t_fus, xtol=1e-3), "", float("nan")
        except ValueError:
            return float("nan"), "openCOSMO-RS did not converge", float("nan")

    def solve(self, n: int = 41, T_lo: float | None = None, progress=None) -> PhaseDiagram:
        T_lo = T_lo or max(50.0, 0.3 * min(self.a.t_fus, self.b.t_fus))
        edge = np.array([1e-4, 1e-3, 1e-2])  # resolve eutectics close to a pure component
        xs = np.unique(np.concatenate([np.linspace(0.0, 1.0, n), edge, 1 - edge]))
        n = len(xs)
        rows, notes = [], set()
        for i, x in enumerate(xs):
            tA, nA, bA = self.liquidus(x, 0, T_lo)
            tB, nB, bB = self.liquidus(x, 1, T_lo)
            below = not (np.isfinite(tA) or np.isfinite(tB)) and "below" in (nA, nB) \
                and all(n in ("below", "") for n in (nA, nB))
            for nt in (nA, nB):
                if nt and nt != "below":
                    notes.add(f"x_A = {x:.2f}: {nt}")
            iA = _t_from_lna(self.ln_a_ideal(x, 0), self.a) if x > 0 else float("nan")
            iB = _t_from_lna(self.ln_a_ideal(x, 1), self.b) if x < 1 else float("nan")
            rows.append({"x_A": x, "T_A": tA, "T_B": tB, "T_ideal_A": iA, "T_ideal_B": iB,
                         "below_range": below, "T_upper": float(np.nanmax([bA, bB])) if below else float("nan")})
            if progress:
                progress((i + 1) / n)
        df = pd.DataFrame(rows)
        df["T_liquidus"] = df[["T_A", "T_B"]].max(axis=1, skipna=True)
        df["T_ideal"] = df[["T_ideal_A", "T_ideal_B"]].max(axis=1, skipna=True)

        x_e, T_e, k = self._intersect(df["T_A"].to_numpy(), df["T_B"].to_numpy(), xs, T_lo, ideal=False)
        if k > 1:
            notes.add("several liquidus crossings (possible compound formation or demixing); lowest reported")
        T_up = float("nan")
        if not np.isfinite(T_e) and df["below_range"].any():
            g = df[df["below_range"]]
            j = g["T_upper"].idxmin()
            x_e, T_up = float(df.at[j, "x_A"]), float(g["T_upper"].min())
            notes.add(f"liquidus lies below the converged range for x_A = {g.x_A.min():.2f}–{g.x_A.max():.2f}: "
                      f"eutectic predicted below {T_up:.0f} K. Ionic HBAs (split ions) are known to be "
                      "over-stabilised by COSMO-RS, so treat this as 'very deep', not as a number.")
        elif not np.isfinite(T_e):
            notes.add("no eutectic found in the computed range")
        x_ei, T_ei, _ = self._intersect(df["T_ideal_A"].to_numpy(), df["T_ideal_B"].to_numpy(), xs, T_lo, ideal=True)
        return PhaseDiagram(self.a, self.b, df, x_e, T_e, x_ei, T_ei, sorted(notes), T_up)

    def _intersect(self, tA, tB, xs, T_lo, ideal: bool) -> tuple[float, float, int]:
        d = tA - tB
        ok = np.isfinite(d)
        idx = [i for i in range(len(xs) - 1) if ok[i] and ok[i + 1] and np.sign(d[i]) != np.sign(d[i + 1])]
        if not idx:
            return float("nan"), float("nan"), 0
        cands = []
        for i in idx:
            if ideal:
                g = lambda x: _t_from_lna(self.ln_a_ideal(x, 0), self.a) - _t_from_lna(self.ln_a_ideal(x, 1), self.b)  # noqa: E731
            else:
                g = lambda x: self.liquidus(x, 0, T_lo)[0] - self.liquidus(x, 1, T_lo)[0]  # noqa: E731
            def g_checked(x, g=g):
                v = g(x)
                if not np.isfinite(v):
                    raise ValueError("liquidus undefined inside bracket")
                return v

            try:
                xe = brentq(g_checked, xs[i], xs[i + 1], xtol=1e-12, rtol=1e-9)
                te = (_t_from_lna(self.ln_a_ideal(xe, 0), self.a) if ideal else self.liquidus(xe, 0, T_lo)[0])
                if not np.isfinite(te):
                    raise ValueError("eutectic temperature undefined")
            except ValueError:  # NaN inside the bracket: fall back to linear interpolation
                w = d[i] / (d[i] - d[i + 1])
                xe, te = xs[i] + w * (xs[i + 1] - xs[i]), tA[i] + w * (tA[i + 1] - tA[i])
            cands.append((te, xe))
        te, xe = min(cands)
        return float(xe), float(te), len(idx)


def phase_diagram(engine: Engine, a: str | Component, b: str | Component, n: int = 41,
                  T_lo: float | None = None, progress=None) -> PhaseDiagram:
    lib = engine.library
    a = a if isinstance(a, Component) else resolve_component(a, lib)
    b = b if isinstance(b, Component) else resolve_component(b, lib)
    return EutecticSolver(engine, a, b).solve(n=n, T_lo=T_lo, progress=progress)


def plotly_figure(pd_: PhaseDiagram, T_op: float | None = None, celsius: bool = True):
    import plotly.graph_objects as go

    from .plots import GRID, INK2, PALETTE

    off = 273.15 if celsius else 0.0
    unit = "°C" if celsius else "K"
    d = pd_.data
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=d.x_A, y=d.T_ideal - off, name="ideal (γ = 1)", mode="lines",
                             line=dict(color=INK2, dash="dash", width=1.5)))
    fig.add_trace(go.Scatter(x=d.x_A, y=d.T_A - off, name=f"liquidus of {pd_.a.name}", mode="lines",
                             line=dict(color=PALETTE[0], width=2.5)))
    fig.add_trace(go.Scatter(x=d.x_A, y=d.T_B - off, name=f"liquidus of {pd_.b.name}", mode="lines",
                             line=dict(color=PALETTE[1], width=2.5)))
    if np.isfinite(pd_.T_e):
        fig.add_trace(go.Scatter(x=[pd_.x_e], y=[pd_.T_e - off], mode="markers+text", name="eutectic",
                                 marker=dict(color=PALETTE[2], size=11, line=dict(color="white", width=1.5)),
                                 text=[f"  x = {pd_.x_e:.2f}, {pd_.T_e - off:.1f}°" + ("C" if celsius else "K")],
                                 textposition="middle right"))
    if d["below_range"].any():
        g = d[d["below_range"]]
        fig.add_vrect(x0=g.x_A.min(), x1=g.x_A.max(), fillcolor=INK2, opacity=0.07, line_width=0,
                      annotation_text=f"liquidus below {g.T_upper.max() - off:.0f}{unit} (model range)",
                      annotation_position="bottom left")
    if T_op is not None:
        fig.add_hline(y=T_op - off, line=dict(color=PALETTE[4], dash="dot"),
                      annotation_text="operating T", annotation_position="top left")
        w = pd_.liquid_window(T_op)
        if w:
            fig.add_vrect(x0=w[0], x1=w[1], fillcolor=PALETTE[4], opacity=0.08, line_width=0)
    unit = "°C" if celsius else "K"
    fig.update_layout(height=500, margin=dict(t=30, l=10, r=10, b=10), plot_bgcolor="rgba(0,0,0,0)",
                      legend=dict(orientation="h", yanchor="top", y=-0.18, x=0),
                      xaxis=dict(title=f"mole fraction of {pd_.a.name}, x_A", range=[0, 1], gridcolor=GRID),
                      yaxis=dict(title=f"Temperature ({unit})", gridcolor=GRID))
    return fig


def matplotlib_png(pd_: PhaseDiagram, T_op: float | None = None,
                   marks: list[tuple[float, str]] | None = None) -> bytes:
    import io

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from .plots import GRID, INK, INK2, PALETTE, SURF

    d, off = pd_.data, 273.15
    fig, ax = plt.subplots(figsize=(8, 5.2), facecolor=SURF)
    ax.set_facecolor(SURF)
    ax.grid(True, color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    ax.plot(d.x_A, d.T_ideal - off, "--", color=INK2, lw=1.5, label="ideal (γ = 1)")
    ax.plot(d.x_A, d.T_A - off, color=PALETTE[0], lw=2.5, label=f"liquidus of {pd_.a.name}")
    ax.plot(d.x_A, d.T_B - off, color=PALETTE[1], lw=2.5, label=f"liquidus of {pd_.b.name}")
    if d["below_range"].any():
        g = d[d["below_range"]]
        ax.axvspan(g.x_A.min(), g.x_A.max(), color=INK2, alpha=0.08, lw=0,
                   label=f"liquidus below {g.T_upper.max() - off:.0f} °C (model range)")
    if np.isfinite(pd_.T_e):
        ax.plot([pd_.x_e], [pd_.T_e - off], "o", ms=9, color=PALETTE[2], mec="white", mew=1.5, zorder=5)
        ax.annotate(f"eutectic  x = {pd_.x_e:.2f}, {pd_.T_e - off:.1f} °C", (pd_.x_e, pd_.T_e - off),
                    xytext=(10, -4), textcoords="offset points", fontsize=9, color=INK)
    if T_op is not None:
        ax.axhline(T_op - off, color=PALETTE[4], ls=":", lw=1.5, label=f"operating T ({T_op - off:.1f} °C)")
    for xm, label in marks or []:
        ax.axvline(xm, color=INK, lw=1.8, ls="-.", label=f"{label} (x = {xm:.3f})")
    ax.set_xlim(0, 1)
    ax.set_xlabel(f"mole fraction of {pd_.a.name}, x_A", color=INK2)
    ax.set_ylabel("Temperature (°C)", color=INK2)
    ax.set_title(f"{pd_.a.name} + {pd_.b.name}: solid–liquid phase diagram (openCOSMO-RS)", loc="left",
                 fontsize=11, fontweight="bold", color=INK)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        ax.spines[sp].set_color(GRID)
    ax.tick_params(colors=INK2)
    ax.legend(frameon=False, fontsize=8.5, loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=2)
    fig.tight_layout()
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=160, facecolor=SURF)
    plt.close(fig)
    return buf.getvalue()


# --------------------------------------------------------------------------- ternary (3 components)
SQ3 = np.sqrt(3.0) / 2


def simplex_grid(step: float, n: int = 3) -> np.ndarray:
    """All compositions (rows summing to 1) on a regular lattice of the 3-component simplex."""
    k = int(round(1 / step))
    pts = [(i, j, k - i - j) for i in range(k + 1) for j in range(k + 1 - i)]
    return np.array(pts, float) / k


def to_xy(Z: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Equilateral-triangle coordinates: A bottom-left, B bottom-right, C top."""
    Z = np.atleast_2d(Z)
    return Z[:, 1] + 0.5 * Z[:, 2], SQ3 * Z[:, 2]


def ratio_label(z, names=None) -> str:
    z = np.asarray(z, float)
    nz = z[z > 1e-6]
    r = nz / nz.min()
    txt = ":".join(f"{v:.2f}".rstrip("0").rstrip(".") for v in r)
    return f"{':'.join(n for n, v in zip(names, z) if v > 1e-6)} = {txt}" if names else txt


@dataclass
class TernaryDiagram:
    comps: tuple[Component, Component, Component]
    data: pd.DataFrame  # z_A, z_B, z_C, T_A, T_B, T_C, T_liquidus, below_range, T_upper, T_ideal, limiting
    z_e: np.ndarray
    T_e: float
    z_e_ideal: np.ndarray
    T_e_ideal: float
    notes: list[str] = field(default_factory=list)
    T_e_upper: float = float("nan")

    @property
    def names(self) -> list[str]:
        return [c.name for c in self.comps]

    @property
    def depression(self) -> float:
        return self.T_e_ideal - self.T_e

    def liquid(self, T: float) -> pd.DataFrame:
        """Computed compositions that are fully liquid at temperature T."""
        d = self.data
        return d[(d["T_liquidus"] < T) | (d["below_range"] & (d["T_upper"] < T))]

    def liquid_fraction(self, T: float) -> float:
        """Share of the regular composition grid (not the refinement points) that is liquid at T."""
        grid = ~self.data["refined"]
        return float(self.liquid(T)["refined"].eq(False).sum() / grid.sum())

    def summary(self) -> dict:
        out = {f"z_{c.name}_eutectic": float(z) for c, z in zip(self.comps, self.z_e)}
        out.update({"T_eutectic_K": self.T_e, "T_eutectic_ideal_K": self.T_e_ideal,
                    "depression_K": self.depression, "T_eutectic_upper_bound_K": self.T_e_upper})
        out.update({f"z_{c.name}_eutectic_ideal": float(z) for c, z in zip(self.comps, self.z_e_ideal)})
        return out


class TernarySolver:
    """Liquidus surface of a 3-component mixture (ideal immiscible solids).

    ln γ is computed exactly on a temperature grid uniform in 1/T for every composition in one
    batched openCOSMO-RS call and interpolated linearly in 1/T (the residual part is ~1/T); ln x and
    the fusion term stay exact, so each liquidus root is solved accurately without one call per root.
    """

    def __init__(self, engine: Engine, comps: list[Component], n_temps: int = 12):
        if len(comps) != 3:
            raise ValueError("A ternary diagram needs exactly three components")
        for i in range(3):
            for j in range(i + 1, 3):
                if set(comps[i].species) & set(comps[j].species):
                    raise ValueError(f"{comps[i].name} and {comps[j].name} share COSMO species")
        self.engine, self.comps, self.n_temps = engine, list(comps), n_temps
        self.species = [s for c in comps for s in c.species]
        self.nu = np.array([[c.species.get(s, 0.0) for s in self.species] for c in comps])
        self.pure = self.nu / self.nu.sum(axis=1, keepdims=True)

    def species_x(self, Z: np.ndarray) -> np.ndarray:
        n = Z @ self.nu
        return n / n.sum(axis=1, keepdims=True)

    def _lng(self, X: np.ndarray, T: np.ndarray) -> np.ndarray:
        try:
            return self.engine.ln_gamma(self.species, X, T, refst="cosmo")
        except Exception:  # noqa: BLE001 - COSMOspace not converged: isolate the failing rows
            if len(X) == 1:
                return np.full((1, len(self.species)), np.nan)
            h = len(X) // 2
            return np.vstack([self._lng(X[:h], T[:h]), self._lng(X[h:], T[h:])])

    def temperatures(self, T_lo: float) -> np.ndarray:
        T_hi = max(c.t_fus for c in self.comps)
        return 1.0 / np.linspace(1.0 / T_hi, 1.0 / T_lo, self.n_temps)  # descending T

    def evaluate(self, Z: np.ndarray, T_lo: float, progress=None) -> pd.DataFrame:
        Z = np.atleast_2d(np.asarray(Z, float))
        temps = self.temperatures(T_lo)
        nZ, nT, nS = len(Z), len(temps), len(self.species)
        Xs = self.species_x(Z)
        X = np.vstack([np.vstack([Xs, self.pure]) for _ in temps])
        Tr = np.repeat(temps, nZ + 3)
        lg = np.empty((len(X), nS))
        chunk = 200
        for a in range(0, len(X), chunk):
            lg[a:a + chunk] = self._lng(X[a:a + chunk], Tr[a:a + chunk])
            if progress:
                progress(min(1.0, (a + chunk) / len(X)))
        lg = lg.reshape(nT, nZ + 3, nS)
        lg_mix, lg_pure = lg[:, :nZ], lg[:, nZ:]
        inv = 1.0 / temps  # ascending
        cols: dict[str, np.ndarray] = {}
        below = np.zeros((nZ, 3), bool)
        bound = np.full((nZ, 3), np.nan)
        notes: set[str] = set()
        for i, c in enumerate(self.comps):
            nu = self.nu[i]
            m = nu > 0
            Ti = np.full(nZ, np.nan)
            for k in range(nZ):
                zi = Z[k, i]
                if zi <= 1e-12:
                    continue
                if zi >= 1 - 1e-12:
                    Ti[k] = c.t_fus
                    continue
                with np.errstate(divide="ignore", invalid="ignore"):
                    dlg = lg_mix[:, k, m] - lg_pure[:, i, m]  # (nT, n_species_of_i)
                    base = np.sum(nu[m] * (np.log(Xs[k, m]) - np.log(self.pure[i, m])))
                fg = base + dlg @ nu[m] - np.array([ideal_rhs(t, c.dh_fus, c.t_fus, c.dcp_fus) for t in temps])
                ok = np.isfinite(fg)
                if not ok.any():
                    notes.add("openCOSMO-RS did not converge for some compositions")
                    continue
                j0 = int(np.argmax(ok))
                if fg[j0] >= 0:
                    notes.add(f"{c.name}: activity ≥ 1 near its melting point (liquid–liquid demixing likely)")
                    continue
                cross = [j for j in range(j0, nT - 1) if ok[j] and ok[j + 1] and fg[j] < 0 <= fg[j + 1]]
                if not cross:
                    last = np.where(ok)[0][-1]
                    below[k, i], bound[k, i] = True, temps[last]
                    continue
                j = cross[0]

                def f(T, j=j, k=k, m=m, nu=nu[m], base=base, c=c):
                    w = (1.0 / T - inv[j]) / (inv[j + 1] - inv[j])
                    d = (1 - w) * (lg_mix[j, k, m] - lg_pure[j, i, m]) + w * (lg_mix[j + 1, k, m] - lg_pure[j + 1, i, m])
                    return base + d @ nu - ideal_rhs(T, c.dh_fus, c.t_fus, c.dcp_fus)

                Ti[k] = brentq(f, temps[j + 1], temps[j], xtol=1e-3)
            cols[f"T_{'ABC'[i]}"] = Ti
        d = pd.DataFrame({"z_A": Z[:, 0], "z_B": Z[:, 1], "z_C": Z[:, 2], **cols})
        T3 = d[["T_A", "T_B", "T_C"]].to_numpy()
        present = Z > 1e-12
        resolved = np.isfinite(T3) | below | ~present
        d["T_liquidus"] = np.where(resolved.all(axis=1) & np.isfinite(T3).any(axis=1), np.nanmax(
            np.where(np.isfinite(T3), T3, -np.inf), axis=1), np.nan)
        d["T_liquidus"] = d["T_liquidus"].replace(-np.inf, np.nan)
        all_below = (below | ~present).all(axis=1)
        d["below_range"] = all_below
        with np.errstate(all="ignore"):
            d["T_upper"] = np.where(all_below, np.nanmax(np.where(below, bound, -np.inf), axis=1), np.nan)
        lim = np.nanargmax(np.where(np.isfinite(T3), T3, -np.inf), axis=1)
        d["limiting"] = [self.comps[j].name if np.isfinite(d.at[k, "T_liquidus"]) else ""
                         for k, j in enumerate(lim)]
        d["T_ideal"] = self.ideal_liquidus(Z)
        d.attrs["notes"] = sorted(notes)
        return d

    def ideal_liquidus(self, Z: np.ndarray) -> np.ndarray:
        Xs = self.species_x(np.atleast_2d(Z))
        out = np.full(len(Xs), -np.inf)
        for i, c in enumerate(self.comps):
            m = self.nu[i] > 0
            with np.errstate(divide="ignore"):
                lna = np.sum(self.nu[i, m] * (np.log(Xs[:, m]) - np.log(self.pure[i, m])), axis=1)
            if c.dcp_fus:
                Ti = np.array([_t_from_lna(v, c) if np.isfinite(v) else np.nan for v in lna])
            else:
                with np.errstate(divide="ignore", invalid="ignore"):
                    Ti = np.where(lna >= 0, c.t_fus, 1.0 / (1.0 / c.t_fus - R * lna / c.dh_fus))
                Ti = np.where(np.isfinite(lna), Ti, np.nan)
            out = np.fmax(out, Ti)
        return np.where(np.isfinite(out), out, np.nan)

    @staticmethod
    def _local(z: np.ndarray, h: float) -> np.ndarray:
        """Small lattice of spacing h/4 around z (inside the simplex)."""
        g = np.arange(-4, 5) * h / 4
        pts = [z + np.array([a, b, -a - b]) for a in g for b in g]
        P = np.array([p for p in pts if p.min() >= -1e-12])
        return np.clip(P, 0, 1) / np.clip(P, 0, 1).sum(axis=1, keepdims=True)

    def solve(self, step: float = 0.1, T_lo: float | None = None, refine: int = 2, progress=None) -> TernaryDiagram:
        # NADES below ~173 K are irrelevant, and that is where COSMO-RS artefacts (ions) dominate
        T_lo = T_lo or max(173.15, 0.5 * min(c.t_fus for c in self.comps))
        Z = simplex_grid(step)
        n_ref = refine * 60
        tot = len(Z) + n_ref
        prog = (lambda f: progress(f * len(Z) / tot)) if progress else None
        d = self.evaluate(Z, T_lo, prog)
        d["refined"] = False
        notes = set(d.attrs.get("notes", []))
        # refine around the lowest computed liquidus point
        h, done = step, len(Z)
        for _ in range(refine):
            liq = d["T_liquidus"]
            if not liq.notna().any():
                break
            z0 = d.loc[liq.idxmin(), ["z_A", "z_B", "z_C"]].to_numpy(float)
            P = self._local(z0, h)
            seen = {tuple(np.round(r, 6)) for r in d[["z_A", "z_B", "z_C"]].to_numpy()}
            P = np.array([p for p in P if tuple(np.round(p, 6)) not in seen])
            if len(P):
                p0 = done
                sub = self.evaluate(P, T_lo, (lambda f, p0=p0, n=len(P): progress((p0 + f * n) / tot))
                                    if progress else None)
                notes |= set(sub.attrs.get("notes", []))
                sub["refined"] = True
                d = pd.concat([d, sub], ignore_index=True)
                done += len(P)
            h /= 4
        if progress:
            progress(1.0)
        d.attrs = {}
        liq = d["T_liquidus"]
        T_up = float("nan")
        if d["below_range"].any() and (not liq.notna().any() or liq.min() < T_lo + 5):
            # the minimum is below the checked range: report the middle of that region as an upper bound
            g = d[d["below_range"]]
            Zg = g[["z_A", "z_B", "z_C"]].to_numpy(float)
            k = g.index[int(np.argmin(((Zg - Zg.mean(axis=0)) ** 2).sum(axis=1)))]
            z_e, T_e, T_up = Zg[list(g.index).index(k)], float("nan"), float(g["T_upper"].max())
            notes.add(f"the liquidus drops below the lowest temperature checked ({T_up:.0f} K) for "
                      f"{len(g)} of {len(d)} compositions (shaded region): the eutectic is only known to be "
                      f"below {T_up - 273.15:.0f} °C.")
        elif liq.notna().any():
            k = liq.idxmin()
            z_e, T_e = d.loc[k, ["z_A", "z_B", "z_C"]].to_numpy(float), float(liq[k])
            if d["below_range"].any() and d.loc[d.below_range, "T_upper"].min() < T_e:
                notes.add("part of the diagram lies below the converged temperature range; the true eutectic "
                          "may be lower than reported")
        elif d["below_range"].any():
            g = d[d["below_range"]]
            k = g["T_upper"].idxmin()
            z_e, T_e, T_up = d.loc[k, ["z_A", "z_B", "z_C"]].to_numpy(float), float("nan"), float(g["T_upper"].min())
            notes.add(f"the liquidus lies below the converged range: eutectic predicted below {T_up:.0f} K. "
                      "Ionic HBAs are over-stabilised by COSMO-RS – treat as 'very deep'.")
        else:
            z_e, T_e = np.full(3, np.nan), float("nan")
            notes.add("no liquidus found in the computed range")
        # ideal eutectic on a fine lattice (cheap, no COSMO calls)
        F = simplex_grid(0.005)
        ti = self.ideal_liquidus(F)
        j = int(np.nanargmin(ti))
        if (d["limiting"].replace("", np.nan).dropna().nunique() < 3) and np.isfinite(T_e):
            notes.add("the minimum lies where fewer than three liquidus surfaces meet (a binary eutectic on an "
                      "edge or a boundary), so the third component does not lower the melting point further")
        if any(len(c.species) > 1 for c in self.comps):
            notes.add("salts are split into ions, which COSMO-RS over-stabilises: salt-rich regions come out far "
                      "too liquid (too deep). Compare compositions relative to each other and check the melting "
                      "point of the chosen mixture experimentally.")
        notes.add(f"lowest temperature checked: {T_lo:.0f} K ({T_lo - 273.15:.0f} °C)")
        notes.add(f"eutectic located to within ±{step / 4 ** refine:.3f} in mole fraction")
        return TernaryDiagram(tuple(self.comps), d.sort_values(["z_A", "z_B"]).reset_index(drop=True),
                              z_e, T_e, F[j], float(ti[j]), sorted(notes), T_up)


def ternary_diagram(engine: Engine, comps: list[str | Component], step: float = 0.1,
                    T_lo: float | None = None, progress=None) -> TernaryDiagram:
    lib = engine.library
    cs = [c if isinstance(c, Component) else resolve_component(c, lib) for c in comps]
    return TernarySolver(engine, cs).solve(step=step, T_lo=T_lo, progress=progress)


def ternary_cost(step: float, has_ions: bool, n_temps: int = 12, refine: int = 2) -> float:
    """Rough run time (s) of a ternary diagram."""
    n = len(simplex_grid(step)) + 60 * refine
    return (n + 3) * n_temps * (0.009 if has_ions else 0.005)


def _ternary_field(res: TernaryDiagram, n: int = 120):
    """Liquidus interpolated onto a fine cartesian grid (NaN outside the triangle)."""
    import matplotlib.tri as mtri

    d = res.data
    T = d["T_liquidus"].to_numpy(float)
    T = np.where(np.isfinite(T), T, np.where(d["below_range"], d["T_upper"], np.nan))
    x, y = to_xy(d[["z_A", "z_B", "z_C"]].to_numpy(float))
    ok = np.isfinite(T)
    tri = mtri.Triangulation(x[ok], y[ok])
    gx, gy = np.linspace(0, 1, n), np.linspace(0, SQ3, n)
    GX, GY = np.meshgrid(gx, gy)
    G = np.asarray(mtri.LinearTriInterpolator(tri, T[ok])(GX, GY).filled(np.nan))
    inside = (GY <= SQ3 * 2 * GX + 1e-9) & (GY <= SQ3 * 2 * (1 - GX) + 1e-9)
    return gx, gy, np.where(inside, G, np.nan)


def ternary_figure(res: TernaryDiagram, T_op: float | None = None, celsius: bool = True,
                   marks: list[tuple[np.ndarray, str]] | None = None):
    import plotly.graph_objects as go

    from .plots import GRID, INK2, PALETTE

    off = 273.15 if celsius else 0.0
    unit = "°C" if celsius else "K"
    gx, gy, G = _ternary_field(res)
    fig = go.Figure()
    fig.add_trace(go.Contour(x=gx, y=gy, z=G - off, colorscale="RdYlBu_r", ncontours=16, connectgaps=False,
                             contours=dict(showlabels=True, labelfont=dict(size=10, color="white")),
                             line=dict(width=0.5, color="rgba(255,255,255,0.6)"), hoverinfo="skip",
                             colorbar=dict(title=f"liquidus ({unit})", thickness=14)))
    if T_op is not None and np.nanmin(G) < T_op < np.nanmax(G):
        fig.add_trace(go.Contour(x=gx, y=gy, z=G - off, showscale=False, hoverinfo="skip",
                                 contours=dict(coloring="none", start=T_op - off, end=T_op - off, size=1,
                                               showlabels=False),
                                 line=dict(width=3, color="black", dash="dot"), name="operating T"))
    d = res.data
    x, y = to_xy(d[["z_A", "z_B", "z_C"]].to_numpy(float))
    a, b, c = res.names
    tl = [("< %.0f" % (u - off)) if br else ("%.1f" % (t - off)) for t, br, u in
          zip(d["T_liquidus"], d["below_range"], d["T_upper"])]
    fig.add_trace(go.Scatter(x=x, y=y, mode="markers", marker=dict(size=4, color="rgba(0,0,0,0.25)"),
                             showlegend=False, customdata=np.c_[d.z_A, d.z_B, d.z_C],
                             text=[f"{tt} {unit} · limiting: {lim or '—'}" for tt, lim in zip(tl, d["limiting"])],
                             hovertemplate=f"{a} %{{customdata[0]:.3f}}<br>{b} %{{customdata[1]:.3f}}<br>"
                                           f"{c} %{{customdata[2]:.3f}}<br>%{{text}}<extra></extra>"))
    # triangle, grid lines and vertex labels
    for t in np.arange(0.2, 1.0, 0.2):
        for p, q in (([0, t, 1 - t], [0, 1 - t, t]), ([t, 0, 1 - t], [1 - t, 0, t]), ([t, 1 - t, 0], [1 - t, t, 0])):
            px, py = to_xy(np.array([p, q]))
            fig.add_trace(go.Scatter(x=px, y=py, mode="lines", line=dict(color=GRID, width=0.8),
                                     hoverinfo="skip", showlegend=False))
    fig.add_trace(go.Scatter(x=[0, 1, 0.5, 0], y=[0, 0, SQ3, 0], mode="lines", line=dict(color=INK2, width=1.5),
                             hoverinfo="skip", showlegend=False))
    for (vx, vy), name, pos in zip(((0, 0), (1, 0), (0.5, SQ3)), res.names,
                                   ("bottom left", "bottom right", "top center")):
        fig.add_annotation(x=vx, y=vy, text=f"<b>{name}</b>", showarrow=False,
                           yshift=14 if pos.startswith("top") else -16, xshift={"bottom left": -10,
                                                                               "bottom right": 10}.get(pos, 0))
    if np.all(np.isfinite(res.z_e)):
        ex, ey = to_xy(res.z_e)
        lbl = f"{res.T_e - off:.1f} {unit}" if np.isfinite(res.T_e) else f"< {res.T_e_upper - off:.0f} {unit}"
        fig.add_trace(go.Scatter(x=ex, y=ey, mode="markers+text", name="eutectic", text=[f"  eutectic {lbl}"],
                                 textposition="middle right",
                                 marker=dict(color=PALETTE[2], size=13, symbol="star",
                                             line=dict(color="white", width=1.5))))
    for z, label in marks or []:
        mx, my = to_xy(np.asarray(z, float))
        fig.add_trace(go.Scatter(x=mx, y=my, mode="markers+text", name=label, text=[label],
                                 textposition="top center",
                                 marker=dict(color=PALETTE[0], size=13, symbol="diamond",
                                             line=dict(color="white", width=1.5))))
    ix, iy = to_xy(res.z_e_ideal)
    fig.add_trace(go.Scatter(x=ix, y=iy, mode="markers", name=f"ideal eutectic {res.T_e_ideal - off:.1f} {unit}",
                             marker=dict(color="white", size=10, symbol="star-open", line=dict(color=INK2, width=2))))
    fig.update_layout(height=560, margin=dict(t=30, l=10, r=10, b=30), plot_bgcolor="rgba(0,0,0,0)",
                      legend=dict(orientation="h", yanchor="top", y=-0.02, x=0),
                      xaxis=dict(visible=False, range=[-0.08, 1.08]),
                      yaxis=dict(visible=False, range=[-0.08, SQ3 + 0.08], scaleanchor="x"))
    return fig


def ternary_png(res: TernaryDiagram, T_op: float | None = None,
                marks: list[tuple[np.ndarray, str]] | None = None) -> bytes:
    import io

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from .plots import INK, INK2, PALETTE, SURF

    gx, gy, G = _ternary_field(res, 200)
    off = 273.15
    fig, ax = plt.subplots(figsize=(7.5, 6.8), facecolor=SURF)
    ax.set_facecolor(SURF)
    cs = ax.contourf(gx, gy, G - off, levels=16, cmap="RdYlBu_r")
    ax.contour(gx, gy, G - off, levels=cs.levels, colors="white", linewidths=0.4)
    fig.colorbar(cs, ax=ax, shrink=0.75, label="liquidus (°C)")
    if T_op is not None and np.nanmin(G) < T_op < np.nanmax(G):
        ax.contour(gx, gy, G, levels=[T_op], colors="k", linewidths=2, linestyles=":")
    ax.plot([0, 1, 0.5, 0], [0, 0, SQ3, 0], color=INK2, lw=1.5)
    for (vx, vy), name, ha, dy in zip(((0, 0), (1, 0), (0.5, SQ3)), res.names, ("right", "left", "center"),
                                      (-0.05, -0.05, 0.03)):
        ax.text(vx, vy + dy, name, ha=ha, fontsize=10, fontweight="bold", color=INK)
    if np.all(np.isfinite(res.z_e)):
        ex, ey = to_xy(res.z_e)
        ax.plot(ex, ey, "*", ms=16, color=PALETTE[2], mec="white", mew=1.2)
        lbl = f"{res.T_e - off:.1f} °C" if np.isfinite(res.T_e) else f"< {res.T_e_upper - off:.0f} °C"
        ax.annotate(f"eutectic {lbl}", (ex[0], ey[0]), xytext=(10, 4), textcoords="offset points", fontsize=9)
    for z, label in marks or []:
        mx, my = to_xy(np.asarray(z, float))
        ax.plot(mx, my, "D", ms=10, color=PALETTE[0], mec="white", mew=1.2, zorder=6)
        ax.annotate(label, (mx[0], my[0]), xytext=(0, 9), textcoords="offset points", ha="center", fontsize=9,
                    color=INK, bbox=dict(boxstyle="round,pad=0.2", fc=SURF, ec="none", alpha=0.8))
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(f"{' + '.join(res.names)}: liquidus (openCOSMO-RS)", loc="left", fontsize=11,
                 fontweight="bold", color=INK, pad=28)
    fig.tight_layout()
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=160, facecolor=SURF)
    plt.close(fig)
    return buf.getvalue()
