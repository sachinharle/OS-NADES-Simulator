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


def matplotlib_png(pd_: PhaseDiagram, T_op: float | None = None) -> bytes:
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
