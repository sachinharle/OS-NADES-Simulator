"""Solubility-vs-temperature curves grouped by solvent class (shared by the web app, CLI and scripts)."""
from __future__ import annotations

import io
from collections import OrderedDict

import numpy as np
import pandas as pd

from .thermo import Engine, SolventSystem, ideal_rhs, presets

# Validated categorical palette (fixed order, max 8 per panel) and chart ink.
PALETTE = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
INK, INK2, GRID, SURF = "#1f1f1e", "#5f5e58", "#e6e5e0", "#fcfcfb"
MAX_PER_PANEL = len(PALETTE)
OTHER = "Custom / other"


def group_systems(names: list[str]) -> "OrderedDict[str, list[str]]":
    """Order systems into the preset solvent classes; unknown ones go to 'Custom / other'.
    Groups larger than the palette are split into numbered panels."""
    groups = presets().get("solvent_groups", {})
    out: OrderedDict[str, list[str]] = OrderedDict()
    placed = set()
    for g, members in groups.items():
        sel = [m for m in members if m in names]
        if sel:
            out[g] = sel
            placed.update(sel)
    rest = [n for n in names if n not in placed]
    if rest:
        out[OTHER] = rest
    final: OrderedDict[str, list[str]] = OrderedDict()
    for g, members in out.items():
        chunks = [members[i : i + MAX_PER_PANEL] for i in range(0, len(members), MAX_PER_PANEL)]
        for k, c in enumerate(chunks):
            final[g if len(chunks) == 1 else f"{g} ({k + 1})"] = c
    return final


def compute_curves(engine: Engine, solute: str, systems: list[SolventSystem], temps, method: str = "iterative",
                   progress=None, **fus) -> pd.DataFrame:
    groups = group_systems([s.name for s in systems])
    by_name = {s.name: s for s in systems}
    frames, n = [], sum(len(v) for v in groups.values())
    i = 0
    for g, names in groups.items():
        for name in names:
            df = engine.solubility_curve(solute, by_name[name], np.asarray(temps, float), method=method, **fus)
            df["group"] = g
            frames.append(df)
            i += 1
            if progress:
                progress(i / n, name)
    return pd.concat(frames, ignore_index=True)


def _use_log(x: pd.Series) -> bool:
    x = x[(x > 0) & np.isfinite(x)]
    return len(x) > 0 and x.max() / x.min() > 4


def plotly_panels(df: pd.DataFrame, dh_fus: float, t_fus: float, dcp_fus: float = 0.0, celsius: bool = True):
    """One plotly figure per solvent group (each with its own legend and y-scale)."""
    import plotly.graph_objects as go

    temps = np.sort(df["T"].unique())
    tx = temps - 273.15 if celsius else temps
    x_ideal = np.exp([ideal_rhs(T, dh_fus, t_fus, dcp_fus) for T in temps])
    figs = []
    for g, d in df.groupby("group", sort=False):
        log = _use_log(d["x"])
        fig = go.Figure()
        lo, hi = d["x"].min(), d["x"].max()
        ideal_visible = log or (x_ideal.max() >= lo * 0.8 and x_ideal.min() <= hi * 1.2)
        if ideal_visible:
            fig.add_trace(go.Scatter(x=tx, y=x_ideal, name="ideal solution", mode="lines",
                                     line=dict(color=INK2, width=1.5, dash="dash"),
                                     hovertemplate="ideal: x = %{y:.4g}<extra></extra>"))
        for i, (name, s) in enumerate(d.groupby("solvent", sort=False)):
            s = s.sort_values("T")
            t = s["T"] - 273.15 if celsius else s["T"]
            fig.add_trace(go.Scatter(
                x=t, y=s["x"], name=name, mode="lines+markers",
                line=dict(color=PALETTE[i % MAX_PER_PANEL], width=2),
                marker=dict(size=7, line=dict(color=SURF, width=1.5)),
                customdata=np.column_stack([s["g_per_kg"], s["ln_gamma"]]),
                hovertemplate=f"<b>{name}</b><br>T = %{{x:.1f}}<br>x = %{{y:.4g}}<br>"
                              "%{customdata[0]:.4g} g/kg · ln γ = %{customdata[1]:.3f}<extra></extra>",
            ))
        fig.update_layout(
            title=dict(text=g, x=0, font=dict(size=14)), height=500, margin=dict(t=50, l=10, r=10, b=10),
            hovermode="x unified", plot_bgcolor="rgba(0,0,0,0)",
            legend=dict(orientation="h", yanchor="top", y=-0.3, x=0, font=dict(size=11)),
            xaxis=dict(title="Temperature (°C)" if celsius else "Temperature (K)", gridcolor=GRID),
            yaxis=dict(title=f"mole fraction x ({'log' if log else 'linear'})", type="log" if log else "linear",
                       gridcolor=GRID, exponentformat="power"),
        )
        figs.append((g, fig, log, ideal_visible))
    return figs


def matplotlib_png(df: pd.DataFrame, solute: str, dh_fus: float, t_fus: float, dcp_fus: float = 0.0,
                   subtitle: str = "") -> bytes:
    """Static small-multiples PNG (2 columns, legends outside each panel)."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    temps = np.sort(df["T"].unique())
    x_ideal = np.exp([ideal_rhs(T, dh_fus, t_fus, dcp_fus) for T in temps])
    groups = list(df.groupby("group", sort=False))
    ncols = 2 if len(groups) > 1 else 1
    nrows = int(np.ceil(len(groups) / ncols))
    plt.rcParams.update({"font.size": 10, "axes.edgecolor": GRID, "axes.labelcolor": INK2,
                         "xtick.color": INK2, "ytick.color": INK2, "text.color": INK})
    fig, axes = plt.subplots(nrows, ncols, figsize=(8.5 * ncols, 4.6 * nrows + 1), squeeze=False, facecolor=SURF)
    for ax, (g, d) in zip(axes.flat, groups):
        ax.set_facecolor(SURF)
        ax.grid(True, color=GRID, lw=0.8)
        ax.set_axisbelow(True)
        log = _use_log(d["x"])
        lo, hi = d["x"].min(), d["x"].max()
        if log or (x_ideal.max() >= lo * 0.8 and x_ideal.min() <= hi * 1.2):
            ax.plot(temps - 273.15, x_ideal, ls="--", lw=1.5, color=INK2, label="ideal solution")
        for i, (name, s) in enumerate(d.groupby("solvent", sort=False)):
            s = s.sort_values("T")
            ax.plot(s["T"] - 273.15, s["x"], lw=2, color=PALETTE[i % MAX_PER_PANEL], label=name)
            ax.plot(s["T"] - 273.15, s["x"], "o", ms=4, color=PALETTE[i % MAX_PER_PANEL], mec=SURF, mew=1)
        if log:
            ax.set_yscale("log")
        else:
            pad = 0.08 * (hi - lo or hi)
            ax.set_ylim(max(0, lo - pad), hi + pad)
        ax.set_title(g, loc="left", fontsize=11, fontweight="bold", color=INK)
        ax.set_xlabel("Temperature (°C)")
        ax.set_ylabel(f"mole fraction x ({'log' if log else 'linear'})")
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.legend(fontsize=8.5, frameon=False, loc="center left", bbox_to_anchor=(1.01, 0.5))
    for ax in axes.flat[len(groups):]:
        ax.set_visible(False)
    fig.suptitle(f"{solute.capitalize()} solubility vs temperature", x=0.03, ha="left", fontsize=14, fontweight="bold")
    sub = (f"openCOSMO-RS · ideal-solid SLE, ΔH_fus = {dh_fus / 1000:.1f} kJ/mol, T_m = {t_fus:.1f} K"
           f" · dashed = ideal solubility · y-axes differ per panel")
    fig.text(0.03, 1 - 0.62 / fig.get_figheight(), sub + (f"\n{subtitle}" if subtitle else ""), fontsize=9, color=INK2)
    fig.tight_layout(rect=(0.01, 0, 1, 1 - 1.0 / fig.get_figheight()))
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=160, facecolor=SURF)
    plt.close(fig)
    return buf.getvalue()
