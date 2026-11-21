"""Statistics, figures and REPORT.md for the accuracy benchmark (see benchmarks/accuracy.py)."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).parent
DATA = HERE / "data"
FIG = HERE / "figures"
CATS = ["pure organic solvent", "water", "organic mixtures", "aqueous mixtures", "NADES / ChCl mixtures",
        "binary SLE (both solid-formers)"]


def load() -> pd.DataFrame:
    df = pd.read_csv(DATA / "predictions.csv")
    df = df[np.isfinite(df["x_pred"]) & (df["x_pred"] > 0) & (df["x_exp"] > 0)].copy()
    df["err"] = np.log10(df["x_pred"] / df["x_exp"])
    df["err_ideal"] = np.log10(df["x_ideal"] / df["x_exp"])
    return df


def stats(g: pd.DataFrame) -> dict:
    e, ei = g["err"], g["err_ideal"]
    return {"points": len(g), "papers": g["doi"].nunique(), "systems": g.groupby(["solute", "solvent"]).ngroups,
            "bias": e.mean(), "AAD": e.abs().mean(), "RMSE": np.sqrt((e ** 2).mean()),
            "within_x2": (e.abs() <= np.log10(2)).mean(), "within_x10": (e.abs() <= 1).mean(),
            "AAD_ideal": ei.abs().mean(), "better_than_ideal": (e.abs() < ei.abs()).mean()}


def solvent_label(comp_json: str, solute: str) -> str:
    comp = {k: v for k, v in json.loads(comp_json).items() if k != solute}
    if len(comp) == 1:
        return next(iter(comp))
    tot = sum(comp.values())
    return " + ".join(f"{k} {v / tot:.2f}" for k, v in sorted(comp.items(), key=lambda kv: -kv[1]))


def solvent_family(comp_json: str, solute: str) -> str:
    return " + ".join(sorted(k for k in json.loads(comp_json) if k != solute))


def ranking(df: pd.DataFrame) -> pd.DataFrame:
    """Spearman rank correlation of predicted vs measured solubility across pure solvents, per solute
    and temperature (only groups with ≥ 4 solvents)."""
    from scipy.stats import spearmanr

    d = df[df["category"].isin(["pure organic solvent", "water"])].copy()
    d["T5"] = (d["T"] / 5).round() * 5
    rows = []
    for (sol, T), g in d.groupby(["solute", "T5"]):
        g = g.groupby("solvent").agg(x_exp=("x_exp", "median"), x_pred=("x_pred", "median"),
                                     x_ideal=("x_ideal", "median")).reset_index()
        if len(g) >= 4:
            rows.append({"solute": sol, "T": T, "solvents": len(g), "rho": spearmanr(g.x_exp, g.x_pred)[0],
                         "span": g.x_exp.max() / g.x_exp.min()})
    return pd.DataFrame(rows)


def fmt_table(t: pd.DataFrame) -> str:
    cols = ["points", "papers", "systems", "bias", "AAD", "RMSE", "within_x2", "within_x10", "AAD_ideal",
            "better_than_ideal"]
    head = ("| group | points | papers | systems | bias (log₁₀) | AAD (log₁₀) | RMSE (log₁₀) | within ×2 | "
            "within ×10 | ideal AAD | better than ideal |\n|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    lines = [head]
    for name, r in t[cols].iterrows():
        lines.append(f"| {name} | {r.points:.0f} | {r.papers:.0f} | {r.systems:.0f} | {r.bias:+.2f} | {r.AAD:.2f} | "
                     f"{r.RMSE:.2f} | {r.within_x2:.0%} | {r.within_x10:.0%} | {r.AAD_ideal:.2f} | "
                     f"{r.better_than_ideal:.0%} |")
    return "\n".join(lines)


def figures(df: pd.DataFrame) -> list[str]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from osnades.plots import GRID, INK, INK2, PALETTE, SURF

    FIG.mkdir(exist_ok=True)
    cats = [c for c in CATS if c in set(df["category"])]
    col = {c: PALETTE[i % len(PALETTE)] for i, c in enumerate(cats)}

    def style(ax):
        ax.set_facecolor(SURF)
        ax.grid(True, color=GRID, lw=0.7)
        ax.set_axisbelow(True)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        for s in ("left", "bottom"):
            ax.spines[s].set_color(GRID)
        ax.tick_params(colors=INK2)

    out = []
    # 1. parity plots
    fig, axes = plt.subplots(1, 2, figsize=(11, 5.4), facecolor=SURF, sharex=True, sharey=True)
    lo = np.floor(np.log10(min(df.x_exp.min(), df.x_pred.min(), df.x_ideal.min())))
    for ax, ycol, title in ((axes[0], "x_pred", "openCOSMO-RS (OS-NADES)"), (axes[1], "x_ideal", "ideal solution (γ = 1)")):
        style(ax)
        for c in cats:
            g = df[df.category == c]
            ax.scatter(np.log10(g.x_exp), np.log10(g[ycol]), s=9, alpha=0.6, color=col[c], label=c, lw=0)
        ax.plot([lo, 0], [lo, 0], color=INK, lw=1)
        for k in (1, -1):
            ax.plot([lo, 0], [lo + k, k], color=INK2, lw=0.8, ls="--")
        ax.set_xlim(lo, 0.05)
        ax.set_ylim(lo, 0.05)
        ax.set_xlabel("log₁₀ x measured", color=INK2)
        ax.set_title(title, loc="left", fontsize=11, fontweight="bold", color=INK)
    axes[0].set_ylabel("log₁₀ x predicted", color=INK2)
    axes[0].legend(frameon=False, fontsize=8, loc="upper left")
    fig.text(0.01, 0.01, "Dashed lines: a factor of 10 off.", fontsize=8, color=INK2)
    fig.tight_layout()
    p = FIG / "parity.png"
    fig.savefig(p, dpi=150, facecolor=SURF)
    plt.close(fig)
    out.append(p.name)

    # 2. error distribution by category
    fig, ax = plt.subplots(figsize=(9, 4.6), facecolor=SURF)
    style(ax)
    data = [df[df.category == c]["err"] for c in cats]
    ideal = [df[df.category == c]["err_ideal"] for c in cats]
    pos = np.arange(len(cats))
    b1 = ax.boxplot(data, positions=pos - 0.18, widths=0.32, patch_artist=True, showfliers=False)
    b2 = ax.boxplot(ideal, positions=pos + 0.18, widths=0.32, patch_artist=True, showfliers=False)
    for b, fc in ((b1, PALETTE[0]), (b2, "#c9c7c0")):
        for patch in b["boxes"]:
            patch.set(facecolor=fc, edgecolor=INK2, alpha=0.9)
        for k in ("medians", "whiskers", "caps"):
            for line in b[k]:
                line.set(color=INK2)
    ax.axhline(0, color=INK, lw=1)
    ax.set_xticks(pos, [c.replace(" (", "\n(").replace(" / ", " /\n") for c in cats], fontsize=8.5)
    ax.set_ylabel("log₁₀ (x predicted / x measured)", color=INK2)
    ax.set_title("Error by solvent type: openCOSMO-RS (blue) vs ideal solution (grey)", loc="left", fontsize=11,
                 fontweight="bold", color=INK)
    fig.tight_layout()
    p = FIG / "error_by_category.png"
    fig.savefig(p, dpi=150, facecolor=SURF)
    plt.close(fig)
    out.append(p.name)

    # 3. bias per solute
    g = df.groupby("solute").agg(bias=("err", "mean"), aad=("err", lambda e: e.abs().mean()),
                                 n=("err", "size")).sort_values("bias")
    fig, ax = plt.subplots(figsize=(8, 0.42 * len(g) + 1.2), facecolor=SURF)
    style(ax)
    ax.barh(g.index, g.bias, color=[PALETTE[0] if b >= 0 else PALETTE[1] for b in g.bias])
    for y, (b, n) in enumerate(zip(g.bias, g.n)):
        ax.text(b + 0.03 if b >= 0 else 0.03, y, f"n={n}", va="center", ha="left", fontsize=8, color=INK2)
    ax.axvline(0, color=INK, lw=1)
    ax.set_xlabel("mean log₁₀ (x predicted / x measured)   ← under-predicts · over-predicts →", color=INK2)
    ax.set_title("Bias per solute", loc="left", fontsize=11, fontweight="bold", color=INK)
    fig.tight_layout()
    p = FIG / "bias_by_solute.png"
    fig.savefig(p, dpi=150, facecolor=SURF)
    plt.close(fig)
    out.append(p.name)
    return out


def cmd_report():
    df = load()
    df["solvent"] = [solvent_label(c, s) for c, s in zip(df.composition, df.solute)]
    df["family"] = [solvent_family(c, s) for c, s in zip(df.composition, df.solute)]
    info = json.loads((DATA / "harvest_info.json").read_text())
    overall = pd.DataFrame([stats(df)], index=["all points"])
    by_cat = pd.DataFrame({c: stats(df[df.category == c]) for c in CATS if (df.category == c).any()}).T
    by_sol = pd.DataFrame({s: stats(g) for s, g in df.groupby("solute")}).T.sort_values("AAD")
    rk = ranking(df)
    # remove each solute's constant offset (mostly fusion-data error, common to all its solvents)
    df["err_centered"] = df["err"] - df.groupby("solute")["err"].transform("mean")
    aad_c = df["err_centered"].abs().mean()
    pure_pts = df[df["category"].isin(["pure organic solvent", "water"])]
    per_solvent = pure_pts.groupby("solvent").agg(solutes=("solute", "nunique"), points=("err", "size"),
                                                  err=("err_centered", "mean")).sort_values("err", ascending=False)
    from benchmarks.accuracy import _lib

    lib = _lib()
    fus = {s: lib.get(s) for s in df.solute.unique()}
    sys_err = df.groupby(["solute", "family"]).agg(points=("err", "size"), bias=("err", "mean"),
                                                   T_min=("T", "min"), T_max=("T", "max"),
                                                   ref=("doi", "first")).reset_index()
    worst = sys_err.reindex(sys_err["bias"].abs().sort_values(ascending=False).index).head(12)
    best = sys_err[sys_err.points >= 3].reindex(sys_err["bias"].abs().sort_values().index).head(8)
    figs = figures(df)
    excluded = pd.read_csv(DATA / "predictions.csv")
    n_fail = int((~np.isfinite(excluded["x_pred"])).sum())
    refs = df.groupby("doi").agg(title=("title", "first"), year=("year", "first"), points=("x_exp", "size"),
                                 solutes=("solute", lambda s: ", ".join(sorted(set(s))))).reset_index()
    o = overall.iloc[0]

    md = [f"# OS-NADES accuracy benchmark\n",
          "Predicted solid solubilities (openCOSMO-RS 24a, library σ-surfaces and fusion data, iterative SLE) "
          "compared with **every** experimental solubility data set in the NIST/TRC ThermoML Archive for the "
          "library's solid compounds whose solvent components are all in the library. Nothing was hand-picked "
          "or tuned.\n",
          "## Headline\n",
          f"* **{int(o.points)} data points** from **{int(o.papers)} papers**, {int(o.systems)} solute–solvent "
          f"systems, {df.solute.nunique()} solutes, {df['T'].min():.0f}–{df['T'].max():.0f} K.",
          f"* Typical error (mean |log₁₀ x_pred/x_exp|): **{o.AAD:.2f}**, i.e. a factor of "
          f"**{10 ** o.AAD:.1f}**; RMSE {o.RMSE:.2f} log units; mean bias {o.bias:+.2f}.",
          f"* **{o.within_x2:.0%}** of points are within a factor of 2, **{o.within_x10:.0%}** within a factor of 10.",
          f"* The ideal-solution model (γ = 1, same fusion data) has AAD {o.AAD_ideal:.2f}; openCOSMO-RS is closer "
          f"than ideal for **{o.better_than_ideal:.0%}** of points.",
          f"* Removing each solute's constant offset (an error in its fusion data shifts all its solvents "
          f"equally) leaves AAD **{aad_c:.2f}** (factor {10 ** aad_c:.1f}): this is the error of the "
          "activity-coefficient model itself, i.e. how well it ranks and scales solvents for one solute.",
          (f"* Ranking of pure solvents for a given solute (what screening needs): median Spearman ρ = "
           f"**{rk.rho.median():.2f}** over {len(rk)} solute/temperature groups with ≥ 4 solvents; "
           f"**{rk[rk.span >= 3].rho.median():.2f}** for the {int((rk.span >= 3).sum())} groups whose measured "
           "solubilities span at least a factor of 3. Solvents within a factor of ~1.5 of each other (e.g. "
           "the C1–C4 alcohols) are often mis-ordered."
           if len(rk) else "* Not enough pure-solvent data to assess solvent ranking."),
          "",
          "",
          "## What is and is not covered\n",
          "* **No deep eutectic solvent was benchmarked.** ThermoML holds no solubility data of the library's "
          "compounds in choline chloride (or other) NADES: the one DES paper it lists for paracetamol "
          "(doi:10.1016/j.fluid.2018.01.017) is archived with its water data only. NADES predictions remain "
          "unvalidated by this benchmark; the known COSMO-RS over-stabilisation of salts (see README) applies.",
          "* No binary eutectic (two solid-formers) data matched: the measured partners of thymol, menthol and "
          "decanoic acid in ThermoML are not in the library.",
          "* Coverage is limited by the library: data sets were skipped when a component was missing. Most "
          "frequent missing components (number of papers): "
          + ", ".join(f"{k} ({v})" for k, v in list(info.get("missing_components", {}).items())[:10]) + ".",
          "* Data per solute are uneven ("
          + ", ".join(f"{k} {v}" for k, v in df.solute.value_counts().items())
          + " points); a few mixed-solvent studies contribute most points for betaine, malic and citric acid.",
          "",
          "## By solvent type\n", fmt_table(by_cat), "",
          "## By solute\n", fmt_table(by_sol), "",
          "Log₁₀ errors: AAD 0.30 = factor 2, 1.0 = factor 10. *Bias* > 0: over-predicts solubility. "
          "*Better than ideal*: share of points where openCOSMO-RS is closer to experiment than γ = 1.\n",
          "### Fusion data used (they set each solute's offset)\n",
          "| solute | ΔH_fus (kJ/mol) | T_m (K) | bias | source |\n|---|---:|---:|---:|---|",
          *[f"| {s} | {fus[s].dh_fus / 1000:.2f} | {fus[s].t_fus:.1f} | {by_sol.loc[s, 'bias']:+.2f} | "
            f"{(fus[s].fusion_source or 'library value').replace('|', '/')[:110]} |" for s in by_sol.index],
          "",
          "![Parity](figures/parity.png)\n", "![Errors by type](figures/error_by_category.png)\n",
          "![Bias per solute](figures/bias_by_solute.png)\n",
          "## Error per pure solvent (solute offsets removed)\n",
          "Mean log₁₀ error per solvent after subtracting each solute's mean error, so that fusion-data "
          "errors cancel: a positive value means the model makes this solvent look like a better solvent "
          "than it is.\n",
          "| solvent | solutes | points | mean centred error | factor |\n|---|---:|---:|---:|---:|",
          *[f"| {r.Index} | {r.solutes} | {r.points} | {r.err:+.2f} | ×{10 ** r.err:.2f} |"
            for r in per_solvent.itertuples()],
          "",
          "## Largest systematic deviations (mean log₁₀ error per solute + solvent family)\n",
          "| solute | solvent | points | T (K) | bias | source |\n|---|---|---:|---|---:|---|"]
    md += [f"| {r.solute} | {r.family} | {r.points} | {r.T_min:.0f}–{r.T_max:.0f} | {r.bias:+.2f} | doi:{r.ref} |"
           for r in worst.itertuples()]
    md += ["", "## Best-predicted systems (≥ 3 points)\n",
           "| solute | solvent | points | T (K) | bias | source |\n|---|---|---:|---|---:|---|"]
    md += [f"| {r.solute} | {r.family} | {r.points} | {r.T_min:.0f}–{r.T_max:.0f} | {r.bias:+.2f} | doi:{r.ref} |"
           for r in best.itertuples()]
    if len(rk):
        md += ["", "## Solvent ranking (pure solvents)\n", "| solute | T (K) | solvents | Spearman ρ |\n|---|---:|---:|---:|"]
        md += [f"| {r.solute} | {r.T:.0f} | {r.solvents} | {r.rho:.2f} |" for r in rk.sort_values(["solute", "T"]).itertuples()]
    md += ["", "## Data and method\n",
           f"* Source: NIST/TRC ThermoML Archive, searched by InChIKey for {df.solute.nunique()} library solids "
           f"(T_m ≥ 285 K); {info['records']} records read. Points were kept when the measured property is the "
           "liquid-phase amount of a component coexisting with its own crystal, all other components are in "
           "the library, and the solvent composition is fully specified. Mass fractions and molalities were "
           "converted to mole fractions; mixed-solvent compositions on a solute-free basis.",
           "* Skipped data sets (counts of points/sets): " + ", ".join(f"{k}: {v}" for k, v in info["skipped"].items()) + ".",
           f"* {n_fail} points could not be predicted (solver did not converge) and are excluded from the statistics.",
           "* Choline chloride is a formula unit in the data and two COSMO species (choline⁺, Cl⁻) in the model; "
           "predictions are converted to formula units before comparison.",
           "* Fusion data are the library values the app uses (see the Compound library page); errors in ΔH_fus "
           "and T_m shift all solvents of a solute equally, which shows up as a per-solute bias.",
           "* Reproduce: `.venv/bin/python -m benchmarks.accuracy all` (downloads are cached in "
           "`benchmarks/cache/`).", "",
           "## References\n", "| DOI | year | points | solutes | title |\n|---|---|---:|---|---|"]
    md += [f"| [{r.doi}](https://doi.org/{r.doi}) | {r.year} | {r.points} | {r.solutes} | {r.title} |"
           for r in refs.sort_values("points", ascending=False).itertuples()]
    (HERE / "REPORT.md").write_text("\n".join(md) + "\n")
    by_cat.to_csv(DATA / "stats_by_category.csv")
    by_sol.to_csv(DATA / "stats_by_solute.csv")
    print(overall.T.to_string())
    print(by_cat[["points", "AAD", "bias", "within_x2", "within_x10", "AAD_ideal"]].to_string())
    print(by_sol[["points", "AAD", "bias", "within_x10", "AAD_ideal"]].to_string())
    print(f"figures: {figs}; report: {HERE / 'REPORT.md'}")
