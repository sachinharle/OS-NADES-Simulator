# EutecticSolvent · OS-NADES Simulator

Computing the solubility and mole fractions of components in complex (deep eutectic) solvents.

Predicts solid solubilities in conventional solvents and **(natural) deep eutectic solvents**
from first principles. It turns the FACCTs
[OPI openCOSMO-RS tutorial](https://www.faccts.de/docs/opi/nightly/docs/contents/notebooks/opencosmors.html)
into an application with a web UI, a CLI, a compound library and support for multi-component solvents.

```
SMILES / XYZ ──RDKit 3D──► GFN2-xTB opt (ORCA) ──► ORCA !COSMORS (BP86/def2-TZVPD gas + CPCM)
      └──► .orcacosmo σ-surface ──► openCOSMO-RS 24a ──► ln γ ──► SLE solubility, screening, Gᴱ
```

## Features

* **NADES workflow** (second in the sidebar, after the Overview; the pages are ordered as you would use them) designs a multicomponent solvent step by step:
  1. **Solute**: pick it from the library, a preset, or enter a SMILES.
  2. **Components**: HBAs (choline chloride, betaine…), HBDs (polyols, amides, acids, sugars,
     terpenes) and optionally water, with min/max mole-fraction bounds. Templates are provided for
     ChCl + polyol, ChCl + acid, betaine and hydrophobic systems.
  3. **Properties**: missing σ-surfaces are queued for ORCA with one click (live progress). Fusion
     data are looked up (Cheméo + NIST) or entered.
  4. **Optimise**: the component mole fractions are searched on a simplex grid, then refined around
     the best points. Each composition is screened with batched openCOSMO-RS calls: the SLE solubility
     of the solute, and whether any component would crystallise (ln a_i vs. its ideal-solid limit), at
     every temperature of the range. Goals: highest solubility over the range, at the lowest
     temperature, in g/kg, or the largest gain for cooling crystallisation.
  5. **Results**: the best composition (mole and mass fractions, molar ratio), exact solubilities
     against the ideal solution, the solubility-vs-T curves of the top candidates, a composition map
     (ternary for 3 components), CSV export, and "Save as solvent". A **Phase diagram** tab computes
     the solid–liquid diagram of the solute-free solvent (binary for 2 components, ternary liquidus
     map for 3) with the chosen composition marked, its predicted melting point, the eutectic, and
     whether it stays liquid down to the coldest temperature of the range. Export the diagram as PNG
     (chosen solvent marked), the liquidus data as CSV, and a one-file summary CSV.

     If no composition on the grid stays liquid over the whole range, **Find a fully liquid
     composition** searches more finely: it looks for the lowest-melting composition within your bounds
     and re-runs the optimisation with every composition liquid at the coldest temperature (the grid
     often misses the narrow region near a ternary eutectic). If none can be liquid, it reports the
     lowest achievable melting point and the bounds it sits on, and offers to start the range just
     above that temperature or to change the components/bounds. `osnades design` does the same
     automatically (`--no-liquid-search` to skip).

  The same is available as `osnades design` (`--phase-diagram` adds the phase diagram).
* **Compound library**: add molecules or ions from presets, SMILES, XYZ, or existing `.orcacosmo`
  files. The library supports several conformers per species and shows a 3D viewer, σ-profiles
  and σ-moment descriptors.
* **Fusion data lookup**: when you add a neutral compound, its measured ΔH_fus and T_m are fetched
  automatically from [Cheméo](https://www.chemeo.com) (NIST and literature DSC values) and the
  [NIST WebBook](https://webbook.nist.gov) (including melting of individual crystal forms). The compound is matched by InChIKey, so only the exact molecule is accepted. The source
  reference (DOI) is stored with the species. If only a melting point is published, ΔH_fus is estimated
  with Walden's rule (ΔS_fus ≈ 56.5 J/mol/K) and flagged as an estimate. The library's **Fusion data**
  tab shows the status of every species (measured / entered by you / estimate / missing), searches
  all missing (or estimated, or all) species at once, and lets you review the results before saving.
  "Look up published values" on a single species lists every value found so you can pick another
  reference. Values you enter yourself are never overwritten unless you tick them.
* **Your own solvents**: switch on "Use as a solvent" for any neutral library species, tick
  "Also use it as a solvent" when adding a new molecule, or press "Save as solvent" on a custom
  mixture. These solvents then show up on the Solubility, Temperature curves and NADES pages (group
  "My solvents / custom") and in the CLI. They are stored in `data/library/solvents.json`.
* **Solubility**: ideal-solid SLE, `ln(xγ) = -ΔH_fus/R(1/T-1/T_m) + ΔCp/R(T_m/T-1-ln(T_m/T))`.
  * The *iterative* mode solves the equation self-consistently with a bracketed root search.
  * The *infinite-dilution* mode is the tutorial's non-iterative variant.
  * You can screen many solvents, scan temperature, add water (wt%) to any solvent, and get
    results as mole fraction, mass fraction and g/kg.
* **Temperature curves**: solubility vs T for many solvents at once, grouped by solvent class
  (ChCl NADES with polyol/amide or acid HBDs, betaine/hydrophobic DES, conventional solvents).
  Each panel picks a linear or log scale on its own. The page includes a temperature-gain table
  for picking cooling-crystallisation solvents, plus CSV and PNG export.
* **NADES**: 26 preset systems.
  * ChCl-based (Reline, Ethaline, Glyceline, organic acids, glucose), betaine-based, and
    hydrophobic (menthol/thymol/decanoic acid) systems.
  * A custom HBA:HBD:… builder with an electroneutrality check.
  * Pick a solute to compare its ideal mole fraction (γ = 1, from ΔH_fus/T_m/ΔCp only) with the
    openCOSMO-RS solubility in the designed mixture (x/x_ideal, ln γ, g/kg, and a T scan).
  * **Eutectic phase diagram** of an HBA + HBD pair: the liquidus of each component, the eutectic
    point, how far it lies below the ideal eutectic (how "deep" it is), and the composition window
    that is liquid at your operating temperature. Salts such as choline chloride are handled as ions,
    referenced to their pure molten salt. Approximate literature fusion data ship for water, ChCl, urea,
    glycerol, ethylene glycol, menthol, thymol, decanoic acid, glucose and DMSO.
  * **Ternary phase diagram** (switch *Components* to 3): the liquidus surface of three components
    (e.g. HBA + two HBDs) drawn as isotherms on a composition triangle, the ternary eutectic (refined
    locally to ±0.006 mole fraction) next to the ideal one, and the region that is liquid at your
    operating temperature (dotted line) with a table of those compositions. Download CSV/PNG or save
    the eutectic composition as a solvent. ln γ is computed on a temperature grid uniform in 1/T for
    all compositions in one batch and interpolated; each edge reproduces the binary diagram within
    ~0.5 K. Takes ~20 s for neutral components and ~15–60 s with a salt at the default 0.1 grid step.
  * Per-species activity coefficients (residual and combinatorial), the mixture σ-profile,
    and binary Gᴱ/RT curves.
* **ORCA backends**:
  * `native`: ORCA runs directly on the host.
  * `docker`: a Linux x86-64 ORCA build runs on macOS/arm64 through an amd64 container.
    ORCA is bind-mounted into the container, never copied into the image.

## Install

```bash
uv venv -p 3.12 .venv
uv pip install -p .venv/bin/python -e ".[dev]"
```

ORCA ≥ 6.0 must be installed separately (free for academic use from the ORCA forum).

```bash
.venv/bin/osnades setup --orca-home ~/orca_6_1_0 --ncores 4   # builds the Docker runtime if needed, runs a test job
```

On this machine, `~/orca_6_1_0` holds the **Linux x86-64** build. The `auto` mode detects
this and uses Docker; Docker Desktop must be running. Another option is to install the native
macOS-arm64 ORCA build and set `--mode native`, which is several times faster.

> ORCA 6.1.0 has no `--version` flag, and OPI 2.0 asks for ORCA ≥ 6.1.1. The app therefore turns
> off OPI's version check (`version_check = false` in `data/osnades.toml`). Upgrading to 6.1.1
> is recommended.

## Use

### Web app

```bash
.venv/bin/osnades app          # http://localhost:8501
```

The app has these pages: Overview (one-click "compute all preset species"), Compound library,
Solubility, Temperature curves, NADES designer, Jobs (live ORCA queue), and Settings.

### CLI

```bash
osnades add paracetamol                                   # preset SMILES + fusion data
osnades add vanillin --smiles "COc1cc(C=O)ccc1O" --dh-fus 22400 --t-fus 355
osnades bootstrap "Reline (ChCl:Urea 1:2)" Water Ethanol  # compute species for these systems
osnades import myspecies conf0.orcacosmo conf1.orcacosmo  # reuse σ-surfaces
osnades list
osnades fusion paraCoumaric_acid                          # list published ΔH_fus / T_m; --apply stores the suggestion
osnades fusion --all --apply                              # fill in every neutral species that has no data
osnades fusion --all --estimates                          # also re-check estimated values
osnades solubility paracetamol -s Water -s Ethanol -s "Reline (ChCl:Urea 1:2)" -T 298.15 318.15
osnades solubility paracetamol -s "choline:chloride:glycerol=1:1:2" --water 10 --csv out.csv
osnades curves paracetamol --tmin 0 --tmax 80 -o results/paracetamol_solubility_vs_T   # CSV + grouped PNG
osnades binary water ethanol -T 298.15
osnades solvent add vanillin                              # use a library species as a solvent
osnades solvent add "Thymol:Decanoic acid 1:2" "thymol:decanoic acid=1:2"
osnades solvent list
osnades design paracetamol -c "choline chloride:0.1:0.5" -c glycerol -c urea -c "water:0:0.3" \
        --tmin 15 --tmax 45 --objective mean_x -o paracetamol_nades   # best multicomponent NADES
osnades design paracetamol -c thymol -c menthol -c "decanoic acid" --phase-diagram -o hyd
        # + binary/ternary phase diagram of the best solvent: its melting point, eutectic, liquid share, CSV + PNG
osnades eutectic thymol menthol -T 298.15 -o thymol_menthol   # eutectic point, liquid window, CSV + PNG
osnades eutectic "choline chloride" urea --dh-b 14600 --tm-b 406
osnades eutectic thymol menthol "decanoic acid" -T 298.15 --step 0.1 -o tmd   # ternary: eutectic, liquid share, CSV + PNG
```

![Paracetamol solubility vs temperature](results/paracetamol_solubility_vs_T.png)

### Reproduce the tutorial

```bash
OPI_ORCA=$PWD/osnades/orca_shim OSNADES_ORCA_HOME=~/orca_6_1_0 .venv/bin/python tutorial/run_tutorial.py
.venv/bin/python -m pytest -q
```

| Paracetamol, 298.15 K | tutorial non-iter. | this app non-iter. | tutorial iter. | this app iter. | experiment |
|---|---|---|---|---|---|
| Water   | 7.424e-4 | 7.424e-4 | 1.949e-4 | **7.670e-4** | 1.773e-3 |
| Ethanol | 1.712e-1 | 1.712e-1 | 9.117e-2 | **9.983e-2** | 7.586e-2 |

The σ-descriptors and non-iterative results match the tutorial exactly. The tutorial's iterative
values come from `fsolve` applied to `|f(x)|`, so they are not converged roots.
`tests/test_thermo.py` checks that this app's iterative answer satisfies the SLE equation to 1e-6.

## Accuracy

[benchmarks/REPORT.md](benchmarks/REPORT.md) compares predictions with every matching experimental
solubility data set in the NIST/TRC ThermoML Archive (1065 points, 24 papers, 8 solutes, water,
C1–C4 alcohols, acetone, acetonitrile, ethyl acetate, DMF, dioxane and their mixtures):

* typical error a factor of 2.3 (mean |log₁₀ x_pred/x_exp| = 0.36); 94% of points within a factor of 10;
  the ideal solution with the same fusion data is off by a factor of 4.9 (AAD 0.69);
* after removing each solute's constant offset, which comes from its fusion data, the error is a factor of 1.6;
* ranking of pure solvents for one solute is only moderate (median Spearman ρ 0.36): solvents within
  a factor of ~1.5 of each other (the alcohols) are often mis-ordered, and acetone, ethyl acetate and
  acetonitrile come out ×1.8–2.6 too good relative to the alcohols and water;
* the largest errors come from estimated fusion data (oxalic acid, betaine), so measured ΔH_fus and T_m matter most;
* the library already covers the solvents of almost all ThermoML data for these solutes (38 data sets
  skipped for a missing component);
* NADES are **not** covered: ThermoML has no DES solubility data for the library compounds.

Reproduce with `.venv/bin/python -m benchmarks.accuracy all`.

## Model notes and limits

* Salts are split into their ions, as is usual in COSMO-RS (choline⁺ + Cl⁻). The mole fraction `x`
  counts every species. Internally, ln γ is computed against the COSMO reference and converted to
  the pure-liquid reference of the solute. This way no unphysical "pure ion" liquid is ever needed.
* The model assumes an ideal solid phase. It ignores polymorphs, solvates, solute
  ionisation/pH, and LLE (a warning is shown if several roots exist).
* ΔH_fus and T_m have a big effect on the result. Paracetamol (NIST, as in the tutorial) and
  p-coumaric acid (34.3 kJ/mol, 494.3 K, doi:10.1021/acs.jced.6b00361) ship with values. Other
  solutes get published values from the automatic lookup (needs internet), or you can enter your own.
  Check the source shown next to the values; published DSC data for the same compound can differ by
  several kJ/mol (thymol: 17.4–22.1 kJ/mol). The shipped library has measured values for every
  neutral species except two estimates: betaine (17.98 kJ/mol, 566.2 K, group contribution, Wang et
  al. 2014; it decomposes on melting) and oxalic acid (Walden's rule; it sublimes and decomposes near
  its melting point, so no DSC fusion enthalpy exists).
* "Citric acid" and "Glucose" as pure solvents are hypothetical sub-cooled liquids (both melt above
  420 K). They show how each HBD interacts with a solute on its own; real extraction media use them
  inside NADES (ChCl:Citric acid, ChCl:Glucose) or with water.
* Eutectic diagrams assume pure, immiscible solids (no co-crystals or solid solutions). For neutral
  pairs (e.g. thymol:menthol, predicted eutectic at x_thymol = 0.39, −1 °C with the library's
  fusion data) results are reasonable. Adding decanoic acid gives a predicted ternary eutectic of
  thymol:menthol:decanoic acid ≈ 1.6:2.1:1 at −13 °C, 12 K below the binary one.
  For choline chloride, COSMO-RS over-stabilises the split ions, so the eutectic comes out far too
  deep (often below the temperatures where the model converges, which is reported as a bound).
  Use these results to compare candidates, not as absolute melting points.
* The NADES optimiser inherits these limits. COSMO-RS strongly favours choline chloride-rich
  mixtures, which in reality are solid, so salts are bounded to 0.1–0.5 by default and the app warns
  when the optimum sits on a bound. Treat the result as a ranking of candidate compositions to test,
  and check the melting point of the chosen mixture.
* Running under Docker emulation, one species takes ~1–3 min (more for big molecules or several
  conformers).

## Precomputed data

`data/library/` holds the σ-surfaces for all 26 preset solvent systems plus paracetamol (BP86/def2-TZVPD,
GFN2-xTB geometries). With them, solubility screening works out of the box, even without ORCA. You
only need ORCA to add new compounds.

## Stand-alone ORCA scripts

The original batch scripts sit in the repository root and work independently of the `osnades` package:

```bash
python runORCA.py input.csv   # each "name","SMILES" row: smiles2xyz.py → createORCAinp.py → orca
```

`smiles2xyz.py` builds a UFF-optimised RDKit geometry. `createORCAinp.py` writes a B3LYP/def2-TZVP
single-point input. `input.csv` lists quercetin, glycerol and choline chloride.

## Layout

```
osnades/        config, orca (OPI + Docker shim), structures (RDKit), library, sigma, thermo,
                plots, jobs, cli, app (Streamlit), data/presets.json, orca_shim/
docker/         amd64 runtime image (Debian + OpenMPI) for Linux ORCA builds
tutorial/       script version of the FACCTs notebook
tests/          regression tests (fixtures = tutorial σ-surfaces)
data/           library (index.json + .orcacosmo); local settings and ORCA job folders are git-ignored
results/        paracetamol screen and temperature curves
```
