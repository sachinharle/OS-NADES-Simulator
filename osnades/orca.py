"""ORCA execution through OPI.

Two backends:

* ``native`` – ORCA binary runs directly on the host.
* ``docker`` – a Linux x86-64 ORCA build is run inside a small amd64 container
  (``docker/Dockerfile``). OPI is pointed at ``orca_shim/`` whose ``orca`` and
  ``orca_2json`` entries forward to ``docker run`` with the ORCA folder and the job
  directory bind-mounted at identical paths.
"""
from __future__ import annotations

import os
import shutil
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

from .config import SHIM_DIR, Settings
from .structures import Geometry, parse_xyz

Log = Callable[[str], None]


class OrcaError(RuntimeError):
    pass


def configure_environment(settings: Settings) -> None:
    """Export the environment variables OPI (and the Docker shim) read."""
    if settings.orca_home is None:
        raise OrcaError("ORCA installation not found. Set orca_home in Settings or OSNADES_ORCA_HOME.")
    home = Path(settings.orca_home).expanduser()
    if settings.resolved_mode == "docker":
        os.environ["OPI_ORCA"] = str(SHIM_DIR)
        os.environ["OSNADES_ORCA_HOME"] = str(home)
        os.environ["OSNADES_DOCKER_IMAGE"] = settings.docker_image
        os.environ.pop("OPI_MPI", None)
    else:
        os.environ["OPI_ORCA"] = str(home)


def build_docker_image(settings: Settings, log: Log = print) -> None:
    dockerfile_dir = Path(__file__).resolve().parent.parent / "docker"
    if not (dockerfile_dir / "Dockerfile").exists():
        raise OrcaError(f"Dockerfile not found in {dockerfile_dir}")
    log(f"Building {settings.docker_image} ...")
    proc = subprocess.run(
        ["docker", "build", "--platform", "linux/amd64", "-t", settings.docker_image, str(dockerfile_dir)],
        capture_output=True, text=True,
    )
    if proc.returncode != 0:
        raise OrcaError(proc.stderr[-2000:])
    log("Docker image ready.")


@dataclass
class BackendStatus:
    ok: bool
    mode: str
    orca_home: str | None
    messages: list[str]


def check_backend(settings: Settings, quick_run: bool = False) -> BackendStatus:
    msgs: list[str] = []
    mode = settings.resolved_mode
    home = settings.orca_home
    ok = True
    if home is None or not (Path(home) / "orca").is_file():
        return BackendStatus(False, mode, str(home) if home else None, ["ORCA binary not found."])
    msgs.append(f"ORCA folder: {home}")
    if mode == "docker":
        if shutil.which("docker") is None:
            return BackendStatus(False, mode, str(home), msgs + ["Docker CLI not found."])
        info = subprocess.run(["docker", "info", "--format", "{{.ServerVersion}}"], capture_output=True, text=True)
        if info.returncode != 0:
            return BackendStatus(False, mode, str(home), msgs + ["Docker daemon is not running."])
        msgs.append(f"Docker {info.stdout.strip()}")
        img = subprocess.run(["docker", "image", "inspect", settings.docker_image], capture_output=True, text=True)
        if img.returncode != 0:
            ok = False
            msgs.append(f"Runtime image {settings.docker_image} missing - build it from Settings or `osnades setup`.")
        else:
            msgs.append(f"Runtime image {settings.docker_image} present")
    if ok and quick_run:
        try:
            t = time.time()
            e = single_point_energy(settings, parse_xyz("3\n\nO 0 0 0.1173\nH 0 0.7572 -0.4692\nH 0 -0.7572 -0.4692\n"))
            msgs.append(f"Test HF/def2-SVP water: E = {e:.6f} Eh ({time.time() - t:.1f} s)")
        except Exception as exc:  # noqa: BLE001
            ok = False
            msgs.append(f"Test run failed: {exc}")
    return BackendStatus(ok, mode, str(home), msgs)


def _run(
    settings: Settings, workdir: Path, basename: str, geom: Geometry, charge: int, mult: int,
    keywords: str, ncores: int | None = None, log: Log = print,
):
    from opi.core import Calculator
    from opi.input.structures.structure import Structure

    configure_environment(settings)
    workdir.mkdir(parents=True, exist_ok=True)
    xyz_path = workdir / f"{basename}.xyz"
    xyz_path.write_text(geom.to_xyz())
    calc = Calculator(basename=basename, working_dir=workdir, version_check=settings.version_check)
    calc.structure = Structure.from_xyz(xyz_path, charge=charge, multiplicity=mult)
    calc.input.add_arbitrary_string(keywords)
    n = ncores or settings.ncores
    if n and n > 1 and geom.symbols.__len__() > 1:
        calc.input.ncores = n
    calc.write_input()
    log(f"ORCA [{basename}] {keywords} (charge {charge}, mult {mult}, {n} cores)")
    t = time.time()
    calc.run()
    out = calc.get_output()
    outfile = workdir / f"{basename}.out"
    if not outfile.exists() or "ORCA TERMINATED NORMALLY" not in outfile.read_text(errors="ignore"):
        tail = outfile.read_text(errors="ignore")[-1500:] if outfile.exists() else "(no output file)"
        raise OrcaError(f"ORCA did not terminate normally ({outfile}):\n{tail}")
    log(f"ORCA [{basename}] finished in {time.time() - t:.0f} s")
    return out


def single_point_energy(settings: Settings, geom: Geometry, workdir: Path | None = None) -> float:
    workdir = workdir or (settings.jobs_dir / f"check_{int(time.time())}")
    _run(settings, workdir, "check", geom, 0, 1, "! HF def2-SVP", ncores=1, log=lambda s: None)
    for line in (workdir / "check.out").read_text().splitlines():
        if line.startswith("FINAL SINGLE POINT ENERGY"):
            return float(line.split()[-1])
    raise OrcaError("Energy not found")


def xtb_optimize(
    settings: Settings, geom: Geometry, workdir: Path, basename: str, charge: int = 0, mult: int = 1,
    solvation: str | None = "water", log: Log = print,
) -> Geometry:
    """GFN2-xTB geometry optimisation (ORCA + bundled xtb), optional ALPB solvation."""
    if len(geom.symbols) == 1:
        return geom
    kw = "! XTB2 Opt" + (f" ALPB({solvation})" if solvation else "")
    _run(settings, workdir, basename, geom, charge, mult, kw, log=log)
    xyz = workdir / f"{basename}.xyz"
    return parse_xyz(xyz.read_text())


def cosmors_orcacosmo(
    settings: Settings, geom: Geometry, workdir: Path, basename: str, charge: int = 0, mult: int = 1,
    log: Log = print,
) -> Path:
    """Run ORCA's COSMO-RS protocol (BP86/def2-TZVPD gas + CPCM) and return the solute .orcacosmo file."""
    _run(settings, workdir, basename, geom, charge, mult, f"!COSMORS({settings.cosmors_solvent})", log=log)
    f = workdir / f"{basename}.solute.orcacosmo"
    if not f.exists():
        raise OrcaError(f"Expected {f} was not produced")
    return f
