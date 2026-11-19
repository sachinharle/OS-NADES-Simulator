"""Runtime configuration.

Values come from (lowest to highest priority): defaults, ``osnades.toml`` in the
data home, environment variables.
"""
from __future__ import annotations

import os
import platform
import shutil
import tomllib
from dataclasses import dataclass, field, asdict
from pathlib import Path

PACKAGE_DIR = Path(__file__).resolve().parent
SHIM_DIR = PACKAGE_DIR / "orca_shim"


def _default_home() -> Path:
    if env := os.environ.get("OSNADES_HOME"):
        return Path(env).expanduser()
    # Project checkout: keep data next to the package; otherwise ~/.osnades
    candidate = PACKAGE_DIR.parent / "data"
    if (PACKAGE_DIR.parent / "pyproject.toml").exists():
        return candidate
    return Path.home() / ".osnades"


def _is_linux_elf(path: Path) -> bool:
    try:
        with open(path, "rb") as f:
            return f.read(4) == b"\x7fELF"
    except OSError:
        return False


def _guess_orca_home() -> Path | None:
    for env in ("OSNADES_ORCA_HOME", "ORCA_HOME"):
        if v := os.environ.get(env):
            return Path(v).expanduser()
    if w := shutil.which("orca"):
        return Path(w).resolve().parent
    for p in sorted(Path.home().glob("orca*")) + sorted(Path("/opt").glob("orca*")):
        if (p / "orca").is_file():
            return p
    return None


@dataclass
class Settings:
    home: Path = field(default_factory=_default_home)
    orca_home: Path | None = field(default_factory=_guess_orca_home)
    # "auto" | "native" | "docker"
    orca_mode: str = "auto"
    docker_image: str = "osnades/orca-runtime:latest"
    ncores: int = 4
    # COSMO-RS reference solvent used for the !COSMORS(<solvent>) keyword.
    # ORCA always needs one; only the *.solute.orcacosmo file is used.
    cosmors_solvent: str = "water"
    # ORCA < 6.1.1 has no `--version`; OPI 2.x refuses it unless the check is disabled.
    version_check: bool = False

    @property
    def library_dir(self) -> Path:
        return self.home / "library"

    @property
    def jobs_dir(self) -> Path:
        return self.home / "jobs"

    @property
    def resolved_mode(self) -> str:
        if self.orca_mode != "auto":
            return self.orca_mode
        if self.orca_home and platform.system() != "Linux" and _is_linux_elf(self.orca_home / "orca"):
            return "docker"
        return "native"

    def to_dict(self) -> dict:
        d = asdict(self)
        return {k: (str(v) if isinstance(v, Path) else v) for k, v in d.items()}

    def save(self) -> Path:
        self.home.mkdir(parents=True, exist_ok=True)
        path = self.home / "osnades.toml"
        lines = []
        for k, v in self.to_dict().items():
            if k == "home" or v is None:
                continue
            lines.append(f"{k} = {v!r}" if not isinstance(v, str) else f'{k} = "{v}"')
        lines = [ln.replace("True", "true").replace("False", "false") for ln in lines]
        path.write_text("\n".join(lines) + "\n")
        return path


def load_settings() -> Settings:
    s = Settings()
    cfg = s.home / "osnades.toml"
    if cfg.exists():
        data = tomllib.loads(cfg.read_text())
        for k, v in data.items():
            if hasattr(s, k):
                setattr(s, k, Path(v).expanduser() if k == "orca_home" else v)
    env_map = {
        "OSNADES_ORCA_MODE": "orca_mode",
        "OSNADES_DOCKER_IMAGE": "docker_image",
        "OSNADES_NCORES": "ncores",
    }
    for env, attr in env_map.items():
        if v := os.environ.get(env):
            setattr(s, attr, int(v) if attr == "ncores" else v)
    if v := os.environ.get("OSNADES_ORCA_HOME"):
        s.orca_home = Path(v).expanduser()
    return s
