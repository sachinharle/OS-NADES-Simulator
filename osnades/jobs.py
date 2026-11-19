"""Background queue for long ORCA calculations (used by the web UI)."""
from __future__ import annotations

import itertools
import threading
import time
import traceback
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from typing import Callable


@dataclass
class Job:
    id: int
    label: str
    status: str = "queued"  # queued | running | done | failed
    log: list[str] = field(default_factory=list)
    error: str = ""
    created: float = field(default_factory=time.time)
    started: float | None = None
    finished: float | None = None

    def write(self, msg: str) -> None:
        self.log.append(f"[{time.strftime('%H:%M:%S')}] {msg}")

    @property
    def elapsed(self) -> float:
        if not self.started:
            return 0.0
        return (self.finished or time.time()) - self.started


class JobManager:
    """Runs jobs one at a time (each ORCA run already uses several cores)."""

    def __init__(self, workers: int = 1):
        self._pool = ThreadPoolExecutor(max_workers=workers, thread_name_prefix="osnades-job")
        self._jobs: dict[int, Job] = {}
        self._ids = itertools.count(1)
        self._lock = threading.Lock()

    def submit(self, label: str, fn: Callable[..., object], *args, **kwargs) -> Job:
        job = Job(next(self._ids), label)
        with self._lock:
            self._jobs[job.id] = job

        def runner():
            job.status, job.started = "running", time.time()
            try:
                fn(*args, log=job.write, **kwargs)
                job.status = "done"
            except Exception as exc:  # noqa: BLE001
                job.status = "failed"
                job.error = f"{exc}"
                job.write(traceback.format_exc(limit=3))
            finally:
                job.finished = time.time()

        self._pool.submit(runner)
        return job

    def jobs(self) -> list[Job]:
        with self._lock:
            return sorted(self._jobs.values(), key=lambda j: -j.id)

    def active(self) -> list[Job]:
        return [j for j in self.jobs() if j.status in ("queued", "running")]

    def clear_finished(self) -> None:
        with self._lock:
            self._jobs = {k: j for k, j in self._jobs.items() if j.status in ("queued", "running")}
