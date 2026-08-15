"""Progress reporting for external energy calculations."""

from __future__ import annotations

import logging
import os
import shutil
import sys
from pathlib import Path
from time import monotonic, sleep

from tqdm import tqdm


TQDM_BAR_FORMAT = (
    '            {l_bar}{bar:100}| {n_fmt}/{total_fmt} '
    '[elapsed: {elapsed} remaining: {remaining}]'
)
PROGRESS_STYLES = ('auto', 'rich', 'classic', 'plain', 'none')
MAX_RICH_WIDTH = 120

try:
    from rich.console import Console
    from rich.progress import (
        BarColumn,
        MofNCompleteColumn,
        Progress,
        ProgressColumn,
        SpinnerColumn,
        TaskProgressColumn,
        TextColumn,
        TimeElapsedColumn,
        TimeRemainingColumn,
    )
    from rich.table import Column
    from rich.text import Text
except ImportError:  # pragma: no cover - exercised only in minimal installations
    Console = Progress = None
    ProgressColumn = object


def resolve_progress_style(style='auto', stream=None):
    """Resolve ``auto`` without assuming that stderr is an interactive terminal."""
    if style not in PROGRESS_STYLES:
        raise ValueError(f'Unknown progress style: {style}')
    if style != 'auto':
        if style == 'rich' and Progress is None:
            logging.warning('Rich progress requested but Rich is unavailable; using classic progress.')
            return 'classic'
        return style

    stream = stream or sys.stderr
    is_terminal = bool(getattr(stream, 'isatty', lambda: False)())
    if is_terminal and os.getenv('TERM', '').lower() != 'dumb' and Progress is not None:
        return 'rich'
    # MPI launchers commonly replace stderr with a forwarding pipe even when the
    # user is watching an interactive terminal. Preserve the live classic bar in
    # that ambiguous case; plain milestone logs remain available explicitly.
    return 'classic'


class FrameCounter:
    """Incrementally count completed frames in rank output files."""

    def __init__(self, output_basename, mpi_size=1, nmode=False):
        self.output_basename = output_basename
        self.mpi_size = mpi_size
        self.marker = 'Total:' if nmode else '                    FINAL RESULTS'
        self._positions = {}
        self._counts = {}

    def count(self):
        if 'gbnsr6' in self.output_basename:
            output_path = Path(self.output_basename)
            output_folder = str(output_path.parent)
            stem = output_path.stem
            return sum(
                len(list(Path(folder).glob(f'{stem}*')))
                for folder in _rank_names(output_folder, self.mpi_size)
            )

        for filename in _rank_names(self.output_basename, self.mpi_size):
            path = Path(filename)
            if not path.exists():
                continue
            size = path.stat().st_size
            position = self._positions.get(path, 0)
            if size < position:
                position = 0
                self._counts[path] = 0
            with path.open(errors='replace') as output:
                output.seek(position)
                new_frames = 0
                while line := output.readline():
                    new_frames += line.startswith(self.marker)
                self._counts[path] = self._counts.get(path, 0) + new_frames
                self._positions[path] = output.tell()
        return sum(self._counts.values())


def _rank_names(template, mpi_size):
    """Expand rank placeholders without returning duplicate paths."""
    return tuple(dict.fromkeys(
        template % rank if '%d' in template else template
        for rank in range(mpi_size)
    ))


class _ClassicReporter:
    def __init__(self, total, stream=None, **_):
        self.progress = tqdm(
            total=total, ascii=True, bar_format=TQDM_BAR_FORMAT,
            file=stream or sys.stderr,
        )

    def update(self, completed):
        self.progress.update(completed - self.progress.n)

    def close(self, completed):
        self.update(completed)
        self.progress.clear()
        self.progress.close()


class _FrameRateColumn(ProgressColumn):
    def render(self, task):
        speed = task.finished_speed or task.speed
        value = '— frame/s' if speed is None else f'{speed:.2f} frame/s'
        return Text(value, style='progress.data.speed', no_wrap=True)


class _RichReporter:
    def __init__(self, total, label, mpi_size, stream=None):
        stream = stream or sys.stderr
        is_terminal = bool(getattr(stream, 'isatty', lambda: False)())
        width = min(
            shutil.get_terminal_size(fallback=(MAX_RICH_WIDTH, 24)).columns,
            MAX_RICH_WIDTH,
        )
        console = Console(
            file=stream,
            force_terminal=None if is_terminal else True,
            width=width,
        )
        self.progress = Progress(
            SpinnerColumn(),
            TextColumn(
                '[bold cyan]{task.description}',
                table_column=Column(no_wrap=True),
            ),
            BarColumn(table_column=Column(ratio=2, min_width=30)),
            TaskProgressColumn(),
            MofNCompleteColumn(),
            _FrameRateColumn(),
            TextColumn('•'),
            TimeElapsedColumn(),
            TextColumn('• ETA'),
            TimeRemainingColumn(),
            console=console,
            expand=True,
        )
        ranks = f' · {mpi_size} rank' + ('s' if mpi_size != 1 else '')
        self.task = self.progress.add_task(f'{label}{ranks}', total=total)
        self.progress.start()

    def update(self, completed):
        self.progress.update(self.task, completed=completed, refresh=True)

    def close(self, completed):
        self.update(completed)
        self.progress.stop()


class _PlainReporter:
    def __init__(self, total, label, mpi_size, log_level=logging.INFO, **_):
        self.total = total
        self.label = label
        self.mpi_size = mpi_size
        self.log_level = log_level
        self.started = monotonic()
        self.last_completed = -1
        self.last_reported_percent = 0

    def update(self, completed):
        percent = int(completed * 100 / self.total) if self.total else 100
        if completed >= self.total:
            self.last_completed = completed
            return
        milestone = percent >= self.last_reported_percent + 10
        if completed != self.last_completed and milestone:
            elapsed = monotonic() - self.started
            rate = completed / elapsed if elapsed else 0
            remaining = _format_duration((self.total - completed) / rate) if rate else '--:--'
            logging.log(
                self.log_level,
                '  %s progress: %d/%d frames (%d%%), %.2f frame/s, elapsed %s, ETA %s [%d MPI rank%s]',
                self.label, completed, self.total, percent, rate,
                _format_duration(elapsed), remaining, self.mpi_size,
                '' if self.mpi_size == 1 else 's',
            )
            self.last_reported_percent = (percent // 10) * 10
        self.last_completed = completed

    def close(self, completed):
        self.update(completed)


def _format_duration(seconds):
    seconds = max(0, int(seconds))
    hours, remainder = divmod(seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    return f'{hours:d}:{minutes:02d}:{seconds:02d}' if hours else f'{minutes:02d}:{seconds:02d}'


def _reporter(style, total, label, mpi_size, stream=None):
    if style == 'rich':
        return _RichReporter(total, label, mpi_size, stream)
    if style == 'classic':
        return _ClassicReporter(total=total, stream=stream)
    if style == 'plain':
        return _PlainReporter(total=total, label=label, mpi_size=mpi_size)
    return None


def monitor_progress(output_basename, nframes=1, mpi_size=1, nmode=False,
                     style='auto', label='Frames', stream=None, poll_interval=1.0):
    """Monitor calculation outputs and render progress until all frames finish."""
    style = resolve_progress_style(style, stream)
    if style == 'none':
        return

    reporter = _reporter(style, nframes, label, mpi_size, stream)
    # DEBUG records are captured by gmx_MMPBSA.log but filtered from the normal
    # INFO-level terminal handler. Plain mode already emits its milestones at INFO.
    log_reporter = None if style == 'plain' else _PlainReporter(
        total=nframes, label=label, mpi_size=mpi_size, log_level=logging.DEBUG
    )
    counter = FrameCounter(output_basename, mpi_size, nmode)
    completed = 0
    started = monotonic()
    try:
        while completed < nframes:
            completed = min(counter.count(), nframes)
            reporter.update(completed)
            if log_reporter:
                log_reporter.update(completed)
            if completed < nframes:
                sleep(poll_interval)
    finally:
        reporter.close(completed)
        if log_reporter:
            log_reporter.close(completed)
        if completed >= nframes:
            elapsed = monotonic() - started
            rate = completed / elapsed if elapsed else 0
            logging.info(
                '  %s completed: %d frames in %s (%.2f frame/s)',
                label, completed, _format_duration(elapsed), rate,
            )
