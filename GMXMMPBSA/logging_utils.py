"""Logging setup shared by the command-line entry points."""

import logging
import shlex
import sys
import threading
import time
from pathlib import Path


class WarningSpacingFormatter(logging.Formatter):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._previous_was_warning = False

    def format(self, record):
        message = super().format(record)
        if record.levelno == logging.WARNING:
            prefix = '' if self._previous_was_warning else '\n'
            self._previous_was_warning = True
            return f'{prefix}{message}\n'
        self._previous_was_warning = False
        return message


class ElapsedProgress:
    """Show a small elapsed-time spinner on an interactive terminal."""

    _frames = ('|', '/', '-', '\\')

    def __init__(self, label='Operation', message='This can take a while; please wait...', stream=None, interval=1):
        self.label = label
        self.message = message
        self.stream = stream or sys.stderr
        self.interval = interval
        self._stop = threading.Event()
        self._thread = None
        self._started_at = None

    def __enter__(self):
        self._started_at = time.monotonic()
        if not getattr(self.stream, 'isatty', lambda: False)():
            return self
        self._thread = threading.Thread(target=self._run, daemon=True)
        self._thread.start()
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        elapsed = time.monotonic() - self._started_at
        if self._thread is not None:
            self._stop.set()
            self._thread.join()
            self.stream.write('\r\033[K')
            self.stream.flush()
        if exc_type is None:
            logging.info('%s completed in %.1f seconds.', self.label, elapsed)
        else:
            logging.info('%s stopped after %.1f seconds.', self.label, elapsed)

    def _run(self):
        tick = 0
        while not self._stop.wait(self.interval):
            elapsed = int(time.monotonic() - self._started_at)
            frame = self._frames[tick % len(self._frames)]
            prefix = f'{self.message} ' if self.message else ''
            self.stream.write(f'\r{prefix}[{frame}] {elapsed}s')
            self.stream.flush()
            tick += 1


def elapsed_progress(label='Operation', message='This can take a while; please wait...', stream=None, interval=1):
    """Return a context manager that reports elapsed time and spins on interactive terminals."""
    return ElapsedProgress(label=label, message=message, stream=stream, interval=interval)


class _MasterConsoleFilter(logging.Filter):
    """Keep informational and warning output on the master MPI rank."""

    def __init__(self, master):
        super().__init__()
        self.master = master

    def filter(self, record):
        return self.master or record.levelno >= logging.ERROR


class RecordCountingHandler(logging.Handler):
    """Count warning and error records without depending on formatted output."""

    def __init__(self):
        super().__init__()
        self.warning_count = 0
        self.error_count = 0

    def emit(self, record):
        if record.levelno == logging.WARNING:
            self.warning_count += 1
        elif record.levelno >= logging.ERROR:
            self.error_count += 1


def get_record_counts():
    """Return warning and error totals from the active CLI logging counter."""
    info = {'warning': 0, 'error': 0}
    for handler in logging.getLogger().handlers:
        if isinstance(handler, RecordCountingHandler):
            info['warning'] += handler.warning_count
            info['error'] += handler.error_count
    return info


def _rank_failure_log_path(log_file, rank):
    """Return the diagnostic log path used when a non-master rank fails."""
    path = Path(log_file)
    suffix = path.suffix
    return path.with_name(f'{path.stem}.rank-{rank}{suffix}')


class _RankFailureFileHandler(logging.Handler):
    """Create a rank-specific log lazily, and only after an error is logged."""

    def __init__(self, log_file, rank):
        super().__init__(level=logging.ERROR)
        self.log_file = _rank_failure_log_path(log_file, rank)
        self._file_handler = None

    def emit(self, record):
        if self._file_handler is None:
            self._file_handler = logging.FileHandler(self.log_file, mode='a', encoding='utf-8')
            self._file_handler.setFormatter(WarningSpacingFormatter("[%(levelname)-7s] %(message)s"))
        self._file_handler.emit(record)

    def close(self):
        if self._file_handler is not None:
            self._file_handler.close()
            self._file_handler = None
        super().close()


def _new_file_handler(log_file, rank, mode='w'):
    file_handler = logging.FileHandler(log_file, mode=mode, encoding='utf-8')
    file_handler.setFormatter(WarningSpacingFormatter("[%(levelname)-7s] %(message)s"))
    file_handler._gmxmmpbsa_master_file = mode == 'w'
    return file_handler


def setup_logging(log_file, master=True, rank=0, *, force=False, file_enabled=True):
    """Configure CLI logging without sharing the master file across MPI ranks."""
    if not force and logging.getLogger().handlers:
        return
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.INFO)
    stream_handler.setFormatter(WarningSpacingFormatter("[%(levelname)-7s] %(message)s"))
    stream_handler.addFilter(_MasterConsoleFilter(master))
    handlers = [stream_handler]
    if master and file_enabled:
        handlers.insert(0, _new_file_handler(log_file, rank))
    elif rank:
        # A non-master rank must never open the master's file. Keep a separate
        # diagnostic file available only if that rank actually reports an error.
        handlers.append(_RankFailureFileHandler(log_file, rank))
    handlers.insert(0, RecordCountingHandler())
    logging.basicConfig(level=logging.DEBUG, handlers=handlers, force=force)


def enable_file_logging(log_file, rank=0):
    """Open the calculation log after parsing confirms a real run is starting."""
    if rank != 0:
        return None
    root = logging.getLogger()
    for handler in root.handlers:
        if getattr(handler, '_gmxmmpbsa_master_file', False):
            return handler
    handler = _new_file_handler(log_file, rank)
    root.addHandler(handler)
    return handler


def format_command_line(args, engine='gmx', mpi_size=1, mpi_requested=False):
    """Format a reproducible command, marking MPI launch details as reconstructed."""
    args = list(args)
    args = [arg for arg in args if arg not in ('mpi', 'MPI')]
    executable = 'amber_MMPBSA' if engine == 'amber' else 'gmx_MMPBSA'
    command = [executable, *args]
    if mpi_requested or mpi_size > 1:
        return shlex.join(['mpirun', '-np', str(mpi_size), *command]) + ' (reconstructed)'
    return shlex.join(command)
