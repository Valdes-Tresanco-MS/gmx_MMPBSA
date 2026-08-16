"""Logging setup shared by the command-line entry points."""

import logging
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
            self._file_handler = logging.FileHandler(self.log_file, mode='a')
            self._file_handler.setFormatter(WarningSpacingFormatter("[%(levelname)-7s] %(message)s"))
        self._file_handler.emit(record)

    def close(self):
        if self._file_handler is not None:
            self._file_handler.close()
            self._file_handler = None
        super().close()


def setup_logging(log_file, master=True, rank=0, *, force=False):
    """Configure CLI logging without sharing the master file across MPI ranks."""
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.INFO)
    stream_handler.setFormatter(WarningSpacingFormatter("[%(levelname)-7s] %(message)s"))
    handlers = [stream_handler]
    if master:
        file_handler = logging.FileHandler(log_file, 'w')
        file_handler.setFormatter(WarningSpacingFormatter("[%(levelname)-7s] %(message)s"))
        handlers.insert(0, file_handler)
    elif rank:
        # A non-master rank must never open the master's file. Keep a separate
        # diagnostic file available only if that rank actually reports an error.
        handlers.append(_RankFailureFileHandler(log_file, rank))
    logging.basicConfig(level=logging.DEBUG, handlers=handlers, force=force)
