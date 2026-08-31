---
template: main.html
title: Logging and progress
---

# Logging and progress

## Calculation log

During a real calculation, `gmx_MMPBSA.log` is created in the working directory. In an MPI run, rank 0 owns this
master log; non-master ranks do not truncate or append to it. A non-master rank creates a rank-specific diagnostic log
only if it emits an error.

The log uses a compact format:

```text
[INFO   ] Complex completed: 10 frames in 00:01 (9.99 frame/s)
```

The file contains `DEBUG` records as well as normal lifecycle, warning, and error records. Terminal output remains at
`INFO` level so it stays readable. `--help`, `--version`, and other command-line parsing actions do not create or
replace a previous calculation log. A real calculation or `--rewrite-output` starts a fresh log, so concurrent runs
should use separate working directories.

Monitor a batch job with:

```bash
tail -f gmx_MMPBSA.log
```

The final summary reports warning and error records counted directly by the logging system. Multiline messages count as
one record. Fatal exceptions terminate the run; a nonfatal error record does not change the exit status.

## Progress display modes

Select a mode with `--progress-style`:

| Mode | Terminal | `gmx_MMPBSA.log` |
|---|---|---|
| `auto` | Rich in a normal interactive terminal; classic when output is forwarded or hidden by MPI | Clean debug checkpoints and one completion record |
| `rich` | Rich live progress display | Clean debug checkpoints and one completion record |
| `classic` | Original-style ASCII progress bar | Clean debug checkpoints and one completion record |
| `plain` | Text milestones | Text milestones and one completion record |
| `none` | No progress monitor | No progress-monitor records |

Rich and classic renderers do not write animated bars or ANSI control sequences to the calculation log. Progress
checkpoints use MPI-rank terminology and are emitted at bounded percentage milestones. If an external calculation stops
producing frames for the default five-minute stall interval, the log emits a rate-limited waiting notice rather than
repeating a warning on every poll.

## Warnings and errors

Warnings are reserved for conditions that may affect interpretation, accuracy, performance, compatibility, or requested
behavior. Expected automatic actions—such as assigning QM charges, preparing GBNSR6 topology copies, or assigning
missing chain IDs—are informational records. Scientific approximations, fallback paths, incomplete convergence, and
user-value mismatches remain warnings.

Fatal input and setup errors are recorded once with a concise message. Unexpected internal exceptions include traceback
context in the log. By default, a failed `gmx_MMPBSA` or `amber_MMPBSA` calculation creates a diagnostic zip bundle and
records its path so the archive can be attached to a bug report. The bundle can contain logs, input and setup files,
generated intermediates, and samples of up to five trajectory frames; review its contents before sharing it.

Use `--no-error-bundle` to disable automatic bundle creation, for example when input coordinates must not be copied into
an archive:

```bash
gmx_MMPBSA --no-error-bundle -O -i mmpbsa.in ...
```

The same option is accepted by `amber_MMPBSA`. It affects only diagnostic bundle creation after a failure; normal
logging and the calculation exit status are unchanged.

Message text is intended for people, not as a stable machine-readable interface. Scripts should use result files and
exit status rather than parsing logging wording.
