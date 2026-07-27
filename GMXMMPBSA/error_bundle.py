"""Best-effort support bundle creation for failed calculations."""

import json
import os
import platform
import re
import subprocess
import sys
import tempfile
import traceback
import zipfile
from datetime import datetime
from pathlib import Path

from GMXMMPBSA import __version__


MAX_REGULAR_FILE_SIZE = 25 * 1024 * 1024
TRAJ_ATTRS = {'complex_trajs', 'receptor_trajs', 'ligand_trajs'}
TRAJ_SUFFIXES = {'.xtc', '.trr', '.mdcrd', '.nc', '.netcdf', '.dcd', '.crd'}


def create_error_bundle(app, exc, tb=None, max_frames=5):
    """Create a zip file with useful inputs and a tiny trajectory sample."""
    files = getattr(app, 'FILES', None)
    prefix = getattr(files, 'prefix', '_GMXMMPBSA_')
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    bundle = Path(f'gmx_MMPBSA_error_bundle_{timestamp}.zip').absolute()
    manifest = _base_manifest(app, exc, tb, max_frames)
    added = set()

    with tempfile.TemporaryDirectory(prefix='gmx_MMPBSA_error_bundle_') as tmpdir:
        tmpdir = Path(tmpdir)
        with zipfile.ZipFile(bundle, 'w', zipfile.ZIP_DEFLATED) as zf:
            _add_existing_file(zf, Path('gmx_MMPBSA.log'), 'logs/gmx_MMPBSA.log', added, manifest)

            if files is not None:
                _add_referenced_files(zf, files, added, manifest)
                _add_topology_includes(zf, files, added, manifest)

            _add_generated_files(zf, prefix, added, manifest)
            _add_trajectory_samples(zf, app, tmpdir, max_frames, added, manifest)
            zf.writestr('manifest.json', json.dumps(manifest, indent=2, sort_keys=True))

    return bundle


def _base_manifest(app, exc, tb, max_frames):
    tb_text = ''.join(traceback.format_exception(type(exc), exc, tb or exc.__traceback__))
    return {
        'bundle_format': 1,
        'created_at': datetime.now().isoformat(timespec='seconds'),
        'cwd': str(Path.cwd()),
        'command_line': sys.argv,
        'gmx_mmpbsa_version': __version__,
        'python': sys.version.replace('\n', ' '),
        'platform': platform.platform(),
        'mpi_rank': getattr(app, 'mpi_rank', 0),
        'mpi_size': getattr(app, 'mpi_size', 1),
        'max_trajectory_frames': max_frames,
        'error': {
            'type': type(exc).__name__,
            'message': str(exc),
            'traceback': tb_text,
        },
        'files': [],
        'skipped_files': [],
        'trajectory_samples': [],
        'notes': [],
    }


def _add_referenced_files(zf, files, added, manifest):
    for attr, value in sorted(vars(files).items()):
        if attr.startswith('_') or attr in TRAJ_ATTRS:
            continue
        for path in _iter_paths(value):
            _add_existing_file(zf, path, f'inputs/{attr}/{path.name}', added, manifest)


def _add_topology_includes(zf, files, added, manifest):
    for attr in ('complex_top', 'receptor_top', 'ligand_top'):
        top = getattr(files, attr, None)
        if not top:
            continue
        for include in _topology_includes(Path(top)):
            _add_existing_file(zf, include, f'inputs/topology_includes/{include.name}', added, manifest)


def _add_generated_files(zf, prefix, added, manifest):
    generated_suffixes = {
        '.pdb', '.prmtop', '.inpcrd', '.ndx', '.in', '.out', '.dat', '.csv', '.info', '.json'
    }
    for path in Path('.').glob(f'{prefix}*'):
        if path.is_file() and path.suffix in generated_suffixes:
            _add_existing_file(zf, path, f'generated/{path.name}', added, manifest)
    for path in Path('.').glob('COMPACT_MMXSA_RESULTS.mmxsa'):
        _add_existing_file(zf, path, f'generated/{path.name}', added, manifest)
    for path in Path('.').glob('*.prmtop'):
        if path.name.startswith(('COM', 'REC', 'LIG', 'MUT_COM', 'MUT_REC', 'MUT_LIG')):
            _add_existing_file(zf, path, f'generated/{path.name}', added, manifest)


def _add_trajectory_samples(zf, app, tmpdir, max_frames, added, manifest):
    files = getattr(app, 'FILES', None)
    if files is None:
        return
    external_progs = getattr(app, 'external_progs', {})
    cpptraj = external_progs.get('cpptraj') if isinstance(external_progs, dict) else None
    trjconv = _command_args(external_progs.get('trjconv')) if isinstance(external_progs, dict) else None
    gmx_check = _gmx_check_command(trjconv)

    specs = [
        ('complex', _first_existing(getattr(files, 'complex_prmtop', None), 'COM.prmtop'),
         getattr(files, 'complex_trajs', None), getattr(files, 'complex_tpr', None)),
        ('receptor', _first_existing(getattr(files, 'receptor_prmtop', None), 'REC.prmtop'),
         getattr(files, 'receptor_trajs', None), getattr(files, 'receptor_tpr', None)),
        ('ligand', _first_existing(getattr(files, 'ligand_prmtop', None), 'LIG.prmtop'),
         getattr(files, 'ligand_trajs', None), getattr(files, 'ligand_tpr', None)),
    ]
    sample_created = False
    for label, prmtop, trajs, tpr in specs:
        if not trajs:
            continue
        for index, traj in enumerate(trajs):
            traj = Path(traj)
            if not traj.is_file():
                continue
            manifest['skipped_files'].append({'path': str(traj), 'reason': 'Full trajectory omitted.'})
            sample = tmpdir / f'{label}_{index}_first_{max_frames}_frames.mdcrd'
            if cpptraj and prmtop and _write_cpptraj_sample(cpptraj, Path(prmtop), traj, sample, max_frames):
                _add_existing_file(zf, sample, f'trajectory_samples/{sample.name}', added, manifest)
                manifest['trajectory_samples'].append({
                    'source': str(traj),
                    'sample': f'trajectory_samples/{sample.name}',
                    'frames': max_frames,
                    'method': 'cpptraj',
                })
                sample_created = True
                continue
            sample = tmpdir / f'{label}_{index}_first_{max_frames}_frames.xtc'
            if trjconv and gmx_check and tpr and _write_gmx_sample(gmx_check, trjconv, Path(tpr), traj, sample, max_frames):
                _add_existing_file(zf, sample, f'trajectory_samples/{sample.name}', added, manifest)
                manifest['trajectory_samples'].append({
                    'source': str(traj),
                    'sample': f'trajectory_samples/{sample.name}',
                    'frames': max_frames,
                    'method': 'gmx trjconv',
                })
                sample_created = True
            else:
                manifest['skipped_files'].append({
                    'path': str(traj),
                    'reason': 'Unable to create trajectory sample with cpptraj or gmx trjconv.',
                })
    if not sample_created and not cpptraj and not trjconv:
        manifest['notes'].append('cpptraj/trjconv were not available; trajectory samples were not created.')


def _write_cpptraj_sample(cpptraj, prmtop, traj, sample, max_frames):
    cpptraj_input = f'trajin {traj.as_posix()} 1 {max_frames} 1\ntrajout {sample.as_posix()}\nrun\n'
    try:
        process = subprocess.run(
            [cpptraj, prmtop.as_posix()],
            input=cpptraj_input.encode(),
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=False,
        )
    except OSError:
        return False
    return process.returncode == 0 and sample.is_file() and sample.stat().st_size > 0


def _write_gmx_sample(gmx_check, trjconv, tpr, traj, sample, max_frames):
    times = _first_gmx_frame_times(gmx_check, traj, max_frames)
    if not times:
        return False
    try:
        process = subprocess.run(
            trjconv + [
                '-s', tpr.as_posix(),
                '-f', traj.as_posix(),
                '-o', sample.as_posix(),
                '-b', str(times[0]),
                '-e', str(times[-1]),
            ],
            input=b'0\n',
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=False,
        )
    except OSError:
        return False
    return process.returncode == 0 and sample.is_file() and sample.stat().st_size > 0


def _first_gmx_frame_times(gmx_check, traj, max_frames):
    try:
        process = subprocess.run(
            gmx_check + ['-f', traj.as_posix()],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=False,
        )
    except OSError:
        return []
    output = process.stdout.decode(errors='replace')
    times = []
    for match in re.finditer(r'Reading frame\s+\d+\s+time\s+([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)', output):
        times.append(float(match.group(1)))
        if len(times) == max_frames:
            break
    return times


def _command_args(command):
    if command is None:
        return None
    if isinstance(command, (list, tuple)):
        return [str(part) for part in command]
    return [str(command)]


def _gmx_check_command(trjconv):
    if not trjconv:
        return None
    if len(trjconv) > 1:
        return [trjconv[0], 'check']
    return None


def _first_existing(*paths):
    for path in paths:
        if path and Path(path).is_file():
            return Path(path)
    return None


def _add_existing_file(zf, path, arcname, added, manifest):
    path = Path(path)
    if not path.is_file():
        return
    resolved = path.resolve()
    if resolved in added:
        return
    if path.suffix.lower() in TRAJ_SUFFIXES and not arcname.startswith('trajectory_samples/'):
        manifest['skipped_files'].append({'path': str(path), 'reason': 'Full trajectory omitted.'})
        return
    size = path.stat().st_size
    if size > MAX_REGULAR_FILE_SIZE:
        manifest['skipped_files'].append({'path': str(path), 'reason': f'File larger than {MAX_REGULAR_FILE_SIZE} bytes.'})
        return
    zf.write(path, arcname)
    added.add(resolved)
    manifest['files'].append({'path': str(path), 'archive_name': arcname, 'size': size})


def _iter_paths(value):
    if value is None or isinstance(value, bool):
        return
    if isinstance(value, (str, os.PathLike)):
        path = Path(value)
        if path.is_file():
            yield path
        return
    if isinstance(value, (list, tuple, set)):
        for item in value:
            yield from _iter_paths(item)


def _topology_includes(topology):
    seen = set()
    pending = [topology]
    include_re = re.compile(r'^\s*#include\s+[<"]([^>"]+)[>"]')
    while pending:
        top = pending.pop()
        if not top.is_file():
            continue
        try:
            lines = top.read_text(errors='replace').splitlines()
        except OSError:
            continue
        for line in lines:
            match = include_re.match(line)
            if not match:
                continue
            include = (top.parent / match.group(1)).resolve()
            if include in seen or not include.is_file():
                continue
            seen.add(include)
            pending.append(include)
            yield include
