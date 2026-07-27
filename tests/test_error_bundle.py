import json
import os
import unittest
import zipfile
from contextlib import contextmanager
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
from unittest.mock import patch

from GMXMMPBSA.commandlineparser import parser
from GMXMMPBSA.error_bundle import create_error_bundle


@contextmanager
def working_directory(path):
    old_cwd = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(old_cwd)


class ErrorBundleTest(unittest.TestCase):
    def test_creates_bundle_with_manifest_and_inputs_without_full_trajectory(self):
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            input_file = root / 'mmpbsa.in'
            topology = root / 'topol.top'
            include = root / 'forcefield.itp'
            traj = root / 'traj.xtc'
            log = root / 'gmx_MMPBSA.log'
            prmtop = root / 'COM.prmtop'
            input_file.write_text('&general\n/\n')
            topology.write_text('#include "forcefield.itp"\n')
            include.write_text('[ defaults ]\n')
            traj.write_text('not a real trajectory')
            log.write_text('log text')
            prmtop.write_text('%VERSION\n')

            app = SimpleNamespace(
                FILES=SimpleNamespace(
                    prefix='_GMXMMPBSA_',
                    input_file=input_file.as_posix(),
                    complex_top=topology.as_posix(),
                    complex_trajs=[traj.as_posix()],
                    no_error_bundle=False,
                ),
                external_progs={},
                mpi_rank=0,
                mpi_size=1,
            )

            with working_directory(root), patch('GMXMMPBSA.error_bundle.shutil.which', return_value=None):
                bundle = create_error_bundle(app, RuntimeError('boom'))

            with zipfile.ZipFile(bundle) as zf:
                names = set(zf.namelist())
                manifest = json.loads(zf.read('manifest.json').decode())

            self.assertIn('manifest.json', names)
            self.assertIn('logs/gmx_MMPBSA.log', names)
            self.assertIn('inputs/input_file/mmpbsa.in', names)
            self.assertIn('inputs/complex_top/topol.top', names)
            self.assertIn('inputs/topology_includes/forcefield.itp', names)
            self.assertIn('generated/COM.prmtop', names)
            self.assertNotIn('inputs/complex_trajs/traj.xtc', names)
            self.assertIn('Trajectory sampling tools were not discovered', ' '.join(manifest['notes']))

    def test_parser_accepts_error_bundle_opt_out(self):
        args = parser.parse_args(['--no-error-bundle', '--create_input', 'gb'])

        self.assertTrue(args.no_error_bundle)


if __name__ == '__main__':
    unittest.main()
