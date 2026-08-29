import logging
import os
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from GMXMMPBSA.exceptions import MMPBSA_Error
from GMXMMPBSA.make_trajs import make_trajectories, warn_concatenated_complex_trajectories


def _input():
    return {
        'general': {
            'netcdf': False,
            'startframe': 1,
            'endframe': 999999,
            'interval': 1,
            'full_traj': False,
            'qh_entropy': False,
            'ligand_mask': ':LIG',
            'receptor_mask': ':REC',
        },
        'gbnsr6': {'gbnsr6run': False},
        'nmode': {'nmoderun': False},
    }


def _files():
    return SimpleNamespace(
        stability=False,
        complex_prmtop='complex.prmtop',
        receptor_prmtop='receptor.prmtop',
        ligand_prmtop='ligand.prmtop',
        complex_trajs=['complex.mdcrd'],
        original_complex_trajs=['full.mdcrd'],
        complex_tpr='full.pdb',
        receptor_trajs=[],
        ligand_trajs=[],
    )


class FakeTrajectory:
    processed_frame_count = 3
    instances = []

    def __init__(self, prmtop, traj_files, cpptraj='cpptraj'):
        self.prmtop = prmtop
        self.traj_files = traj_files
        self.outputs = []
        self.processed_frames = 0
        self.total_frames = self.processed_frame_count
        self.setup_args = None
        self.membrane_output_dir = None
        FakeTrajectory.instances.append(self)

    def Setup(self, *args):
        self.setup_args = args
        self.processed_frames = self.processed_frame_count

    def rms(self, mask):
        pass

    def ExtractMembraneAtoms(self, atom_names, output_dir):
        self.membrane_output_dir = Path(output_dir)
        self.membrane_atom_names = atom_names

    def Strip(self, mask):
        pass

    def Unstrip(self, restrip_solvent=True):
        pass

    def Outtraj(self, fname, **kwargs):
        self.outputs.append((fname, kwargs))

    def Run(self, output):
        if self.membrane_output_dir is not None:
            for frame, z in enumerate(((-20.0, 20.0), (-21.0, 21.0)), start=1):
                (self.membrane_output_dir / f'P.pdb.{frame}').write_text(
                    ''.join(
                        f'HETATM    1  P   POPC A   1       1.000   2.000  {value:6.3f}  1.00  0.00           P\n'
                        for value in z
                    )
                )


class MakeTrajectoriesMPIFrameTest(unittest.TestCase):
    def setUp(self):
        FakeTrajectory.instances = []
        FakeTrajectory.processed_frame_count = 3

    def test_limits_active_ranks_to_selected_frames(self):
        with patch('GMXMMPBSA.make_trajs.Trajectory', FakeTrajectory):
            logging.disable(logging.CRITICAL)
            try:
                com_frames, rec_frames, lig_frames, nmode_frames, mpi_size = make_trajectories(
                    _input(), _files(), 4, 'cpptraj', '_GMXMMPBSA_'
                )
            finally:
                logging.disable(logging.NOTSET)

        outputs = [fname for fname, _ in FakeTrajectory.instances[0].outputs]

        self.assertEqual((com_frames, rec_frames, lig_frames, nmode_frames, mpi_size), (3, 3, 3, 0, 3))
        self.assertIn('_GMXMMPBSA_complex.mdcrd.2', outputs)
        self.assertNotIn('_GMXMMPBSA_complex.mdcrd.3', outputs)

    def test_rejects_zero_selected_frames(self):
        FakeTrajectory.processed_frame_count = 0

        with patch('GMXMMPBSA.make_trajs.Trajectory', FakeTrajectory):
            with self.assertRaises(MMPBSA_Error):
                make_trajectories(_input(), _files(), 4, 'cpptraj', '_GMXMMPBSA_')

    def test_preselected_explicit_water_trajs_are_not_subsampled_again(self):
        files = _files()
        files.explicit_waters_preselected = True
        inp = _input()
        inp['general']['startframe'] = 2
        inp['general']['endframe'] = 20
        inp['general']['interval'] = 3

        with patch('GMXMMPBSA.make_trajs.Trajectory', FakeTrajectory):
            make_trajectories(inp, files, 1, 'cpptraj', '_GMXMMPBSA_')

        self.assertEqual(FakeTrajectory.instances[0].setup_args, (1, 3, 1))

    def test_automatic_membrane_parameters_use_the_original_full_trajectory(self):
        inp = _input()
        inp['pb'] = {
            'memopt': 1,
            'mctrdz': 'automatic',
            'mthick': 'automatic',
            'membrane_atoms': 'P',
        }
        with tempfile.TemporaryDirectory() as tmpdir:
            old_cwd = os.getcwd()
            os.chdir(tmpdir)
            try:
                with patch('GMXMMPBSA.make_trajs.Trajectory', FakeTrajectory):
                    make_trajectories(inp, _files(), 1, 'cpptraj', '_GMXMMPBSA_')
            finally:
                os.chdir(old_cwd)

            self.assertEqual(inp['pb']['mctrdz'], 0.0)
            self.assertEqual(inp['pb']['mthick'], 41.0)
            self.assertTrue(Path(tmpdir, 'GMXMMPBSA_membrane_parameters.csv').exists())
            self.assertTrue(Path(tmpdir, 'GMXMMPBSA_membrane_parameters.png').exists())

        self.assertEqual(len(FakeTrajectory.instances), 2)
        self.assertEqual(FakeTrajectory.instances[0].prmtop, 'complex.prmtop')
        self.assertEqual(FakeTrajectory.instances[0].traj_files, ['complex.mdcrd'])
        self.assertEqual(FakeTrajectory.instances[1].prmtop, 'full.pdb')
        self.assertEqual(FakeTrajectory.instances[1].traj_files, ['full.mdcrd'])

    def test_multiple_complex_trajectories_warn_about_pooled_statistics(self):
        with self.assertLogs(level=logging.WARNING) as captured:
            warn_concatenated_complex_trajectories(['complex_0.mdcrd', 'complex_1.mdcrd'])

        self.assertEqual(len(captured.records), 1)
        self.assertIn('one concatenated trajectory', captured.records[0].getMessage())
        self.assertIn('not independent-replica statistics', captured.records[0].getMessage())

    def test_multiple_unbound_trajectories_use_their_cli_option_in_warning(self):
        with self.assertLogs(level=logging.WARNING) as captured:
            warn_concatenated_complex_trajectories(
                ['ligand_0.mdcrd', 'ligand_1.mdcrd'], option='-lt', label='ligand'
            )

        self.assertEqual(len(captured.records), 1)
        self.assertIn('Multiple ligand trajectories were supplied with -lt',
                      captured.records[0].getMessage())
