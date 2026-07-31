import logging
import unittest
from types import SimpleNamespace
from unittest.mock import patch

from GMXMMPBSA.exceptions import MMPBSA_Error
from GMXMMPBSA.make_trajs import make_trajectories


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
        FakeTrajectory.instances.append(self)

    def Setup(self, *args):
        self.setup_args = args
        self.processed_frames = self.processed_frame_count

    def rms(self, mask):
        pass

    def Strip(self, mask):
        pass

    def Unstrip(self, restrip_solvent=True):
        pass

    def Outtraj(self, fname, **kwargs):
        self.outputs.append((fname, kwargs))

    def Run(self, output):
        pass


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
