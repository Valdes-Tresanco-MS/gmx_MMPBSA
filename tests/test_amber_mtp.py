import unittest
from types import SimpleNamespace
from unittest.mock import patch

from GMXMMPBSA.exceptions import MMPBSA_Error
from GMXMMPBSA.make_top_amber import CheckAmberTop


def _files(**overrides):
    values = {
        'complex_top': 'complex.prmtop',
        'complex_trajs': ['complex.mdcrd'],
        'complex_mask': (':1-10', ':11-12'),
        'receptor_top': None,
        'receptor_mask': None,
        'receptor_trajs': None,
        'ligand_top': None,
        'ligand_mask': None,
        'ligand_trajs': None,
    }
    values.update(overrides)
    return SimpleNamespace(**values)


class AmberMTPInputValidationTest(unittest.TestCase):
    def test_ligand_mtp_requires_topology(self):
        checker = object.__new__(CheckAmberTop)
        checker.FILES = _files(ligand_trajs=['ligand.mdcrd'])

        with self.assertRaisesRegex(
                MMPBSA_Error,
                r'ligand topology \(-lp\) is required'):
            checker.checkFiles()

    def test_ligand_mtp_requires_mask(self):
        checker = object.__new__(CheckAmberTop)
        checker.FILES = _files(
            ligand_top='ligand.prmtop', ligand_trajs=['ligand.mdcrd']
        )

        with self.assertRaisesRegex(
                MMPBSA_Error,
                r'ligand mask \(-lm\) is required'):
            checker.checkFiles()

    def test_receptor_mtp_requires_topology(self):
        checker = object.__new__(CheckAmberTop)
        checker.FILES = _files(receptor_trajs=['receptor.mdcrd'])

        with self.assertRaisesRegex(
                MMPBSA_Error,
                r'receptor topology \(-rp\) is required'):
            checker.checkFiles()

    def test_receptor_mtp_requires_mask(self):
        checker = object.__new__(CheckAmberTop)
        checker.FILES = _files(
            receptor_top='receptor.prmtop', receptor_trajs=['receptor.mdcrd']
        )

        with self.assertRaisesRegex(
                MMPBSA_Error,
                r'receptor mask \(-rm\) is required'):
            checker.checkFiles()

    def test_complete_mtp_inputs_pass_validation(self):
        checker = object.__new__(CheckAmberTop)
        checker.FILES = _files(
            receptor_top='receptor.prmtop', receptor_mask=':1-10',
            receptor_trajs=['receptor.mdcrd'],
            ligand_top='ligand.prmtop', ligand_mask=':1-2',
            ligand_trajs=['ligand.mdcrd'],
        )

        checker.checkFiles()


class AmberMTPStructureExtractionTest(unittest.TestCase):
    def test_extraction_passes_all_unbound_trajectories_to_cpptraj(self):
        class FakeTrajectory:
            instance = None

            def __init__(self, topology, trajectories):
                self.topology = topology
                self.trajectories = trajectories
                self.actions = []
                FakeTrajectory.instance = self

            def Setup(self):
                self.actions.append(('setup',))

            def Strip(self, mask):
                self.actions.append(('strip', mask))

            def Outtraj(self, output, frames=None, filetype=''):
                self.actions.append(('outtraj', output, frames, filetype))

            def Run(self, output):
                self.actions.append(('run', output))

        with patch('GMXMMPBSA.make_top_amber.Trajectory', FakeTrajectory):
            CheckAmberTop._extract_mtp_structure(
                'ligand.prmtop', ['ligand_0.mdcrd', 'ligand_1.mdcrd'],
                ':1-2', 'LIG.pdb', 'ligand'
            )

        self.assertEqual(FakeTrajectory.instance.topology, 'ligand.prmtop')
        self.assertEqual(
            FakeTrajectory.instance.trajectories,
            ['ligand_0.mdcrd', 'ligand_1.mdcrd'],
        )
        self.assertEqual(
            FakeTrajectory.instance.actions,
            [
                ('setup',),
                ('strip', '!:1-2'),
                ('outtraj', 'LIG.pdb', '1', 'pdb'),
                ('run', 'ligand_pdb.out'),
            ],
        )


if __name__ == '__main__':
    unittest.main()
