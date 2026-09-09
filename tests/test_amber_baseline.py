import hashlib
import shutil
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
from unittest.mock import patch

from GMXMMPBSA.exceptions import MMPBSA_Error
from GMXMMPBSA.make_trajs import Trajectory

try:
    import parmed
except ImportError:  # pragma: no cover - exercised by dependency-light test runners
    parmed = None


REPO_ROOT = Path(__file__).resolve().parents[1]
EXAMPLE_DIR = REPO_ROOT / 'examples' / 'AMBER'
CPPTRAJ = shutil.which('cpptraj')
AMBER_RUNTIME = parmed is not None and CPPTRAJ is not None


def _atom_signature(structure):
    return [
        (atom.name.strip(), atom.residue.name.strip())
        for atom in structure.atoms
    ]


def _assert_atom_signatures_match(testcase, actual, expected):
    actual_signature = _atom_signature(actual)
    expected_signature = _atom_signature(expected)
    testcase.assertEqual(len(actual_signature), len(expected_signature))
    for index, (actual_atom, expected_atom) in enumerate(
            zip(actual_signature, expected_signature)):
        with testcase.subTest(atom=index):
            testcase.assertEqual(actual_atom, expected_atom)


@unittest.skipUnless(AMBER_RUNTIME, 'AMBER topology/trajectory tools are unavailable')
class AmberDryBaselineTest(unittest.TestCase):
    def test_checked_in_example_is_a_dry_zero_water_system(self):
        from GMXMMPBSA.make_top_amber import CheckAmberTop

        complex_structure = parmed.load_file(str(EXAMPLE_DIR / '_GMXMMPBSA_COM.pdb'))
        checker = object.__new__(CheckAmberTop)
        checker.complex_str = complex_structure

        # The dry AMBER workflow must continue to accept a structure without
        # solvent or ions alongside the new explicit-water path.
        checker.check4water()

        self.assertFalse(
            [res.name for res in complex_structure.residues
             if res.name.upper() in {'WAT', 'SOL', 'TIP3P', 'NA', 'CL', 'K'}]
        )

    def test_water_residue_remains_rejected_by_current_dry_builder(self):
        from GMXMMPBSA.make_top_amber import CheckAmberTop

        checker = object.__new__(CheckAmberTop)
        checker.complex_str = SimpleNamespace(
            residues=[SimpleNamespace(name='WAT')]
        )

        with self.assertRaises(MMPBSA_Error):
            checker.check4water()

    def test_hoh_residue_is_rejected_by_dry_builder(self):
        from GMXMMPBSA.make_top_amber import CheckAmberTop

        checker = object.__new__(CheckAmberTop)
        checker.complex_str = SimpleNamespace(
            residues=[SimpleNamespace(name='HOH')]
        )

        with self.assertRaises(MMPBSA_Error):
            checker.check4water()

    def test_warns_for_nonconventional_gb_radius_pairing(self):
        from GMXMMPBSA.make_top_amber import CheckAmberTop

        checker = object.__new__(CheckAmberTop)
        checker.INPUT = {'gb': {'gbrun': True, 'igb': 8}}
        parm = SimpleNamespace(
            parm_data={'RADIUS_SET': ['modified Bondi radii (mbondi2)']}
        )

        with self.assertLogs(level='WARNING') as messages:
            checker._warn_gb_radius_compatibility(parm, 'complex')

        self.assertIn("uses 'mbondi2' radii, while igb=8", messages.output[0])
        self.assertIn("with 'mbondi3' radii", messages.output[0])

    def test_does_not_warn_for_conventional_gb_radius_pairing(self):
        from GMXMMPBSA.make_top_amber import CheckAmberTop

        checker = object.__new__(CheckAmberTop)
        checker.INPUT = {'gb': {'gbrun': True, 'igb': 8}}
        parm = SimpleNamespace(
            parm_data={'RADIUS_SET': ['modified Bondi radii (mbondi3)']}
        )

        with self.assertNoLogs(level='WARNING'):
            checker._warn_gb_radius_compatibility(parm, 'complex')

    def test_native_builder_preserves_topology_radii_over_input_pbradii(self):
        from GMXMMPBSA.make_top_amber import CheckAmberTop

        checker = object.__new__(CheckAmberTop)
        checker.INPUT = {'general': {'PBRadii': 4}}
        source = SimpleNamespace(
            atoms=[SimpleNamespace(), SimpleNamespace()],
            parm_data={
                'RADII': [1.1, 1.2],
                'SCREEN': [0.8, 0.9],
                'RADIUS_SET': ['modified Bondi radii (mbondi2)'],
            },
        )
        target = SimpleNamespace(
            atoms=[SimpleNamespace(), SimpleNamespace()],
            parm_data={},
        )

        checker._copy_implicit_radii(source, target, 'Complex')

        self.assertEqual(target.parm_data['RADII'], [1.1, 1.2])
        self.assertEqual(target.parm_data['SCREEN'], [0.8, 0.9])
        self.assertEqual(target.parm_data['RADIUS_SET'], ['modified Bondi radii (mbondi2)'])

    def test_amber_residue_numbers_are_one_based_for_masks(self):
        from GMXMMPBSA.make_top_amber import CheckAmberTop

        residues = [SimpleNamespace(idx=0), SimpleNamespace(idx=1), SimpleNamespace(idx=4)]

        self.assertEqual(CheckAmberTop._amber_residue_numbers(residues), [1, 2, 5])

    def test_solvated_masks_classify_complete_receptor_ligand_water_and_ions(self):
        from GMXMMPBSA.make_top_amber import CheckAmberTop

        def residue(index, name, atom_indexes):
            return SimpleNamespace(
                idx=index,
                name=name,
                atoms=[SimpleNamespace(idx=atom_index) for atom_index in atom_indexes],
            )

        parm = SimpleNamespace(residues=[
            residue(0, 'ALA', [0, 1]),
            residue(1, 'LIG', [2, 3]),
            residue(2, 'WAT', [4, 5, 6]),
            residue(3, 'Na+', [7]),
        ])
        selections = {
            ':1': [1, 1, 0, 0, 0, 0, 0, 0],
            ':2': [0, 0, 1, 1, 0, 0, 0, 0],
        }
        checker = object.__new__(CheckAmberTop)
        checker.FILES = SimpleNamespace(complex_mask=(':1', ':2'))
        checker.explicit_waters = 1
        checker.explicit_waters_group = ''

        with patch('GMXMMPBSA.make_top_amber.parmed.amber.AmberMask',
                   side_effect=lambda parm, mask: SimpleNamespace(Selection=lambda: selections[mask])):
            rec, lig, waters, ions = checker._source_residue_classes(parm)

        self.assertEqual([res.idx for res in rec], [0])
        self.assertEqual([res.idx for res in lig], [1])
        self.assertEqual([res.idx for res in waters], [2])
        self.assertEqual([res.idx for res in ions], [3])

    def test_solvated_mask_must_select_complete_residues(self):
        from GMXMMPBSA.make_top_amber import CheckAmberTop

        residue = SimpleNamespace(
            idx=0,
            name='ALA',
            atoms=[SimpleNamespace(idx=0), SimpleNamespace(idx=1)],
        )
        parm = SimpleNamespace(residues=[residue])
        selections = {':1': [1, 0], ':2': [0, 0]}
        checker = object.__new__(CheckAmberTop)
        checker.FILES = SimpleNamespace(complex_mask=(':1', ':2'))
        checker.explicit_waters = 0
        checker.explicit_waters_group = ''

        with patch('GMXMMPBSA.make_top_amber.parmed.amber.AmberMask',
                   side_effect=lambda parm, mask: SimpleNamespace(Selection=lambda: selections[mask])):
            with self.assertRaises(MMPBSA_Error):
                checker._source_residue_classes(parm)


@unittest.skipUnless(AMBER_RUNTIME, 'AMBER topology/trajectory tools are unavailable')
class AmberExplicitWaterPreprocessingTest(unittest.TestCase):
    def test_cleanup_uses_amber_solvent_masks_and_global_frame_selection(self):
        from GMXMMPBSA.make_top_amber import CheckAmberTop

        class FakeTrajectory:
            def __init__(self, *args):
                self.traj_sizes = [11, 11]

        class FakeProcess:
            cpptraj_input = b''

            def communicate(self, data=None):
                self.__class__.cpptraj_input = data or b''

            def wait(self):
                return 0

        with TemporaryDirectory() as tmpdir:
            prefix = f'{tmpdir}/_GMXMMPBSA_'
            checker = object.__new__(CheckAmberTop)
            checker.FILES = SimpleNamespace(
                prefix=prefix,
                complex_trajs=['first.mdcrd', 'second.mdcrd'],
            )
            checker.INPUT = {'general': {'startframe': 10, 'endframe': 15, 'interval': 2}}
            checker.external_progs = {'cpptraj': 'cpptraj'}
            checker.explicit_waters = 10
            checker.explicit_waters_mask = ':1-166'
            checker.explicit_water_prmtop = f'{prefix}COM_FULL_SOLVENT.prmtop'
            checker.explicit_water_source_ion_mask = ':NA,CL'
            checker.explicit_water_source_all_mask = ':243-11638'
            checker.explicit_water_source_extra_points = ''

            with patch('GMXMMPBSA.make_top_amber.Trajectory', FakeTrajectory):
                with patch('GMXMMPBSA.make_top_amber.subprocess.Popen', return_value=FakeProcess()):
                    checker._cleanup_explicit_water_trajs()

        cpptraj_input = FakeProcess.cpptraj_input.decode()
        self.assertIn('trajin first.mdcrd 10 10 2', cpptraj_input)
        self.assertIn('trajin second.mdcrd 1 3 2', cpptraj_input)
        self.assertIn('strip :NA,CL', cpptraj_input)
        self.assertIn('closest 10 (:1-166)&(!:243-11638) solventmask :243-11638 noimage', cpptraj_input)
        self.assertTrue(checker.FILES.explicit_waters_preselected)

    def test_checked_in_topologies_match_generated_structure_atom_order(self):
        pairs = (
            ('ras-raf_complex.prmtop', '_GMXMMPBSA_COM.pdb'),
            ('ras.prmtop', '_GMXMMPBSA_REC.pdb'),
            ('raf.prmtop', '_GMXMMPBSA_LIG.pdb'),
        )

        for topology_name, structure_name in pairs:
            with self.subTest(system=topology_name):
                topology = parmed.load_file(str(EXAMPLE_DIR / topology_name))
                structure = parmed.load_file(str(EXAMPLE_DIR / structure_name))
                self.assertEqual(len(topology.atoms), len(structure.atoms))
                self.assertEqual(len(topology.residues), len(structure.residues))
                _assert_atom_signatures_match(self, topology, structure)

    def test_first_trajectory_frame_matches_complex_topology_and_structure(self):
        topology_path = EXAMPLE_DIR / 'ras-raf_complex.prmtop'
        trajectory_path = EXAMPLE_DIR / 'prod_complex.mdcrd'
        expected_path = EXAMPLE_DIR / '_GMXMMPBSA_COM.pdb'

        expected = parmed.load_file(str(expected_path))
        trajectory = Trajectory(str(topology_path), str(trajectory_path), CPPTRAJ)
        self.assertEqual(trajectory.total_frames, 5)
        trajectory.Setup(1, 1, 1)

        with self.subTest(selection='first frame'):
            self.assertEqual(trajectory.processed_frames, 1)

        with self.subTest(selection='all five frames'):
            trajectory.Setup(1, 5, 1)
            self.assertEqual(trajectory.processed_frames, 5)

        # Re-query frame 1 in a temporary output and compare the complete
        # atom/residue sequence, not only the atom count.
        import tempfile
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / 'frame1.pdb'
            cpptraj_log = Path(tmpdir) / 'cpptraj.out'
            trajectory.Setup(1, 1, 1)
            trajectory.Outtraj(str(output_path), frames='1', filetype='pdb')
            trajectory.Run(str(cpptraj_log))
            extracted = parmed.load_file(str(output_path))

        self.assertEqual(len(extracted.atoms), len(expected.atoms))
        _assert_atom_signatures_match(self, extracted, expected)

    def test_checked_in_result_is_the_current_five_frame_baseline(self):
        log_text = (EXAMPLE_DIR / 'gmx_MMPBSA.log').read_text()
        csv_bytes = (EXAMPLE_DIR / 'FINAL_RESULTS_MMPBSA.csv').read_bytes()

        self.assertIn('5 frames were prepared by cpptraj', log_text)
        self.assertIn('Run completed with 0 errors and 0 warnings.', log_text)
        self.assertNotIn(' -cs ', log_text)
        # Reproduced byte-for-byte with local master 5b9f2f95 and the repaired branch
        # using the bundled five frames and native topology radii (igb=1).
        self.assertEqual(
            hashlib.sha256(csv_bytes).hexdigest(),
            'd95c50d03a418e49dd0dea2955aab60fb957e9689ed0417837740c82a0b28801',
        )


if __name__ == '__main__':
    unittest.main()
