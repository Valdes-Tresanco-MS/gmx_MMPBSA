import hashlib
import shutil
import unittest
from pathlib import Path
from types import SimpleNamespace

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

        # The current dry AMBER workflow must continue to accept a structure
        # without solvent or ions before explicit-water support is added.
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

        self.assertIn('5 frames were processed by cpptraj', log_text)
        self.assertIn('[ERROR  ] = 0', log_text)
        self.assertNotIn(' -cs ', log_text)
        self.assertEqual(
            hashlib.sha256(csv_bytes).hexdigest(),
            '51f945296ddcecc3ddcaac49d24f27ad6e99ee0de12f44c8637f20e74e2267c7',
        )


if __name__ == '__main__':
    unittest.main()
