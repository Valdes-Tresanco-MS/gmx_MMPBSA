import unittest
from pathlib import Path

from GMXMMPBSA.input_parser import SUPPORTED_QM_THEORIES


ROOT = Path(__file__).resolve().parents[1]


class QMMMGBSASupportDocumentationTest(unittest.TestCase):
    def test_input_documentation_lists_exact_canonical_theories(self):
        documentation = (ROOT / 'docs' / 'input_file.md').read_text()
        for theory in SUPPORTED_QM_THEORIES:
            with self.subTest(theory=theory):
                self.assertIn(f'`{theory}`', documentation)

    def test_input_documentation_does_not_advertise_rejected_aliases(self):
        documentation = (ROOT / 'docs' / 'input_file.md').read_text()
        accepted_section = documentation.split('The dispersion correction', 1)[0]
        for theory in ('PDDG-PM3', 'PM3PDDG', 'PDDG-MNDO', 'PDDGMNDO', 'SCC-DFTB',
                       'PM3-ZnB', 'PM3ZNB', 'MNDOD'):
            with self.subTest(theory=theory):
                self.assertNotIn(f'`{theory}`', accepted_section)

    def test_qmmm_example_matches_its_checked_in_input(self):
        documentation = (ROOT / 'docs' / 'examples' / 'QM_MMGBSA' / 'README.md').read_text()
        example = (ROOT / 'examples' / 'QM_MMGBSA' / 'mmpbsa.in').read_text()
        self.assertIn('ifqnt=1, qm_theory=PM6-DH+,', documentation)
        self.assertIn('ifqnt=1, qm_theory=PM6-DH+,', example)
        self.assertIn('If `qm_theory` is omitted', documentation)

    def test_qmmm_documentation_describes_scf_iteration_limit(self):
        documentation = (ROOT / 'docs' / 'input_file.md').read_text()
        self.assertIn('`itrmax` (Default = 1000)', documentation)
        self.assertIn('`ndiis_attempts` (Default = None)', documentation)
        self.assertIn('calculation stops rather than including the unconverged energy', documentation)

    def test_qmmm_documentation_distinguishes_validation_from_runtime_support(self):
        documentation = (ROOT / 'docs' / 'input_file.md').read_text()
        self.assertIn('one-frame QM/MMGBSA run', documentation)
        self.assertIn('not a guarantee of universal support', documentation)
        self.assertIn('PM3-MAIS', documentation)
        self.assertIn('missing PM3-MAIS parameters for nitrogen', documentation)
        self.assertIn('gmx_MMPBSA` validates the method name and writes the QM/MM input', documentation)

    def test_qmmm_documentation_describes_runtime_diagnostics(self):
        documentation = (ROOT / 'docs' / 'input_file.md').read_text()
        for phrase in (
            'Fatal diagnostics such as SCF',
            'nonconvergence, missing Hamiltonian parameters',
            'missing Hamiltonian parameters',
            'missing dispersion-correction parameters',
            'missing DFTB',
            'numerical-derivative message for d orbitals',
            'final binding result is not considered valid',
        ):
            with self.subTest(phrase=phrase):
                self.assertIn(phrase, documentation)


if __name__ == '__main__':
    unittest.main()
