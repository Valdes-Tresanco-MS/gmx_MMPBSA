import unittest

from GMXMMPBSA.utils import residue_names_match


class ReferenceResidueMatchingTest(unittest.TestCase):
    def test_terminal_prefixes_match_base_residue_names(self):
        self.assertTrue(residue_names_match('NVAL', 'VAL'))
        self.assertTrue(residue_names_match('ALA', 'CALA'))

    def test_unrelated_residue_names_do_not_match(self):
        self.assertFalse(residue_names_match('VAL', 'ASP'))
        self.assertFalse(residue_names_match('NVAL', 'NALA'))


if __name__ == '__main__':
    unittest.main()
