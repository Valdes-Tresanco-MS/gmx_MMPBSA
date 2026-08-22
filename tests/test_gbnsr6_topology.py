import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from GMXMMPBSA.gbnsr6_topology import prepare_gbnsr6_topology


class GBNSR6TopologyTest(unittest.TestCase):
    def test_strips_dihedral_pointers_in_prmtop_copy(self):
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            source = root / 'COM.prmtop'
            output = root / 'COM_gbnsr6.prmtop'

            source.write_text(
                '%VERSION  VERSION_STAMP = V0001.000\n'
                '%FLAG POINTERS\n'
                '%FORMAT(10I8)\n'
                '       4       2       1       2       3       4       5       6       0       0\n'
                '       7       1       2       3       4       5       6       7       1       0\n'
                '       0       0       0       0       0       0       0       0       1       0\n'
                '       0\n'
                '%FLAG ATOM_NAME\n'
                '%FORMAT(20a4)\n'
                'C   H   H   H   \n'
                '%FLAG DIHEDRAL_FORCE_CONSTANT\n'
                '%FORMAT(5E16.8)\n'
                '  1.00000000E+00  2.00000000E+00\n'
                '%FLAG SCEE_SCALE_FACTOR\n'
                '%FORMAT(5E16.8)\n'
                '  1.20000000E+00  1.20000000E+00\n'
                '%FLAG DIHEDRALS_INC_HYDROGEN\n'
                '%FORMAT(10I8)\n'
                '       1       2       3       4       1\n'
                '%FLAG CHARGE\n'
                '%FORMAT(5E16.8)\n'
                '  0.00000000E+00  0.00000000E+00\n'
            )

            prepared = Path(prepare_gbnsr6_topology(source, output))
            text = prepared.read_text()

        self.assertIn('%FLAG ATOM_NAME\n%FORMAT(20a4)\nC   H   H   H   \n', text)
        self.assertIn('%FLAG CHARGE\n%FORMAT(5E16.8)\n  0.00000000E+00', text)
        self.assertIn(
            '       4       2       1       2       3       4       0       0       0       0\n'
            '       7       1       2       3       0       5       6       0       1       0\n',
            text,
        )
        self.assertIn('%FLAG DIHEDRAL_FORCE_CONSTANT\n%FORMAT(5E16.8)\n', text)
        self.assertIn('%FLAG SCEE_SCALE_FACTOR\n%FORMAT(5E16.8)\n', text)
        self.assertIn('%FLAG DIHEDRALS_INC_HYDROGEN\n%FORMAT(10I8)\n', text)
        self.assertIn('  1.00000000E+00  2.00000000E+00', text)
        self.assertIn('       1       2       3       4       1', text)


if __name__ == '__main__':
    unittest.main()
