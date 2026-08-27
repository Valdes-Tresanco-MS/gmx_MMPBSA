import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

from GMXMMPBSA.gbnsr6_topology import _compact_equivalent_lj_types, prepare_gbnsr6_topology


class GBNSR6TopologyTest(unittest.TestCase):
    def test_compaction_preserves_solty_length(self):
        class FakeAmberFormat:
            def __init__(self, _path):
                type(self).last_instance = self
                self.parm_data = {
                    'POINTERS': [4, 3] + [0] * 29,
                    'ATOM_TYPE_INDEX': [1, 2, 1, 2],
                    'NONBONDED_PARM_INDEX': [1] * 9,
                    'LENNARD_JONES_ACOEF': [1.0] * 6,
                    'LENNARD_JONES_BCOEF': [1.0] * 6,
                    'SOLTY': [0.0] * 4,
                }

            def write_parm(self, path):
                Path(path).write_text('prepared')

        with TemporaryDirectory() as tmpdir:
            source = Path(tmpdir) / 'LIG.prmtop'
            source.write_text('source')
            with patch('parmed.amber.AmberFormat', FakeAmberFormat):
                _compact_equivalent_lj_types(source)

        self.assertEqual(FakeAmberFormat.last_instance.parm_data['SOLTY'], [0.0] * 4)

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
