import io
import sys
import types
import unittest
from unittest.mock import patch


if 'pandas' not in sys.modules:
    pandas = types.ModuleType('pandas')
    pandas.MultiIndex = type('MultiIndex', (), {})
    pandas.Index = type('Index', (), {})
    sys.modules['pandas'] = pandas
if 'parmed' not in sys.modules:
    parmed = types.ModuleType('parmed')
    parmed.__version__ = 'stub'
    sys.modules['parmed'] = parmed

from GMXMMPBSA.calculation import ListEnergyCalculation, MergeGBNSR6Output


class ListEnergyCalculationTest(unittest.TestCase):
    def test_gbnsr6_postprocess_uses_original_topology(self):
        calc = ListEnergyCalculation(
            'gbnsr6', 'stripped.prmtop', 'gbnsr6.mdin',
            ['frame.inpcrd'], ['frame.mdout'],
            postprocess_prmtop='original.prmtop', keep_mdouts=True
        )
        calc.setup()

        class FakeProcess:
            def __init__(self, command_args, **kwargs):
                self.command_args = command_args

            def wait(self):
                return 0

        with patch('subprocess.Popen', FakeProcess), patch('GMXMMPBSA.calculation.mdout2json') as mdout2json:
            calc.run(0, stdout=io.StringIO(), stderr=io.StringIO())

        command_args = calc.list_calc[0]
        self.assertEqual(command_args[command_args.index('-p') + 1], 'stripped.prmtop')
        postprocess_args = mdout2json.call_args.args[0]
        self.assertEqual(postprocess_args[postprocess_args.index('-p') + 1], 'original.prmtop')
        self.assertTrue(mdout2json.call_args.kwargs['keep_mdout'])

    def test_gbnsr6_merge_parses_chamber_energy_row(self):
        results = [
            'minimizing coord set #       1\n',
            ' BOND = 1.0 ANGLE = 2.0 DIHED = 3.0\n',
            ' UB = 4.0 IMP = 5.0 CMAP = 6.0\n',
            ' VDWAALS = 7.0 EEL = 8.0 EGB = 9.0\n',
            ' 1-4 VDW = 10.0 1-4 EEL = 11.0 RESTRAINT = 0.0\n',
            'minimization completed\n',
        ]

        parsed = MergeGBNSR6Output._get_energy_decomp(results)['energy'][1]

        self.assertEqual(parsed['UB'], 4.0)
        self.assertEqual(parsed['IMP'], 5.0)
        self.assertEqual(parsed['CMAP'], 6.0)
        self.assertEqual(parsed['1-4 VDW'], 10.0)
        self.assertEqual(parsed['1-4 EEL'], 11.0)

    def test_gbnsr6_merge_parses_all_coordinate_sets(self):
        results = [
            'minimizing coord set #       1\n',
            ' BOND = 1.0 ANGLE = 2.0 DIHED = 3.0\n',
            ' VDWAALS = 7.0 EEL = 8.0 EGB = 9.0\n',
            ' 1-4 VDW = 10.0 1-4 EEL = 11.0 RESTRAINT = 0.0\n',
            'minimizing coord set #       2\n',
            ' BOND = 11.0 ANGLE = 12.0 DIHED = 13.0\n',
            ' VDWAALS = 17.0 EEL = 18.0 EGB = 19.0\n',
            ' 1-4 VDW = 20.0 1-4 EEL = 21.0 RESTRAINT = 0.0\n',
        ]

        parsed = MergeGBNSR6Output._get_energy_decomp(results)['energy']

        self.assertEqual(sorted(parsed), [1, 2])
        self.assertEqual(parsed[1]['VDWAALS'], 7.0)
        self.assertEqual(parsed[2]['VDWAALS'], 17.0)


if __name__ == '__main__':
    unittest.main()
