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

from GMXMMPBSA.calculation import ListEnergyCalculation


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


if __name__ == '__main__':
    unittest.main()
