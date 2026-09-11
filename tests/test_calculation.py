import io
import math
import sys
import tempfile
import types
import unittest
from unittest.mock import patch

import numpy as np

if 'pandas' not in sys.modules:
    pandas = types.ModuleType('pandas')
    pandas.MultiIndex = type('MultiIndex', (), {})
    pandas.Index = type('Index', (), {})
    sys.modules['pandas'] = pandas
if 'parmed' not in sys.modules:
    parmed = types.ModuleType('parmed')
    parmed.__version__ = 'stub'
    sys.modules['parmed'] = parmed

from GMXMMPBSA.amber_outputs import C2out, IEout
from GMXMMPBSA.calculation import (
    C2EntropyCalc, InteractionEntropyCalc, ListEnergyCalculation, MergeGBNSR6Output,
    MolsurfCalc,
)


class MolsurfCalcTest(unittest.TestCase):
    def test_cpptraj_input_terminates_molsurf_action(self):
        calculation = MolsurfCalc(
            'cpptraj', 'COM.prmtop', 'complex.mdcrd.%d',
            'complex_surf.dat.%d', probe=1.6, offset=0.2,
        )
        script = calculation._get_instring(0)
        self.assertTrue(script.endswith('\n'))
        self.assertIn('molsurf :* out complex_surf.dat.0 probe 1.6 offset 0.2\n', script)


class InteractionEntropyCalcTest(unittest.TestCase):
    @staticmethod
    def _input(temperature=300.0, ie_segment=100):
        return {
            'general': {
                'temperature': temperature,
                'ie_segment': ie_segment,
                'interaction_entropy': 1,
                'startframe': 1,
                'interval': 1,
            }
        }

    @staticmethod
    def _paper_curve(energies, temperature=300.0):
        energies = np.asarray(energies, dtype=float)
        kT = InteractionEntropyCalc.GAS_CONSTANT * temperature
        result = []
        for nframes in range(1, energies.size + 1):
            prefix = energies[:nframes]
            centered = (prefix - prefix.mean()) / kT
            max_centered = centered.max()
            result.append(kT * (
                max_centered + math.log(np.exp(centered - max_centered).mean())
            ))
        return np.asarray(result)

    def _calculate(self, energies, temperature=300.0):
        with patch('GMXMMPBSA.calculation.tqdm', side_effect=lambda values, **kwargs: values):
            return InteractionEntropyCalc(
                np.asarray(energies, dtype=float), self._input(temperature), 'gb'
            )

    def test_matches_paper_equations_for_every_prefix(self):
        energies = np.array([-2.0, -1.0, 0.0, 1.0, 2.0])

        calculated = self._calculate(energies).data
        expected = self._paper_curve(energies)

        np.testing.assert_allclose(calculated, expected, rtol=0.0, atol=1e-12)

    def test_final_value_is_order_invariant_and_nonnegative(self):
        energies = np.array([0.0, 1.0, 2.0])

        forward = self._calculate(energies).data
        reverse = self._calculate(energies[::-1]).data

        self.assertTrue(np.all(forward >= 0.0))
        self.assertTrue(np.all(reverse >= 0.0))
        self.assertAlmostEqual(forward[-1], reverse[-1], places=12)

    def test_is_invariant_to_constant_energy_offset(self):
        energies = np.array([-3.2, 0.5, 1.7, 4.1])

        original = self._calculate(energies).data
        shifted = self._calculate(energies + 10000.0).data

        np.testing.assert_allclose(original, shifted, rtol=0.0, atol=1e-10)

    def test_large_energy_fluctuation_remains_finite(self):
        calculated = self._calculate([0.0, 1000.0]).data

        self.assertTrue(np.all(np.isfinite(calculated)))
        np.testing.assert_allclose(
            calculated, self._paper_curve([0.0, 1000.0]), rtol=0.0, atol=1e-12
        )

    def test_full_ensemble_is_primary_and_tail_is_only_diagnostic(self):
        energies = np.linspace(-2.0, 3.0, 64) ** 2

        calculated = self._calculate(energies)

        self.assertEqual(calculated.ie_value, calculated.data[-1])
        self.assertAlmostEqual(calculated.tail_mean, calculated.iedata.mean())
        self.assertEqual(calculated.block_nblocks, 8)
        self.assertGreater(calculated.block_size, 1)
        self.assertTrue(np.isfinite(calculated.block_std))

    def test_ie_segment_zero_keeps_empty_tail_diagnostic(self):
        energies = np.array([-2.0, -1.0, 0.0, 1.0, 2.0])
        with patch('GMXMMPBSA.calculation.tqdm', side_effect=lambda values, **kwargs: values):
            calculated = InteractionEntropyCalc(
                energies, self._input(ie_segment=25), 'gb', iesegment=0
            )

        self.assertEqual(calculated.isegment, 0)
        self.assertEqual(calculated.ieframes, 0)
        self.assertEqual(calculated.iedata.size, 0)
        self.assertTrue(np.isnan(calculated.tail_mean))
        self.assertAlmostEqual(calculated.ie_value, calculated.data[-1])
        # Falsy ``or`` would have fallen back to INPUT ie_segment=25 (2 frames).
        self.assertNotEqual(calculated.ieframes, math.ceil(energies.size * 0.25))

    def test_ie_output_round_trip_preserves_block_diagnostics(self):
        energies = np.sin(np.arange(64, dtype=float) / 5.0)
        calculated = self._calculate(energies)

        with tempfile.NamedTemporaryFile(mode='w+', suffix='.dat') as output:
            calculated.save_output(output.name)
            parsed = IEout(self._input(), 'gb')
            parsed.parse_from_file(output.name, energies.size)

        self.assertAlmostEqual(parsed['ie_value'], calculated.ie_value, places=4)
        self.assertAlmostEqual(parsed['block_std'], calculated.block_std, places=4)
        self.assertEqual(parsed['block_size'], calculated.block_size)
        self.assertEqual(len(parsed.block_analysis), len(calculated.block_analysis))
        self.assertAlmostEqual(parsed.summary()[1][3], calculated.ie_value, places=4)

    def test_legacy_ie_output_remains_readable(self):
        # 1.6.x reported the tail mean as "| Interaction Entropy (-TΔS):".
        # The last cumulative curve point must not become the primary value.
        legacy_output = """| Interaction Entropy results for gb calculations
IE-frames: last 2
Internal Energy SD (sigma):      1.00
| Interaction Entropy (-TΔS):      0.20 +/-    0.14

| Interaction Entropy per-frame:
Frame # | IE value
1  0.00
2  0.10
3  0.30
"""

        with tempfile.NamedTemporaryFile(mode='w+', suffix='.dat') as output:
            output.write(legacy_output)
            output.flush()
            parsed = IEout(self._input(), 'gb')
            parsed.parse_from_file(output.name, 3)

        self.assertAlmostEqual(parsed['ie_value'], 0.20)
        self.assertAlmostEqual(parsed['tail_mean'], 0.20)
        self.assertEqual(parsed['block_size'], 0)
        self.assertAlmostEqual(parsed['block_std'], parsed['tail_std'])

    def test_legacy_ie_without_header_uses_tail_mean(self):
        legacy_output = """| Interaction Entropy results for gb calculations
IE-frames: last 2
Internal Energy SD (sigma):      1.00

| Interaction Entropy per-frame:
Frame # | IE value
1  0.00
2  0.10
3  0.30
"""

        with tempfile.NamedTemporaryFile(mode='w+', suffix='.dat') as output:
            output.write(legacy_output)
            output.flush()
            parsed = IEout(self._input(), 'gb')
            parsed.parse_from_file(output.name, 3)

        self.assertAlmostEqual(parsed['ie_value'], 0.20)
        self.assertNotAlmostEqual(parsed['ie_value'], float(parsed['data'][-1]))


class C2EntropyCalcTest(unittest.TestCase):
    @staticmethod
    def _input(temperature=300.0):
        return {'general': {'temperature': temperature}}

    def test_matches_second_order_cumulant_and_uses_blocks(self):
        energies = np.sin(np.arange(64, dtype=float) / 5.0)
        calculated = C2EntropyCalc(energies, self._input(), 'gb')
        kT = InteractionEntropyCalc.GAS_CONSTANT * 300.0

        self.assertAlmostEqual(calculated.c2data, energies.std() ** 2 / (2 * kT), places=12)
        self.assertEqual(calculated.block_nblocks, 8)
        self.assertTrue(np.isfinite(calculated.c2_std))
        self.assertTrue(np.isfinite(calculated.c2_sem))

    def test_c2_output_round_trip_preserves_block_diagnostics(self):
        energies = np.cos(np.arange(64, dtype=float) / 7.0)
        calculated = C2EntropyCalc(energies, self._input(), 'gb')

        with tempfile.NamedTemporaryFile(mode='w+', suffix='.dat') as output:
            calculated.save_output(output.name)
            parsed = C2out('gb')
            parsed.parse_from_file(output.name)

        self.assertAlmostEqual(parsed['c2data'], calculated.c2data, places=4)
        self.assertAlmostEqual(parsed['c2_std'], calculated.c2_std, places=4)
        self.assertAlmostEqual(parsed['c2_sem'], calculated.c2_sem, places=4)
        self.assertEqual(parsed['block_size'], calculated.block_size)
        self.assertEqual(len(parsed.block_analysis), len(calculated.block_analysis))

    def test_legacy_c2_output_remains_readable(self):
        legacy_output = """| C2 Entropy results for gb calculations
C2 Entropy (-TΔS): 1.2500
C2 Entropy SD: 0.4000
Internal Energy SD (sigma):      2.00
C2 Entropy CI: 0.5000 2.0000
"""

        with tempfile.NamedTemporaryFile(mode='w+', suffix='.dat') as output:
            output.write(legacy_output)
            output.flush()
            parsed = C2out('gb')
            parsed.parse_from_file(output.name)

        self.assertAlmostEqual(parsed['c2data'], 1.25)
        self.assertAlmostEqual(parsed['c2_sem'], 0.40)
        self.assertEqual(parsed['block_size'], 0)
        self.assertEqual(parsed['c2_ci'], [0.5, 2.0])


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
