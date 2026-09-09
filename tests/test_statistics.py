import sys
import types
import unittest

import numpy as np

if 'parmed' not in sys.modules:
    parmed = types.ModuleType('parmed')
    parmed.__version__ = 'stub'
    sys.modules['parmed'] = parmed

from GMXMMPBSA.amber_outputs import DecompOut
from GMXMMPBSA.output_file import _combine_entropy
from GMXMMPBSA.utils import EnergyVector, block_statistics, calc_sum


class EnergyVectorStatisticsTest(unittest.TestCase):
    def test_frame_statistics_use_explicit_population_convention(self):
        values = EnergyVector([1.0, 2.0, 3.0, 4.0])

        self.assertAlmostEqual(float(values.mean()), 2.5)
        self.assertAlmostEqual(float(values.std()), np.std(values, ddof=0))
        self.assertAlmostEqual(float(values.sem()), float(values.std()) / np.sqrt(4))
        self.assertNotAlmostEqual(float(values.std()), float(np.asarray(values).std(ddof=1)))

    def test_zero_propagated_uncertainty_is_not_replaced_by_frame_sd(self):
        values = EnergyVector([1.0, 2.0, 3.0, 4.0], com_std=0.0)

        self.assertEqual(values.stdev(), 0.0)
        self.assertEqual(values.semp(), 0.0)

    def test_propagated_statistics_are_distinct_from_frame_statistics(self):
        values = EnergyVector([1.0, 2.0, 3.0, 4.0], com_std=2.0)

        self.assertEqual(values.stdev(), 2.0)
        self.assertEqual(values.semp(), 1.0)
        self.assertAlmostEqual(values.std(), np.std(values, ddof=0))

    def test_scalar_combinations_use_primary_block_sem(self):
        values = EnergyVector([1.0, 2.0, 3.0, 4.0])
        expected_uncertainty = np.std([1.5, 3.5], ddof=1) / np.sqrt(2)

        self.assertEqual(calc_sum(values, 2.0), (4.5, expected_uncertainty))
        self.assertEqual(_combine_entropy(values, 2.0, float('nan')), (4.5, expected_uncertainty))

    def test_block_statistics_use_sample_sd_of_nonoverlapping_block_means(self):
        values = EnergyVector([1.0, 2.0, 3.0, 4.0])

        block_size, nblocks, block_sd, block_sem = block_statistics(values)

        self.assertEqual((block_size, nblocks), (2, 2))
        self.assertAlmostEqual(block_sd, np.std([1.5, 3.5], ddof=1))
        self.assertAlmostEqual(block_sem, block_sd / np.sqrt(2))


class DecompositionStatisticsTest(unittest.TestCase):
    def test_csv_summary_prints_raw_and_propagated_statistics(self):
        decomp = DecompOut('delta')
        decomp.numframes = 4
        decomp['TDC'] = {
            'R:A:ALA:1': {
                'int': EnergyVector([0.0, 0.0, 0.0, 0.0]),
                'vdw': EnergyVector([1.0, 2.0, 3.0, 4.0], com_std=2.0),
                'eel': EnergyVector([0.0, 0.0, 0.0, 0.0]),
                'pol': EnergyVector([0.0, 0.0, 0.0, 0.0]),
                'sas': EnergyVector([0.0, 0.0, 0.0, 0.0]),
                'tot': EnergyVector([1.0, 2.0, 3.0, 4.0], com_std=3.0),
            }
        }

        summary = decomp.summary('csv')
        headers = summary[3]
        row = summary[4]

        self.assertEqual(headers[1:8], ['Avg.', 'SD(Prop.)', 'SD', 'SEM(Prop.)', 'SEM', 'Block SD', 'Block SEM'])
        self.assertEqual(row[0], 'R:A:ALA:1')
        np.testing.assert_allclose(row[8:15], [2.5, 2.0, np.std([1.0, 2.0, 3.0, 4.0]), 1.0,
                                               np.std([1.0, 2.0, 3.0, 4.0]) / 2,
                                               np.std([1.5, 3.5], ddof=1), np.std([1.5, 3.5], ddof=1) / np.sqrt(2)])
        np.testing.assert_allclose(row[-7:], [2.5, 3.0, np.std([1.0, 2.0, 3.0, 4.0]), 1.5,
                                              np.std([1.0, 2.0, 3.0, 4.0]) / 2,
                                              np.std([1.5, 3.5], ddof=1), np.std([1.5, 3.5], ddof=1) / np.sqrt(2)])

        ascii_summary = decomp.summary('ascii')
        self.assertIn('Avg +/- Block SEM', ascii_summary)
        self.assertIn('/ Block SD', ascii_summary)

if __name__ == '__main__':
    unittest.main()
