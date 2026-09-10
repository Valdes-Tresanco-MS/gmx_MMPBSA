import unittest
from types import SimpleNamespace

import pandas as pd

from GMXMMPBSA.API import MMPBSA_API


class StabilityAPITest(unittest.TestCase):
    def test_qh_entropy_summary_uses_complex_for_stability(self):
        api = MMPBSA_API()
        api.app_namespace = SimpleNamespace(FILES=SimpleNamespace(stability=True))
        api.data = {
            'normal': {
                'qh': {
                    'complex': {'TOTAL': pd.Series([1.0, 2.0])},
                },
            },
        }

        result = api.get_qh_entropy()

        self.assertEqual(result['summary'].columns.tolist(), [('normal', 'qh')])
        self.assertEqual(result['summary'].iloc[:, 0].tolist(), [1.0, 2.0])


class BindingSummaryCompatibilityTest(unittest.TestCase):
    def test_legacy_extended_and_mixed_summaries(self):
        import math
        for stability in (False, True):
            for energy_extended, entropy_extended in ((False, False), (True, True), (False, True), (True, False)):
                with self.subTest(stability=stability, energy=energy_extended, entropy=entropy_extended):
                    api = MMPBSA_API()
                    api.app_namespace = SimpleNamespace(FILES=SimpleNamespace(stability=stability))
                    mol = 'complex' if stability else 'delta'
                    energy = pd.Series({'Average': -10., 'SD': 4., 'SEM': 1.})
                    entropy = pd.Series({'Average': 3., 'SD': 3., 'SEM': .5})
                    if energy_extended:
                        energy['Block SD'], energy['Block SEM'] = 2., .8
                    if entropy_extended:
                        entropy['Block SD'], entropy['Block SEM'] = 1., .6
                    result = api.get_binding({'normal': {'gb': {mol: {'TOTAL': energy}}}},
                                             {'normal': {'nmode': {mol: {'TOTAL': entropy}}}})
                    total = result['data']['normal']['gb']['nmode']['ΔG']
                    self.assertEqual(total['Average'], -7.)
                    self.assertEqual(total['SD'], 5.)
                    self.assertAlmostEqual(total['SEM'], math.sqrt(1.25))
                    if energy_extended and entropy_extended:
                        self.assertAlmostEqual(total['Block SD'], math.sqrt(5.))
                        self.assertAlmostEqual(total['Block SEM'], 1.)
                    else:
                        self.assertTrue(math.isnan(total['Block SD']))
                        self.assertTrue(math.isnan(total['Block SEM']))

    def test_ie_and_c2_dataframe_summaries_preserve_block_uncertainties(self):
        import math

        energy = pd.Series({'Average': -10., 'SD': 4., 'SEM': 1., 'Block SD': 2., 'Block SEM': .8})
        for entropy_name in ('ie', 'c2'):
            with self.subTest(entropy=entropy_name):
                api = MMPBSA_API()
                api.app_namespace = SimpleNamespace(FILES=SimpleNamespace(stability=False))
                entropy = pd.DataFrame(
                    {entropy_name: [3., 3., .5, 1., .6], 'sigma': [0., 0., 0., 0., 0.]},
                    index=['Average', 'SD', 'SEM', 'Block SD', 'Block SEM'],
                )

                result = api.get_binding(
                    {'normal': {'gb': {'delta': {'TOTAL': energy}}}},
                    {'normal': {entropy_name: {'gb': {entropy_name: entropy}}}},
                )
                total = result['data']['normal']['gb'][entropy_name]['ΔG']

                self.assertEqual(total['Average'], -7.)
                self.assertAlmostEqual(total['Block SD'], math.sqrt(5.))
                self.assertAlmostEqual(total['Block SEM'], 1.)


class ReferenceStatisticsTest(unittest.TestCase):
    def test_reference_adjusts_each_uncertainty_without_mutating_source(self):
        from GMXMMPBSA.utils import adjust_reference_statistics
        original = pd.DataFrame({('GB', 'Average'): [-10., -7.], ('GB', 'SD'): [8., 12.],
                                 ('GB', 'SEM'): [2., 3.], ('GB', 'Block SD'): [2., 6.],
                                 ('GB', 'Block SEM'): [.5, 1.5]})
        result = adjust_reference_statistics(original, original.iloc[0])
        self.assertEqual(result.iloc[0].tolist(), [0.] * 5)
        self.assertEqual(result.iloc[1].tolist(), [3., 4., 1., 4., 1.])
        self.assertEqual(original[('GB', 'Block SEM')].tolist(), [.5, 1.5])

    def test_missing_block_uncertainties_remain_unavailable(self):
        from GMXMMPBSA.utils import adjust_reference_statistics
        original = pd.DataFrame({('GB', 'Block SEM'): [float('nan'), 1.5]})
        result = adjust_reference_statistics(original, original.iloc[0])
        self.assertTrue(result[('GB', 'Block SEM')].isna().all())


if __name__ == '__main__':
    unittest.main()
