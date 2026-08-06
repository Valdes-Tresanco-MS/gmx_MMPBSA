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


if __name__ == '__main__':
    unittest.main()
