import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import numpy as np

from GMXMMPBSA import output_file


class _Value:
    def __array__(self, dtype=None):
        return np.asarray([1.0], dtype=dtype)

    def mean(self):
        return 1.0

    def std(self):
        return 0.1

    def stdev(self):
        return 0.2


class _Stats:
    def summary_output(self):
        return 'summary'

    def mean(self):
        return 1.0

    def std(self):
        return 0.1

    def __getitem__(self, key):
        if key == 'TOTAL':
            return _Value()
        raise KeyError(key)


class StabilityOutputTest(unittest.TestCase):
    def test_output_separator_covers_extended_statistics_columns(self):
        with TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / 'output.dat'
            output = output_file.OutputFile(output_path, 'w')
            output.separate()
            output._handle.close()

            self.assertEqual(output_path.read_text().splitlines(), ['-' * 101, '-' * 101])

    def test_stability_alanine_scan_uses_complex_mutant_delta(self):
        normal = _Stats()
        mutant = _Stats()
        mutant_normal = {'complex': _Stats()}
        app = SimpleNamespace(
            FILES=SimpleNamespace(output_file='ignored.dat', energyout=None),
            normal_system=None,
            numframes=1,
            INPUT={
                'general': {
                    'sys_name': 'stability',
                    'interaction_entropy': 0,
                    'c2_entropy': 0,
                    'qh_entropy': 0,
                    'temperature': 298.15,
                },
                'ala': {'alarun': True, 'mutant_only': False},
                'nmode': {'nmoderun': False},
                'gb': {'gbrun': True, 'molsurf': False, 'ifqnt': 0},
                'pb': {'pbrun': False},
            },
            calc_types=SimpleNamespace(
                normal={'gb': {'complex': normal}},
                mutant={'gb': {'complex': mutant}},
                mut_norm={'gb': mutant_normal},
            ),
            mut_str='ALA',
            stability=True,
        )

        with patch.object(output_file, 'OutputFile', return_value=MagicMock()):
            output_file.write_outputs(app)

    def test_stability_alanine_scan_uses_complex_nmode_delta(self):
        app = SimpleNamespace(
            FILES=SimpleNamespace(output_file='ignored.dat', energyout=None),
            normal_system=None,
            numframes=1,
            numframes_nmode=1,
            INPUT={
                'general': {
                    'sys_name': 'stability',
                    'interaction_entropy': 0,
                    'c2_entropy': 0,
                    'qh_entropy': 0,
                    'temperature': 298.15,
                },
                'ala': {'alarun': True, 'mutant_only': False},
                'nmode': {'nmoderun': True},
                'gb': {'gbrun': False, 'molsurf': False, 'ifqnt': 0},
                'pb': {'pbrun': False},
            },
            calc_types=SimpleNamespace(
                normal={'nmode': {'complex': _Stats()}},
                mutant={'nmode': {'complex': _Stats()}},
                mut_norm={'nmode': {'complex': _Stats()}},
            ),
            mut_str='ALA',
            stability=True,
        )

        with patch.object(output_file, 'OutputFile', return_value=MagicMock()):
            output_file.write_outputs(app)


class InteractionEntropyCompatibilityTest(unittest.TestCase):
    def test_legacy_compact_result_uses_primary_fallback_uncertainty(self):
        legacy = {
            'data': np.asarray([0.0, 0.1, 0.3]),
            'iedata': np.asarray([0.1, 0.3]),
        }

        value, uncertainty = output_file._ie_result(legacy)

        self.assertAlmostEqual(value, 0.2)
        self.assertAlmostEqual(uncertainty, np.std([0.1, 0.3], ddof=0) / np.sqrt(2))

    def test_new_compact_result_uses_explicit_primary_and_block_values(self):
        current = {
            'data': np.asarray([0.0, 0.1, 0.3]),
            'iedata': np.asarray([0.1, 0.3]),
            'ie_value': 0.25,
            'block_std': 0.05,
            'block_sem': 0.025,
        }

        value, uncertainty = output_file._ie_result(current)

        self.assertAlmostEqual(value, 0.25)
        self.assertAlmostEqual(uncertainty, 0.025)


if __name__ == '__main__':
    unittest.main()
