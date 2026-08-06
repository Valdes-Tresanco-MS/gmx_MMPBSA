import unittest
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

from GMXMMPBSA import output_file


class _Value:
    def mean(self):
        return 1.0

    def std(self):
        return 0.1


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


if __name__ == '__main__':
    unittest.main()
