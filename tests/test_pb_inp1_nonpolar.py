import copy
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest import mock

from GMXMMPBSA.utils import PB_INP1_DEFAULTS, PB_INP1_NONPOLAR, sync_pb_nonpolar_for_inp

try:
    from GMXMMPBSA.createinput import SanderPBSAInput
except ModuleNotFoundError as exc:
    if exc.name and exc.name.startswith('parmed'):
        SanderPBSAInput = None
    else:
        raise


def _pb_input(inp=1, **overrides):
    pb = {
        'pbrun': True,
        'inp': inp,
        'sprob': 0.557,
        'cavity_surften': 0.0378,
        'cavity_offset': -0.5692,
        'radiopt': 1,
    }
    pb.update(overrides)
    return {'pb': pb}


class SyncPbNonpolarForInpTest(unittest.TestCase):
    def test_disabled_pb_preserves_parameters_without_warning(self):
        INPUT = _pb_input(pbrun=False)
        original = copy.deepcopy(INPUT)

        with self.assertNoLogs(level='WARNING'):
            changed = sync_pb_nonpolar_for_inp(INPUT)

        self.assertFalse(changed)
        self.assertEqual(INPUT, original)

    def test_missing_pb_run_flag_preserves_parameters_without_warning(self):
        INPUT = _pb_input()
        del INPUT['pb']['pbrun']
        original = copy.deepcopy(INPUT)

        with self.assertNoLogs(level='WARNING'):
            changed = sync_pb_nonpolar_for_inp(INPUT)

        self.assertFalse(changed)
        self.assertEqual(INPUT, original)

    def test_inp1_resets_nonpolar_and_radiopt_to_amber_guidance(self):
        INPUT = _pb_input(inp=1)

        with self.assertLogs(level='WARNING') as logged:
            changed = sync_pb_nonpolar_for_inp(INPUT)

        self.assertTrue(changed)
        self.assertEqual(INPUT['pb']['sprob'], PB_INP1_DEFAULTS['sprob'])
        self.assertEqual(INPUT['pb']['cavity_surften'], PB_INP1_DEFAULTS['cavity_surften'])
        self.assertEqual(INPUT['pb']['cavity_offset'], PB_INP1_DEFAULTS['cavity_offset'])
        self.assertEqual(INPUT['pb']['radiopt'], 0)
        self.assertTrue(any('cavity_surften: 0.0378 -> 0.005' in msg for msg in logged.output))
        self.assertTrue(any('radiopt: 1 -> 0' in msg for msg in logged.output))

    def test_inp1_already_aligned_is_noop(self):
        INPUT = _pb_input(inp=1, **PB_INP1_NONPOLAR)

        with mock.patch('GMXMMPBSA.utils.logging.warning') as warning:
            changed = sync_pb_nonpolar_for_inp(INPUT)

        self.assertFalse(changed)
        warning.assert_not_called()
        self.assertEqual(INPUT['pb']['cavity_surften'], 0.005)
        self.assertEqual(INPUT['pb']['radiopt'], 0)

    def test_inp2_preserves_cavity_and_radiopt(self):
        INPUT = _pb_input(inp=2, radiopt=1)

        with mock.patch('GMXMMPBSA.utils.logging.warning') as warning:
            changed = sync_pb_nonpolar_for_inp(INPUT)

        self.assertFalse(changed)
        warning.assert_not_called()
        self.assertEqual(INPUT['pb']['cavity_surften'], 0.0378)
        self.assertEqual(INPUT['pb']['cavity_offset'], -0.5692)
        self.assertEqual(INPUT['pb']['sprob'], 0.557)
        self.assertEqual(INPUT['pb']['radiopt'], 1)

    def test_process_input_applies_sync(self):
        from GMXMMPBSA.main import MMPBSA_App
        from GMXMMPBSA.membrane import AUTOMATIC

        app = MMPBSA_App.__new__(MMPBSA_App)
        app.INPUT = {
            'pb': {
                'pbrun': True,
                'inp': 1,
                'scale': 2.0,
                'memopt': 0,
                'mthick': AUTOMATIC,
                'mctrdz': AUTOMATIC,
                'sprob': 0.557,
                'cavity_surften': 0.0378,
                'cavity_offset': -0.5692,
                'radiopt': 1,
            },
            'general': {'netcdf': 0},
            'decomp': {'decomprun': False},
            'rism': {
                'solvcut': None,
                'buffer': 14.0,
                'rismrun': False,
                'gfcorrection': 0,
                'pcpluscorrection': 0,
            },
        }

        with self.assertLogs(level='WARNING'):
            app.process_input()

        self.assertEqual(app.INPUT['pb']['cavity_surften'], 0.005)
        self.assertEqual(app.INPUT['pb']['cavity_offset'], 0.0)
        self.assertEqual(app.INPUT['pb']['sprob'], 1.4)
        self.assertEqual(app.INPUT['pb']['radiopt'], 0)
        self.assertEqual(app.INPUT['pb']['scale'], 0.5)

    @staticmethod
    def _parsed_app(*namelists):
        from GMXMMPBSA.input_parser import input_file
        from GMXMMPBSA.main import MMPBSA_App

        with TemporaryDirectory() as directory:
            path = Path(directory) / 'mmpbsa.in'
            path.write_text('&general\n/\n' + ''.join(f'&{name}\n/\n' for name in namelists))
            inputs = copy.deepcopy(input_file).Parse(str(path))
        app = MMPBSA_App.__new__(MMPBSA_App)
        app.INPUT = inputs
        return app

    def test_process_input_non_pb_runs_do_not_adjust_pb_or_warn(self):
        for namelist in ('gb', 'gbnsr6', 'rism'):
            with self.subTest(namelist=namelist):
                app = self._parsed_app(namelist)
                self.assertFalse(app.INPUT['pb']['pbrun'])
                self.assertEqual(app.INPUT['pb']['inp'], 1)
                original = {key: app.INPUT['pb'][key] for key in PB_INP1_DEFAULTS}

                with self.assertNoLogs(level='WARNING'):
                    app.process_input()

                self.assertEqual({key: app.INPUT['pb'][key] for key in PB_INP1_DEFAULTS}, original)

    def test_process_input_combined_gb_pb_still_adjusts_pb_and_warns(self):
        app = self._parsed_app('gb', 'pb')
        self.assertTrue(app.INPUT['gb']['gbrun'])
        self.assertTrue(app.INPUT['pb']['pbrun'])

        with self.assertLogs(level='WARNING') as logged:
            app.process_input()

        self.assertEqual({key: app.INPUT['pb'][key] for key in PB_INP1_DEFAULTS}, PB_INP1_DEFAULTS)
        self.assertEqual(len(logged.output), 1)
        self.assertIn('PB inp=1', logged.output[0])


@unittest.skipIf(SanderPBSAInput is None, 'ParmEd is required to generate PB input files')
class PBMdinFallbackTest(unittest.TestCase):
    def test_pb_mdin_fallback_uses_inp1_aligned_values(self):
        pb = SanderPBSAInput({
            'general': {'netcdf': False},
            'pb': {'istrng': 0.0, 'inp': 1},
            'decomp': {'idecomp': 0, 'dec_verbose': 0},
        })
        pb.make_mdin()

        # Incomplete INPUT should not resurrect inp=2 / Tan-Luo defaults when inp=1.
        self.assertEqual(pb.input_items['cavity_surften'], 0.005)
        self.assertEqual(pb.input_items['cavity_offset'], 0.0)
        self.assertEqual(pb.input_items['sprob'], 1.4)
        self.assertEqual(pb.input_items['radiopt'], 0)


class PBRadioptParserDefaultTest(unittest.TestCase):
    def test_parser_default_radiopt_is_zero(self):
        from GMXMMPBSA.input_parser import input_file

        self.assertEqual(input_file.namelists['pb'].variables['radiopt'].value, 0)


if __name__ == '__main__':
    unittest.main()
