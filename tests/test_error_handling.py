import io
import logging
import sys
import types
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

from GMXMMPBSA.exceptions import GMXMMPBSA_ERROR, InputError, MMPBSA_Error
from GMXMMPBSA.infofile import _determine_type
from GMXMMPBSA.input_parser import input_file


def _install_import_stubs():
    if 'pandas' not in sys.modules:
        pandas = types.ModuleType('pandas')
        pandas.MultiIndex = type('MultiIndex', (), {})
        pandas.Index = type('Index', (), {})
        sys.modules['pandas'] = pandas

    if 'parmed' not in sys.modules:
        parmed = types.ModuleType('parmed')
        parmed.__version__ = 'stub'
        sys.modules['parmed'] = parmed


def _stub_module(name, **attrs):
    module = types.ModuleType(name)
    for attr, value in attrs.items():
        setattr(module, attr, value)
    sys.modules[name] = module
    return module


def _import_main_with_stubs():
    _install_import_stubs()
    dummy = type('Dummy', (), {})
    for module_name, names in {
        'GMXMMPBSA.amber_outputs': [
            'QHout', 'NMODEout', 'QMMMout', 'GBout', 'PBout', 'PolarRISM_std_Out',
            'RISM_std_Out', 'PolarRISM_gf_Out', 'RISM_gf_Out', 'PolarRISM_pcplus_Out',
            'RISM_pcplus_Out', 'BindingStatistics', 'IEout', 'C2out',
            'DeltaDeltaStatistics', 'DeltaIEC2Statistic', 'DeltaDeltaQH',
            'GBNSR6out', 'MMout'
        ],
        'GMXMMPBSA.calculation': [
            'CalculationList', 'EnergyCalculation', 'PBEnergyCalculation', 'NmodeCalc',
            'QuasiHarmCalc', 'CopyCalc', 'PrintCalc', 'LcpoCalc', 'MolsurfCalc',
            'InteractionEntropyCalc', 'C2EntropyCalc', 'MergeOut', 'ListEnergyCalculation'
        ],
    }.items():
        _stub_module(module_name, **{name: dummy for name in names})

    _stub_module('GMXMMPBSA.createinput', create_inputs=lambda *args, **kwargs: None,
                 SanderRISMInput=dummy)
    _stub_module('GMXMMPBSA.make_top_amber', CheckAmberTop=dummy)
    _stub_module('GMXMMPBSA.output_file', write_outputs=lambda *args, **kwargs: None,
                 write_decomp_output=lambda *args, **kwargs: None,
                 data2pkl=lambda *args, **kwargs: None)
    _stub_module('GMXMMPBSA.parm_setup', MMPBSA_System=dummy)
    _stub_module('GMXMMPBSA.make_top', CheckMakeTop=dummy)

    sys.modules.pop('GMXMMPBSA.main', None)
    from GMXMMPBSA import main
    return main


class ErrorHandlingTest(unittest.TestCase):
    def test_malformed_namelist_field_raises_input_error(self):
        with TemporaryDirectory() as tmpdir:
            infile = Path(tmpdir) / 'mmpbsa.in'
            infile.write_text('&general\norphan_value,\n/\n')

            with self.assertRaises(InputError):
                input_file.Parse(infile.as_posix())

    def test_gmxmmpbsa_error_handles_non_exception_class_argument(self):
        logging.disable(logging.CRITICAL)
        try:
            with self.assertRaises(MMPBSA_Error) as exc:
                GMXMMPBSA_ERROR('', 'original failure')
        finally:
            logging.disable(logging.NOTSET)

        self.assertIn('original failure', str(exc.exception))

    def test_infofile_list_parsing_does_not_execute_code(self):
        with TemporaryDirectory() as tmpdir:
            marker = Path(tmpdir) / 'executed'
            payload = f"[__import__('os').system('touch {marker}')]"

            with self.assertWarns(Warning):
                parsed = _determine_type(payload)

            self.assertEqual(parsed, payload)
            self.assertFalse(marker.exists())

    def test_excepthook_does_not_traceback_project_errors(self):
        main = _import_main_with_stubs()

        class FakeComm:
            def Abort(self, status):
                raise SystemExit(status)

        class FakeMPI:
            COMM_WORLD = FakeComm()

        stderr = io.StringIO()
        old_mpi, old_stderr = main._MPI, main._stderr
        old_mpi_size, old_rank = main._mpi_size, main._rank
        try:
            main._MPI = FakeMPI()
            main._stderr = stderr
            main._mpi_size = 1
            main._rank = 0
            with patch('traceback.print_tb') as print_tb:
                with self.assertRaises(SystemExit):
                    main.excepthook(InputError, InputError('bad input'), None)
        finally:
            main._MPI = old_mpi
            main._stderr = old_stderr
            main._mpi_size = old_mpi_size
            main._rank = old_rank

        print_tb.assert_not_called()
        self.assertIn('InputError: bad input', stderr.getvalue())

    def test_selector_missing_within_distance_reports_project_error(self):
        _install_import_stubs()
        from GMXMMPBSA.utils import selector

        logging.disable(logging.CRITICAL)
        try:
            with self.assertRaises(MMPBSA_Error) as exc:
                selector('within')
        finally:
            logging.disable(logging.NOTSET)

        self.assertIn('Invalid distance value', str(exc.exception))

    def test_stability_skips_interaction_and_c2_entropy(self):
        main = _import_main_with_stubs()
        app = object.__new__(main.MMPBSA_App)
        app.stability = True
        app.INPUT = {'general': {'interaction_entropy': 1, 'c2_entropy': 1}}

        with self.assertLogs(level='WARNING') as logs:
            app.get_iec2entropy(from_calc=True)

        self.assertEqual(app.INPUT['general']['interaction_entropy'], 0)
        self.assertEqual(app.INPUT['general']['c2_entropy'], 0)
        self.assertIn('doesn\'t support stability calculations', '\n'.join(logs.output))

    def test_qh_entropy_is_rejected_for_new_calculations(self):
        main = _import_main_with_stubs()
        app = object.__new__(main.MMPBSA_App)
        app.master = True
        app.INPUT = {'general': {'qh_entropy': 1}}
        app.FILES = types.SimpleNamespace(input_file='mmpbsa.in')

        with self.assertRaises(InputError) as exc:
            app.check_for_bad_input()

        self.assertIn('qh_entropy=1', str(exc.exception))
        self.assertIn('not supported', str(exc.exception))


if __name__ == '__main__':
    unittest.main()
