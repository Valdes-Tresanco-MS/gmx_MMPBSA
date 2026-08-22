import copy
import logging
import sys
import types
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np

from GMXMMPBSA.amber_outputs import NMODEout
from GMXMMPBSA.exceptions import InputError
from GMXMMPBSA.input_parser import DEFAULT_QM_THEORY, SUPPORTED_QM_THEORIES, input_file
from GMXMMPBSA.utils import EnergyVector


def _parse_base_input():
    with TemporaryDirectory() as tmpdir:
        infile = Path(tmpdir) / 'mmpbsa.in'
        infile.write_text('&general\n/\n&gb\n/\n')
        parser = copy.deepcopy(input_file)
        for namelist in parser.namelists.values():
            namelist.open = False
        parsed = parser.Parse(infile.as_posix())
    parsed['gb']['gbrun'] = True
    return parsed


def _import_main_with_stubs():
    if 'pandas' not in sys.modules:
        pandas = types.ModuleType('pandas')
        pandas.MultiIndex = type('MultiIndex', (), {})
        pandas.Index = type('Index', (), {})
        sys.modules['pandas'] = pandas

    if 'parmed' not in sys.modules:
        parmed = types.ModuleType('parmed')
        parmed.__version__ = 'stub'
        sys.modules['parmed'] = parmed

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
        module = types.ModuleType(module_name)
        for name in names:
            setattr(module, name, dummy)
        sys.modules[module_name] = module

    createinput = types.ModuleType('GMXMMPBSA.createinput')
    createinput.create_inputs = lambda *args, **kwargs: None
    createinput.SanderRISMInput = dummy
    sys.modules['GMXMMPBSA.createinput'] = createinput
    sys.modules['GMXMMPBSA.make_top_amber'] = types.SimpleNamespace(CheckAmberTop=dummy)
    sys.modules['GMXMMPBSA.make_top'] = types.SimpleNamespace(CheckMakeTop=dummy)
    output_file = types.ModuleType('GMXMMPBSA.output_file')
    output_file.write_outputs = lambda *args, **kwargs: None
    output_file.write_decomp_output = lambda *args, **kwargs: None
    output_file.data2pkl = lambda *args, **kwargs: None
    sys.modules['GMXMMPBSA.output_file'] = output_file
    sys.modules['GMXMMPBSA.parm_setup'] = types.SimpleNamespace(MMPBSA_System=dummy)

    sys.modules.pop('GMXMMPBSA.main', None)
    from GMXMMPBSA import main
    return main


class Phase4ValidationTest(unittest.TestCase):
    def _app(self):
        main = _import_main_with_stubs()
        app = main.MMPBSA_App.__new__(main.MMPBSA_App)
        app.master = True
        app.stability = False
        app.FILES = SimpleNamespace(
            input_file='mmpbsa.in', ligand_mol2=None, complex_top='topol.top',
            receptor_tpr=None, receptor_trajs=[], receptor_top=None,
            ligand_tpr=None, ligand_trajs=[], ligand_top=None,
        )
        app.INPUT = _parse_base_input()
        return app

    def test_pb_validation_reports_invalid_value_from_pb_section(self):
        for name, value, expected in (
            ('prbrad', 1.5, 'PRBRAD (1.5)'),
            ('inp', 3, 'INP/NPOPT (3)'),
            ('radiopt', 2, 'RADIOPT (2)'),
        ):
            with self.subTest(name=name):
                app = self._app()
                app.INPUT['pb'][name] = value
                with self.assertRaises(InputError) as exc:
                    app.check_for_bad_input()
                self.assertIn(expected, str(exc.exception))

    def test_gbnsr6_coordinate_rank_accepts_unsuffixed_and_ranked_files(self):
        main = _import_main_with_stubs()

        self.assertEqual(
            main._gbnsr6_coordinate_rank(Path('_GMXMMPBSA_complex.inpcrd')),
            0,
        )
        self.assertEqual(
            main._gbnsr6_coordinate_rank(Path('_GMXMMPBSA_complex.3.inpcrd')),
            3,
        )

    def test_linit_zero_is_valid(self):
        app = self._app()
        app.INPUT['pb']['linit'] = 0
        app.check_for_bad_input()

    def test_idecomp_zero_is_valid_when_decomposition_is_disabled(self):
        app = self._app()
        app.INPUT['decomp']['idecomp'] = 0
        app.INPUT['decomp']['decomprun'] = False
        app.check_for_bad_input()

    def test_ligand_mol2_requires_gaff_forcefield(self):
        app = self._app()
        app.FILES.ligand_mol2 = 'ligand.mol2'
        app.INPUT['general']['forcefields'] = ['oldff/leaprc.ff99SB']
        with self.assertRaises(InputError) as exc:
            app.check_for_bad_input()
        self.assertIn('ligand MOL2 file (-lm)', str(exc.exception))

    def test_qmmm_default_theory_is_pm6_dh_plus(self):
        app = self._app()
        self.assertEqual(app.INPUT['gb']['qm_theory'], DEFAULT_QM_THEORY)

    def test_qmmm_validation_accepts_all_canonical_theories(self):
        for theory in SUPPORTED_QM_THEORIES:
            with self.subTest(theory=theory):
                app = self._app()
                app.INPUT['gb']['ifqnt'] = 1
                app.INPUT['gb']['qm_residues'] = ':1'
                app.INPUT['gb']['qm_theory'] = theory
                app.check_for_bad_input()

    def test_qmmm_validation_rejects_undocumented_legacy_aliases(self):
        for theory in ('PDDG-PM3', 'PM3PDDG', 'PDDG-MNDO', 'PDDGMNDO', 'SCC-DFTB',
                       'PM3-ZnB', 'PM3ZNB', 'MNDOD'):
            with self.subTest(theory=theory):
                app = self._app()
                app.INPUT['gb']['ifqnt'] = 1
                app.INPUT['gb']['qm_residues'] = ':1'
                app.INPUT['gb']['qm_theory'] = theory
                with self.assertRaises(InputError):
                    app.check_for_bad_input()

    def test_qmmm_validation_rejects_ndiis_attempts_outside_sander_range(self):
        for value in (-1, 1001):
            with self.subTest(value=value):
                app = self._app()
                app.INPUT['gb']['ifqnt'] = 1
                app.INPUT['gb']['qm_residues'] = ':1'
                app.INPUT['gb']['ndiis_attempts'] = value
                with self.assertRaises(InputError) as exc:
                    app.check_for_bad_input()
                self.assertIn('NDIIS_ATTEMPTS must be between 0 and 1000', str(exc.exception))

    def test_qmmm_validation_rejects_ndiis_attempts_for_dftb(self):
        app = self._app()
        app.INPUT['gb']['ifqnt'] = 1
        app.INPUT['gb']['qm_residues'] = ':1'
        app.INPUT['gb']['qm_theory'] = 'DFTB'
        app.INPUT['gb']['ndiis_attempts'] = 700
        with self.assertRaises(InputError) as exc:
            app.check_for_bad_input()
        self.assertIn('NDIIS_ATTEMPTS is not available for DFTB', str(exc.exception))

    def test_analyzer_startup_failure_is_reported_as_nonfatal(self):
        main = _import_main_with_stubs()
        with self.assertLogs(level='WARNING') as logs:
            main.MMPBSA_App._report_analyzer_startup_failure('results.info')
        self.assertIn('Calculation completed', logs.output[0])
        self.assertIn('gmx_MMPBSA_ana -f results.info', logs.output[0])

    def test_missing_analyzer_executable_is_handled_as_startup_failure(self):
        main = _import_main_with_stubs()
        app = object.__new__(main.MMPBSA_App)
        app.master = True
        app.FILES = SimpleNamespace(
            rewrite_output=False, gui=True, prefix='results-',
        )
        app.INPUT = {'general': {'keep_files': 0}}

        class FakeTimer:
            def done(self):
                pass

            def print_(self, *args, **kwargs):
                pass

        class FakeMPI:
            def Finalize(self):
                pass

        app.timer = FakeTimer()
        app.MPI = FakeMPI()
        with patch.object(app, 'remove'), patch.object(app, '_finalize_timers'), patch.object(main.utils, 'get_warnings', return_value={
            'error': 0, 'warning': 0,
        }), patch('subprocess.Popen', side_effect=FileNotFoundError('missing analyzer')):
            with self.assertRaises(SystemExit) as exc:
                app.finalize()

        self.assertEqual(exc.exception.code, 0)


class Phase4ScientificWarningTest(unittest.TestCase):
    @staticmethod
    def _nmode(values):
        output = NMODEout.__new__(NMODEout)
        dict.__init__(output)
        output.mol = 'complex'
        output.numframes = len(values)
        output.data_keys = ['TRANSLATIONAL', 'ROTATIONAL', 'VIBRATIONAL']
        output.no_nmode_convergence = False
        for key in output.data_keys + ['TOTAL']:
            output[key] = EnergyVector(values)
        return output

    def test_nmode_warning_counts_frames_and_describes_substitution(self):
        output = self._nmode([1.0, np.nan, 3.0])
        with self.assertLogs(level='WARNING') as logs:
            output._fill_nmode_values()
        self.assertIn('1 of 3 NMODE frames', logs.output[0])
        self.assertIn('replaced with the mean', logs.output[0])
        self.assertAlmostEqual(output['TRANSLATIONAL'][1], 2.0)

    def test_nmode_all_unconverged_frames_are_not_substituted(self):
        output = self._nmode([np.nan, np.nan])
        with self.assertLogs(level='WARNING') as logs:
            output._fill_nmode_values()
        self.assertTrue(output.no_nmode_convergence)
        self.assertIn('2 of 2 NMODE frames', logs.output[0])
        self.assertIn('no substitution was performed', logs.output[0])


if __name__ == '__main__':
    unittest.main()
