import copy
import io
import logging
import sys
import types
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
from unittest.mock import patch, mock_open

from GMXMMPBSA.exceptions import InputError, MMPBSA_Error
from GMXMMPBSA.input_parser import input_file


def _base_input():
    with TemporaryDirectory() as tmpdir:
        infile = Path(tmpdir) / 'mmpbsa.in'
        infile.write_text('&general\n/\n&gb\n/\n')
        parsed = copy.deepcopy(_parse_input(infile))
        parsed['gb']['gbrun'] = True
        return parsed


def _parse_input(infile):
    parser = copy.deepcopy(input_file)
    for namelist in parser.namelists.values():
        namelist.open = False
    return parser.Parse(infile.as_posix())


def _base_files():
    return SimpleNamespace(
        input_file='mmpbsa.in',
        ligand_mol2=None,
        complex_top='topol.top',
        receptor_tpr=None,
        receptor_trajs=[],
        receptor_top=None,
        ligand_tpr=None,
        ligand_trajs=[],
        ligand_top=None,
    )


def _stub_module(name, **attrs):
    module = types.ModuleType(name)
    for attr, value in attrs.items():
        setattr(module, attr, value)
    sys.modules[name] = module
    return module


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


class ExplicitWaterInputTest(unittest.TestCase):
    def test_input_parser_sets_explicit_water_defaults(self):
        with TemporaryDirectory() as tmpdir:
            infile = Path(tmpdir) / 'mmpbsa.in'
            infile.write_text('&general\n/\n&gb\n/\n')
            parsed = _parse_input(infile)

        self.assertEqual(parsed['general']['explicit_waters'], 0)
        self.assertEqual(parsed['general']['explicit_waters_mask'], '')
        self.assertEqual(parsed['general']['explicit_waters_group'], '')
        self.assertEqual(parsed['general']['explicit_waters_pymol_cutoff'], 0.5)
        self.assertEqual(parsed['general']['explicit_waters_as'], 'receptor')
        self.assertEqual(parsed['general']['explicit_waters_extra_points'], 'error')

    def test_rejects_missing_explicit_water_mask(self):
        main = _import_main_with_stubs()
        app = main.MMPBSA_App.__new__(main.MMPBSA_App)
        app.master = True
        app.FILES = _base_files()
        app.INPUT = _base_input()
        app.INPUT['general']['explicit_waters'] = 10

        logging.disable(logging.CRITICAL)
        try:
            with self.assertRaises(InputError):
                app.check_for_bad_input()
        finally:
            logging.disable(logging.NOTSET)

    def test_rejects_invalid_explicit_water_extra_point_mode(self):
        main = _import_main_with_stubs()
        app = main.MMPBSA_App.__new__(main.MMPBSA_App)
        app.master = True
        app.FILES = _base_files()
        app.INPUT = _base_input()
        app.INPUT['general']['explicit_waters_extra_points'] = 'ignore'

        logging.disable(logging.CRITICAL)
        try:
            with self.assertRaises(InputError):
                app.check_for_bad_input()
        finally:
            logging.disable(logging.NOTSET)

    def test_accepts_pb_explicit_water_calculation(self):
        main = _import_main_with_stubs()
        app = main.MMPBSA_App.__new__(main.MMPBSA_App)
        app.master = True
        app.FILES = _base_files()
        app.INPUT = _base_input()
        app.INPUT['gb']['gbrun'] = False
        app.INPUT['pb']['pbrun'] = True
        app.INPUT['general']['explicit_waters'] = 10
        app.INPUT['general']['explicit_waters_mask'] = ':1-608'

        app.check_for_bad_input()

    def test_rejects_explicit_waters_without_gb_or_pb(self):
        main = _import_main_with_stubs()
        app = main.MMPBSA_App.__new__(main.MMPBSA_App)
        app.master = True
        app.FILES = _base_files()
        app.INPUT = _base_input()
        app.INPUT['gb']['gbrun'] = False
        app.INPUT['pb']['pbrun'] = False
        app.INPUT['general']['explicit_waters'] = 10
        app.INPUT['general']['explicit_waters_mask'] = ':1-608'

        logging.disable(logging.CRITICAL)
        try:
            with self.assertRaises(InputError):
                app.check_for_bad_input()
        finally:
            logging.disable(logging.NOTSET)

    def test_accepts_explicit_waters_with_decomposition(self):
        main = _import_main_with_stubs()
        app = main.MMPBSA_App.__new__(main.MMPBSA_App)
        app.master = True
        app.FILES = _base_files()
        app.INPUT = _base_input()
        app.INPUT['general']['explicit_waters'] = 10
        app.INPUT['general']['explicit_waters_mask'] = ':1-608'
        app.INPUT['decomp']['decomprun'] = True
        app.INPUT['decomp']['idecomp'] = 1

        app.check_for_bad_input()

    def test_accepts_explicit_waters_with_alanine_scanning(self):
        main = _import_main_with_stubs()
        app = main.MMPBSA_App.__new__(main.MMPBSA_App)
        app.master = True
        app.FILES = _base_files()
        app.INPUT = _base_input()
        app.INPUT['general']['explicit_waters'] = 10
        app.INPUT['general']['explicit_waters_mask'] = ':1-608'
        app.INPUT['general']['netcdf'] = ''
        app.INPUT['ala']['alarun'] = True
        app.INPUT['ala']['mutant_res'] = 'A/1'

        app.check_for_bad_input()

    def test_rejects_explicit_waters_outside_st(self):
        main = _import_main_with_stubs()
        app = main.MMPBSA_App.__new__(main.MMPBSA_App)
        app.master = True
        app.FILES = _base_files()
        app.FILES.receptor_tpr = 'rec.tpr'
        app.INPUT = _base_input()
        app.INPUT['general']['explicit_waters'] = 10
        app.INPUT['general']['explicit_waters_mask'] = ':1-608'

        logging.disable(logging.CRITICAL)
        try:
            with self.assertRaises(InputError):
                app.check_for_bad_input()
        finally:
            logging.disable(logging.NOTSET)

    def test_requires_pymol_when_pymol_interface_selection_is_requested(self):
        from GMXMMPBSA.utils import find_progs

        parsed = _base_input()
        parsed['general']['explicit_waters'] = 10
        parsed['general']['explicit_waters_mask'] = 'pymol'

        def fake_which(program, path=None):
            return None if program == 'pymol' else f'/usr/bin/{program}'

        logging.disable(logging.CRITICAL)
        try:
            with patch('GMXMMPBSA.utils.shutil.which', side_effect=fake_which):
                with self.assertRaises(MMPBSA_Error):
                    find_progs(parsed, engine='amber')
        finally:
            logging.disable(logging.NOTSET)


class ExplicitWaterCleanupTest(unittest.TestCase):
    @staticmethod
    def _import_make_top_with_stubs():
        sys.modules.pop('GMXMMPBSA.make_top', None)
        parmed = sys.modules.get('parmed') or types.ModuleType('parmed')
        parmed.__version__ = getattr(parmed, '__version__', 'stub')
        parmed.__path__ = []
        sys.modules['parmed'] = parmed
        tools = types.ModuleType('parmed.tools')
        tools.__path__ = []
        sys.modules['parmed.tools'] = tools
        changeradii = types.ModuleType('parmed.tools.changeradii')
        changeradii.ChRad = type('ChRad', (), {})
        sys.modules['parmed.tools.changeradii'] = changeradii
        from GMXMMPBSA.make_top import CheckMakeTop
        return CheckMakeTop

    def test_within_explicit_water_mask_resolves_to_interface_residue_mask(self):
        CheckMakeTop = self._import_make_top_with_stubs()
        from GMXMMPBSA.utils import Residue

        maketop = CheckMakeTop.__new__(CheckMakeTop)
        maketop.explicit_waters = 10
        maketop.explicit_waters_mask = 'within 4'
        maketop.INPUT = {'general': {'explicit_waters_mask': 'within 4'}}
        selected = [
            Residue(1, 1, 'A', 'R', 1, 'ALA'),
            Residue(3, 3, 'A', 'R', 3, 'ASP'),
            Residue(4, 4, 'B', 'L', 1, 'LYS'),
        ]
        maketop.get_selected_residues = lambda selection: selected

        maketop._resolve_explicit_waters_mask()

        self.assertEqual(maketop.explicit_waters_mask, ':1,3,4')
        self.assertEqual(maketop.INPUT['general']['explicit_waters_mask'], ':1,3,4')

    def test_pymol_explicit_water_mask_resolves_to_interface_residue_mask(self):
        CheckMakeTop = self._import_make_top_with_stubs()
        from GMXMMPBSA.utils import Residue

        maketop = CheckMakeTop.__new__(CheckMakeTop)
        maketop.explicit_waters = 10
        maketop.explicit_waters_mask = 'pymol'
        maketop.external_progs = {'pymol': 'pymol'}
        maketop.complex_str_file = '_GMXMMPBSA_COM.pdb'
        maketop.FILES = SimpleNamespace(prefix='_GMXMMPBSA_')
        maketop.INPUT = {'general': {'explicit_waters_mask': 'pymol', 'explicit_waters_pymol_cutoff': 0.5}}
        maketop.resl = [
            Residue(1, 1, 'A', 'R', 1, 'ALA'),
            Residue(2, 2, 'A', 'R', 2, 'ASP'),
            Residue(3, 1, 'B', 'L', 1, 'LYS'),
        ]

        class FakeProcess:
            def wait(self):
                return 0

        def fake_open(name, mode='r', *args, **kwargs):
            if name == '_GMXMMPBSA_explicit_waters_interface.dat' and 'r' in mode:
                return io.StringIO('A\t2\tASP\nB\t1\tLYS\n')
            return mock_open()()

        with patch('subprocess.Popen', return_value=FakeProcess()) as popen:
            with patch('builtins.open', side_effect=fake_open):
                maketop._resolve_explicit_waters_mask()

        popen.assert_called_once()
        self.assertEqual(maketop.explicit_waters_mask, ':2,3')
        self.assertEqual(maketop.INPUT['general']['explicit_waters_mask'], ':2,3')

    def test_cleanup_trajs_emits_closest_cpptraj_action(self):
        CheckMakeTop = self._import_make_top_with_stubs()

        class FakeProcess:
            def __init__(self, *args, **kwargs):
                self.args = args
                self.kwargs = kwargs
                self.stdout = None

            def wait(self):
                return 0

            def communicate(self, data=None):
                self.__class__.cpptraj_input = data or b''
                if data and b'closest 10 (:1-608)&(!:NA,CL,K' in data and b'solventmask' in data:
                    self.__class__.saw_closest = True
                return b'', b''

        FakeProcess.saw_closest = False
        FakeProcess.cpptraj_input = b''
        maketop = CheckMakeTop.__new__(CheckMakeTop)
        maketop.INPUT = {'general': {'solvated_trajectory': 1, 'startframe': 2, 'endframe': 20, 'interval': 3}}
        maketop.FILES = SimpleNamespace(
            complex_trajs=['com.xtc'],
            complex_tpr='com.tpr',
            complex_index='index.ndx',
            receptor_tpr=None,
            ligand_tpr=None,
            prefix='_GMXMMPBSA_',
        )
        maketop.trjconv = ['gmx', 'trjconv']
        maketop.external_progs = {'cpptraj': 'cpptraj'}
        maketop.explicit_waters = 10
        maketop.explicit_waters_mask = ':1-608'
        maketop.explicit_water_extra_point_mask = '@100'
        maketop.explicit_water_prmtop = '_GMXMMPBSA_COM_FULL_SOLVENT.prmtop'

        with patch('subprocess.Popen', side_effect=lambda *args, **kwargs: FakeProcess(*args, **kwargs)):
            with patch('GMXMMPBSA.make_top.log_subprocess_output', lambda proc: None):
                with patch('builtins.open', mock_open()):
                    maketop.cleanup_trajs()

        self.assertTrue(FakeProcess.saw_closest)
        self.assertIn(b'trajin _GMXMMPBSA_COM_full_traj_0.xtc 2 20 3', FakeProcess.cpptraj_input)
        self.assertIn(b'closestout _GMXMMPBSA_explicit_waters_closest_0.dat\nstrip @100\ntrajout',
                      FakeProcess.cpptraj_input)
        self.assertEqual(maketop.FILES.complex_trajs, ['_GMXMMPBSA_COM_traj_0.mdcrd'])
        self.assertTrue(maketop.FILES.explicit_waters_preselected)

    def test_rejects_explicit_waters_with_extra_points_by_default(self):
        CheckMakeTop = self._import_make_top_with_stubs()

        maketop = CheckMakeTop.__new__(CheckMakeTop)
        maketop.explicit_waters = 10
        maketop.explicit_waters_extra_points = 'error'
        maketop.INPUT = {'gb': {'gbrun': True, 'igb': 8}}
        parm = SimpleNamespace(atoms=[
            SimpleNamespace(name='O', type='OW', atomic_number=8, mass=16.0),
            SimpleNamespace(name='EPW', type='EP', atomic_number=0, mass=0.0),
        ])

        logging.disable(logging.CRITICAL)
        try:
            with self.assertRaises(MMPBSA_Error):
                maketop._check_explicit_waters_supported_by_energy_model(parm)
        finally:
            logging.disable(logging.NOTSET)

    def test_strips_explicit_water_extra_points_when_requested(self):
        CheckMakeTop = self._import_make_top_with_stubs()

        maketop = CheckMakeTop.__new__(CheckMakeTop)
        maketop.explicit_waters = 10
        maketop.explicit_waters_extra_points = 'strip'
        maketop.explicit_water_extra_point_mask = ''
        maketop.INPUT = {'gb': {'gbrun': False, 'igb': 8}}
        stripped = []
        parm = SimpleNamespace(atoms=[
            SimpleNamespace(name='O', type='OW', atomic_number=8, mass=16.0),
            SimpleNamespace(name='EPW', type='EP', atomic_number=0, mass=0.0),
        ], strip=lambda mask: stripped.append(mask))

        maketop._check_explicit_waters_supported_by_energy_model(parm)

        self.assertEqual(stripped, ['@2'])
        self.assertEqual(maketop.explicit_water_extra_point_mask, '@2')

    def test_explicit_water_residues_are_mapped_as_receptor_residues(self):
        CheckMakeTop = self._import_make_top_with_stubs()

        maketop = CheckMakeTop.__new__(CheckMakeTop)
        maketop.explicit_waters = 2
        maketop.explicit_water_range = '4-5'
        maketop.resi = {'REC': {'num': [[1, 2]]}}
        maketop.resl = []
        maketop.complex_pmrtop = 'COM.prmtop'
        complex_prmtop = SimpleNamespace(residues=[
            SimpleNamespace(number=1, chain='A', name='ALA', insertion_code=''),
            SimpleNamespace(number=2, chain='A', name='ASP', insertion_code=''),
            SimpleNamespace(number=1, chain='B', name='LYS', insertion_code=''),
            SimpleNamespace(number=101, chain='', name='WAT', insertion_code=''),
            SimpleNamespace(number=102, chain='', name='WAT', insertion_code=''),
        ])

        with patch('GMXMMPBSA.make_top.parmed.load_file', return_value=complex_prmtop, create=True):
            water_res = maketop._ensure_explicit_water_residues_mapped()

        self.assertEqual([res.index for res in water_res], [4, 5])
        self.assertTrue(all(res.is_receptor() for res in water_res))
        self.assertEqual([res.id_index for res in water_res], [3, 4])
        self.assertEqual([res.name for res in water_res], ['WAT', 'WAT'])
        self.assertEqual([res.index for res in maketop.resl], [4, 5])

    def test_explicit_water_group_candidates_include_common_water_names(self):
        CheckMakeTop = self._import_make_top_with_stubs()

        maketop = CheckMakeTop.__new__(CheckMakeTop)
        maketop.explicit_waters_group = ''

        candidates = maketop._explicit_water_group_candidates()

        self.assertIn('SOLV', candidates)
        self.assertIn('OPC', candidates)
        self.assertIn('TP3', candidates)

    def test_explicit_water_group_override_is_used_as_only_candidate(self):
        CheckMakeTop = self._import_make_top_with_stubs()

        maketop = CheckMakeTop.__new__(CheckMakeTop)
        maketop.explicit_waters_group = 'Water_and_ions'

        self.assertEqual(maketop._explicit_water_group_candidates(), ['Water_and_ions'])

    def test_index_group_parser_reads_bracket_names(self):
        CheckMakeTop = self._import_make_top_with_stubs()

        with TemporaryDirectory() as tmpdir:
            ndx = Path(tmpdir) / 'index.ndx'
            ndx.write_text('[ System ]\n1 2 3\n[ OPC ]\n4 5 6\n')

            groups = CheckMakeTop._get_index_group_names(ndx)

        self.assertEqual(groups, ['System', 'OPC'])

    def test_within_decomposition_selection_includes_explicit_water_residues(self):
        CheckMakeTop = self._import_make_top_with_stubs()
        from GMXMMPBSA.utils import Residue

        decomp_res = [
            Residue(1, 1, 'A', 'R', 1, 'ALA'),
            Residue(3, 1, 'B', 'L', 1, 'LYS'),
        ]
        water_res = [
            Residue(4, 101, '', 'R', 2, 'WAT'),
            Residue(5, 102, '', 'R', 3, 'WAT'),
        ]

        selected = CheckMakeTop._include_explicit_waters_in_decomp(decomp_res, water_res)

        self.assertEqual([res.index for res in selected], [1, 3, 4, 5])
        self.assertTrue(selected[-1].is_receptor())


if __name__ == '__main__':
    unittest.main()
