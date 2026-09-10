import csv
import json
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
from unittest.mock import patch

from GMXMMPBSA.exceptions import PrmtopError
from GMXMMPBSA.parm_setup import MMPBSA_System, Residue
from GMXMMPBSA.radii import (_log_advisories, _record, collect_radii_provenance, radius_checksum,
                              source_force_field_family)


def _parm(radius_set='modified Bondi radii (mbondi3)', radii=None, screen=None):
    residue = SimpleNamespace(name='ALA', number=1, idx=0)
    atoms = [
        SimpleNamespace(name='CA', type='CT', element='C', atomic_number=6, charge=-0.1,
                        sigma=1.9, residue=residue),
        SimpleNamespace(name='O', type='O', element='O', atomic_number=8, charge=-0.5,
                        sigma=1.7, residue=residue),
    ]
    return SimpleNamespace(
        atoms=atoms,
        residues=[residue],
        chamber=False,
        parm_data={
            'RADIUS_SET': [radius_set],
            'RADII': radii or [1.7, 1.5],
            'SCREEN': screen or [0.8, 0.85],
        },
    )


class RadiusProvenanceTest(unittest.TestCase):
    def test_checksum_is_stable(self):
        self.assertEqual(radius_checksum([1, 2.0]), radius_checksum([1.0, 2]))
        self.assertNotEqual(radius_checksum([1, 2]), radius_checksum([1, 2.01]))

    def test_force_field_override_wins_over_auto_detection(self):
        input_data = {'general': {'source_force_field': 'opls', 'forcefields': ['leaprc.protein.ff19SB']}}
        self.assertEqual(source_force_field_family(input_data, _parm()), 'opls')

    def test_audit_flags_metals_extra_points_and_unknown_elements(self):
        parm = _parm()
        residue = parm.atoms[0].residue
        parm.atoms.extend([
            SimpleNamespace(name='ZN', type='ZN', element='Zn', atomic_number=30, charge=2.0,
                            sigma=1.2, residue=residue),
            SimpleNamespace(name='EP', type='EP', element='', atomic_number=0, charge=0.0,
                            sigma=0.0, residue=residue),
        ])
        parm.parm_data['RADII'].extend([1.1, 0.0])
        parm.parm_data['SCREEN'].extend([0.5, 0.0])
        record = _record(parm, 'complex', 'mbondi3', 'native', {'gb': {'gbrun': False}}, 'amber')
        self.assertEqual(record['flag_counts']['metal'], 1)
        self.assertEqual(record['flag_counts']['dummy_or_extra_point'], 1)
        self.assertEqual(record['flag_counts']['unknown_element'], 1)

    def test_real_parmed_elements_in_records_and_csv(self):
        import parmed
        from GMXMMPBSA.radii import _write_audit_csv
        parm = _parm()
        structure = parmed.Structure()
        for name, number in [('ZN', 30), ('FE', 26), ('C', 6), ('EP', 0)]:
            structure.add_atom(parmed.Atom(name=name, type=name, atomic_number=number), 'ALA', 1)
        parm.atoms = structure.atoms
        parm.parm_data['RADII'] = [1.1, 1.2, 1.7, 0.]
        parm.parm_data['SCREEN'] = [.5, .5, .8, 0.]
        record = _record(parm, 'complex', 'mbondi3', 'native', {}, 'amber')
        self.assertEqual(record['flag_counts']['metal'], 2)
        self.assertEqual(record['flag_counts']['unknown_element'], 1)
        self.assertEqual(record['flag_counts']['dummy_or_extra_point'], 1)
        self.assertEqual(json.loads(json.dumps(record))['flag_counts']['metal'], 2)
        with TemporaryDirectory() as directory:
            path = Path(directory) / 'audit.csv'
            _write_audit_csv(parm, 'complex', path)
            with path.open() as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual([r['element'] for r in rows], ['ZN', 'FE', 'C', ''])
            self.assertEqual([r['assignment_category'] for r in rows],
                             ['metal_unvalidated', 'metal_unvalidated', 'radius_set_rule', 'dummy_or_extra_point'])

    def test_native_mutants_reapply_inherited_radii_but_aliases_are_preserved(self):
        files = SimpleNamespace(complex_prmtop='COM.prmtop', receptor_prmtop='REC.prmtop',
                                ligand_prmtop='LIG.prmtop', mutant_complex_prmtop='MUT_COM.prmtop',
                                mutant_receptor_prmtop='MUT_REC.prmtop', mutant_ligand_prmtop='LIG.prmtop')
        data = {'general': {'PBRadii': 2}}
        import os
        with TemporaryDirectory() as directory, patch('GMXMMPBSA.radii.parmed.load_file', return_value=_parm()):
            old = os.getcwd()
            try:
                os.chdir(directory)
                for filename in vars(files).values():
                    Path(filename).touch()
                records = collect_radii_provenance(files, data, 'amber')['components']
                for name in ('mutant_complex', 'mutant_receptor'):
                    self.assertEqual(records[name]['assignment_route'], 'native_amber_mutant_inherited_ChRad')
                    self.assertFalse(records[name]['input_pbradii_applied'])
                    self.assertIn('reapplied', records[name]['effective_radius_source'])
                for name in ('complex', 'receptor', 'ligand', 'mutant_ligand'):
                    self.assertEqual(records[name]['assignment_route'], 'native_amber_topology_preserved')
                    self.assertFalse(records[name]['input_pbradii_applied'])
                gmx = collect_radii_provenance(files, data, 'gmx')['components']
                self.assertTrue(all(r['assignment_route'] == 'parmed_ChRad' and r['input_pbradii_applied']
                                    for r in gmx.values()))
            finally:
                os.chdir(old)

    def test_force_specific_advisories_are_emitted(self):
        base = _record(_parm(), 'complex', 'mbondi3', 'native',
                       {'gb': {'gbrun': True, 'igb': 8}, 'pb': {'pbrun': True},
                        'gbnsr6': {'gbnsr6run': True}}, 'charmm')
        with self.assertLogs(level='INFO') as messages:
            _log_advisories(base, {'gb': {'gbrun': True, 'igb': 8}, 'pb': {'pbrun': True},
                                    'gbnsr6': {'gbnsr6run': True}})
        output = '\n'.join(messages.output)
        self.assertIn('cross-parameterization', output)
        self.assertIn('charmm_radii', output)
        self.assertIn('GBNSR6 radius provenance is reported independently', output)

    def test_writes_compact_json_and_optional_complete_csv(self):
        files = SimpleNamespace(
            complex_prmtop='COM.prmtop', receptor_prmtop='REC.prmtop', ligand_prmtop='LIG.prmtop',
            mutant_complex_prmtop='MUT_COM.prmtop', mutant_receptor_prmtop='MUT_REC.prmtop',
            mutant_ligand_prmtop=None,
        )
        input_data = {
            'general': {'PBRadii': 4, 'radii_audit': 1, 'source_force_field': 'amber'},
            'gb': {'gbrun': True, 'igb': 8}, 'pb': {'pbrun': False},
            'gbnsr6': {'gbnsr6run': False},
        }
        with TemporaryDirectory() as directory, patch('GMXMMPBSA.radii.parmed.load_file', return_value=_parm()):
            # The collector's JSON path is relative to the current directory;
            # run it in an isolated directory without changing the test process.
            import os
            old_cwd = os.getcwd()
            os.chdir(directory)
            try:
                for filename in ('COM.prmtop', 'REC.prmtop', 'LIG.prmtop'):
                    Path(filename).touch()
                result = collect_radii_provenance(files, input_data, 'amber')
            finally:
                os.chdir(old_cwd)
            self.assertEqual(set(result['components']), {'complex', 'receptor', 'ligand'})
            self.assertEqual(result['components']['complex']['effective_radius_set'], 'mbondi3')
            self.assertFalse(result['components']['complex']['input_pbradii_applied'])
            self.assertIn('native topology RADII/SCREEN', result['input_pbradii_semantics'])
            self.assertTrue((Path(directory) / 'GMXMMPBSA_radii.json').exists())
            with open(Path(directory) / 'GMXMMPBSA_radii_complex.csv', newline='', encoding='utf-8') as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 2)
            self.assertEqual(rows[0]['assigned_continuum_radius'], '1.7')
            with open(Path(directory) / 'GMXMMPBSA_radii.json', encoding='utf-8') as handle:
                payload = json.load(handle)
            self.assertEqual(payload['components']['complex']['screen_sha256'], radius_checksum([0.8, 0.85]))

    def test_records_prepared_gbnsr6_topologies_separately(self):
        files = SimpleNamespace(
            complex_prmtop='COM.prmtop', receptor_prmtop='REC.prmtop', ligand_prmtop='LIG.prmtop',
            mutant_complex_prmtop=None, mutant_receptor_prmtop=None, mutant_ligand_prmtop=None,
        )
        input_data = {
            'general': {'PBRadii': 4, 'radii_audit': 0},
            'gb': {'gbrun': False}, 'pb': {'pbrun': False},
            'gbnsr6': {'gbnsr6run': True},
        }
        with TemporaryDirectory() as directory, patch('GMXMMPBSA.radii.parmed.load_file', return_value=_parm()):
            old_cwd = __import__('os').getcwd()
            __import__('os').chdir(directory)
            try:
                for filename in ('COM.prmtop', 'REC.prmtop', 'LIG.prmtop', 'prepared.prmtop'):
                    Path(filename).touch()
                result = collect_radii_provenance(
                    files, input_data, 'gmx',
                    additional_topologies={'gbnsr6_complex': 'prepared.prmtop'},
                )
            finally:
                __import__('os').chdir(old_cwd)

        prepared = result['prepared_topologies']['gbnsr6_complex']
        self.assertEqual(prepared['assignment_route'], 'gbnsr6_prepared_copy')
        self.assertEqual(prepared['effective_radius_set'], 'mbondi3')


class MappedRadiusConsistencyTest(unittest.TestCase):
    @staticmethod
    def _system(receptor_radii):
        complex_top = SimpleNamespace(
            atoms=[SimpleNamespace(residue=SimpleNamespace(idx=0)),
                   SimpleNamespace(residue=SimpleNamespace(idx=0))],
            parm_data={'RESIDUE_POINTER': [1], 'CHARGE': [0.1, -0.2],
                       'RADIUS_SET': ['complex-label'],
                       'RADII': [1.7, 1.5], 'SCREEN': [0.8, 0.85]},
            ptr=lambda key: {'natom': 2, 'nres': 1}[key],
        )
        receptor_top = SimpleNamespace(
            parm_data={'RESIDUE_POINTER': [1], 'CHARGE': [0.1, -0.2],
                       'RADIUS_SET': ['component-label'],
                       'RADII': receptor_radii, 'SCREEN': [0.8, 0.85]},
            ptr=lambda key: {'natom': 2, 'nres': 1}[key],
        )
        ligand_top = SimpleNamespace(
            parm_data={'RADII': [], 'SCREEN': []},
            ptr=lambda key: {'natom': 0, 'nres': 0}[key],
        )
        system = object.__new__(MMPBSA_System)
        system.stability = False
        system.mapped = True
        mapped = Residue(1, 'ALA')
        mapped.receptor_number = 1
        system.res_list = [mapped]
        system.complex_prmtop = complex_top
        system.receptor_prmtop = receptor_top
        system.ligand_prmtop = ligand_top
        return system

    def test_different_radius_labels_are_allowed_when_mapped_arrays_match(self):
        # Labels are not consulted by CheckConsistency; the mapped arrays are.
        self._system([1.7, 1.5]).CheckConsistency()

    def test_mapped_radius_difference_is_rejected(self):
        with self.assertRaisesRegex(PrmtopError, 'Inconsistent RADII definition'):
            self._system([1.7, 1.6]).CheckConsistency()


if __name__ == '__main__':
    unittest.main()
