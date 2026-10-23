import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace

import parmed

from GMXMMPBSA.make_top import CheckMakeTop
from GMXMMPBSA.make_top_amber import CheckAmberTop
from GMXMMPBSA.mutation import MutationForceField, missing_lj14_pairs, fill_lj14_pairs


class MutationLennardJonesTest(unittest.TestCase):
    builder = CheckMakeTop

    @staticmethod
    def _topology(residues, chamber=False):
        """Small parameterized residues, independent of installed force-field files."""
        structure = parmed.Structure()
        types = {}
        for name, number, mass, radius, epsilon in (
                ('CT', 6, 12.01, 1.908, 0.1094),
                ('CT3', 6, 12.01, 2.060, 0.0800),
                ('OH', 8, 16.00, 1.721, 0.2104),
                ('SH', 16, 32.06, 2.000, 0.2500),
                ('HC', 1, 1.008, 1.487, 0.0157),
                ('H1', 1, 1.008, 1.387, 0.0157),
                ('H0', 1, 1.008, 1.387, 0.0157),
                ('HA3', 1, 1.008, 1.340, 0.0240),
                ('HB2', 1, 1.008, 1.340, 0.0280)):
            atom_type = parmed.AtomType(name, None, mass, number)
            atom_type.set_lj_params(eps=epsilon, rmin=radius,
                                    eps14=epsilon / 2, rmin14=radius - 0.1)
            types[name] = atom_type
        for ri, (resname, atoms) in enumerate(residues):
            for ai, (name, typename) in enumerate(atoms):
                atom_type = types[typename]
                atom = parmed.Atom(name=name, type=typename, mass=atom_type.mass,
                                   atomic_number=atom_type.atomic_number)
                atom.atom_type = atom_type
                atom.xx, atom.xy, atom.xz = ri * 10.0 + ai, float(ai % 2), 0.0
                structure.add_atom(atom, resname, ri + 1, chain='A')
        parm_class = parmed.amber.ChamberParm if chamber else parmed.amber.AmberParm
        return parm_class.from_structure(structure)

    def _mutate(self, topology, indices=0, target='ALA', pdb=False, checker=None):
        checker = checker or object.__new__(self.builder)
        checker.INPUT = {'ala': {'mutant': target, 'cas_intdiel': 0}}
        return checker.makeMutTop(topology, indices, pdb=pdb)

    @staticmethod
    def _reference(target='ALA', hydrogen_type=None):
        if target == 'GLY':
            return ('GLY', [('CA', 'CT'), ('HA2', hydrogen_type or 'H1'),
                            ('HA3', hydrogen_type or 'H1')])
        return ('ALA', [('CA', 'CT'), ('HA', 'H1'), ('CB', 'CT'),
                        ('HB1', hydrogen_type or 'HC'), ('HB2', hydrogen_type or 'HC'),
                        ('HB3', hydrogen_type or 'HC')])

    @staticmethod
    def _roundtrip(topology):
        with TemporaryDirectory() as directory:
            path = Path(directory) / 'mutant.prmtop'
            topology.write_parm(str(path))
            return parmed.load_file(str(path))

    def _assert_lj(self, topology, residue_index, names, reference, chamber=False):
        atoms = {a.name: a for a in topology.residues[residue_index].atoms}
        for name in names:
            with self.subTest(atom=name):
                atom = atoms[name]
                self.assertEqual(atom.type, reference.type)
                self.assertEqual(atom.nb_idx, reference.nb_idx)
                self.assertAlmostEqual(atom.rmin, reference.rmin, places=7)
                self.assertAlmostEqual(atom.epsilon, reference.epsilon, places=7)
                if chamber:
                    self.assertAlmostEqual(atom.rmin_14, reference.rmin_14, places=7)
                    self.assertAlmostEqual(atom.epsilon_14, reference.epsilon_14, places=7)

    def test_alanine_hydrogens_use_hc_in_written_prmtop(self):
        for resname, sidechain in (
                ('THR', [('HB', 'H1'), ('CG2', 'CT'), ('OG1', 'OH')]),
                ('ILE', [('HB', 'HC'), ('CG2', 'CT'), ('CG1', 'CT')]),
                ('VAL', [('HB', 'HC'), ('CG2', 'CT'), ('CG1', 'CT')]),
                ('SER', [('HB2', 'H1'), ('HB3', 'H1'), ('OG', 'OH')]),
                ('CYS', [('HB2', 'HC'), ('HB3', 'HC'), ('SG', 'SH')]),
                ('LEU', [('HB2', 'HC'), ('HB3', 'HC'), ('CG', 'CT')])):
            with self.subTest(residue=resname):
                wt = self._topology([
                    (resname, [('CA', 'CT'), ('HA', 'H1'), ('CB', 'CT')] + sidechain),
                    self._reference()])
                original = [(a.name, a.type, a.nb_idx) for a in wt.atoms]
                wt_disk = self._roundtrip(wt)
                mutant = self._roundtrip(self._mutate(wt))
                reference = mutant.residues[1].atoms[-1]
                self._assert_lj(mutant, 0, ('HB1', 'HB2', 'HB3'), reference)
                self.assertEqual([(a.name, a.type, a.nb_idx) for a in wt.atoms], original)
                for flag in ('NONBONDED_PARM_INDEX', 'LENNARD_JONES_ACOEF', 'LENNARD_JONES_BCOEF'):
                    self.assertEqual(mutant.parm_data[flag], wt_disk.parm_data[flag])

    def test_glycine_hydrogens_use_the_native_type_in_written_prmtop(self):
        wt = self._topology([('SER', [('CA', 'CT'), ('HA', 'H1'), ('CB', 'CT'), ('OG', 'OH')]),
                             self._reference('GLY', 'H0')])
        reference = self._roundtrip(wt).residues[1].atoms[-1]
        mutant = self._roundtrip(self._mutate(wt, target='GLY'))
        self._assert_lj(mutant, 0, ('HA1', 'HA2'), reference)

    def test_reference_removed_by_stripping_is_still_available(self):
        wt = self._topology([('THR', [('CA', 'CT'), ('CB', 'CT'), ('HB', 'H1'),
                                      ('HA', 'H1'), ('CG2', 'CT'), ('HG21', 'HC'), ('OG1', 'OH')]),
                             self._reference()])
        checker = object.__new__(self.builder)
        checker._remember_native_mutant_parameters(wt, 'ALA')
        wt.strip(':2')
        reference = self._roundtrip(wt).atoms[5]
        mutant = self._roundtrip(self._mutate(wt, checker=checker))
        self._assert_lj(mutant, 0, ('HB1', 'HB2', 'HB3'), reference)

    def test_composite_mutations_keep_the_same_hc_type(self):
        ser = ('SER', [('CA', 'CT'), ('HA', 'H1'), ('CB', 'CT'),
                       ('HB2', 'H1'), ('HB3', 'H1'), ('OG', 'OH')])
        wt = self._topology([ser, ser, self._reference()])
        mutant = self._roundtrip(self._mutate(wt, indices=[0, 1]))
        reference = mutant.residues[2].atoms[-1]
        for index in (0, 1):
            self._assert_lj(mutant, index, ('HB1', 'HB2', 'HB3'), reference)

    def test_chamber_preserves_target_normal_and_14_lj_parameters(self):
        for target, typename, names in (('ALA', 'HA3', ('HB1', 'HB2', 'HB3')),
                                        ('GLY', 'HB2', ('HA1', 'HA2'))):
            with self.subTest(target=target):
                wt = self._topology([
                    ('SER', [('CA', 'CT'), ('HA', 'H1'), ('CB', 'CT'),
                             ('HB2', 'H1'), ('HB3', 'H1'), ('OG', 'OH')]),
                    self._reference(target, typename)], chamber=True)
                wt_disk = self._roundtrip(wt)
                mutant = self._roundtrip(self._mutate(wt, target=target))
                self._assert_lj(mutant, 0, names, mutant.residues[1].atoms[-1], chamber=True)
                for flag in ('LENNARD_JONES_14_ACOEF', 'LENNARD_JONES_14_BCOEF'):
                    self.assertEqual(mutant.parm_data[flag], wt_disk.parm_data[flag])

    def test_pdb_mutation_does_not_require_lj_parameters(self):
        wt = self._topology([('SER', [('CA', 'CT'), ('HA', 'H1'), ('CB', 'CT'), ('OG', 'OH')])])
        mutant = self._mutate(wt, pdb=True)
        self.assertEqual([a.name for a in mutant.atoms], ['CA', 'HA', 'CB'])

    def test_missing_target_type_is_added_from_the_cached_native_parameters(self):
        wt = self._topology([('SER', [('CA', 'CT'), ('HA', 'H1'), ('CB', 'CT'),
                                    ('HB2', 'H1'), ('HB3', 'H1'), ('OG', 'OH')]), self._reference()])
        checker = object.__new__(self.builder)
        checker._remember_native_mutant_parameters(wt, 'ALA')
        wt = self._topology([('SER', [('CA', 'CT'), ('HA', 'H1'), ('CB', 'CT'),
                                    ('HB2', 'H1'), ('HB3', 'H1'), ('OG', 'OH')])])
        mutant = self._roundtrip(self._mutate(wt, checker=checker))
        for atom in mutant.residues[0].atoms:
            if atom.name.startswith('HB'):
                self.assertEqual(atom.type, 'HC')
                self.assertAlmostEqual(atom.rmin, 1.487, places=7)
                self.assertAlmostEqual(atom.epsilon, 0.0157, places=7)

    def test_target_carbon_lj_is_assigned_as_well_as_hydrogen_lj(self):
        native = self._reference()
        native[1][2] = ('CB', 'CT3')
        wt = self._topology([
            ('SER', [('CA', 'CT'), ('HA', 'H1'), ('CB', 'CT'),
                     ('HB2', 'H1'), ('HB3', 'H1'), ('OG', 'OH')]), native], chamber=True)
        mutant = self._roundtrip(self._mutate(wt))
        self._assert_lj(mutant, 0, ('CB',), mutant.residues[1].atoms[2], chamber=True)

    def test_ambiguous_lj_can_select_the_hydrogen_index_with_a_warning(self):
        wt = self._topology([
            ('SER', [('CA', 'CT'), ('HA', 'H1'), ('CB', 'CT'),
                     ('HB2', 'H1'), ('HB3', 'H1'), ('OG', 'OH')]), self._reference()])
        wt.residues[1].atoms[-1].nb_idx = wt.residues[1].atoms[0].nb_idx
        wt.remake_parm()
        with self.assertLogs(level='WARNING') as records:
            mutant = self._roundtrip(self._mutate(wt))
        self.assertEqual(len(records.output), 1)
        self.assertIn('Estimated mutation LJ', records.output[0])
        for atom in mutant.residues[0].atoms:
            if atom.name.startswith('HB'):
                self.assertAlmostEqual(atom.rmin,1.487,places=7)
                self.assertAlmostEqual(atom.epsilon,.0157,places=7)

    def test_chamber_pair_specific_normal_and_14_coefficients_are_preserved(self):
        wt = self._topology([
            ('SER', [('CA', 'CT'), ('HA', 'H1'), ('CB', 'CT'),
                     ('HB2', 'H1'), ('HB3', 'H1'), ('OG', 'OH')]),
            self._reference(hydrogen_type='HA3')], chamber=True)
        ntypes = wt.ptr('ntypes')
        target_index = wt.residues[1].atoms[-1].nb_idx
        oxygen_index = wt.residues[0].atoms[-1].nb_idx
        pair = wt.parm_data['NONBONDED_PARM_INDEX'][ntypes * (target_index - 1) + oxygen_index - 1] - 1
        flags = ('LENNARD_JONES_ACOEF', 'LENNARD_JONES_BCOEF',
                 'LENNARD_JONES_14_ACOEF', 'LENNARD_JONES_14_BCOEF')
        for flag in flags:
            # A non-combining-rule interaction must survive assignment/copying.
            wt.parm_data[flag][pair] *= 1.234
        reference = self._roundtrip(wt)
        mutant = self._roundtrip(self._mutate(wt))
        for flag in flags:
            self.assertEqual(mutant.parm_data[flag], reference.parm_data[flag])
        for atom in mutant.residues[0].atoms:
            if atom.name in ('HB1', 'HB2', 'HB3'):
                self.assertEqual(atom.nb_idx, target_index)

    def test_cached_template_resolves_indices_in_the_current_component_table(self):
        ser = ('SER', [('CA', 'CT'), ('HA', 'H1'), ('CB', 'CT'),
                       ('HB2', 'H1'), ('HB3', 'H1'), ('OG', 'OH')])
        complex_top = self._topology([self._reference(), ser])
        checker = object.__new__(self.builder)
        checker._remember_native_mutant_parameters(complex_top, 'ALA')
        component = self._topology([ser, ('LIG', [('X', 'HC')])])
        reference = self._roundtrip(component).residues[1].atoms[0]
        self.assertNotEqual(reference.nb_idx, complex_top.residues[0].atoms[-1].nb_idx)
        mutant = self._roundtrip(self._mutate(component, checker=checker))
        self._assert_lj(mutant, 0, ('HB1', 'HB2', 'HB3'), reference)


class AmberMutationLennardJonesTest(MutationLennardJonesTest):
    builder = CheckAmberTop


class ChamberMissingPairsTest(unittest.TestCase):
    def test_native_builder_preserves_chamber_tables_even_without_ureys_or_impropers(self):
        for chamber in (False, True):
            with self.subTest(chamber=chamber), TemporaryDirectory() as tmp:
                root = Path(tmp)
                source = MutationLennardJonesTest._topology([
                    ('SER', [('CA', 'CT'), ('CB', 'CT3')]),
                    ('ALA', [('CA', 'CT'), ('HB1', 'H1')]),
                    ('GLY', [('CA', 'CT'), ('HA1', 'H1')])], chamber=chamber)
                for atom in source.atoms:
                    atom.solvent_radius, atom.screen = 1.5, .8
                source.parm_data['RADIUS_SET'] = ['ArgH and AspGluO modified Bondi2 radii (mbondi3)']
                source.write_parm(str(root / 'input.prmtop'))
                checker = object.__new__(CheckAmberTop)
                checker.INPUT = {'ala': {'alarun': False}, 'gb': {'gbrun': False},
                                 'pb': {'pbrun': False}, 'general': {'PBRadii': 4}}
                checker.FILES = SimpleNamespace(complex_top=str(root / 'input.prmtop'),
                    complex_mask=[':1-2', ':3'], receptor_top=None, ligand_top=None,
                    prefix=str(root / '_test_'))
                checker.explicit_waters = False
                checker.complex_str = source.__copy__()
                checker.complex_str_file = str(root / 'reference.pdb')
                checker.resi = {'REC': {'string': ['1-2']}, 'LIG': {'string': ['3']}}
                for prefix in ('complex', 'receptor', 'ligand', 'mutant_complex',
                               'mutant_receptor', 'mutant_ligand'):
                    setattr(checker, prefix + '_pmrtop', str(root / (prefix + '.prmtop')))
                checker.ambertop2prmtop()
                for component in ('complex', 'receptor', 'ligand'):
                    actual = parmed.amber.LoadParm(str(root / (component + '.prmtop')))
                    self.assertEqual(actual.chamber, chamber)
                    if chamber:
                        for flag in ('LENNARD_JONES_14_ACOEF', 'LENNARD_JONES_14_BCOEF'):
                            self.assertEqual(actual.parm_data[flag], source.parm_data[flag])

    @staticmethod
    def _pair_fixture():
        parm = MutationLennardJonesTest._topology([
            ('SER', [('CA', 'CT'), ('CB', 'CT3'), ('OG', 'OH'), ('HB', 'H1')])], chamber=True)
        dtype = parmed.DihedralType(1, 1, 0)
        parm.dihedral_types.append(dtype)
        parm.dihedrals.append(parmed.Dihedral(
            parm.atoms[0], parm.atoms[3], parm.atoms[1], parm.atoms[2], type=dtype))
        parm.remake_parm()
        original = parm.__copy__()
        parm.atoms[0].type = parm.atoms[1].type
        parm.atoms[0].nb_idx = parm.atoms[1].nb_idx
        pair = parm.ptr('ntypes') * (parm.atoms[0].nb_idx - 1) + parm.atoms[2].nb_idx - 1
        slot = parm.parm_data['NONBONDED_PARM_INDEX'][pair] - 1
        for flag in ('LENNARD_JONES_14_ACOEF', 'LENNARD_JONES_14_BCOEF'):
            parm.parm_data[flag][slot] = 0
        return parm, original, slot

    def test_newly_active_sparse_pair_uses_explicit_source_14_parameters(self):
        parm, original, slot = self._pair_fixture()
        saved = {flag: list(parm.parm_data[flag]) for flag in
                 ('LENNARD_JONES_14_ACOEF', 'LENNARD_JONES_14_BCOEF')}
        library = MutationForceField(pair_types={('CT3', 'OH'): parmed.NonbondedExceptionType(3.3, .035)})
        missing = missing_lj14_pairs(parm, original, parm.residues[0])
        self.assertEqual([pair[0] for pair in missing], [slot])
        self.assertFalse(fill_lj14_pairs(parm, missing, library))
        disk = MutationLennardJonesTest._roundtrip(parm)
        for flag, expected in (('LENNARD_JONES_14_ACOEF', .035 * 3.3**12),
                               ('LENNARD_JONES_14_BCOEF', 2 * .035 * 3.3**6)):
            self.assertAlmostEqual(disk.parm_data[flag][slot], expected, places=4)
            self.assertEqual(parm.parm_data[flag][:slot], saved[flag][:slot])
            self.assertEqual(parm.parm_data[flag][slot+1:], saved[flag][slot+1:])

    def test_previously_active_zero_pair_is_preserved(self):
        parm, _, slot = self._pair_fixture()
        original = parm.__copy__()
        original.atoms[0].nb_idx = parm.atoms[0].nb_idx
        self.assertEqual(missing_lj14_pairs(parm, original, parm.residues[0]), [])
        self.assertEqual(parm.parm_data['LENNARD_JONES_14_ACOEF'][slot], 0)
        self.assertEqual(parm.parm_data['LENNARD_JONES_14_BCOEF'][slot], 0)

    def test_explicit_zero_14_depth_is_not_replaced_by_normal_depth(self):
        parm, original, slot = self._pair_fixture()
        types = {}
        for name in ('CT3', 'OH'):
            atom_type = parmed.AtomType(name, None, 12, 6)
            atom_type.set_lj_params(.2, 2, eps14=0, rmin14=2)
            types[name] = atom_type
        library = MutationForceField(types)
        missing = missing_lj14_pairs(parm, original, parm.residues[0])
        self.assertFalse(fill_lj14_pairs(parm, missing, library))
        self.assertEqual(parm.parm_data['LENNARD_JONES_14_ACOEF'][slot], 0)
        self.assertEqual(parm.parm_data['LENNARD_JONES_14_BCOEF'][slot], 0)


if __name__ == '__main__':
    unittest.main()
