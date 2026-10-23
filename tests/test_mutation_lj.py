import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import parmed

from GMXMMPBSA.exceptions import MutateError
from GMXMMPBSA.make_top import CheckMakeTop


class MutationLennardJonesTest(unittest.TestCase):
    @staticmethod
    def _topology(residues, chamber=False):
        """Small parameterized residues, independent of installed force-field files."""
        structure = parmed.Structure()
        types = {}
        for name, number, mass, radius, epsilon in (
                ('CT', 6, 12.01, 1.908, 0.1094),
                ('OH', 8, 16.00, 1.721, 0.2104),
                ('SH', 16, 32.06, 2.000, 0.2500),
                ('HC', 1, 1.008, 1.487, 0.0157),
                ('H1', 1, 1.008, 1.387, 0.0157),
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

    @staticmethod
    def _mutate(topology, indices=0, target='ALA', pdb=False):
        checker = object.__new__(CheckMakeTop)
        checker.INPUT = {'ala': {'mutant': target, 'cas_intdiel': 0}}
        return checker.makeMutTop(topology, indices, pdb=pdb)

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
                    ('ALA', [('CB', 'CT'), ('HB1', 'HC')])])
                original = [(a.name, a.type, a.nb_idx) for a in wt.atoms]
                wt_disk = self._roundtrip(wt)
                mutant = self._roundtrip(self._mutate(wt))
                reference = mutant.residues[1].atoms[-1]
                self._assert_lj(mutant, 0, ('HB1', 'HB2', 'HB3'), reference)
                self.assertEqual([(a.name, a.type, a.nb_idx) for a in wt.atoms], original)
                for flag in ('NONBONDED_PARM_INDEX', 'LENNARD_JONES_ACOEF', 'LENNARD_JONES_BCOEF'):
                    self.assertEqual(mutant.parm_data[flag], wt_disk.parm_data[flag])

    def test_glycine_uses_h1_including_reference_in_mutated_residue(self):
        wt = self._topology([('SER', [('CA', 'CT'), ('HA', 'H1'), ('CB', 'CT'), ('OG', 'OH')])])
        reference = self._roundtrip(wt).atoms[1]
        mutant = self._roundtrip(self._mutate(wt, target='GLY'))
        self._assert_lj(mutant, 0, ('HA1', 'HA2'), reference)

    def test_reference_removed_by_stripping_is_still_available(self):
        wt = self._topology([('THR', [('CA', 'CT'), ('CB', 'CT'), ('HB', 'H1'),
                                      ('CG2', 'CT'), ('HG21', 'HC'), ('OG1', 'OH')])])
        reference = self._roundtrip(wt).atoms[4]
        mutant = self._roundtrip(self._mutate(wt))
        self._assert_lj(mutant, 0, ('HB1', 'HB2', 'HB3'), reference)

    def test_composite_mutations_keep_the_same_hc_type(self):
        ser = ('SER', [('CA', 'CT'), ('HA', 'H1'), ('CB', 'CT'),
                       ('HB2', 'H1'), ('HB3', 'H1'), ('OG', 'OH')])
        wt = self._topology([ser, ser, ('ALA', [('CB', 'CT'), ('HB1', 'HC')])])
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
                    ('ALA', [('CB', 'CT'), ('HB1', typename)])], chamber=True)
                wt_disk = self._roundtrip(wt)
                mutant = self._roundtrip(self._mutate(wt, target=target))
                self._assert_lj(mutant, 0, names, mutant.residues[1].atoms[-1], chamber=True)
                for flag in ('LENNARD_JONES_14_ACOEF', 'LENNARD_JONES_14_BCOEF'):
                    self.assertEqual(mutant.parm_data[flag], wt_disk.parm_data[flag])

    def test_pdb_mutation_does_not_require_lj_parameters(self):
        wt = self._topology([('SER', [('CA', 'CT'), ('HA', 'H1'), ('CB', 'CT'), ('OG', 'OH')])])
        mutant = self._mutate(wt, pdb=True)
        self.assertEqual([a.name for a in mutant.atoms], ['CA', 'HA', 'CB'])

    def test_missing_target_type_fails_instead_of_writing_heavy_atom_lj(self):
        wt = self._topology([('SER', [('CA', 'CT'), ('CB', 'CT'), ('OG', 'OH')])])
        with self.assertRaisesRegex(MutateError, 'Lennard-Jones.*HC'):
            self._mutate(wt)

    def test_ambiguous_target_type_fails_instead_of_choosing_arbitrarily(self):
        wt = self._topology([
            ('SER', [('CA', 'CT'), ('CB', 'CT'), ('OG', 'OH')]),
            ('ALA', [('CB', 'CT'), ('HB1', 'HC'), ('HB2', 'HC')])])
        wt.residues[1].atoms[-1].nb_idx = wt.residues[1].atoms[0].nb_idx
        wt.remake_parm()
        with self.assertRaisesRegex(MutateError, 'Lennard-Jones.*HC'):
            self._mutate(wt)


if __name__ == '__main__':
    unittest.main()
