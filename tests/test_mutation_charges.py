import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import parmed

from GMXMMPBSA.exceptions import MutateError
from GMXMMPBSA.make_top import CheckMakeTop
from GMXMMPBSA.make_top_amber import CheckAmberTop
from GMXMMPBSA.mutation import MutationForceField


# Native ALA entries from Amber amino12/aminont12/aminoct12 and all_amino03 libraries.
ALA14 = dict(N=-0.4157, H=0.2719, CA=0.0337, HA=0.0823, CB=-0.1825,
             HB1=0.0603, HB2=0.0603, HB3=0.0603, C=0.5973, O=-0.5679)
ALA03 = dict(N=-0.404773, H=0.294276, CA=-0.027733, HA=0.120802, CB=-0.229951,
             HB1=0.077428, HB2=0.077428, HB3=0.077428, C=0.570224, O=-0.555129)
NALA = dict(N=0.1414, H1=0.1997, H2=0.1997, H3=0.1997, CA=0.0962, HA=0.0889,
            CB=-0.0597, HB1=0.0300, HB2=0.0300, HB3=0.0300, C=0.6163, O=-0.5722)
CALA = dict(N=-0.3821, H=0.2681, CA=-0.1747, HA=0.1067, CB=-0.2093,
            HB1=0.0764, HB2=0.0764, HB3=0.0764, C=0.7731, O=-0.8055, OXT=-0.8055)
THR = dict(N=-0.4157, H=0.2719, CA=-0.0389, HA=0.1007, CB=0.3654,
           HB=0.0043, CG2=-0.2438, HG21=0.0642, OG1=-0.6761, C=0.5973, O=-0.5679)
GLY14 = dict(N=-0.4157, H=0.2719, CA=-0.0252, HA2=0.0698, HA3=0.0698, C=0.5973, O=-0.5679)
GLY03 = dict(N=-0.374282, H=0.253981, CA=-0.128844, HA2=0.088859, HA3=0.088859,
             C=0.580584, O=-0.509157)
NGLY = dict(N=0.2943, H1=0.1642, H2=0.1642, H3=0.1642, CA=-0.01,
            HA2=0.0895, HA3=0.0895, C=0.6163, O=-0.5722)
CGLY = dict(N=-0.3821, H=0.2681, CA=-0.2493, HA2=0.1056, HA3=0.1056,
            C=0.7231, O=-0.7855, OXT=-0.7855)


class MutationChargeTest(unittest.TestCase):
    builders = (CheckMakeTop, CheckAmberTop)

    @staticmethod
    def _topology(residues, chamber=False, gly_type='H1'):
        structure = parmed.Structure()
        types = {}
        for name, number, mass, radius, epsilon in (
                ('CT', 6, 12.01, 1.908, 0.1094), ('N', 7, 14.01, 1.824, 0.17),
                ('OH', 8, 16.00, 1.721, 0.2104), ('H1', 1, 1.008, 1.387, 0.0157),
                ('HC', 1, 1.008, 1.487, 0.0157), ('HA3', 1, 1.008, 1.340, 0.024),
                ('HB2', 1, 1.008, 1.340, 0.028), ('H0', 1, 1.008, 1.387, 0.0157)):
            atom_type = parmed.AtomType(name, None, mass, number)
            atom_type.set_lj_params(eps=epsilon, rmin=radius)
            types[name] = atom_type
        for ri, (resname, charges) in enumerate(residues):
            for ai, (name, charge) in enumerate(charges.items()):
                if name.startswith('H'):
                    typename = ('HA3' if chamber else 'HC') if name.startswith(('HB1', 'HB2', 'HB3', 'HG')) else 'H1'
                    if resname in ('GLY', 'NGLY', 'CGLY') and name.startswith('HA'):
                        typename = 'HB2' if chamber else gly_type
                elif name.startswith('N'):
                    typename = 'N'
                elif name.startswith('O'):
                    typename = 'OH'
                else:
                    typename = 'CT'
                atom_type = types[typename]
                atom = parmed.Atom(name=name, type=typename, charge=charge,
                                   atomic_number=atom_type.atomic_number, mass=atom_type.mass)
                atom.atom_type = atom_type
                atom.xx, atom.xy, atom.xz = ri * 15.0 + ai, float(ai % 2), 0.0
                structure.add_atom(atom, resname, ri + 1, chain='A')
        parm_class = parmed.amber.ChamberParm if chamber else parmed.amber.AmberParm
        return parm_class.from_structure(structure)

    @staticmethod
    def _checker(builder, target='ALA'):
        checker = object.__new__(builder)
        checker.INPUT = {'ala': {'mutant': target, 'cas_intdiel': 0}}
        return checker

    @staticmethod
    def _roundtrip(topology):
        with TemporaryDirectory() as directory:
            path = Path(directory) / 'mutant.prmtop'
            topology.write_parm(str(path))
            return parmed.load_file(str(path))

    def _assert_charges(self, residue, expected, net=0):
        actual = {a.name: a.charge for a in residue.atoms}
        self.assertEqual(set(actual), set(expected))
        for name, charge in expected.items():
            self.assertAlmostEqual(actual[name], charge, places=7, msg=name)
        self.assertAlmostEqual(sum(actual.values()), net, places=6)

    @staticmethod
    def _gly_mutant_charges(reference):
        return {('HA1' if name == 'HA3' else name): charge for name, charge in reference.items()}

    def test_all_glycine_atoms_use_native_force_field_charges_and_types(self):
        for builder in self.builders:
            for reference, typename in ((GLY14, 'H1'), (GLY03, 'H0')):
                with self.subTest(builder=builder.__name__, force_field=typename):
                    wt = self._topology([('THR', THR), ('GLY', reference)], gly_type=typename)
                    original = [(a.name, a.type, a.charge) for a in wt.atoms]
                    mutant = self._roundtrip(self._checker(builder, 'GLY').makeMutTop(wt, 0))
                    self._assert_charges(mutant.residues[0], self._gly_mutant_charges(reference))
                    self._assert_charges(mutant.residues[1], reference)
                    self.assertEqual([a.type for a in mutant.residues[0].atoms if a.name.startswith('HA')],
                                     [typename, typename])
                    self.assertEqual([(a.name, a.type, a.charge) for a in wt.atoms], original)

    def test_glycine_terminal_templates_keep_terminal_charges(self):
        nterm = {name: THR.get(name, 0.2) for name in NALA if not name.startswith('HB')}
        nterm.update({name: THR[name] for name in ('HB', 'CG2', 'HG21', 'OG1')})
        for builder in self.builders:
            for source, reference, resname, net in ((nterm, NGLY, 'NGLY', 1),
                                                   (dict(THR, OXT=-0.8), CGLY, 'CGLY', -1)):
                with self.subTest(builder=builder.__name__, terminal=resname):
                    wt = self._topology([('THR', source), ('GLY', GLY14), (resname, reference)])
                    mutant = self._roundtrip(self._checker(builder, 'GLY').makeMutTop(wt, 0))
                    self._assert_charges(mutant.residues[0], self._gly_mutant_charges(reference), net)

    def test_c_terminal_oxygen_aliases_survive_mutation_in_both_builders(self):
        for builder in self.builders:
            for target,reference in (('ALA',CALA),('GLY',CGLY)):
                for first, second in (('O1', 'O2'), ('OC1', 'OC2')):
                    with self.subTest(builder=builder.__name__,target=target,oxygen=first):
                        terminal={(first if name=='O' else name):q for name,q in THR.items()}
                        terminal[second]=-.8
                        wanted={{'O':first,'OXT':second}.get(name,name):q for name,q in reference.items()}
                        wt=self._topology([('THR',terminal),('C'+target,wanted)])
                        mutant=self._roundtrip(self._checker(builder,target).makeMutTop(wt,0))
                        expected=self._gly_mutant_charges(wanted) if target=='GLY' else wanted
                        self._assert_charges(mutant.residues[0],expected,net=-1)

    def test_glycine_component_and_composite_use_the_source_template(self):
        for builder in self.builders:
            with self.subTest(builder=builder.__name__):
                wt = self._topology([('THR', THR), ('THR', THR), ('GLY', GLY03)], gly_type='H0')
                checker = self._checker(builder, 'GLY')
                mutant = self._roundtrip(checker.makeMutTop(wt, [0, 1]))
                for residue in mutant.residues[:2]:
                    self._assert_charges(residue, self._gly_mutant_charges(GLY03))
                component = wt.__copy__()
                component.strip(':3')
                # A separately supplied component can have a different LJ table;
                # it must still contain the reference atom types.
                mutant = self._roundtrip(checker.makeMutTop(component, [0, 1]))
                for residue in mutant.residues:
                    self._assert_charges(residue, self._gly_mutant_charges(GLY03))

    def test_glycine_unavailable_terminal_templates_report_an_error(self):
        for builder in self.builders:
            for residues, message in (([('THR', THR)], 'native GLY'),
                                      ([('THR', dict(THR, OXT=-0.8)), ('GLY', GLY14)], 'terminal state')):
                with self.subTest(builder=builder.__name__, message=message):
                    with patch.object(MutationForceField, 'candidate_sources', return_value=[]), \
                            self.assertLogs(level='ERROR'), self.assertRaisesRegex(MutateError, message):
                        self._checker(builder, 'GLY').makeMutTop(self._topology(residues), 0)

    def test_glycine_ambiguous_templates_use_closest_match_with_a_warning(self):
        other = dict(GLY14, CA=GLY14['CA'] + 0.01, N=GLY14['N'] - 0.01)
        for builder in self.builders:
            with self.subTest(builder=builder.__name__):
                wt = self._topology([('THR', THR), ('GLY', GLY14), ('GLY', other)])
                with patch.object(MutationForceField, 'candidate_sources', return_value=[]), \
                        self.assertLogs(level='WARNING') as records:
                    mutant = self._roundtrip(self._checker(builder, 'GLY').makeMutTop(wt, 0))
                self._assert_charges(mutant.residues[0], self._gly_mutant_charges(GLY14))
                self.assertEqual(len(records.output), 1)
                self.assertIn('ambiguous native references', records.output[0])

    def test_charmm_glycine_uses_native_charges_and_ha1_ha2_names(self):
        reference = dict(N=-0.47, HN=0.31, CA=-0.02, HA1=0.09, HA2=0.09, C=0.51, O=-0.51)
        source = {('HN' if name == 'H' else name): charge for name, charge in THR.items()}
        for builder in self.builders:
            with self.subTest(builder=builder.__name__):
                wt = self._topology([('THR', source), ('GLY', reference)], chamber=True)
                mutant = self._roundtrip(self._checker(builder, 'GLY').makeMutTop(wt, 0))
                self._assert_charges(mutant.residues[0], reference)

    def test_one_letter_mutant_names_use_the_same_templates(self):
        for builder in self.builders:
            for target, reference, expected in (('A', ALA14, ALA14),
                                                 ('G', GLY14, self._gly_mutant_charges(GLY14))):
                with self.subTest(builder=builder.__name__, target=target):
                    name = 'ALA' if target == 'A' else 'GLY'
                    wt = self._topology([('THR', THR), (name, reference)])
                    mutant = self._roundtrip(self._checker(builder, target).makeMutTop(wt, 0))
                    self.assertEqual(mutant.residues[0].name, name)
                    self._assert_charges(mutant.residues[0], expected)

    def test_all_alanine_atoms_use_the_source_force_field_charges(self):
        for builder in self.builders:
            for reference in (ALA14, ALA03):
                with self.subTest(builder=builder.__name__, hb=reference['HB1']):
                    wt = self._topology([('THR', THR), ('ALA', reference)])
                    original = [a.charge for a in wt.atoms]
                    mutant = self._roundtrip(self._checker(builder).makeMutTop(wt, 0))
                    self._assert_charges(mutant.residues[0], reference)
                    self._assert_charges(mutant.residues[1], reference)
                    self.assertEqual([a.charge for a in wt.atoms], original)

    def test_terminal_atom_names_select_the_correct_charge_template(self):
        nterm = {name: THR.get(name, 0.2) for name in NALA if not name.startswith('HB')}
        nterm.update({name: THR[name] for name in ('HB', 'CG2', 'HG21', 'OG1')})
        cterm = dict(THR, OXT=-0.8)
        for builder in self.builders:
            for source, reference, resname, net in ((nterm, NALA, 'NALA', 1), (cterm, CALA, 'CALA', -1)):
                with self.subTest(builder=builder.__name__, terminal=resname):
                    wt = self._topology([('THR', source), ('ALA', ALA14), (resname, reference)])
                    mutant = self._roundtrip(self._checker(builder).makeMutTop(wt, 0))
                    self._assert_charges(mutant.residues[0], reference, net)

    def test_separate_component_without_ala_uses_the_complex_template(self):
        for builder in self.builders:
            with self.subTest(builder=builder.__name__):
                wt = self._topology([('THR', THR), ('ALA', ALA14)])
                checker = self._checker(builder)
                checker.makeMutTop(wt, 0)
                component = wt.__copy__()
                component.strip(':2')
                mutant = self._roundtrip(checker.makeMutTop(component, 0))
                self._assert_charges(mutant.residues[0], ALA14)

    def test_composite_mutations_do_not_use_uncorrected_mutants_as_references(self):
        for builder in self.builders:
            with self.subTest(builder=builder.__name__):
                wt = self._topology([('THR', THR), ('THR', THR), ('ALA', ALA03)])
                mutant = self._roundtrip(self._checker(builder).makeMutTop(wt, [0, 1]))
                for residue in mutant.residues:
                    self._assert_charges(residue, ALA03)

    def test_missing_native_ala_and_unavailable_libraries_report_an_error(self):
        for builder in self.builders:
            with self.subTest(builder=builder.__name__):
                wt = self._topology([('THR', THR)])
                with patch.object(MutationForceField, 'candidate_sources', return_value=[]), \
                        self.assertLogs(level='ERROR'), self.assertRaisesRegex(MutateError, 'native ALA'):
                    self._checker(builder).makeMutTop(wt, 0)

    def test_internal_ala_cannot_supply_c_terminal_charges(self):
        for builder in self.builders:
            with self.subTest(builder=builder.__name__):
                wt = self._topology([('THR', dict(THR, OXT=-0.8)), ('ALA', ALA14)])
                with patch.object(MutationForceField, 'candidate_sources', return_value=[]), \
                        self.assertLogs(level='ERROR'), self.assertRaisesRegex(MutateError, 'terminal state'):
                    self._checker(builder).makeMutTop(wt, 0)

    def test_conflicting_native_ala_charge_patterns_use_closest_match_with_a_warning(self):
        other = dict(ALA14, CA=ALA14['CA'] - 0.01, CB=ALA14['CB'] + 0.01)
        for builder in self.builders:
            with self.subTest(builder=builder.__name__):
                wt = self._topology([('THR', THR), ('ALA', ALA14), ('ALA', other)])
                with patch.object(MutationForceField, 'candidate_sources', return_value=[]), \
                        self.assertLogs(level='WARNING') as records:
                    mutant = self._roundtrip(self._checker(builder).makeMutTop(wt, 0))
                self._assert_charges(mutant.residues[0], other)
                self.assertEqual(len(records.output), 1)
                self.assertIn('ambiguous native references', records.output[0])

    def test_current_topology_reference_takes_precedence_over_cached_charges(self):
        for builder in self.builders:
            with self.subTest(builder=builder.__name__):
                checker = self._checker(builder)
                checker.makeMutTop(self._topology([('THR', THR), ('ALA', ALA14)]), 0)
                mutant = checker.makeMutTop(self._topology([('THR', THR), ('ALA', ALA03)]), 0)
                self._assert_charges(self._roundtrip(mutant).residues[0], ALA03)

    def test_charmm_uses_its_native_ala_charges(self):
        reference = dict(N=-0.47, HN=0.31, CA=0.07, HA=0.09, CB=-0.27,
                         HB1=0.09, HB2=0.09, HB3=0.09, C=0.51, O=-0.51)
        source = {('HN' if name == 'H' else name): charge for name, charge in THR.items()}
        for builder in self.builders:
            with self.subTest(builder=builder.__name__):
                wt = self._topology([('THR', source), ('ALA', reference)], chamber=True)
                mutant = self._roundtrip(self._checker(builder).makeMutTop(wt, 0))
                self._assert_charges(mutant.residues[0], reference)

    def test_pdb_route_does_not_require_a_charge_template(self):
        for builder in self.builders:
            with self.subTest(builder=builder.__name__):
                wt = self._topology([('THR', THR)])
                mutant = self._checker(builder).makeMutTop(wt, 0, pdb=True)
                for atom in mutant.atoms:
                    self.assertEqual(atom.charge, THR[atom.name])


if __name__ == '__main__':
    unittest.main()
