import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace

import parmed

from GMXMMPBSA.alamdcrd import GlyMutantMdcrd, MutantMdcrd, _getnumatms
from GMXMMPBSA.exceptions import MMPBSA_Error
from GMXMMPBSA.make_top import CheckMakeTop
from GMXMMPBSA.make_top_amber import CheckAmberTop
from GMXMMPBSA.utils import Residue, selector


class CompositeMutationSelectionTest(unittest.TestCase):
    def _checker(self, checker_type, selection, cas_intdiel=0):
        checker = object.__new__(checker_type)
        checker.INPUT = {
            'ala': {
                'mutant_res': selection,
                'mutant': 'ALA',
                'cas_intdiel': cas_intdiel,
            }
        }
        checker.complex_str = SimpleNamespace(residues=[
            SimpleNamespace(name='TYR', chain='A', number=13, insertion_code=''),
            SimpleNamespace(name='VAL', chain='A', number=25, insertion_code=''),
        ])
        refs = [
            Residue(1, 13, 'A', 'R', 1, 'TYR'),
            Residue(2, 25, 'A', 'R', 2, 'VAL'),
        ]
        checker.get_selected_residues = lambda value: refs
        return checker

    def test_selector_accepts_multiple_residues_in_existing_syntax(self):
        self.assertEqual(selector('A/13,25')[1], [['A', 13, ''], ['A', 25, '']])

    def test_selector_requires_direct_insertion_code_suffix_and_excludes_it_from_ranges(self):
        self.assertEqual(selector('A/27B')[1], [['A', 27, 'B']])
        self.assertEqual(
            selector('A/5-7')[1],
            [['A', 5, ''], ['A', 6, ''], ['A', 7, '']],
        )
        with self.assertRaises(MMPBSA_Error):
            selector('A/27:B')

    def test_both_topology_builders_return_composite_indices(self):
        for checker_type in (CheckMakeTop, CheckAmberTop):
            with self.subTest(checker=checker_type.__name__):
                checker = self._checker(checker_type, 'A/13,25')
                self.assertEqual(checker.getMutationInfo(), ([0, 1], 'REC', [0, 1]))

    def test_mixed_receptor_ligand_selection_is_rejected(self):
        for checker_type in (CheckMakeTop, CheckAmberTop):
            with self.subTest(checker=checker_type.__name__):
                checker = self._checker(checker_type, 'A/13 B/25')
                refs = [
                    Residue(1, 13, 'A', 'R', 1, 'TYR'),
                    Residue(2, 25, 'B', 'L', 1, 'VAL'),
                ]
                checker.get_selected_residues = lambda value, refs=refs: refs
                with self.assertRaisesRegex(MMPBSA_Error, 'cannot mix receptor and ligand'):
                    checker.getMutationInfo()

    def test_cas_intdiel_is_rejected_for_composite_selection(self):
        for checker_type in (CheckMakeTop, CheckAmberTop):
            with self.subTest(checker=checker_type.__name__):
                checker = self._checker(checker_type, 'A/13,25', cas_intdiel=1)
                with self.assertRaisesRegex(MMPBSA_Error, 'cas_intdiel=1 is ambiguous'):
                    checker.getMutationInfo()


class ReferenceInsertionCodePropagationTest(unittest.TestCase):
    @staticmethod
    def _structure(residues):
        structure = parmed.Structure()
        for index, (name, chain, number, insertion_code) in enumerate(residues, start=1):
            atom = parmed.Atom(name='CA', atomic_number=6, mass=12.011)
            atom.xx, atom.xy, atom.xz = float(index), 0.0, 0.0
            structure.add_atom(atom, name, number, chain=chain, inscode=insertion_code)
        return structure

    def test_reference_insertion_code_reaches_all_structure_maps(self):
        reference_pdb = (
            'ATOM      1  CA  SER A  27B      1.000   0.000   0.000  1.00  0.00           C\n'
            'ATOM      2  CA  LYS B   1      2.000   0.000   0.000  1.00  0.00           C\n'
            'END\n'
        )
        for builder in (CheckMakeTop, CheckAmberTop):
            with self.subTest(builder=builder.__name__), TemporaryDirectory() as directory:
                reference_path = Path(directory) / 'reference.pdb'
                reference_path.write_text(reference_pdb)
                complex_str = self._structure([
                    ('SER', 'X', 27, ''),
                    ('LYS', 'Y', 1, ''),
                ])
                receptor_str = self._structure([('SER', 'X', 27, '')])
                ligand_str = self._structure([('LYS', 'Y', 1, '')])
                checker = object.__new__(builder)
                checker.FILES = SimpleNamespace(
                    reference_structure=str(reference_path),
                    prefix=str(Path(directory) / 'fixed_'),
                )
                checker.INPUT = {'ala': {'alarun': True}}
                checker.complex_str = complex_str
                checker.resl = [
                    Residue(1, 27, 'X', 'R', 1, 'SER'),
                    Residue(2, 1, 'Y', 'L', 1, 'LYS'),
                ]

                checker.check_structures(complex_str, receptor_str, ligand_str)

                self.assertEqual(complex_str.residues[0].chain, 'A')
                self.assertEqual(complex_str.residues[0].insertion_code, 'B')
                self.assertEqual(receptor_str.residues[0].insertion_code, 'B')
                self.assertEqual(checker.resl[0].chain, 'A')
                self.assertEqual(checker.resl[0].icode, 'B')
                self.assertEqual(checker.get_selected_residues('A/27B'), [checker.resl[0]])

class CompositeTopologyMutationDispatchTest(unittest.TestCase):
    def test_single_and_multiple_indices_are_dispatched_for_both_builders(self):
        for checker_type in (CheckMakeTop, CheckAmberTop):
            with self.subTest(checker=checker_type.__name__):
                checker = object.__new__(checker_type)
                calls = []
                checker.molstr = lambda topology: ('copy', topology)
                checker._makeMutTopSingle = lambda topology, index, pdb=False, copy=True: calls.append(
                    (topology, index, pdb, copy)
                ) or ('mutated', index)

                self.assertEqual(checker.makeMutTop('WT', [3, 7]), ('mutated', 7))
                self.assertEqual(checker.makeMutTop('WT', 3, True), ('mutated', 3))
                self.assertEqual(calls, [
                    (('copy', 'WT'), 3, False, False),
                    (('mutated', 3), 7, False, False),
                    ('WT', 3, True, True),
                ])


class CompositeTrajectoryMutationTest(unittest.TestCase):
    @staticmethod
    def _topology(labels, pointers, natom):
        return SimpleNamespace(
            parm_data={'RESIDUE_LABEL': labels, 'RESIDUE_POINTER': pointers},
            prm_name='topology.prmtop',
            ptr=lambda name: {'natom': natom, 'ifbox': 0}[name],
        )

    @staticmethod
    def _write_trajectory(path, values, frames=2):
        with path.open('w') as handle:
            handle.write('composite test\n')
            for frame in range(frames):
                for index, value in enumerate(values):
                    if index % 10 == 0:
                        handle.write('\n')
                    handle.write('%8.3f' % (value + frame))
            handle.write('\n')

    def _run(self, mutant_type):
        labels = ['TYR', 'VAL']
        pointers = [1, _getnumatms('TYR') + 1]
        original_atoms = _getnumatms('TYR') + _getnumatms('VAL')
        mutant_atoms = _getnumatms('ALA') * 2 if mutant_type is MutantMdcrd else _getnumatms('GLY') * 2
        original = self._topology(labels, pointers, original_atoms)
        target = 'ALA' if mutant_type is MutantMdcrd else 'GLY'
        mutant = self._topology([target, target], [1, _getnumatms(target) + 1], mutant_atoms)
        values = []
        for residue_name in labels:
            for atom in range(_getnumatms(residue_name)):
                values.extend((float(atom + 1), float(atom + 2), float(atom + 3)))

        with TemporaryDirectory() as tmpdir:
            source = Path(tmpdir) / 'normal.mdcrd'
            output = Path(tmpdir) / 'mutant.mdcrd'
            self._write_trajectory(source, values)
            trajectory = mutant_type(str(source), original, mutant)
            trajectory.MutateTraj(str(output))
            result_values = []
            with output.open() as handle:
                next(handle)
                for line in handle:
                    for start in range(0, len(line.rstrip('\n')), 8):
                        field = line.rstrip('\n')[start:start + 8]
                        if field.strip():
                            result_values.append(float(field))

        expected_values = mutant_atoms * 3 * 2
        self.assertEqual(len(result_values), expected_values)
        self.assertEqual(trajectory.mutres, [1, 2])
        self.assertEqual(str(trajectory), f'Y1{target[0]}; V2{target[0]}')

    def test_double_alanine_trajectory_mutation(self):
        self._run(MutantMdcrd)

    def test_double_glycine_trajectory_mutation(self):
        self._run(GlyMutantMdcrd)


class SegmentedMutationTest(unittest.TestCase):
    def test_three_segments_keep_global_indices_and_residue_identity(self):
        import os
        import parmed
        for builder in (CheckMakeTop, CheckAmberTop):
            for component in ('REC', 'LIG'):
                with self.subTest(builder=builder.__name__, component=component), TemporaryDirectory() as directory:
                    structure = parmed.Structure()
                    for residue in range(9):
                        for atom_index, name in enumerate(('N', 'CA', 'C', 'O', 'CB', 'OG')):
                            atom = parmed.Atom(name=name, atomic_number=8 if name.startswith('O') else 6)
                            atom.xx, atom.xy, atom.xz = residue * 4 + atom_index, 1., 2.
                            structure.add_atom(atom, 'SER', residue + 10, chain='A', inscode='B' if residue == 4 else '')
                    checker = object.__new__(builder)
                    checker.FILES = SimpleNamespace(prefix=str(Path(directory) / ''))
                    checker.INPUT = {'general': {'PBRadii': 4}, 'ala': {'alarun': True, 'mutant': 'ALA', 'cas_intdiel': 0}}
                    checker.complex_str = structure
                    checker.receptor_str = structure
                    checker.ligand_str = structure
                    checker.resi = {key: {'num': [[1, 3], [7, 9], [13, 15]]} for key in ('REC', 'LIG')}
                    checker.fixparm2amber = lambda *args: None
                    checker._warn_gmx_gb_radius_compatibility = lambda: None
                    selected = [0, 2, 3, 4, 8]
                    checker.getMutationInfo = lambda: (selected, component, selected)
                    checker.pdb2prmtop()
                    paths = checker.mut_receptor_list if component == 'REC' else checker.mut_ligand_list
                    residues = [r for path in paths.values() for r in parmed.load_file(path).residues]
                    self.assertEqual(len(residues), 9)
                    self.assertEqual([r.name for r in residues], ['ALA' if i in selected else 'SER' for i in range(9)])
                    self.assertEqual([(r.number, r.chain, r.insertion_code) for r in residues],
                                     [(r.number, r.chain, r.insertion_code) for r in structure.residues])
                    self.assertEqual([r.name for r in structure.residues], ['SER'] * 9)


class StreamingMutationTest(CompositeTrajectoryMutationTest):
    def test_streaming_reader_does_not_read_ahead_past_a_frame(self):
        obj = object.__new__(MutantMdcrd)
        obj.traj = 'generated'
        consumed = []
        def lines():
            for i in range(10000):
                consumed.append(i)
                yield ''.join('%8.3f' % value for value in (1, 2, 3)) + '\n'
        frames = obj._iter_frames(lines(), 3)
        self.assertEqual(next(frames), [1, 2, 3])
        self.assertEqual(len(consumed), 1)
        self.assertEqual(sum(1 for _ in frames), 9999)

    def test_empty_truncated_and_wrong_atom_count_remove_partial_output(self):
        from GMXMMPBSA.exceptions import MutateError
        for kind in ('empty', 'truncated', 'count'):
            with self.subTest(kind=kind), TemporaryDirectory() as directory:
                src, dst = Path(directory) / 'in.mdcrd', Path(directory) / 'out.mdcrd'
                orig = self._topology(['VAL'], [1], _getnumatms('VAL'))
                target = self._topology(['ALA'], [1], _getnumatms('ALA') + (kind == 'count'))
                coords = [float(i + 1) for i in range(_getnumatms('VAL') * 3)]
                if kind == 'empty':
                    src.write_text('empty\n')
                else:
                    self._write_trajectory(src, coords, frames=2)
                    if kind == 'truncated':
                        with src.open('a') as out:
                            out.write('%8.3f' % 1.)
                with self.assertRaises(MutateError):
                    MutantMdcrd(str(src), orig, target).MutateTraj(str(dst))
                self.assertFalse(dst.exists())

    def test_box_coordinates_are_preserved(self):
        for mutant_type, name in ((MutantMdcrd, 'ALA'), (GlyMutantMdcrd, 'GLY')):
            with self.subTest(name=name), TemporaryDirectory() as directory:
                src, dst = Path(directory) / 'in.mdcrd', Path(directory) / 'out.mdcrd'
                orig = self._topology(['VAL'], [1], _getnumatms('VAL'))
                target = self._topology([name], [1], _getnumatms(name))
                orig.ptr = lambda key: {'natom': _getnumatms('VAL'), 'ifbox': 1}[key]
                coords = [float(i + 1) for i in range(_getnumatms('VAL') * 3)] + [30., 31., 32.]
                self._write_trajectory(src, coords)
                obj = mutant_type(str(src), orig, target)
                obj.MutateTraj(str(dst))
                with dst.open() as out:
                    next(out)
                    frames = list(obj._iter_frames(out, _getnumatms(name) * 3 + 3))
                self.assertEqual([f[-3:] for f in frames], [[30., 31., 32.], [31., 32., 33.]])


if __name__ == '__main__':
    unittest.main()
