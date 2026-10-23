import math
import os
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import parmed
import test_mutation_charges as charges
import test_mutation_lj as lj
from GMXMMPBSA.make_top import CheckMakeTop
from GMXMMPBSA.make_top_amber import CheckAmberTop
from GMXMMPBSA.mutation import MutationForceField, ensure_lj_type


class MutationForceFieldTest(unittest.TestCase):
    def test_truncated_opls_type_names_still_select_the_opls_family(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            for family in ('oplsaa', 'amber99sb', 'charmm27', 'gromos54a7'):
                directory = root / (family + '.ff')
                directory.mkdir()
                (directory / 'forcefield.itp').touch()
            topology = SimpleNamespace(atoms=[SimpleNamespace(type='opls')])
            self.assertEqual(MutationForceField.candidate_sources(topology, [root]), [root / 'oplsaa.ff'])

    def test_selected_gromacs_reported_data_prefix_handles_wrappers_and_output_streams(self):
        for stream, suffix, relative in [('stdout', '', 'share/gromacs/top'),
                                          ('stderr', '', 'share/gromacs/top'),
                                          ('stdout', ' (source tree)', 'share/top')]:
            with self.subTest(stream=stream, suffix=suffix), TemporaryDirectory() as tmp:
                prefix = Path(tmp) / 'external installation'
                library = prefix / relative
                library.mkdir(parents=True)
                result = SimpleNamespace(stdout='', stderr='')
                setattr(result, stream, f'Executable: /another/bin/gmx\nData prefix: {prefix}{suffix}\n')
                with patch.dict(os.environ, {'GMXLIB': ''}), \
                        patch('GMXMMPBSA.mutation.subprocess.run', return_value=result):
                    roots = MutationForceField.gromacs_library_dirs([str(Path(tmp)/'wrapper/gmx'), 'trjconv'])
                self.assertEqual(roots[0], library.resolve())

    def test_selected_gromacs_symlink_can_locate_data_without_version_output(self):
        for error in (OSError('unavailable'), subprocess.TimeoutExpired('gmx', 5)):
            with self.subTest(error=type(error).__name__), TemporaryDirectory() as tmp:
                installation = Path(tmp) / 'external'
                library = installation / 'share/gromacs/top'
                library.mkdir(parents=True)
                binary = installation / 'bin/gmx'
                binary.parent.mkdir()
                binary.touch()
                link = Path(tmp) / 'links/gmx'
                link.parent.mkdir()
                link.symlink_to(binary)
                with patch.dict(os.environ, {'GMXLIB': ''}), \
                        patch('GMXMMPBSA.mutation.subprocess.run', side_effect=error):
                    roots = MutationForceField.gromacs_library_dirs([str(link), 'trjconv'])
                self.assertEqual(roots[0], library.resolve())

    def test_gmxlib_overrides_the_selected_installation_data_directory(self):
        with TemporaryDirectory() as tmp:
            override = Path(tmp) / 'custom-data'
            installed = Path(tmp) / 'external/share/gromacs/top'
            override.mkdir()
            installed.mkdir(parents=True)
            result = SimpleNamespace(stdout=f'Data prefix: {Path(tmp)/"external"}\n', stderr='')
            with patch.dict(os.environ, {'GMXLIB': str(override)}), \
                    patch('GMXMMPBSA.mutation.subprocess.run', return_value=result):
                roots = MutationForceField.gromacs_library_dirs(['/external/bin/gmx', 'trjconv'])
            self.assertEqual(roots[:2], (override.resolve(), installed.resolve()))

    def test_selected_include_directory_is_restored_after_a_parse_error(self):
        original = parmed.gromacs.GROMACS_TOPDIR
        with self.assertRaisesRegex(ValueError, 'invalid topology'):
            with MutationForceField.gromacs_library(['/external/share/gromacs/top']):
                self.assertEqual(parmed.gromacs.GROMACS_TOPDIR, '/external/share/gromacs/top')
                raise ValueError('invalid topology')
        self.assertEqual(parmed.gromacs.GROMACS_TOPDIR, original)

    @staticmethod
    def _forcefield(directory, reference, target='ALA', gly_type='H1'):
        path = Path(directory) / 'protein.ff'
        path.mkdir()
        template = charges.MutationChargeTest._topology([(target, reference)], gly_type=gly_type)
        entries = {a.type: a for a in template.atoms}
        # Include the parent types as a real force-field directory would.
        parent = charges.MutationChargeTest._topology([('THR', charges.THR)])
        entries.update({a.type:a for a in parent.atoms if a.type not in entries})
        text = '[ defaults ]\n1 2 yes 0.5 0.83333333\n[ atomtypes ]\n'
        for typename, atom in entries.items():
            sigma = atom.rmin * 2 / 2**(1/6) / 10
            text += f'{typename} {atom.atomic_number} {atom.mass} 0 A {sigma:.12g} {atom.epsilon*4.184:.12g}\n'
        (path/'forcefield.itp').write_text(text)
        atoms = '\n'.join(f'{a.name} {a.type} {a.charge:.8f} 1' for a in template.atoms)
        (path/'aminoacids.rtp').write_text(f'[ {target} ]\n[ atoms ]\n{atoms}\n[ bonds ]\n')
        return path

    @staticmethod
    def _attach_source(checker, path):
        """Metadata retained from the active GROMACS topology include tree."""
        checker._mutation_source_files = (path/'forcefield.itp',)
        checker._mutation_gromacs_parameters = parmed.gromacs.GromacsTopologyFile(
            str(path/'forcefield.itp'), parametrize=False).parameterset

    def test_no_native_target_uses_referenced_gromacs_forcefield_templates(self):
        for builder in (CheckMakeTop, CheckAmberTop):
            for target, reference in [('ALA',charges.ALA14),('GLY',charges.GLY03)]:
                with self.subTest(builder=builder.__name__, target=target), TemporaryDirectory() as tmp:
                    path = self._forcefield(tmp, reference, target, gly_type='H0')
                    wt = charges.MutationChargeTest._topology([('THR',charges.THR)])
                    checker = charges.MutationChargeTest._checker(builder,target)
                    self._attach_source(checker, path)
                    with self.assertNoLogs(level='WARNING'):
                        mut = charges.MutationChargeTest._roundtrip(checker.makeMutTop(wt,0))
                    wanted = charges.MutationChargeTest._gly_mutant_charges(reference) if target=='GLY' else reference
                    for atom in mut.residues[0].atoms:
                        self.assertAlmostEqual(atom.charge,wanted[atom.name],places=7)
                    self.assertAlmostEqual(sum(a.charge for a in mut.atoms),0,places=6)
                    if target=='GLY':
                        self.assertTrue(all(a.type=='H0' for a in mut.atoms if a.name.startswith('HA')))

    def test_missing_source_is_estimated_from_charge_patterns_with_a_warning(self):
        for builder in (CheckMakeTop,CheckAmberTop):
            for target,reference in [('ALA',charges.ALA14),('GLY',charges.GLY03)]:
                with self.subTest(builder=builder.__name__,target=target),TemporaryDirectory() as tmp:
                    source=self._forcefield(tmp,reference,target,gly_type='H0')
                    wt=charges.MutationChargeTest._topology([('THR',charges.THR)])
                    checker=charges.MutationChargeTest._checker(builder,target)
                    with patch.object(MutationForceField,'candidate_sources',return_value=[source]), \
                            self.assertLogs(level='WARNING') as logs:
                        mutated=charges.MutationChargeTest._roundtrip(checker.makeMutTop(wt,0))
                    wanted=charges.MutationChargeTest._gly_mutant_charges(reference) if target=='GLY' else reference
                    for atom in mutated.residues[0].atoms:
                        self.assertAlmostEqual(atom.charge,wanted[atom.name],places=7)
                    self.assertAlmostEqual(sum(a.charge for a in mutated.atoms),0,places=6)
                    self.assertTrue(any('Estimated' in line and 'inferred' in line for line in logs.output))

    def test_both_builders_can_infer_from_the_selected_external_gromacs_libraries(self):
        for builder in (CheckMakeTop, CheckAmberTop):
            for target, reference in [('ALA', charges.ALA14), ('GLY', charges.GLY03)]:
                with self.subTest(builder=builder.__name__, target=target), TemporaryDirectory() as tmp:
                    library = Path(tmp) / 'external/share/gromacs/top'
                    library.mkdir(parents=True)
                    field = self._forcefield(library, reference, target, gly_type='H0')
                    wt = charges.MutationChargeTest._topology([('THR', charges.THR)])
                    with (field/'aminoacids.rtp').open('a') as stream:
                        stream.write('[ THR ]\n[ atoms ]\n')
                        for atom in wt.atoms:
                            stream.write(f'{atom.name} {atom.type} {atom.charge} 1\n')
                    binary = Path(tmp)/'wrapper/gmx'
                    binary.parent.mkdir()
                    binary.touch()
                    binary.chmod(0o755)
                    checker = charges.MutationChargeTest._checker(builder, target)
                    if builder is CheckMakeTop:
                        checker.external_progs = {'trjconv': [str(binary), 'trjconv']}
                    else:
                        checker.INPUT['general'] = {'gmx_path': str(binary.parent)}
                    result = SimpleNamespace(stdout=f'Data prefix: {Path(tmp)/"external"}\n', stderr='')
                    with patch.dict(os.environ, {'GMXLIB': ''}), \
                            patch('GMXMMPBSA.mutation.subprocess.run', return_value=result), \
                            self.assertLogs(level='WARNING') as records:
                        mutant = charges.MutationChargeTest._roundtrip(checker.makeMutTop(wt, 0))
                    wanted = charges.MutationChargeTest._gly_mutant_charges(reference) if target=='GLY' else reference
                    for atom in mutant.residues[0].atoms:
                        self.assertAlmostEqual(atom.charge, wanted[atom.name], places=7)
                    self.assertEqual(len(records.output), 1)
                    self.assertIn(str(field), records.output[0])

    def test_topology_includes_use_the_selected_external_gromacs_installation(self):
        with TemporaryDirectory() as tmp:
            library = Path(tmp)/'external/share/gromacs/top'
            library.mkdir(parents=True)
            field = self._forcefield(library, charges.ALA14)
            wt = charges.MutationChargeTest._topology([('THR', charges.THR)])
            topology = Path(tmp)/'input.top'
            atoms = '\n'.join(f'{a.idx+1} {a.type} 1 THR {a.name} 1 {a.charge} {a.mass}' for a in wt.atoms)
            topology.write_text(f'#include "{field.name}/forcefield.itp"\n[ moleculetype ]\nProtein 3\n'
                                f'[ atoms ]\n{atoms}\n[ system ]\nexternal-library test\n'
                                '[ molecules ]\nProtein 1\n')
            checker = charges.MutationChargeTest._checker(CheckMakeTop)
            checker.external_progs = {'trjconv': [str(Path(tmp)/'wrapper/gmx'), 'trjconv']}
            result = SimpleNamespace(stdout=f'Data prefix: {Path(tmp)/"external"}\n', stderr='')
            original = parmed.gromacs.GROMACS_TOPDIR
            with patch.dict(os.environ, {'GMXLIB': ''}), \
                    patch('GMXMMPBSA.mutation.subprocess.run', return_value=result):
                selected, error = checker._cleantop_with_retry(topology, list(range(1, len(wt.atoms)+1)), wt)
            self.assertFalse(error)
            self.assertEqual([a.charge for a in selected.atoms], [a.charge for a in wt.atoms])
            self.assertIn(field/'forcefield.itp', selected._mutation_source_files)
            self.assertEqual(parmed.gromacs.GROMACS_TOPDIR, original)

    def test_known_source_resolves_ambiguous_native_references_before_inference(self):
        for builder in (CheckMakeTop, CheckAmberTop):
            with self.subTest(builder=builder.__name__), TemporaryDirectory() as tmp:
                field = self._forcefield(tmp, charges.ALA14)
                wt = charges.MutationChargeTest._topology([
                    ('THR', charges.THR), ('ALA', charges.ALA14), ('ALA', charges.ALA03)])
                checker = charges.MutationChargeTest._checker(builder)
                self._attach_source(checker, field)
                with patch.object(MutationForceField, 'candidate_sources', side_effect=AssertionError), \
                        self.assertLogs(level='WARNING') as records:
                    mutant = charges.MutationChargeTest._roundtrip(checker.makeMutTop(wt, 0))
                for atom in mutant.residues[0].atoms:
                    self.assertAlmostEqual(atom.charge, charges.ALA14[atom.name], places=7)
                self.assertEqual(len(records.output), 1)
                self.assertIn('discovered force-field source', records.output[0])

    def test_estimation_prefers_the_closer_residue_and_hydrogen_charge_profile(self):
        with TemporaryDirectory() as tmp:
            correct=self._forcefield(tmp,charges.ALA14)
            wrong=Path(tmp)/'wrong';wrong.mkdir()
            incorrect=self._forcefield(wrong,charges.ALA03)
            wt=charges.MutationChargeTest._topology([('THR',charges.THR)])
            for directory,offset in [(correct,0),(incorrect,.05)]:
                with (directory/'aminoacids.rtp').open('a') as stream:
                    stream.write('[ THR ]\n[ atoms ]\n')
                    for atom in wt.atoms:
                        stream.write(f'{atom.name} {atom.type} {atom.charge+offset} 1\n')
            names=frozenset(charges.ALA14)
            with patch.object(MutationForceField,'candidate_sources',return_value=[incorrect,correct]):
                library,fit=MutationForceField.infer(wt,'ALA',names)
            self.assertEqual(fit['source'],str(correct))
            self.assertAlmostEqual(fit['hydrogen_rms'],0)
            self.assertEqual({a.name:a.charge for a in library.topology('ALA').residues[0].atoms},charges.ALA14)

    def test_composite_source_and_estimated_matching_use_unmodified_parent_charges(self):
        for builder in (CheckMakeTop,CheckAmberTop):
            for inferred in (False,True):
                with self.subTest(builder=builder.__name__,inferred=inferred),TemporaryDirectory() as tmp:
                    field=self._forcefield(tmp,charges.ALA14)
                    wt=charges.MutationChargeTest._topology([('THR',charges.THR),('THR',charges.THR)])
                    with (field/'aminoacids.rtp').open('a') as stream:
                        stream.write('[ THR ]\n[ atoms ]\n')
                        for atom in wt.residues[0].atoms:
                            stream.write(f'{atom.name} {atom.type} {atom.charge} 1\n')
                    checker=charges.MutationChargeTest._checker(builder)
                    if not inferred:
                        self._attach_source(checker, field)
                    with patch.object(MutationForceField,'candidate_sources',return_value=[field]):
                        if inferred:
                            with self.assertLogs(level='WARNING') as records:
                                mutant=checker.makeMutTop(wt,[0,1])
                            self.assertEqual(sum('inferred assignment' in line for line in records.output),1)
                        else:
                            mutant=checker.makeMutTop(wt,[0,1])
                    back=charges.MutationChargeTest._roundtrip(mutant)
                    for residue in back.residues:
                        for atom in residue.atoms:
                            self.assertAlmostEqual(atom.charge,charges.ALA14[atom.name],places=7)
                    self.assertTrue(all(r.name=='THR' for r in wt.residues))

    def test_referenced_rtp_is_discovered_without_explicit_option(self):
        with TemporaryDirectory() as tmp:
            path = self._forcefield(tmp,charges.ALA14)
            lib = MutationForceField.load(includes=[path/'forcefield.itp'],
                                          gromacs_parameters=parmed.gromacs.GromacsTopologyFile(
                                              str(path/'forcefield.itp'),parametrize=False).parameterset)
            self.assertEqual(len(lib.topology('ALA').residues),1)
            self.assertEqual(lib.topology('ALA').residues[0].atoms[0].charge,charges.ALA14['N'])

    def test_charmm_gui_sibling_toppar_is_discovered_automatically(self):
        with TemporaryDirectory() as tmp:
            ff=self._forcefield(tmp,charges.ALA14)
            top=Path(tmp)/'bundle/gromacs/topol.top';top.parent.mkdir(parents=True);top.write_text('; input\n')
            toppar=top.parent.parent/'toppar';toppar.mkdir()
            native=charges.MutationChargeTest._topology([('ALA',charges.ALA14)])
            types={a.type:a for a in native.atoms}
            text='* mutation test\n*\n36 1\n'
            text+=''.join(f'MASS -1 {name} {a.mass}\n' for name,a in types.items())
            text+='RESI ALA 0.0\nGROUP\n'
            text+=''.join(f'ATOM {a.name} {a.type} {a.charge}\n' for a in native.atoms)
            text+='END\n';(toppar/'top_all36_prot.rtf').write_text(text)
            gmx=parmed.gromacs.GromacsTopologyFile(str(ff/'forcefield.itp'),parametrize=False)
            library=MutationForceField.load(topology_file=top,gromacs_parameters=gmx.parameterset)
            residue=library.topology('ALA').residues[0]
            self.assertEqual({a.name:a.charge for a in residue.atoms},charges.ALA14)

    def test_missing_reference_and_no_source_still_reports_actionable_error(self):
        wt = charges.MutationChargeTest._topology([('THR',charges.THR)])
        checker = charges.MutationChargeTest._checker(CheckMakeTop)
        with patch.object(MutationForceField, 'candidate_sources', return_value=[]), \
                self.assertLogs(level='ERROR'), self.assertRaisesRegex(Exception,'force-field libraries'):
            checker.makeMutTop(wt,0)

    def test_a_forcefield_with_conflicting_parent_charges_is_rejected(self):
        wt=charges.MutationChargeTest._topology([('THR',charges.THR)])
        library=MutationForceField()
        library.records.append(('THR', {name:('CT', charge + (0.1 if name=='CA' else 0))
                                        for name,charge in charges.THR.items()}, 'wrong.ff'))
        with self.assertRaisesRegex(ValueError, 'does not match.*charges'):
            library.validate_topology(wt)

    def test_small_source_charge_difference_can_continue_with_a_warning(self):
        with TemporaryDirectory() as tmp:
            field=self._forcefield(tmp,charges.ALA14)
            wt=charges.MutationChargeTest._topology([('THR',charges.THR)])
            with (field/'aminoacids.rtp').open('a') as stream:
                stream.write('[ THR ]\n[ atoms ]\n')
                for atom in wt.atoms:
                    stream.write(f'{atom.name} {atom.type} {atom.charge+(0.002 if atom.name=="CA" else 0)} 1\n')
            checker=charges.MutationChargeTest._checker(CheckMakeTop)
            self._attach_source(checker, field)
            with self.assertLogs(level='WARNING') as records:
                mutated=charges.MutationChargeTest._roundtrip(checker.makeMutTop(wt,0))
            self.assertTrue(any('Continuing with estimated' in line for line in records.output))
            for atom in mutated.residues[0].atoms:
                self.assertAlmostEqual(atom.charge,charges.ALA14[atom.name],places=7)

    def test_nested_amber_source_cycle_is_rejected(self):
        with TemporaryDirectory() as tmp:
            a=Path(tmp)/'leaprc.a';b=Path(tmp)/'leaprc.b'
            a.write_text('source leaprc.b\n');b.write_text('source leaprc.a\n')
            with self.assertRaisesRegex(ValueError, 'Cyclic'):
                MutationForceField.load(str(a))

    def test_tdb_terminal_changes_are_read_from_source(self):
        with TemporaryDirectory() as tmp:
            path = self._forcefield(tmp,charges.ALA14)
            (path/'aminoacids.c.tdb').write_text('[ COO- ]\n[ replace ]\nC C CT 12.01 0.7731\n'
                'O O OH 16.00 -0.8055\nN N N 14.01 -0.3821\nH H H1 1.008 0.2681\n'
                'CA CA CT 12.01 -0.1747\nHA HA H1 1.008 0.1067\nCB CB CT 12.01 -0.2093\n'
                'HB1 HB1 HC 1.008 0.0764\nHB2 HB2 HC 1.008 0.0764\nHB3 HB3 HC 1.008 0.0764\n'
                '[ add ]\n1 1 OXT C CA N\nOH 16.00 -0.8055\n')
            lib=MutationForceField.load(str(path))
            terminal=next(r for r in lib.topology('ALA').residues if r.name=='CALA')
            self.assertEqual({a.name:a.charge for a in terminal.atoms},charges.CALA)

    def test_tdb_terminal_oxygen_aliases_and_mass_overrides_are_preserved(self):
        with TemporaryDirectory() as tmp:
            path=self._forcefield(tmp,charges.ALA14)
            rows=[]
            for name,charge in charges.CALA.items():
                renamed={'O':'O1','OXT':'O2'}.get(name,name)
                typename='N' if name=='N' else 'H1' if name.startswith('H') else 'OH' if name.startswith('O') else 'CT'
                mass=14.0027 if name=='N' else 1.008 if name.startswith('H') else 16.0 if name.startswith('O') else 12.01
                rows.append(f'{name} {renamed} {typename} {mass} {charge}')
            (path/'aminoacids.c.tdb').write_text('[ COO- ]\n[ replace ]\n'+'\n'.join(rows)+
                '\n[ add ]\n2 8 O C CA N\nOH 16.0 -0.8055\n')
            library=MutationForceField.load(str(path))
            terminal=next(r for r in library.topology('ALA').residues if r.name=='CALA')
            expected={{'O':'O1','OXT':'O2'}.get(name,name):q for name,q in charges.CALA.items()}
            self.assertEqual({a.name:a.charge for a in terminal.atoms},expected)
            self.assertAlmostEqual(next(a.mass for a in terminal.atoms if a.name=='N'),14.0027)

    def test_rtp_uses_atomtypes_atp_masses(self):
        with TemporaryDirectory() as tmp:
            path=self._forcefield(tmp,charges.ALA14)
            (path/'atomtypes.atp').write_text('N 14.007\nOH 15.999\n')
            residue=MutationForceField.load(str(path)).topology('ALA').residues[0]
            self.assertAlmostEqual(next(a.mass for a in residue.atoms if a.name=='N'),14.007)
            self.assertAlmostEqual(next(a.mass for a in residue.atoms if a.name=='O'),15.999)

    def test_new_lj_type_keeps_old_pairs_and_source_nbfix_and_14_pairs(self):
        wt=lj.MutationLennardJonesTest._topology([
            ('SER',[('CA','CT'),('HA','H1'),('CB','CT'),('OG','OH')])],chamber=True)
        new=parmed.AtomType('NEW',None,1.008,1)
        new.set_lj_params(eps=.024,rmin=1.34,eps14=.012,rmin14=1.24)
        new.add_nbfix('OH',3.45,.078,3.12,.045)
        indices={name:{i} for name,i in wt.LJ_types.items()}
        flags=['LENNARD_JONES_ACOEF','LENNARD_JONES_BCOEF','LENNARD_JONES_14_ACOEF','LENNARD_JONES_14_BCOEF']
        old={flag:list(wt.parm_data[flag]) for flag in flags}
        idx=ensure_lj_type(wt,dict(type='NEW',atom_type=new),wt.atoms[2],indices)
        self.assertEqual(idx,wt.ptr('ntypes'))
        for flag in flags:
            self.assertEqual(wt.parm_data[flag][:len(old[flag])],old[flag])
        other=next(iter(indices['OH']))
        pair=wt.parm_data['NONBONDED_PARM_INDEX'][wt.ptr('ntypes')*(idx-1)+other-1]-1
        expected=[.078*3.45**12,2*.078*3.45**6,.045*3.12**12,2*.045*3.12**6]
        for flag,value in zip(flags,expected):
            self.assertAlmostEqual(wt.parm_data[flag][pair],value,places=7)

    def test_geometric_combining_rule_is_used_when_adding_a_type(self):
        wt=lj.MutationLennardJonesTest._topology([('LIG',[('C','CT'),('H','H1')])])
        new=parmed.AtomType('NEW',None,1.008,1);new.set_lj_params(eps=.024,rmin=1.34)
        indices={name:{i} for name,i in wt.LJ_types.items()}
        carbon=next(iter(indices['CT']))
        idx=ensure_lj_type(wt,dict(type='NEW',atom_type=new),wt.atoms[1],indices,combining_rule='geometric')
        pair=wt.parm_data['NONBONDED_PARM_INDEX'][wt.ptr('ntypes')*(idx-1)+carbon-1]-1
        radius=2*math.sqrt(1.34*wt.LJ_radius[carbon-1]);epsilon=math.sqrt(.024*wt.LJ_depth[carbon-1])
        self.assertAlmostEqual(wt.parm_data['LENNARD_JONES_ACOEF'][pair],epsilon*radius**12,places=7)

    def test_new_lj_type_preserves_self_nbfix_in_normal_and_14_tables(self):
        wt=lj.MutationLennardJonesTest._topology([('LIG',[('C','CT'),('H','H1')])],chamber=True)
        new=parmed.AtomType('NEW',None,1.008,1)
        new.set_lj_params(eps=.024,rmin=1.34,eps14=.012,rmin14=1.24)
        new.add_nbfix('NEW',2.8,.033,2.6,.022)
        indices={name:{i} for name,i in wt.LJ_types.items()}
        idx=ensure_lj_type(wt,dict(type='NEW',atom_type=new),wt.atoms[1],indices)
        pair=wt.parm_data['NONBONDED_PARM_INDEX'][wt.ptr('ntypes')*(idx-1)+idx-1]-1
        expected=[.033*2.8**12,2*.033*2.8**6,.022*2.6**12,2*.022*2.6**6]
        flags=['LENNARD_JONES_ACOEF','LENNARD_JONES_BCOEF','LENNARD_JONES_14_ACOEF','LENNARD_JONES_14_BCOEF']
        for flag,value in zip(flags,expected):
            self.assertAlmostEqual(wt.parm_data[flag][pair],value,places=7)

    def test_new_lj_type_preserves_cached_native_self_pair(self):
        wt=lj.MutationLennardJonesTest._topology([('LIG',[('C','CT'),('H','H1')])],chamber=True)
        new=parmed.AtomType('NEW',None,1.008,1)
        new.set_lj_params(eps=.024,rmin=1.34,eps14=.012,rmin14=1.24)
        indices={name:{i} for name,i in wt.LJ_types.items()}
        expected=(.033*2.8**12,2*.033*2.8**6,.022*2.6**12,2*.022*2.6**6)
        idx=ensure_lj_type(wt,dict(type='NEW',atom_type=new,pairs={'NEW':expected}),wt.atoms[1],indices)
        pair=wt.parm_data['NONBONDED_PARM_INDEX'][wt.ptr('ntypes')*(idx-1)+idx-1]-1
        flags=['LENNARD_JONES_ACOEF','LENNARD_JONES_BCOEF','LENNARD_JONES_14_ACOEF','LENNARD_JONES_14_BCOEF']
        for flag,value in zip(flags,expected):
            self.assertAlmostEqual(wt.parm_data[flag][pair],value,places=7)

    def test_explicit_gromacs_pairtypes_override_combined_14_parameters(self):
        wt=lj.MutationLennardJonesTest._topology([('LIG',[('C','CT'),('H','H1')])],chamber=True)
        new=parmed.AtomType('NEW',None,1.008,1);new.set_lj_params(eps=.024,rmin=1.34)
        indices={name:{i} for name,i in wt.LJ_types.items()}
        pairtype=parmed.NonbondedExceptionType(3.456,.0789)
        idx=ensure_lj_type(wt,dict(type='NEW',atom_type=new),wt.atoms[1],indices,
                           pair_types={('NEW','CT'):pairtype})
        pair=wt.parm_data['NONBONDED_PARM_INDEX'][wt.ptr('ntypes')*(idx-1)+next(iter(indices['CT']))-1]-1
        self.assertAlmostEqual(wt.parm_data['LENNARD_JONES_14_ACOEF'][pair],.0789*3.456**12,places=7)

    def test_missing_1264_pair_data_is_rejected_before_lj_table_changes(self):
        wt=lj.MutationLennardJonesTest._topology([('LIG',[('C','CT'),('H','H1')])])
        wt.add_flag('LENNARD_JONES_CCOEF','5E16.8',data=[1.0]*len(wt.parm_data['LENNARD_JONES_ACOEF']))
        new=parmed.AtomType('NEW',None,1.008,1);new.set_lj_params(eps=.024,rmin=1.34)
        indices={name:{i} for name,i in wt.LJ_types.items()}
        before=wt.ptr('ntypes')
        with self.assertRaisesRegex(ValueError,'C-coefficients'):
            ensure_lj_type(wt,dict(type='NEW',atom_type=new),wt.atoms[1],indices)
        self.assertEqual(wt.ptr('ntypes'),before)


if __name__=='__main__':
    unittest.main()
