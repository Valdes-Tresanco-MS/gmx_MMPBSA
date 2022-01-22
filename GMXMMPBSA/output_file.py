"""
This module is responsible for creating the final output file for MM/PBSA
statistics printing.
"""

# ##############################################################################
#                           GPLv3 LICENSE INFO                                 #
#                                                                              #
#  Copyright (C) 2020  Mario S. Valdes-Tresanco and Mario E. Valdes-Tresanco   #
#  Copyright (C) 2014  Jason Swails, Bill Miller III, and Dwight McGee         #
#                                                                              #
#   Project: https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA                  #
#                                                                              #
#   This program is free software; you can redistribute it and/or modify it    #
#  under the terms of the GNU General Public License version 3 as published    #
#  by the Free Software Foundation.                                            #
#                                                                              #
#  This program is distributed in the hope that it will be useful, but         #
#  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY  #
#  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License    #
#  for more details.                                                           #
# ##############################################################################

import io
import numpy as np
from GMXMMPBSA.amber_outputs import EnergyVector
from GMXMMPBSA.exceptions import LengthError
from GMXMMPBSA import utils
from math import sqrt, ceil
from os import linesep as ls
import h5py

h5py.get_config().track_order = True

class Data2h5:
    def __init__(self, app):
        self.app = app
        self.h5f = h5py.File('RESULTS_gmx_MMPBSA.h5', 'w')
        grp = self.h5f.create_group('normal')
        self._e2h5(app.calc_types.normal, grp)
        if app.calc_types.mutant:
            grp = self.h5f.create_group('mutant')
            self._e2h5(app.calc_types.mutant, grp)
        if app.calc_types.decomp_normal:
            grp = self.h5f.create_group('decomp_normal')
            self._decomp2h5(app.calc_types.decomp_normal, grp)
        if app.calc_types.decomp_mutant:
            grp = self.h5f.create_group('decomp_mutant')
            self._decomp2h5(app.calc_types.decomp_mutant, grp)
        self._info2h5()
        self.h5f.close()

    def _info2h5(self):
        grp = self.h5f.create_group('INPUT')
        for x in self.app.INPUT:
            data =  np.nan if self.app.INPUT[x] is None else self.app.INPUT[x]
            dset = grp.create_dataset(x, data=data)

        grp = self.h5f.create_group('FILES')
        for x in dir(self.app.FILES):
            # this must be equal to the info saved in the info file
            if x.startswith('_') or x in ('rewrite_output', 'energyout', 'dec_energies', 'overwrite'):
                continue
            d = getattr(self.app.FILES, x)
            data = np.nan if d is None else d
            dset = grp.create_dataset(x, data=data)

        grp = self.h5f.create_group('INFO')
        dset = grp.create_dataset('size', data=self.app.mpi_size)
        dset = grp.create_dataset('numframes', data=self.app.numframes)
        dset = grp.create_dataset('numframes_nmode', data=self.app.numframes_nmode)
        dset = grp.create_dataset('mut_str', data=self.app.mut_str)
        dset = grp.create_dataset('using_chamber', data=self.app.using_chamber)
        dset = grp.create_dataset('input_file', data=self.app.input_file_text)

        # save the complex fixed structure
        com_fixed = ''.join(open(self.app.FILES.complex_fixed).readlines())
        dset = grp.create_dataset('COM_PDB', data=com_fixed)

        # get output files
        outfile = ''.join(open(self.app.FILES.output_file).readlines())
        dset = grp.create_dataset('output_file', data=outfile)
        if self.app.INPUT['decomprun']:
            doutfile = ''.join(open(self.app.FILES.decompout).readlines())
            dset = grp.create_dataset('decomp_output_file', data=doutfile)

    @staticmethod
    def _e2h5(d, f):
        # key  Energy: [gb, pb, rism std, rism gf], Decomp: [gb, pb], Entropy: [nmode, qh, ie, c2]
        for key in d:
            if key in ['gb', 'pb', 'rism std', 'rism gf']:
                grp = f.create_group(key)
                # key2 is complex, receptor, ligand, delta
                for key2 in d[key]:
                    grp2 = grp.create_group(key2)
                    # complex, receptor, etc., is a class and the data is contained in the attribute data
                    for key3 in d[key][key2]:
                        dset = grp2.create_dataset(key3, data=d[key][key2][key3])
            elif key in ['nmode', 'qh']:
                grp = f.create_group(key)
                # key2 is complex, receptor, ligand, delta
                for key2 in d[key]:
                    grp2 = grp.create_group(key2)
                    # vibrational, translational, rotational, total
                    for key3 in d[key][key2]:
                        dset = grp2.create_dataset(key3, data=d[key][key2][key3])
            elif key in ['ie', 'c2']:
                grp = f.create_group(key)
                # key2 is PB, GB or RISM?
                for key2 in d[key]:
                    grp2 = grp.create_group(key2)
                    for key3 in d[key][key2]:
                        dset = grp2.create_dataset(key3, data=d[key][key2][key3])

    @staticmethod
    def _decomp2h5(d, g):
        for key in d:
            # model
            grp = g.create_group(key)
            # key2 is complex, receptor, ligand, delta
            for key2 in d[key]:
                grp2 = grp.create_group(key2)
                # TDC, SDC, BDC
                for key3 in d[key][key2]:
                    grp3 = grp2.create_group(key3)
                    # residue first level
                    for key4 in d[key][key2][key3]:
                        grp4 = grp3.create_group(key4)
                        # residue sec level
                        for key5 in d[key][key2][key3][key4]:
                            if isinstance(d[key][key2][key3][key4][key5], dict):
                                grp5 = grp4.create_group(key5)
                                # energy terms
                                for key6 in d[key][key2][key3][key4][key5]:
                                    dset = grp5.create_dataset(key6, data=d[key][key2][key3][key4][key5][key6])
                            else:
                                dset = grp4.create_dataset(key5, data=d[key][key2][key3][key4][key5])

# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

def _get_diff(arr2, arr1):
    if isinstance(arr2, EnergyVector) and isinstance(arr1, EnergyVector):
        dmean = arr2.mean() - arr1.mean()
        dstd = utils.get_std(arr2.stdev(), arr1.stdev())
    elif isinstance(arr2, EnergyVector) and isinstance(arr1, (int, float)):
        dmean = arr2.mean() - arr1
        dstd = utils.get_std(arr2.stdev(), 0)
    elif isinstance(arr2, (int, float)) and isinstance(arr1, EnergyVector):
        dmean = arr2 - arr1.mean()
        dstd = utils.get_std(0, arr1.stdev())
    else:
        dmean = arr2 - arr1
        dstd = utils.get_std(0, 0)
    return dmean, dstd


def _get_sum(arr2: EnergyVector, arr1: EnergyVector):
    if isinstance(arr2, EnergyVector) and isinstance(arr1, EnergyVector):
        dmean = arr2.mean() + arr1.mean()
        dstd = utils.get_std(arr2.stdev(), arr1.stdev())
    elif isinstance(arr2, EnergyVector) and isinstance(arr1, (int, float)):
        dmean = arr2.mean() + arr1
        dstd = utils.get_std(arr2.stdev(), 0)
    elif isinstance(arr2, (int, float)) and isinstance(arr1, EnergyVector):
        dmean = arr2 + arr1.mean()
        dstd = utils.get_std(0, arr1.stdev())
    else:
        dmean = arr2 + arr1
        dstd = utils.get_std(0, 0)
    return dmean, dstd


# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

def write_stability_output(app):
    """ Writes output files based on stability calculations """
    import csv
    # Load some objects into top-level name space
    FILES = app.FILES
    INPUT = app.INPUT
    mut_str = app.mut_str

    if FILES.energyout:
        ene_csv = open(FILES.energyout, 'w')
        energyvectors = csv.writer(ene_csv, dialect='excel')

    final_output = OutputFile(FILES.output_file, 'w')
    final_output.write_date()
    final_output.add_comment('')
    final_output.write(app.input_file_text)
    final_output.print_file_info(FILES, INPUT)
    final_output.add_comment('')
    final_output.add_comment('Calculations performed using %s complex frames.' % app.numframes)
    if INPUT['nmoderun']:
        final_output.add_comment('NMODE calculations performed using %s frames.' % app.numframes_nmode)

    if INPUT['pbrun']:
        if INPUT['sander_apbs']:
            final_output.add_comment('Poisson Boltzmann calculations performed ' +
                                     'using iAPBS interface to sander (sander.APBS)')
        else:
            final_output.add_comment('Poisson Boltzmann calculations performed ' +
                                     'using internal PBSA solver in sander.')

    final_output.add_comment('')

    if INPUT['gbrun']:
        if INPUT['molsurf']:
            final_output.add_comment('Generalized Born ESURF calculated using \'molsurf\' surface areas')
        else:
            final_output.add_comment('Generalized Born ESURF calculated using \'LCPO\' surface areas')
        final_output.add_comment('')

    final_output.add_comment('Using temperature = %.2f K)\n' % INPUT['temperature'])
    final_output.add_comment('All units are reported in kcal/mol.')
    # if INPUT['nmoderun'] or INPUT['qh_entropy'] or INPUT['interaction_entropy' or INPUT['c2_entropy']]:
    #     final_output.add_comment('All entropy results have units kcal/mol')
    #     final_output.add_comment('(Temperature for NMODE and QH is %.2f K)\n' % INPUT['temperature'])
    #     final_output.add_comment('(Temperature for IE and C2 entropy is %.2f K)\n' % INPUT['temperature'])
    if INPUT['ifqnt']:
        final_output.add_comment(('QM/MM: Residues %s are treated with the ' +
                                  'Quantum Hamiltonian %s') % (INPUT['qm_residues'], INPUT['qm_theory']))
    final_output.separate()

    # Start with entropies
    if INPUT['qh_entropy']:
        if not INPUT['mutant_only']:
            qhnorm = app.calc_types['qh']
            final_output.writeline('ENTROPY RESULTS (QUASI-HARMONIC APPROXIMATION) CALCULATED WITH PTRAJ:')
            final_output.add_section(qhnorm.print_summary())
        if INPUT['alarun']:
            qhmutant = app.calc_types.mutant['qh']
            final_output.writeline(mut_str + ' MUTANT')
            final_output.writeline('ENTROPY RESULTS (QUASI-HARMONIC APPROXIMATION) CALCULATED WITH PTRAJ:')
            final_output.add_section(qhmutant.print_summary())
        if INPUT['alarun'] and not INPUT['mutant_only']:
            davg, dstd = _get_diff(qhmutant.total_avg(), qhnorm.total_avg())
            final_output.add_section(f'\nRESULT OF ALANINE SCANNING ({mut_str}):\n'
                                     f'ΔΔS = {davg:9.4f} +/ {dstd:7.4f}\n')
    # end if INPUT['entropy']

    # Now print out the normal mode results
    if INPUT['nmoderun']:
        if not INPUT['mutant_only']:
            nm_norm = app.calc_types['nmode']['complex']
            final_output.write('ENTROPY RESULTS (HARMONIC APPROXIMATION) CALCULATED WITH NMODE:\n')
            final_output.add_section('Complex:\n' + nm_norm.print_summary())
            # Now dump the energy vectors in CSV format
            if FILES.energyout:
                energyvectors.writerow(['NMODE entropy results'])
                nm_norm.print_vectors(energyvectors)
                energyvectors.writerow([])

        if INPUT['alarun']:
            nm_mut = app.calc_types.mutant['nmode']['complex']
            final_output.write(mut_str + ' MUTANT\nENTROPY RESULTS (HARMONIC APPROXIMATION) CALCULATED WITH NMODE:\n')
            final_output.add_section('Complex:\n' + nm_mut.print_summary())
            # Now dump the energy vectors in CSV format
            if FILES.energyout:
                energyvectors.writerow([mut_str + ' Mutant NMODE entropy results'])
                nm_mut.print_vectors(energyvectors)
                energyvectors.writerow([])

        # Now calculate the effect of alanine scanning
        if INPUT['alarun'] and not INPUT['mutant_only']:
            davg, dstd = _get_diff(nm_mut.data['Total'], nm_norm.data['Total'])
            final_output.add_section(f'\nRESULT OF ALANINE SCANNING ({mut_str}):\n'
                                     f'ΔΔS = {davg:9.4f} +/- {dstd:7.4f}\n')

    # end if INPUT['nmoderun']

    # Now print out the Free Energy results
    triggers = ('gbrun', 'pbrun', 'rismrun_std', 'rismrun_gf')
    outkeys = ('gb', 'pb', 'rism std', 'rism gf')
    headers = ('\nGENERALIZED BORN:\n\n', '\nPOISSON BOLTZMANN:\n\n',
               '\n3D-RISM:\n\n', '\n3D-RISM (Gauss. Fluct.):\n\n')
    for i, key in enumerate(outkeys):
        if not INPUT[triggers[i]]:
            continue
        if not INPUT['mutant_only']:
            com_norm = app.calc_types[key]['complex']
            final_output.write(headers[i])
            final_output.add_section('Complex:\n' + com_norm.print_summary())
            # Dump energy vectors to a CSV
            if FILES.energyout:
                energyvectors.writerow([headers[i].strip()])
                com_norm.print_vectors(energyvectors)
                energyvectors.writerow([])

            # Combine with the entropy(ies)
            if INPUT['qh_entropy']:
                qh_davg, qh_dstd = _get_diff(com_norm.data['TOTAL'], qhnorm.total_avg())
                final_output.add_section('Using Quasi-harmonic Entropy Approximation: '
                                         f'ΔG = {qh_davg:9.4f} +/- {qh_dstd:7.4f}\n')
            if INPUT['nmoderun']:
                nm_davg, nm_dstd = _get_diff(com_norm.data['TOTAL'], nm_norm.data['Total'])
                final_output.add_section('Using Normal Mode Entropy Approximation: '
                                         f'ΔG = {nm_davg:9.4f} +/- {nm_dstd:7.4f}\n')

        if INPUT['alarun']:
            com_mut = app.calc_types.mutant[key]['complex']
            final_output.write('%s MUTANT:%s' % (mut_str, headers[i]))
            final_output.add_section('Complex:\n' + com_mut.print_summary())
            # Dump energy vectors to a CSV
            if FILES.energyout:
                energyvectors.writerow([mut_str + ' Mutant ' + headers[i]])
                com_mut.print_vectors(energyvectors)
                energyvectors.writerow([])

            # Combine with the entropy(ies)
            if INPUT['qh_entropy']:
                mqh_davg, mqh_dstd = _get_diff(com_mut.data['TOTAL'], qhmutant.total_avg())
                final_output.add_section('Using Quasi-harmonic Entropy Approximation: '
                                         f'ΔG = {mqh_davg:9.4f} +/- {mqh_dstd:7.4f}\n')
            if INPUT['nmoderun']:
                mnm_davg, mnm_dstd = _get_diff(com_mut.data['TOTAL'], nm_norm.data['Total'])
                final_output.add_section('Using Normal Mode Entropy Approximation: '
                                         f'ΔG = {mnm_davg:9.4f} +/- {mnm_dstd:7.4f}\n')

        if INPUT['alarun'] and not INPUT['mutant_only']:
            davg, dstd = _get_diff(com_mut.data['TOTAL'], com_norm.data['TOTAL'])
            final_output.write(f'\nRESULT OF ALANINE SCANNING ({mut_str}):\n'
                               f'ΔΔH = {davg:9.4f} +/- {dstd:7.4f}\n')

            if INPUT['qh_entropy']:
                davg, dstd = _get_diff(qhmutant.total_avg(), qhnorm.total_avg())
                final_output.write('\n   (quasi-harmonic entropy)\n'
                                   f'({mut_str}) ΔΔG = {davg:9.4f} +/- {dstd:7.4f}\n')

            if INPUT['nmoderun']:
                davg, _ = _get_diff(mnm_davg, nm_davg)
                dstd = utils.get_std(mnm_dstd, nm_dstd)
                final_output.write('\n   (normal mode entropy)\n'
                                   f'({mut_str}) ΔΔG = {davg:9.4f} +/- {dstd:7.4f}\n')
            final_output.separate()

    # end for solv in ['gbrun', 'pbrun', ...]

    if FILES.energyout: ene_csv.close()


# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

def write_binding_output(app):
    """ Writes a binding output file """
    import csv

    FILES = app.FILES
    INPUT = app.INPUT
    mut_str = app.mut_str
    prmtop_system = app.normal_system

    # Open the energy vector CSV output file if we are writing one
    if FILES.energyout:
        ene_csv = open(FILES.energyout, 'w')
        energyvectors = csv.writer(ene_csv, dialect='excel')

    # Open the file and write some initial data to it
    final_output = OutputFile(FILES.output_file, 'w')
    final_output.write_date()
    final_output.add_comment('')
    final_output.write(app.input_file_text)
    final_output.print_file_info(FILES, INPUT)
    final_output.add_comment('')
    final_output.add_comment('Receptor mask:                  "%s"' % INPUT['receptor_mask'])
    final_output.add_comment('Ligand mask:                    "%s"' % INPUT['ligand_mask'])
    if prmtop_system.ligand_prmtop.ptr('nres') == 1:
        final_output.add_comment('Ligand residue name is "%s"' %
                                 prmtop_system.ligand_prmtop.parm_data['RESIDUE_LABEL'][0])
    final_output.add_comment('')
    final_output.add_comment('Calculations performed using %s complex frames.' % app.numframes)
    if INPUT['nmoderun']:
        final_output.add_comment('NMODE calculations performed using %s frames.' % app.numframes_nmode)
    if INPUT['interaction_entropy']:
        final_output.add_comment('Interaction Entropy calculations performed using last %s frames.' %
                                 ceil(app.numframes * (INPUT['ie_segment'] / 100)))
    if INPUT['c2_entropy']:
        final_output.add_comment('C2 Entropy calculations performed using last %s frames.' %
                                 ceil(app.numframes * (INPUT['c2_segment'] / 100)))
        final_output.add_comment('C2 Entropy Std. Dev. and Conf. Interv. (95%) have been obtained by '
                                 'bootstrapping with number_of_resamplings = 2000')


    if INPUT['pbrun']:
        if INPUT['sander_apbs']:
            final_output.add_comment('Poisson Boltzmann calculations performed ' +
                                     'using iAPBS interface to sander (sander.APBS)')
        else:
            final_output.add_comment('Poisson Boltzmann calculations performed ' +
                                     'using internal PBSA solver in sander.')
    final_output.add_comment('')

    if INPUT['gbrun']:
        if INPUT['molsurf']:
            final_output.add_comment('Generalized Born ESURF calculated using \'molsurf\' surface areas')
        else:
            final_output.add_comment('Generalized Born ESURF calculated using \'LCPO\' surface areas')
        final_output.add_comment('')

    final_output.add_comment('Using temperature = %.2f K)\n' % INPUT['temperature'])
    final_output.add_comment('All units are reported in kcal/mol.')
    # if INPUT['nmoderun'] or INPUT['qh_entropy'] or INPUT['interaction_entropy' or INPUT['c2_entropy']]:
    #     final_output.add_comment('All entropy results have units kcal/mol\n' +
    #                              '(Temperature for NMODE and QH is %.2f K)\n' % INPUT['temperature'] +
    #                              '(Temperature for IE and C2 is %.2f K)\n' % INPUT['temperature'])
    if INPUT['ifqnt']:
        final_output.add_comment(('QM/MM: Residues %s are treated with the ' +
                                  'Quantum Hamiltonian %s') % (INPUT['qm_residues'], INPUT['qm_theory']))
    final_output.separate()

    # First do the entropies
    if INPUT['qh_entropy']:
        if not INPUT['mutant_only']:
            qhnorm = app.calc_types['qh']
            final_output.writeline('ENTROPY RESULTS (QUASI-HARMONIC APPROXIMATION) CALCULATED WITH PTRAJ:')
            final_output.add_section(qhnorm.print_summary())
        if INPUT['alarun']:
            qhmutant = app.calc_types.mutant['qh']
            final_output.writeline(mut_str + ' MUTANT')
            final_output.writeline('ENTROPY RESULTS (QUASI-HARMONIC APPROXIMATION) CALCULATED WITH PTRAJ:')
            final_output.add_section(qhmutant.print_summary())
        if INPUT['alarun'] and not INPUT['mutant_only']:
            davg, dstd = _get_diff(qhmutant.total_avg(), qhnorm.total_avg())
            final_output.add_section(f'\nRESULT OF ALANINE SCANNING ({mut_str}):\n'
                                     f'ΔΔS binding = {davg:9.4f} +/ {dstd:7.4f}\n')

    # end if INPUT['entropy']
    if INPUT['interaction_entropy']:
        ie_inconsistent = False
        if not INPUT['mutant_only']:
            ienorm = app.calc_types['ie']
            final_output.writeline('ENTROPY RESULTS (INTERACTION ENTROPY):')
            final_output.add_section(ienorm.print_summary())
            for m in ienorm.data:
                if ienorm.data[m]['sigma'] > 3.6:
                    ie_inconsistent = True
        if INPUT['alarun']:
            iemutant = app.calc_types.mutant['ie']
            final_output.writeline(mut_str + ' MUTANT')
            final_output.writeline('ENTROPY RESULTS (INTERACTION ENTROPY):')
            final_output.add_section(iemutant.print_summary())
            for m in iemutant.data:
                if iemutant.data[m]['sigma'] > 3.6:
                    ie_inconsistent = True
        if INPUT['alarun'] and not INPUT['mutant_only']:
            text = f'\nRESULT OF ALANINE SCANNING ({mut_str}):\n'
            for model in ienorm.data:
                davg, dstd = _get_diff(iemutant.data[model]['iedata'], ienorm.data[model]['iedata'])
                text += f'ΔΔS binding ({model.upper()}) = {davg:9.4f} +/- {dstd:9.4f}\n'
            final_output.add_section(text)
        if ie_inconsistent:
            final_output.writeline(
                'WARNING: THE INTERACTION ENERGY STANDARD DEVIATION [ σ(Int. Energy) ]\n'
                'IS GREATER THAN 3.6 kcal/mol (~15 kJ/mol). THUS, THE INTERACTION ENTROPY VALUES ARE\n'
                'NOT RELIABLE. CHECK THIS PAPER FOR MORE INFO (https://doi.org/10.1021/acs.jctc.1c00374)\n\n')

    if INPUT['c2_entropy']:
        c2_inconsistent = False
        if not INPUT['mutant_only']:
            c2norm = app.calc_types['c2']
            final_output.writeline('ENTROPY RESULTS (C2 ENTROPY):')
            final_output.add_section(c2norm.print_summary())
            for m in c2norm.data:
                if c2norm.data[m]['sigma'] > 3.6:
                    c2_inconsistent = True
        if INPUT['alarun']:
            c2mutant = app.calc_types.mutant['c2']
            final_output.writeline(mut_str + ' MUTANT')
            final_output.writeline('ENTROPY RESULTS (C2 ENTROPY):')
            final_output.add_section(c2mutant.print_summary())
            for m in c2mutant.data:
                if c2mutant.data[m]['sigma'] > 3.6:
                    c2_inconsistent = True
        if INPUT['alarun'] and not INPUT['mutant_only']:
            text = '\nRESULT OF ALANINE SCANNING (%s):\n' % mut_str
            for model in c2norm.data:
                davg, dstd = _get_diff(c2mutant.data[model]['c2data'], c2norm.data[model]['c2data'])
                text += f'ΔΔS binding ({model.upper()}) = {davg:9.4f} +/- {dstd:9.4f}\n'
            final_output.add_section(text)
        if c2_inconsistent:
            final_output.writeline(
                'WARNING: THE INTERACTION ENERGY STANDARD DEVIATION [ σ(Int. Energy)]\n'
                'IS GREATER THAN 3.6 kcal/mol (~15 kJ/mol). THUS, THE C2 ENTROPY VALUES ARE NOT\n'
                'RELIABLE. CHECK THIS PAPER FOR MORE INFO (https://doi.org/10.1021/acs.jctc.1c00374)\n')

    # Now print out the normal mode results
    if INPUT['nmoderun']:
        if not INPUT['mutant_only']:
            nm_sys_norm = app.calc_types['nmode']['delta']
            final_output.write('ENTROPY RESULTS (HARMONIC APPROXIMATION) CALCULATED WITH NMODE:\n\n')
            final_output.add_section(nm_sys_norm.print_summary())
            # Now dump the energy vectors in CSV format
            if FILES.energyout:
                energyvectors.writerow(['NMODE entropy results'])
                nm_sys_norm.print_vectors(energyvectors)
                energyvectors.writerow([])

        if INPUT['alarun']:
            nm_sys_mut = app.calc_types.mutant['nmode']['delta']
            final_output.write(mut_str + ' MUTANT')
            final_output.write('ENTROPY RESULTS (HARMONIC APPROXIMATION) CALCULATED WITH NMODE:\n')
            final_output.add_section(nm_sys_mut.print_summary())
            # Now dump the energy vectors in CSV format
            if FILES.energyout:
                energyvectors.writerow([mut_str + ' Mutant NMODE entropy results'])
                nm_sys_mut.print_vectors(energyvectors)
                energyvectors.writerow([])

        # Now calculate the effect of alanine scanning
        if INPUT['alarun'] and not INPUT['mutant_only']:
            davg, dstd = _get_diff(nm_sys_mut.data['Total'], nm_sys_norm.data['Total'])
            final_output.add_section(f'\nRESULT OF ALANINE SCANNING ({mut_str}):\n'
                                      f'ΔΔS binding = {davg:9.4f} +/- {dstd:9.4f}\n')

    # end if INPUT['nmoderun']

    triggers = ('gbrun', 'pbrun', 'rismrun_std', 'rismrun_gf')
    outkeys = ('gb', 'pb', 'rism std', 'rism gf')
    headers = ('\nGENERALIZED BORN:\n\n', '\nPOISSON BOLTZMANN:\n\n',
               '\n3D-RISM:\n\n', '\n3D-RISM (Gauss. Fluct.):\n\n')
    # Now print out the Free Energy results
    for i, key in enumerate(outkeys):
        if not INPUT[triggers[i]]:
            continue

        if not INPUT['mutant_only']:
            sys_norm = app.calc_types[key]['delta']
            final_output.write(headers[i])
            final_output.add_section(sys_norm.print_summary())
            # Dump energy vectors to a CSV
            if FILES.energyout:
                energyvectors.writerow([headers[i].strip()])
                sys_norm.print_vectors(energyvectors)
                energyvectors.writerow([])

            # Combine with the entropy(ies)
            if INPUT['qh_entropy']:
                qh_davg, qh_dstd = _get_diff(sys_norm.data['DELTA TOTAL'], qhnorm.total_avg())
                final_output.add_section('Using Quasi-harmonic Entropy Approximation:\n'
                                         f'ΔG binding = {qh_davg:9.4f} +/- {qh_dstd:7.4f}\n')
            if INPUT['interaction_entropy']:
                ie_davg, ie_dstd = _get_sum(sys_norm.data['DELTA TOTAL'], ienorm.data[key]['iedata'])
                final_output.add_section(f"Using Interaction Entropy Approximation:\n"
                                         f"ΔG binding = {ie_davg:9.4f} +/- {ie_dstd:7.4f}\n")
            if INPUT['c2_entropy']:
                c2_davg, c2_dstd = _get_sum(sys_norm.data['DELTA TOTAL'], c2norm.data[key]['c2data'])
                final_output.add_section(f"Using C2 Entropy Approximation:\n"
                                         f"ΔG binding = {c2_davg:9.4f} +/- {c2_dstd:7.4f}\n")
            if INPUT['nmoderun']:
                nm_davg, nm_dstd = _get_diff(sys_norm.data['DELTA TOTAL'], nm_sys_norm.data['Total'])
                final_output.add_section('Using Normal Mode Entropy Approximation:\n'
                                         f"ΔG binding = {nm_davg:9.4f} +/- {nm_dstd:7.4f}\n")

        if INPUT['alarun']:
            sys_mut = app.calc_types.mutant[key]['delta']
            final_output.write('%s MUTANT:%s' % (mut_str, headers[i]))
            final_output.add_section(sys_mut.print_summary())
            # Dump energy vectors to a CSV
            if FILES.energyout:
                energyvectors.writerow([mut_str + ' Mutant ' + headers[i]])
                sys_mut.print_vectors(energyvectors)
                energyvectors.writerow([])

            if INPUT['qh_entropy']:
                mqh_davg, mqh_dstd = _get_diff(sys_mut.data['DELTA TOTAL'], qhmutant.total_avg())
                final_output.add_section('Using Quasi-harmonic Entropy Approximation:\n'
                                         f'ΔG binding = {mqh_davg:9.4f} +/- {mqh_dstd:7.4f}\n')
            if INPUT['interaction_entropy']:
                mie_davg, mie_dstd = _get_sum(sys_mut.data['DELTA TOTAL'], iemutant.data[key]['iedata'])
                final_output.add_section(f"Using Interaction Entropy Approximation:\n"
                                         f"ΔG binding = {mie_davg:9.4f} +/- {mie_dstd:7.4f}\n")
            if INPUT['c2_entropy']:
                mc2_davg, mc2_dstd = _get_sum(sys_mut.data['DELTA TOTAL'], c2mutant.data[key]['c2data'])
                final_output.add_section(f"Using C2 Entropy Approximation:\n"
                                         f"ΔG binding = {mc2_davg:9.4f} +/- {mc2_dstd:7.4f}\n")
            if INPUT['nmoderun']:
                mnm_davg, mnm_dstd = _get_diff(sys_mut.data['DELTA TOTAL'], nm_sys_mut.data['Total'])
                final_output.add_section('Using Normal Mode Entropy Approximation:\n'
                                         f"ΔG binding = {mnm_davg:9.4f} +/- {mnm_dstd:7.4f}\n")

        if INPUT['alarun'] and not INPUT['mutant_only']:
            davg, dstd = _get_diff(sys_mut.data['DELTA TOTAL'], sys_norm.data['DELTA TOTAL'])
            final_output.write(f'\nRESULT OF ALANINE SCANNING ({mut_str}):\n' 
                               f'ΔΔH binding = {davg:9.4f} +/- {dstd:7.4f}\n')

            if INPUT['qh_entropy']:
                davg, _ =_get_diff(mqh_davg, qh_davg)
                dstd = utils.get_std(mqh_dstd, qh_dstd)
                final_output.write('\n   (quasi-harmonic entropy)\n'
                                   f'ΔΔG binding = {davg:9.4f} +/- {dstd:7.4f}\n')
            if INPUT['interaction_entropy']:
                davg, _ = _get_diff(mie_davg, ie_davg)
                dstd = utils.get_std(mie_dstd, ie_dstd)
                final_output.write('\n   (interaction entropy)\n'
                                   f'ΔΔG binding = {davg:9.4f} +/- {dstd:7.4f}\n')
            if INPUT['c2_entropy']:
                davg, _ = _get_diff(mc2_davg, c2_davg)
                dstd = utils.get_std(mc2_dstd, c2_dstd)
                final_output.write('\n   (C2 entropy)\n'
                                   f'ΔΔG binding = {davg:9.4f} +/- {dstd:7.4f}\n')
            if INPUT['nmoderun']:
                davg, _ = _get_diff(mnm_davg, nm_davg)
                dstd = utils.get_std(mnm_dstd, nm_dstd)
                final_output.write('\n   (normal mode entropy)\n'
                                   f'ΔΔG binding = {davg:9.4f} +/- {dstd:7.4f}\n')
            final_output.separate()

    # end for solv in ['gbrun', 'pbrun', ...]

    if FILES.energyout:
        ene_csv.close()


# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

def write_decomp_stability_output(FILES, INPUT, size, prmtop_system,
                                  mutant_system, mutstr, pre):
    """ Write output file for stability decomposition calculations """
    from csv import writer
    from datetime import datetime
    from GMXMMPBSA.amber_outputs import DecompOut, PairDecompOut, idecompString

    DecompClass = DecompOut if INPUT['idecomp'] in [1, 2] else PairDecompOut

    # Open up the files
    if INPUT['csv_format']:
        dec_out_file = OutputFile(FILES.decompout, 'wb', 0)
        decompout = writer(dec_out_file)
        decompout.writerow(['| Run on %s' % datetime.now().ctime()])
        if INPUT['gbrun']:
            decompout.writerow(['| GB non-polar solvation energies calculated with gbsa=2'])
        decompout.writerow([idecompString[INPUT['idecomp']]])
    else:
        dec_out_file = OutputFile(FILES.decompout, 'w', 0)
        decompout = dec_out_file
        decompout.write_date()
        if INPUT['gbrun']:
            decompout.add_comment('| GB non-polar solvation energies calculated with gbsa=2')
        decompout.writeline(idecompString[INPUT['idecomp']])

    # Open up the energyvector CSV file
    if FILES.dec_energies:
        dec_energies = open(FILES.dec_energies, 'w')

    # GB first
    if INPUT['gbrun']:
        # See if we're dumping all energy terms
        csv_prefix = ''
        if FILES.dec_energies: csv_prefix = pre + 'gb'

        # See if we do normal system
        if not INPUT['mutant_only']:
            gb_com = DecompClass(pre + 'complex_gb.mdout', prmtop_system.complex_prmtop, INPUT['surften'], csv_prefix,
                                 size, INPUT['dec_verbose'])
            gb_com.fill_all_terms()
            if INPUT['csv_format']:
                decompout.writerow(['Energy Decomposition Analysis (All units kcal/mol): Generalized Born Solvent'])
                decompout.writerow([])
                gb_com.write_summary_csv(gb_com.numframes, decompout)
                decompout.writerow([])
            else:
                decompout.writeline('Energy Decomposition Analysis (All units kcal/mol): Generalized Born Solvent')
                decompout.writeline('')
                gb_com.write_summary(gb_com.numframes, decompout)
                decompout.writeline('')
            # Copy over all temporary .csv files created by gb_com to the total
            # utils.concatenate removes the temporary files after they are copied
            if FILES.dec_energies:
                # Close out the CSV files
                del gb_com.csvwriter
                dec_energies.write('Generalized Born Decomposition Energies' + ls)
                for token in gb_com.allowed_tokens:
                    utils.concatenate(dec_energies, csv_prefix + '.' + token + '.csv')
                dec_energies.write(ls)

        if INPUT['alarun']:
            gb_com = DecompClass(pre + 'mutant_complex_gb.mdout', mutant_system.complex_prmtop, INPUT['surften'],
                                 csv_prefix, size, INPUT['dec_verbose'])
            gb_com.fill_all_terms()
            if INPUT['csv_format']:
                decompout.writerow(['Energy Decomposition Analysis (All units kcal/mol): Generalized Born Solvent '
                                    f'({mutstr})'])
                decompout.writerow([])
                gb_com.write_summary_csv(gb_com.numframes, decompout)
                decompout.writerow([])
            else:
                decompout.writeline('Energy Decomposition Analysis (All units kcal/mol): Generalized Born Solvent '
                                    f'({mutstr})')
                decompout.writeline('')
                gb_com.write_summary(gb_com.numframes, decompout)
                decompout.writeline('')
            # Copy over all temporary .csv files created by gb_com to the total
            # utils.concatenate removes the temporary files after they are copied
            if FILES.dec_energies:
                # Close out the CSV files
                del gb_com.csvwriter
                dec_energies.write('GB Decomposition Energies (%s mutant)' %
                                   mutstr + ls)
                for token in gb_com.allowed_tokens:
                    utils.concatenate(dec_energies, csv_prefix + '.' + token + '.csv')
                dec_energies.write(ls)

    # PB next
    if INPUT['pbrun']:
        csv_prefix = pre + 'pb' if FILES.dec_energies else ''
        # See if we do normal system
        if not INPUT['mutant_only']:
            pb_com = DecompClass(pre + 'complex_pb.mdout', prmtop_system.complex_prmtop, INPUT['surften'], csv_prefix,
                                 size, INPUT['dec_verbose'])
            pb_com.fill_all_terms()
            if INPUT['csv_format']:
                decompout.writerow(['Energy Decomposition Analysis (All units kcal/mol): Generalized Born Solvent'])
                decompout.writerow([])
                pb_com.write_summary_csv(pb_com.numframes, decompout)
                decompout.writerow([])
            else:
                decompout.writeline('Energy Decomposition Analysis (All units kcal/mol): Poisson Boltzmann Solvent')
                decompout.writeline('')
                pb_com.write_summary(pb_com.numframes, decompout)
                decompout.writeline('')
            # Dump the energy vectors if requested
            if FILES.dec_energies:
                # Close out the CSV files
                del pb_com.csvwriter
                dec_energies.write('Poisson Boltzmann Decomposition Energies' + ls)
                for token in pb_com.allowed_tokens:
                    utils.concatenate(dec_energies, csv_prefix + '.' + token + '.csv')
                dec_energies.write(ls)
        if INPUT['alarun']:
            pb_com = DecompClass(pre + 'mutant_complex_pb.mdout', mutant_system.complex_prmtop, INPUT['surften'],
                                 csv_prefix, size, INPUT['dec_verbose'])
            pb_com.fill_all_terms()
            if INPUT['csv_format']:
                decompout.writerow(['Energy Decomposition Analysis (All units kcal/mol): Poisson Boltzmann Solvent '
                                    f'({mutstr})'])
                decompout.writerow([])
                pb_com.write_summary_csv(pb_com.numframes, decompout)
                decompout.writerow([])
            else:
                decompout.writeline('Energy Decomposition Analysis (All units kcal/mol): Generalized Born Solvent '
                                    f'({mutstr})')
                decompout.writeline('')
                pb_com.write_summary(pb_com.numframes, decompout)
                decompout.writeline('')
            # Dump the energy vectors if requested
            if FILES.dec_energies:
                del pb_com.csvwriter
                dec_energies.write('PB Decomposition Energies (%s mutant)' % mutstr + ls)
                for token in gb_com.allowed_tokens:
                    utils.concatenate(dec_energies, csv_prefix + '.' + token + '.csv')
                dec_energies.write(ls)

    if FILES.dec_energies: dec_energies.close()
    dec_out_file.close()


# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

def write_decomp_binding_output(FILES, INPUT, size, prmtop_system,
                                mutant_system, mutstr, pre):
    """ Write output file for binding free energy decomposition calculations """
    from csv import writer
    from datetime import datetime
    from GMXMMPBSA.amber_outputs import (DecompOut, PairDecompOut, DecompBinding,
                               PairDecompBinding, MultiTrajDecompBinding,
                               MultiTrajPairDecompBinding)

    multitraj = bool(FILES.receptor_trajs or FILES.ligand_trajs)
    # Single trajectory
    if not multitraj:
        # Per-residue
        if INPUT['idecomp'] in [1, 2]:
            BindingClass = DecompBinding
            SingleClass = DecompOut
        # Pairwise
        else:
            BindingClass = PairDecompBinding
            SingleClass = PairDecompOut
    # Multiple Traj
    # per-residue
    elif INPUT['idecomp'] in [1, 2]:
        BindingClass = MultiTrajDecompBinding
        SingleClass = DecompOut
    # Pairwise
    else:
        BindingClass = MultiTrajPairDecompBinding
        SingleClass = PairDecompOut

    # First open up our output file and turn it into a CSV writer if necessary
    if INPUT['csv_format']:
        dec_out_file = OutputFile(FILES.decompout, 'wb', 0)
        decompout = writer(dec_out_file)
        decompout.writerow(['| Run on %s' % datetime.now().ctime()])
        if INPUT['gbrun']:
            decompout.writerow(['| GB non-polar solvation energies calculated with gbsa=2'])
    else:
        dec_out_file = OutputFile(FILES.decompout, 'w', 0)
        decompout = dec_out_file
        decompout.write_date()
        if INPUT['gbrun']:
            decompout.add_comment('| GB non-polar solvation energies calculated with gbsa=2')

    # Open up the CSV energy vector file
    if FILES.dec_energies: dec_energies = open(FILES.dec_energies, 'w')

    # First we do GB
    if INPUT['gbrun']:
        # Setup for energy dump
        if FILES.dec_energies:
            csv_prefix = pre + 'gb_%s'
        else:
            csv_prefix = csv_pre = ''
        # Normal system
        if not INPUT['mutant_only']:
            if csv_prefix:
                csv_pre = csv_prefix % 'com'
            gb_com = SingleClass(pre + 'complex_gb.mdout', prmtop_system.complex_prmtop, INPUT['surften'], csv_pre,
                                 size, INPUT['dec_verbose'])
            if csv_prefix:
                csv_pre = csv_prefix % 'rec'
            gb_rec = SingleClass(pre + 'receptor_gb.mdout', prmtop_system.receptor_prmtop, INPUT['surften'], csv_pre,
                                 size, INPUT['dec_verbose'])
            if csv_prefix:
                csv_pre = csv_prefix % 'lig'
            gb_lig = SingleClass(pre + 'ligand_gb.mdout', prmtop_system.ligand_prmtop, INPUT['surften'], csv_pre,
                                 size, INPUT['dec_verbose'])
            if csv_prefix:
                csv_pre = pre + 'gb_bind'
            gb_bind = BindingClass(gb_com, gb_rec, gb_lig, prmtop_system, INPUT['idecomp'], INPUT['dec_verbose'],
                                   decompout, csv_pre, 'Energy Decomposition Analysis (All units kcal/mol): '
                                                       'Generalized Born solvent')
            # Write the data to the output file
            gb_bind.parse_all()
            # Now it's time to dump everything to the CSV file
            if FILES.dec_energies:
                del gb_com.csvwriter, gb_rec.csvwriter, gb_lig.csvwriter
                del gb_bind.csvwriter
                dec_energies.write('Generalized Born Decomposition Energies' + ls)
                dec_energies.write(ls + 'Complex:' + ls)
                for token in gb_com.allowed_tokens:
                    utils.concatenate(dec_energies, pre + 'gb_com.' + token + '.csv')
                dec_energies.write(ls + 'Receptor:' + ls)
                for token in gb_rec.allowed_tokens:
                    utils.concatenate(dec_energies, pre + 'gb_rec.' + token + '.csv')
                dec_energies.write(ls + 'Ligand:' + ls)
                for token in gb_lig.allowed_tokens:
                    utils.concatenate(dec_energies, pre + 'gb_lig.' + token + '.csv')
                # DELTAs for each frame only computed if this is a single trajectory
                if not multitraj:
                    dec_energies.write(ls + 'DELTAS:' + ls)
                    for token in gb_com.allowed_tokens:
                        utils.concatenate(dec_energies, pre + 'gb_bind.' + token + '.csv')
                    dec_energies.write(ls)
        # Mutant system
        if INPUT['alarun']:
            if csv_prefix:
                csv_pre = csv_prefix % 'com'
            gb_com = SingleClass(pre + 'mutant_complex_gb.mdout', mutant_system.complex_prmtop, INPUT['surften'],
                                 csv_pre, size, INPUT['dec_verbose'])
            if csv_prefix:
                csv_pre = csv_prefix % 'rec'
            gb_rec = SingleClass(pre + 'mutant_receptor_gb.mdout', mutant_system.receptor_prmtop, INPUT['surften'],
                                 csv_pre, size, INPUT['dec_verbose'])
            if csv_prefix:
                csv_pre = csv_prefix % 'lig'
            gb_lig = SingleClass(pre + 'mutant_ligand_gb.mdout', mutant_system.ligand_prmtop, INPUT['surften'],
                                 csv_pre, size, INPUT['dec_verbose'])
            if csv_prefix:
                csv_pre = pre + 'gb_bind'
            gb_bind = BindingClass(gb_com, gb_rec, gb_lig, mutant_system, INPUT['idecomp'], INPUT['dec_verbose'],
                                   decompout, csv_pre, 'Energy Decomposition Analysis (All units kcal/mol): '
                                                       f'Generalized Born solvent ({mutstr})')
            # Write the data to the output file
            gb_bind.parse_all()
            # Now it's time to dump everything to the CSV file
            if FILES.dec_energies:
                del gb_com.csvwriter, gb_rec.csvwriter, gb_lig.csvwriter
                del gb_bind.csvwriter
                dec_energies.write('GB Decomposition Energies (%s mutant)' %
                                   mutstr + ls)
                dec_energies.write(ls + 'Complex:' + ls)
                for token in gb_com.allowed_tokens:
                    utils.concatenate(dec_energies, pre + 'gb_com.' + token + '.csv')
                dec_energies.write(ls + 'Receptor:' + ls)
                for token in gb_rec.allowed_tokens:
                    utils.concatenate(dec_energies, pre + 'gb_rec.' + token + '.csv')
                dec_energies.write(ls + 'Ligand:' + ls)
                for token in gb_lig.allowed_tokens:
                    utils.concatenate(dec_energies, pre + 'gb_lig.' + token + '.csv')
                # DELTAs for each frame only computed if this is a single trajectory
                if not multitraj:
                    dec_energies.write(ls + 'DELTAS:' + ls)
                    for token in gb_com.allowed_tokens:
                        utils.concatenate(dec_energies, pre + 'gb_bind.' + token + '.csv')
                    dec_energies.write(ls)
    # Next do PB
    if INPUT['pbrun']:
        # Setup for energy dump
        if FILES.dec_energies:
            csv_prefix = pre + 'pb_%s'
        else:
            csv_prefix = csv_pre = ''
        # Normal system
        if not INPUT['mutant_only']:
            if csv_prefix:
                csv_pre = csv_prefix % 'com'
            pb_com = SingleClass(pre + 'complex_pb.mdout', prmtop_system.complex_prmtop, INPUT['cavity_surften'],
                                 csv_pre, size, INPUT['dec_verbose'])
            if csv_prefix:
                csv_pre = csv_prefix % 'rec'
            pb_rec = SingleClass(pre + 'receptor_pb.mdout', prmtop_system.receptor_prmtop, INPUT['cavity_surften'],
                                 csv_pre, size, INPUT['dec_verbose'])
            if csv_prefix:
                csv_pre = csv_prefix % 'lig'
            pb_lig = SingleClass(pre + 'ligand_pb.mdout', prmtop_system.ligand_prmtop, INPUT['cavity_surften'],
                                 csv_pre, size, INPUT['dec_verbose'])
            if csv_prefix:
                csv_pre = pre + 'pb_bind'
            pb_bind = BindingClass(pb_com, pb_rec, pb_lig, prmtop_system, INPUT['idecomp'], INPUT['dec_verbose'],
                                   decompout, csv_pre, 'Energy Decomposition Analysis (All units kcal/mol): Poisson '
                                                       'Boltzmann solvent')
            # Write the data to the output file
            pb_bind.parse_all()
            # Now it's time to dump everything to the CSV file
            if FILES.dec_energies:
                del pb_com.csvwriter, pb_rec.csvwriter, pb_lig.csvwriter
                del pb_bind.csvwriter
                dec_energies.write('Poisson Boltzmann Decomposition Energies' + ls)
                dec_energies.write(ls + 'Complex:' + ls)
                for token in pb_com.allowed_tokens:
                    utils.concatenate(dec_energies, pre + 'pb_com.' + token + '.csv')
                dec_energies.write(ls + 'Receptor:' + ls)
                for token in pb_rec.allowed_tokens:
                    utils.concatenate(dec_energies, pre + 'pb_rec.' + token + '.csv')
                dec_energies.write(ls + 'Ligand:' + ls)
                for token in pb_lig.allowed_tokens:
                    utils.concatenate(dec_energies, pre + 'pb_lig.' + token + '.csv')
                # DELTAs for each frame only computed if this is a single trajectory
                if not multitraj:
                    dec_energies.write(ls + 'DELTAS:' + ls)
                    for token in pb_com.allowed_tokens:
                        utils.concatenate(dec_energies, pre + 'pb_bind.' + token + '.csv')
                    dec_energies.write(ls)
        # Mutant system
        if INPUT['alarun']:
            if csv_prefix:
                csv_pre = csv_prefix % 'com'
            pb_com = SingleClass(pre + 'mutant_complex_pb.mdout', mutant_system.complex_prmtop, INPUT['cavity_surften'],
                                 csv_pre, size, INPUT['dec_verbose'])
            if csv_prefix:
                csv_pre = csv_prefix % 'rec'
            pb_rec = SingleClass(pre + 'mutant_receptor_pb.mdout', mutant_system.receptor_prmtop,
                                 INPUT['cavity_surften'],
                                 csv_pre, size, INPUT['dec_verbose'])
            if csv_prefix:
                csv_pre = csv_prefix % 'lig'
            pb_lig = SingleClass(pre + 'mutant_ligand_pb.mdout', mutant_system.ligand_prmtop, INPUT['cavity_surften'],
                                 csv_pre, size, INPUT['dec_verbose'])
            if csv_prefix:
                csv_pre = pre + 'pb_bind'
            pb_bind = BindingClass(pb_com, pb_rec, pb_lig, mutant_system, INPUT['idecomp'], INPUT['dec_verbose'],
                                   decompout, csv_pre, 'Energy Decomposition Analysis (All units kcal/mol): Poisson '
                                                       f'Boltzmann solvent ({mutstr})')
            # Write the data to the output file
            pb_bind.parse_all()
            if FILES.dec_energies:
                del pb_com.csvwriter, pb_rec.csvwriter, pb_lig.csvwriter
                del pb_bind.csvwriter
                dec_energies.write('PB Decomposition Energies (%s mutant)' % mutstr + ls)
                dec_energies.write(ls + 'Complex:' + ls)
                for token in pb_com.allowed_tokens:
                    utils.concatenate(dec_energies, pre + 'pb_com.' + token + '.csv')
                dec_energies.write(ls + 'Receptor:' + ls)
                for token in pb_rec.allowed_tokens:
                    utils.concatenate(dec_energies, pre + 'pb_rec.' + token + '.csv')
                dec_energies.write(ls + 'Ligand:' + ls)
                for token in pb_lig.allowed_tokens:
                    utils.concatenate(dec_energies, pre + 'pb_lig.' + token + '.csv')
                # DELTAs for each frame only computed if this is a single trajectory
                if not multitraj:
                    dec_energies.write(ls + 'DELTAS:' + ls)
                    for token in gb_com.allowed_tokens:
                        utils.concatenate(dec_energies, pre + 'pb_bind.' + token + '.csv')
                    dec_energies.write(ls)

    # Close the file(s)
    if FILES.dec_energies: dec_energies.close()
    dec_out_file.close()


# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class OutputFile(object):
    """ Main output file """

    # ==================================================

    def __init__(self, fname, mode='r', buffer=0):
        # Buffer is ignored
        self._binary = 'b' in mode.lower()
        self._handle = open(fname, mode)

    # ==================================================

    def write(self, stuff):
        # Handle bytes/str division in Py3 transparently
        if self._binary and isinstance(stuff, str):
            self._handle.write(stuff.encode())
        elif self._binary or isinstance(stuff, str):
            self._handle.write(stuff)
        else:
            self._handle.write(stuff.decode())
        self._handle.flush()

    # ==================================================

    def separate(self):
        """ Delimiter between fields in output file """
        for _ in range(2):
            self.write('-' * 79)
            self.write(ls)

    # ==================================================

    def add_section(self, input_str):
        """ Adds a section to the output file """
        self.write(input_str)
        self.separate()

    # ==================================================

    def write_date(self):
        """ Writes a date string to the output file """
        from datetime import datetime
        self.writeline('| Run on %s' % datetime.now().ctime())

    # ==================================================

    def print_file_info(self, FILES, INPUT):
        """ Prints the summary information to a file """
        from GMXMMPBSA import __version__, __mmpbsa_version__
        stability = not FILES.receptor_prmtop

        self.writeline('|gmx_MMPBSA Version=%s based on MMPBSA.py v.%s' % (__version__, __mmpbsa_version__))

        # if FILES.solvated_prmtop:
        #    self.writeline('|Solvated complex topology file:  %s' %
        #               FILES.solvated_prmtop)
        self.writeline(f'{"|Complex Structure file:":40}{FILES.complex_tpr:>40}')
        if FILES.complex_top:
            self.writeline(f'{"|Complex (GROMACS) topology file:":40}{FILES.complex_top:>40}')
        self.writeline(f'{"|Complex (AMBER) topology file:":40}{FILES.complex_prmtop:>40}')

        if not stability:
            # if FILES.receptor_trajs:
            #    if FILES.solvated_receptor_prmtop:
            #       self.writeline('|Solvated receptor topology file: %s' %
            #                  FILES.solvated_receptor_prmtop)
            if FILES.receptor_tpr:
                self.writeline(f'{"|Receptor Structure file:":40}{FILES.receptor_tpr:>40}')
            if FILES.receptor_top:
                self.writeline(f'{"|Receptor (GROMACS) topology file:":40}{FILES.receptor_top:>40}')
            self.writeline(f'{"|Receptor (AMBER) topology file:":40}{FILES.receptor_prmtop:>40}')

            # if FILES.ligand_trajs:
            # if FILES.solvated_ligand_prmtop:
            #    self.writeline('|Solvated ligand topology file:   %s' %
            #               FILES.solvated_ligand_prmtop)
            if FILES.ligand_mol2:
                self.writeline(f'{"|Ligand Structure file:":40}{FILES.ligand_mol2:>40}')
            if FILES.ligand_tpr:
                self.writeline(f'{"|Complex Structure file:":40}{FILES.ligand_tpr:>40}')
            if FILES.ligand_top:
                self.writeline(f'{"|Complex (GROMACS) topology file:":40}{FILES.ligand_top:>40}')
            self.writeline(f'{"|Complex (AMBER) topology file:":40}{FILES.ligand_prmtop:>40}')

        if INPUT['alarun']:
            self.writeline(f'{"|Mutant (AMBER) complex topology file:":40}{FILES.mutant_complex_prmtop:>40}')
            if not stability:
                self.writeline(f'{"|Mutant (AMBER) receptor topology file:":40}{FILES.mutant_receptor_prmtop:>40}')
                self.writeline(f'{"|Mutant (AMBER) ligand topology file:":40}{FILES.mutant_ligand_prmtop:>40}')

        self.write(f'{"|Initial trajectories:":40}')
        for i in range(len(FILES.complex_trajs)):
            if i == 0:
                self.writeline(f'{FILES.complex_trajs[i]:>40}')
            else:
                self.writeline(f'{"|":40}{FILES.complex_trajs[i]:>40}')

        if FILES.receptor_trajs:
            self.write(f'{"|Initial Receptor mdcrd(s):":40}')
            for i in range(len(FILES.receptor_trajs)):
                if i == 0:
                    self.writeline(f'{FILES.receptor_trajs[i]:>40}')
                else:
                    self.writeline(f'{"|":40}{FILES.receptor_trajs[i]:>40}')

        if FILES.ligand_trajs:
            self.write(f'{"|Initial Ligand mdcrd(s):":40}')
            for i in range(len(FILES.ligand_trajs)):
                if i == 0:
                    self.writeline(f'{FILES.ligand_trajs[i]:>40}')
                else:
                    self.writeline(f'{"|":40}{FILES.ligand_trajs[i]:>40}')

    # ==================================================

    def add_comment(self, comment):
        """ Adds a comment string to the output file """
        self.writeline('|' + comment)

    # ==================================================

    def writeline(self, line):
        """ Adds a newline to the end of the passed line and writes it """
        self.write(line + ls)

    # ==================================================

    def __getattr__(self, attr):
        return getattr(self._handle, attr)

    # ==================================================

    def __del__(self):
        try:
            self._handle.close()
        except:
            pass

# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
