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
            data = np.nan if self.app.INPUT[x] is None else self.app.INPUT[x]
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
        dset = grp.create_dataset('mutant_index', data= np.nan if self.app.mutant_index is None else self.app.mutant_index)
        dset = grp.create_dataset('mut_str', data=self.app.mut_str)
        dset = grp.create_dataset('using_chamber', data=self.app.using_chamber)
        dset = grp.create_dataset('input_file', data=self.app.input_file_text)

        # save the complex fixed structure
        com_fixed = ''.join(open(self.app.FILES.complex_fixed).readlines())
        dset = grp.create_dataset('COM_PDB', data=com_fixed)

        # get output files
        outfile = ''.join(open(self.app.FILES.output_file).readlines())
        dset = grp.create_dataset('output_file', data=outfile)
        if self.app.INPUT['decomp']['decomprun']:
            doutfile = ''.join(open(self.app.FILES.decompout).readlines())
            dset = grp.create_dataset('decomp_output_file', data=doutfile)

    @staticmethod
    def _e2h5(d, f):
        # key  Energy: [gb, pb, rism std, rism gf], Entropy: [nmode, qh, ie, c2]
        for key in d:
            grp = f.create_group(key)
            # key2 is complex, receptor, ligand, delta or model for ie and c2
            for key2 in d[key]:
                grp2 = grp.create_group(key2)
                for key3 in d[key][key2]:
                    dset = grp2.create_dataset(key3, data=d[key][key2][key3])

    @staticmethod
    def _decomp2h5(d, g):
        for key in d:
            # model (GB or PB)
            grp = g.create_group(key)
            # complex, receptor, ligand, delta
            for key2 in d[key]:
                grp2 = grp.create_group(key2)
                # TDC, SDC, BDC
                for key3 in d[key][key2]:
                    grp3 = grp2.create_group(key3)
                    # residue first level
                    for key4 in d[key][key2][key3]:
                        # we need to convert the res number in str since h5 only admit str as keys
                        grp4 = grp3.create_group(str(key4))
                        # residue sec level for per-wise or energy terms for per-residue
                        for key5 in d[key][key2][key3][key4]:
                            if isinstance(d[key][key2][key3][key4][key5], dict):
                                grp5 = grp4.create_group(str(key5))
                                # energy terms
                                for key6 in d[key][key2][key3][key4][key5]:
                                    dset = grp5.create_dataset(key6, data=d[key][key2][key3][key4][key5][key6])
                            else:
                                # energy terms
                                dset = grp4.create_dataset(key5, data=d[key][key2][key3][key4][key5])


def write_outputs(app):
    """ Writes stability or binding output file """
    import csv
    # Load some objects into top-level name space
    FILES = app.FILES
    INPUT = app.INPUT
    mut_str = app.mut_str
    prmtop_system = app.normal_system
    stability = app.stability

    # Open the energy vector CSV output file if we are writing one
    if FILES.energyout:
        ene_csv = open(FILES.energyout, 'w')
        energyvectors = csv.writer(ene_csv, dialect='excel')

    # Open the file and write some initial data to it
    final_output = OutputFile(FILES.output_file, 'w')
    final_output.write_date()
    final_output.add_comment('')
    final_output.print_file_info(FILES, INPUT)
    if not stability:
        final_output.add_comment('')
        final_output.add_comment('Receptor mask:                  "%s"' % INPUT['general']['receptor_mask'])
        final_output.add_comment('Ligand mask:                    "%s"' % INPUT['general']['ligand_mask'])
        if prmtop_system.ligand_prmtop.ptr('nres') == 1:
            final_output.add_comment('Ligand residue name is:         "%s"' %
                                     prmtop_system.ligand_prmtop.parm_data['RESIDUE_LABEL'][0])
    final_output.add_comment('')
    final_output.add_comment('Calculations performed using %s complex frames' % app.numframes)
    if INPUT['nmode']['nmoderun']:
        final_output.add_comment('NMODE calculations performed using %s frames' % app.numframes_nmode)
    if not stability:
        if INPUT['general']['interaction_entropy']:
            final_output.add_comment('Interaction Entropy calculations performed using last %s frames' %
                                     ceil(app.numframes * (INPUT['general']['ie_segment'] / 100)))
        if INPUT['general']['c2_entropy']:
            final_output.add_comment('C2 Entropy Std. Dev. and Conf. Interv. (95%) have been obtained by '
                                     'bootstrapping with number of re-samplings = 2000')
    if INPUT['pb']['pbrun']:
        if INPUT['pb']['sander_apbs']:
            final_output.add_comment('Poisson Boltzmann calculations performed using iAPBS interface to sander '
                                     '(sander.APBS)')
        else:
            final_output.add_comment('Poisson Boltzmann calculations performed using internal PBSA solver in sander')
    final_output.add_comment('')

    if INPUT['gb']['gbrun']:
        if INPUT['gb']['molsurf']:
            final_output.add_comment('Generalized Born ESURF calculated using \'molsurf\' surface areas')
        else:
            final_output.add_comment('Generalized Born ESURF calculated using \'LCPO\' surface areas')
        final_output.add_comment('')

    final_output.add_comment('Using temperature = %.2f K' % INPUT['general']['temperature'])
    final_output.add_comment('All units are reported in kcal/mol')
    final_output.add_comment('')
    final_output.add_comment('SD - Sample standard deviation, SEM - Sample standard error of the mean')
    final_output.add_comment('SD(Prop.), SEM(Prop.) - SD and SEM obtained with propagation of uncertainty formula')
    final_output.add_comment('https://en.wikipedia.org/wiki/Propagation_of_uncertainty#Example_formulae')
    final_output.add_comment('')

    if INPUT['gb']['ifqnt']:
        final_output.add_comment(('QM/MM: Residues %s are treated with the ' +
                                  'Quantum Hamiltonian %s') % (INPUT['gb']['qm_residues'], INPUT['gb']['qm_theory']))
    final_output.separate()

    # First do the entropies
    if INPUT['general']['qh_entropy']:
        if not INPUT['ala']['mutant_only']:
            if stability:
                qhnorm = app.calc_types.normal['qh']['complex']
            else:
                qhnorm = app.calc_types.normal['qh']['delta']
            final_output.writeline('ENTROPY RESULTS (QUASI-HARMONIC APPROXIMATION) CALCULATED WITH CPPTRAJ:')
            final_output.add_section(qhnorm.summary_output())
        if INPUT['ala']['alarun']:
            if stability:
                qhnorm = app.calc_types.mutant['qh']['complex']
            else:
                qhnorm = app.calc_types.mutant['qh']['delta']
            qhmutant = app.calc_types.mutant['qh']
            final_output.writeline(mut_str + ' MUTANT')
            final_output.writeline('ENTROPY RESULTS (QUASI-HARMONIC APPROXIMATION) CALCULATED WITH CPPTRAJ:')
            final_output.add_section(qhmutant.summary_output())
        if INPUT['ala']['alarun'] and not INPUT['ala']['mutant_only']:
            ddqh_davg, ddqh_dstd = utils.calc_sub(qhmutant['TOTAL'], qhnorm['TOTAL'], mut=True)
            final_output.add_section(f'\nRESULT OF ALANINE SCANNING ({mut_str}):\n'
                                     f'-TΔΔS{"" if stability else " binding"} = {ddqh_davg:9.2f} +/ {ddqh_dstd:7.2f}\n')

    if not stability:
        # end if INPUT['entropy']
        if INPUT['general']['interaction_entropy']:
            ie_inconsistent = False
            if not INPUT['ala']['mutant_only']:
                ienorm = app.calc_types.normal['ie']
                final_output.writeline('ENTROPY RESULTS (INTERACTION ENTROPY):')
                final_output.add_section(ienorm.summary_output())
                # Now dump the energy vectors in CSV format
                if FILES.energyout:
                    energyvectors.writerow(['Interaction entropy results'])
                    ienorm._print_vectors(energyvectors)
                    energyvectors.writerow([])

                for m in ienorm:
                    if ienorm[m]['sigma'] > 3.6:
                        ie_inconsistent = True
            if INPUT['ala']['alarun']:
                iemutant = app.calc_types.mutant['ie']
                final_output.writeline(mut_str + ' MUTANT')
                final_output.writeline('ENTROPY RESULTS (INTERACTION ENTROPY):')
                final_output.add_section(iemutant.summary_output())
                for m in iemutant:
                    if iemutant[m]['sigma'] > 3.6:
                        ie_inconsistent = True
            if INPUT['ala']['alarun'] and not INPUT['ala']['mutant_only']:
                text = f'\nRESULT OF ALANINE SCANNING ({mut_str}):\n'
                for model in ienorm:
                    davg, dstd = utils.calc_sub(iemutant[model]['iedata'], ienorm[model]['iedata'], mut=True)
                    text += f'-TΔΔS binding ({model.upper()}) = {davg:9.2f} +/- {dstd:9.2f}\n'
                final_output.add_section(text)
            if ie_inconsistent:
                final_output.writeline(
                    'WARNING: THE INTERACTION ENERGY STANDARD DEVIATION [ σ(Int. Energy) ]\n'
                    'IS GREATER THAN 3.6 kcal/mol (~15 kJ/mol). THUS, THE INTERACTION ENTROPY VALUES ARE\n'
                    'NOT RELIABLE. CHECK THIS PAPER FOR MORE INFO (https://doi.org/10.1021/acs.jctc.1c00374)\n\n')

        if INPUT['general']['c2_entropy']:
            c2_inconsistent = False
            if not INPUT['ala']['mutant_only']:
                c2norm = app.calc_types.normal['c2']
                final_output.writeline('ENTROPY RESULTS (C2 ENTROPY):')
                final_output.add_section(c2norm.summary_output())
                for m in c2norm:
                    if c2norm[m]['sigma'] > 3.6:
                        c2_inconsistent = True
            if INPUT['ala']['alarun']:
                c2mutant = app.calc_types.mutant['c2']
                final_output.writeline(mut_str + ' MUTANT')
                final_output.writeline('ENTROPY RESULTS (C2 ENTROPY):')
                final_output.add_section(c2mutant.summary_output())
                for m in c2mutant:
                    if c2mutant[m]['sigma'] > 3.6:
                        c2_inconsistent = True
            if INPUT['ala']['alarun'] and not INPUT['ala']['mutant_only']:
                text = '\nRESULT OF ALANINE SCANNING (%s):\n' % mut_str
                for model in c2norm:
                    davg, dstd = utils.calc_sub(c2mutant[model]['c2data'], c2norm[model]['c2data'], mut=True)
                    text += f'-TΔΔS binding ({model.upper()}) = {davg:9.2f} +/- {dstd:9.2f}\n'
                final_output.add_section(text)
            if c2_inconsistent:
                final_output.writeline(
                    'WARNING: THE INTERACTION ENERGY STANDARD DEVIATION [ σ(Int. Energy)]\n'
                    'IS GREATER THAN 3.6 kcal/mol (~15 kJ/mol). THUS, THE C2 ENTROPY VALUES ARE NOT\n'
                    'RELIABLE. CHECK THIS PAPER FOR MORE INFO (https://doi.org/10.1021/acs.jctc.1c00374)\n')

    # Now print out the normal mode results
    if INPUT['nmode']['nmoderun']:
        if not INPUT['ala']['mutant_only']:
            if stability:
                nm_sys_norm = app.calc_types.normal['nmode']['complex']
            else:
                nm_sys_norm = app.calc_types.normal['nmode']['delta']

            final_output.write('ENTROPY RESULTS (HARMONIC APPROXIMATION) CALCULATED WITH NMODE:\n\n')
            final_output.add_section(nm_sys_norm.summary_output())
            # Now dump the energy vectors in CSV format
            if FILES.energyout:
                energyvectors.writerow(['NMODE entropy results'])
                nm_sys_norm._print_vectors(energyvectors)
                energyvectors.writerow([])

        if INPUT['ala']['alarun']:
            if stability:
                nm_sys_mut = app.calc_types.mutant['nmode']['complex']
            else:
                nm_sys_mut = app.calc_types.mutant['nmode']['delta']
            final_output.writeline(mut_str + ' MUTANT')
            final_output.writeline('ENTROPY RESULTS (HARMONIC APPROXIMATION) CALCULATED WITH NMODE:')
            final_output.add_section(nm_sys_mut.summary_output())
            # Now dump the energy vectors in CSV format
            if FILES.energyout:
                energyvectors.writerow([mut_str + ' Mutant NMODE entropy results'])
                nm_sys_mut._print_vectors(energyvectors)
                energyvectors.writerow([])

        # Now calculate the effect of alanine scanning
        if INPUT['ala']['alarun'] and not INPUT['ala']['mutant_only']:
            nm_sys_norm = app.calc_types.mut_norm['nmode']
            mut_norm = app.calc_types.mut_norm['nmode']['delta']
            final_output.add_section(mut_norm.summary_output())

            final_output.write('ENTROPY RESULTS (HARMONIC APPROXIMATION) CALCULATED WITH NMODE:\n\n')
            final_output.add_section(nm_sys_norm.summary_output())

            davg, dstd = utils.calc_sub(nm_sys_mut['TOTAL'], nm_sys_norm['TOTAL'], mut=True)
            final_output.add_section(f'\nRESULT OF ALANINE SCANNING ({mut_str}):\n'
                                     f'-TΔΔS{"" if stability else " binding"} = {davg:9.2f} +/- {dstd:9.2f}\n')

    # end if INPUT['nmode']['nmoderun']
    nmls = ('gb', 'pb', 'rism', 'rism', 'rism', 'gbnsr6')
    triggers = ('gbrun', 'pbrun', 'rismrun_std', 'rismrun_gf', 'rismrun_pcplus', 'gbnsr6run')
    outkeys = ('gb', 'pb', 'rism std', 'rism gf', 'rism pcplus', 'gbnsr6')
    headers = ('\nGENERALIZED BORN:\n\n', '\nPOISSON BOLTZMANN:\n\n',
               '\n3D-RISM:\n\n', '\n3D-RISM (Gauss. Fluct.):\n\n', '\n3D-RISM (PC+):\n\n',
               '\nGENERALIZED BORN (GBNSR6):\n\n')
    # Now print out the Free Energy results
    for i, key in enumerate(outkeys):
        if not INPUT[nmls[i]][triggers[i]]:
            continue
        if not INPUT['ala']['mutant_only']:
            final_output.write(headers[i])
            if stability:
                sys_norm = app.calc_types.normal[key]['complex']
            else:
                sys_norm = app.calc_types.normal[key]['delta']
                if sys_norm.inconsistent:
                    final_output.writeline(sys_norm.report_inconsistency())
                final_output.add_section(app.calc_types.normal[key]['complex'].summary_output())
                final_output.add_section(app.calc_types.normal[key]['receptor'].summary_output())
                final_output.add_section(app.calc_types.normal[key]['ligand'].summary_output())
            final_output.add_section(sys_norm.summary_output())
            # Dump energy vectors to a CSV
            if FILES.energyout:
                energyvectors.writerow([headers[i].strip()])
                sys_norm._print_vectors(energyvectors)
                energyvectors.writerow([])

            # Combine with the entropy(ies)
            if INPUT['general']['qh_entropy']:
                dg_qh_davg, dg_qh_dstd = utils.calc_sum(sys_norm['TOTAL'], qhnorm['TOTAL'])
                final_output.add_section('Using Quasi-harmonic Entropy Approximation:\n'
                                         f'ΔG{"" if stability else " binding"} = {dg_qh_davg:9.2f} +'
                                         f'/- {dg_qh_dstd:7.2f}\n')
            if not stability:
                if INPUT['general']['interaction_entropy']:
                    dg_ie_davg, dg_ie_dstd = utils.calc_sum(sys_norm['TOTAL'], ienorm[key]['iedata'])
                    final_output.add_section(f"Using Interaction Entropy Approximation:\n"
                                             f"ΔG binding = {dg_ie_davg:9.2f} +/- {dg_ie_dstd:7.2f}\n")
                if INPUT['general']['c2_entropy']:
                    dg_c2_davg, dh_dstd = utils.calc_sum(sys_norm['TOTAL'], c2norm[key]['c2data'])
                    dg_c2_dstd = utils.get_std(dh_dstd, c2norm[key]['c2_std'])
                    final_output.add_section(f"Using C2 Entropy Approximation:\n"
                                             f"ΔG binding = {dg_c2_davg:9.2f} +/- {dg_c2_dstd:7.2f}\n")
            if INPUT['nmode']['nmoderun']:
                dg_nm_davg, dg_nm_dstd = utils.calc_sum(sys_norm['TOTAL'], nm_sys_norm['TOTAL'])
                final_output.add_section('Using Normal Mode Entropy Approximation:\n'
                                         f'ΔG{"" if stability else " binding"} = {dg_nm_davg:9.2f} +/-'
                                         f' {dg_nm_dstd:7.2f}\n')

        if INPUT['ala']['alarun']:
            final_output.write('%s MUTANT:%s' % (mut_str, headers[i]))
            if stability:
                sys_mut = app.calc_types.mutant[key]['complex']
            else:
                sys_mut = app.calc_types.mutant[key]['delta']
                if sys_mut.inconsistent:
                    final_output.writeline(sys_mut.report_inconsistency())
                final_output.add_section(app.calc_types.mutant[key]['complex'].summary_output())
                final_output.add_section(app.calc_types.mutant[key]['receptor'].summary_output())
                final_output.add_section(app.calc_types.mutant[key]['ligand'].summary_output())
            final_output.add_section(sys_mut.summary_output())
            # Dump energy vectors to a CSV
            if FILES.energyout:
                energyvectors.writerow([mut_str + ' Mutant ' + headers[i]])
                sys_mut._print_vectors(energyvectors)
                energyvectors.writerow([])

            if INPUT['general']['qh_entropy']:
                mqh_davg, mqh_dstd = utils.calc_sum(sys_mut['TOTAL'], qhmutant['TOTAL'])
                final_output.add_section('Using Quasi-harmonic Entropy Approximation:\n'
                                         f'ΔG{"" if stability else " binding"} = {mqh_davg:9.2f} +/- {mqh_dstd:7.2f}\n')
            if not stability:
                if INPUT['general']['interaction_entropy']:
                    mie_davg, mie_dstd = utils.calc_sum(sys_mut['TOTAL'], iemutant[key]['iedata'])
                    final_output.add_section(f"Using Interaction Entropy Approximation:\n"
                                             f"ΔG binding = {mie_davg:9.2f} +/- {mie_dstd:7.2f}\n")
                if INPUT['general']['c2_entropy']:
                    mc2_davg, mc2_dstd = utils.calc_sum(sys_mut['TOTAL'], c2mutant[key]['c2data'])
                    final_output.add_section(f"Using C2 Entropy Approximation:\n"
                                             f"ΔG binding = {mc2_davg:9.2f} +/- {mc2_dstd:7.2f}\n")
            if INPUT['nmode']['nmoderun']:
                mnm_davg, mnm_dstd = utils.calc_sum(sys_mut['TOTAL'], nm_sys_mut['TOTAL'])
                final_output.add_section('Using Normal Mode Entropy Approximation:\n'
                                         f'ΔG{"" if stability else " binding"} = {mnm_davg:9.2f} +/- {mnm_dstd:7.2f}\n')

        if INPUT['ala']['alarun'] and not INPUT['ala']['mutant_only']:
            mut_norm = app.calc_types.mut_norm[key]['delta']
            final_output.add_section(mut_norm.summary_output())
            ddh_davg = mut_norm['TOTAL'].mean()
            ddh_dstd = mut_norm['TOTAL'].std()

            final_output.write(f'\nRESULT OF ALANINE SCANNING ({mut_str}):\n' 
                               f"ΔΔH binding = {ddh_davg:9.2f} +/- {ddh_dstd:7.2f}\n")

            if INPUT['general']['qh_entropy']:
                ddgqh_davg = ddh_davg + ddqh_davg
                # this std is the same of ΔΔH
                final_output.write('\n   (quasi-harmonic entropy)\n'
                                   f'ΔΔG{"" if stability else " binding"} = {ddgqh_davg:9.2f} +/- {ddh_dstd:7.2f}\n')
            if not stability:
                if INPUT['general']['interaction_entropy']:
                    ddie_davg, ddie_dstd = utils.calc_sub(iemutant[key]['iedata'], ienorm[key]['iedata'], mut=True)
                    ddgie_davg = ddh_davg + ddie_davg
                    ddgie_dstd = utils.get_std(ddh_dstd, ddie_dstd)
                    final_output.write('\n   (interaction entropy)\n'
                                       f'ΔΔG binding = {ddgie_davg:9.2f} +/- {ddgie_dstd:7.2f}\n')
                if INPUT['general']['c2_entropy']:
                    ddc2_davg = c2mutant[key]['c2data'] - c2norm[key]['c2data']
                    ddc2_dstd = c2mutant[key]['c2_std'] - c2norm[key]['c2_std']
                    ddgc2_davg = ddh_davg + ddc2_davg
                    ddgc2_dstd = utils.get_std(ddh_dstd, ddc2_dstd)
                    final_output.write('\n   (C2 entropy)\n'
                                       f'ΔΔG binding = {ddgc2_davg:9.2f} +/- {ddgc2_dstd:7.2f}\n')
            if INPUT['nmode']['nmoderun']:
                ddnm_davg, ddnm_dstd = utils.calc_sub(nm_sys_mut['TOTAL'], nm_sys_norm['TOTAL'], mut=True)
                ddgnm_davg = ddh_davg + ddnm_davg
                ddgnm_dstd = utils.get_std(ddh_dstd, ddnm_dstd)
                final_output.write('\n   (normal mode entropy)\n'
                                   f'ΔΔG{"" if stability else " binding"} = {ddgnm_davg:9.2f} +/- {ddgnm_dstd:7.2f}\n')
            final_output.separate()

    # end for solv in ['gbrun', 'pbrun', ...]

    if FILES.energyout:
        ene_csv.close()


def write_decomp_output(app):
    """ Write output file for binding free energy decomposition calculations """
    from csv import writer
    from datetime import datetime

    FILES = app.FILES
    INPUT = app.INPUT

    stability = app.stability
    outkeys = ('gb', 'pb')
    triggers = ('gbrun', 'pbrun')
    headers = ('Generalized Born Decomposition Energies', 'Poisson Boltzmann Decomposition Energies')

# First open up our output file and turn it into a CSV writer if necessary
    if INPUT['csv_format']:
        dec_out_file = OutputFile(FILES.decompout, 'wb', 0)
        decompout = writer(dec_out_file)
        decompout.writerow(['| Run on %s' % datetime.now().ctime()])
        if INPUT['gb']['gbrun']:
            decompout.writerow(['| GB non-polar solvation energies calculated with gbsa=2'])
    else:
        dec_out_file = OutputFile(FILES.decompout, 'w', 0)
        decompout = dec_out_file
        decompout.write_date()
        if INPUT['gb']['gbrun']:
            decompout.add_comment('| GB non-polar solvation energies calculated with gbsa=2')

    # Open up the CSV energy vector file
    if FILES.dec_energies:
        dec_energies_f = open(FILES.dec_energies, 'w')
        dec_energies = writer(dec_energies_f)

    for i, key in enumerate(outkeys):
        if not INPUT[triggers[i]]:
            continue
        if not INPUT['ala']['mutant_only']:
            if stability:
                decomp_norm = app.calc_types.decomp_normal[key]['complex']
            else:
                decomp_norm = app.calc_types.decomp_normal[key]['delta']

            if INPUT['csv_format']:
                decompout.writerows(decomp_norm.summary('csv'))
            else:
                decompout.writeline(decomp_norm.summary())

            # Now it's time to dump everything to the CSV file
            if FILES.dec_energies:
                dec_energies.writerow([headers[i] + '\n'])
                dec_energies.writerow(['Complex:'])
                app.calc_types.decomp_normal[key]['complex']._print_vectors(dec_energies)
                if not stability:
                    dec_energies.writerow(['\n'])
                    dec_energies.writerow(['Receptor:'])
                    app.calc_types.decomp_normal[key]['receptor']._print_vectors(dec_energies)
                    dec_energies.writerow(['Ligand:'])
                    app.calc_types.decomp_normal[key]['ligand']._print_vectors(dec_energies)
                    dec_energies.writerow(['DELTAS:'])
                    app.calc_types.decomp_normal[key]['delta']._print_vectors(dec_energies)
                dec_energies.writerow('\n')
        # Mutant system
        if INPUT['ala']['alarun']:
            if stability:
                decomp_mut = app.calc_types.decomp_mutant[key]['complex']
            else:
                decomp_mut = app.calc_types.decomp_mutant[key]['delta']

            # Write the data to the output file
            if INPUT['csv_format']:
                decompout.writerows(decomp_mut.summary('csv'))
            else:
                decompout.writeline(decomp_mut.summary())
            # Now it's time to dump everything to the CSV file
            if FILES.dec_energies:
                dec_energies.writerow([headers[i] + '(%s mutant)\n' % app.mut_str])
                dec_energies.writerow(['Complex:'])
                app.calc_types.decomp_mutant[key]['complex']._print_vectors(dec_energies)
                if not stability:
                    dec_energies.writerow(['\n'])
                    dec_energies.writerow(['Receptor:'])
                    app.calc_types.decomp_mutant[key]['receptor']._print_vectors(dec_energies)
                    dec_energies.writerow(['\n'])
                    dec_energies.writerow(['Ligand:'])
                    app.calc_types.decomp_mutant[key]['ligand']._print_vectors(dec_energies)
                    dec_energies.writerow(['\n'])
                    dec_energies.writerow(['DELTAS:'])
                    app.calc_types.decomp_mutant[key]['delta']._print_vectors(dec_energies)
                dec_energies.writerow(['\n'])

    # Close the file(s)
    if FILES.dec_energies:
        dec_energies_f.close()
    dec_out_file.close()


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
                self.writeline(f'{"|Ligand Structure file:":40}{FILES.ligand_tpr:>40}')
            if FILES.ligand_top:
                self.writeline(f'{"|Ligand (GROMACS) topology file:":40}{FILES.ligand_top:>40}')
            self.writeline(f'{"|Ligand (AMBER) topology file:":40}{FILES.ligand_prmtop:>40}')

        if INPUT['ala']['alarun']:
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
