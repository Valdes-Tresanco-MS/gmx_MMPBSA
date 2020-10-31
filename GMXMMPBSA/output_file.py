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
from GMXMMPBSA.amber_outputs import EnergyVector
from GMXMMPBSA.exceptions import LengthError
from GMXMMPBSA import utils
from math import sqrt, ceil
from os import linesep as ls


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
    final_output.add_comment('Calculations performed using %s complex frames.' %
                             app.numframes)
    if INPUT['nmoderun']:
        final_output.add_comment('NMODE calculations performed using %s frames.' %
                                 app.numframes_nmode)

    if INPUT['pbrun']:
        if INPUT['sander_apbs']:
            final_output.add_comment('Poisson Boltzmann calculations performed ' +
                                     'using iAPBS interface to sander (sander.APBS)')
        elif INPUT['use_sander'] or INPUT['decomprun']:
            final_output.add_comment('Poisson Boltzmann calculations performed ' +
                                     'using internal PBSA solver in sander.')
        else:
            final_output.add_comment('Poisson Boltzmann calculations performed ' +
                                     'using internal PBSA solver in mmpbsa_py_energy')

    final_output.add_comment('')

    if INPUT['gbrun']:
        if INPUT['molsurf']:
            final_output.add_comment('Generalized Born ESURF calculated using '
                                     '\'molsurf\' surface areas')
        else:
            final_output.add_comment('Generalized Born ESURF calculated using '
                                     '\'LCPO\' surface areas')
        final_output.add_comment('')

    final_output.add_comment('All units are reported in kcal/mole.')
    if INPUT['nmoderun'] or INPUT['entropy']:
        if INPUT['entropy'] == 2:
            final_output.add_comment('All entropy results have units kcal/mol ' +
                                 '(Temperature is %.2f K).' % INPUT['entropy_temp'])
        else:
            final_output.add_comment('All entropy results have units kcal/mol ' +
                                     '(Temperature is %.2f K).' % INPUT['temp'])
    if INPUT['ifqnt']:
        final_output.add_comment(('QM/MM: Residues %s are treated with the ' +
                                  'Quantum Hamiltonian %s') % (INPUT['qm_residues'], INPUT['qm_theory']))
    final_output.separate()

    # Start with entropies
    if INPUT['entropy'] == 1:
        if not INPUT['mutant_only']:
            qhnorm = app.calc_types['qh']
            final_output.writeline('ENTROPY RESULTS (QUASI-HARMONIC ' +
                                   'APPROXIMATION) CALCULATED WITH PTRAJ:')
            final_output.add_section(qhnorm.print_summary())
        if INPUT['alarun']:
            qhmutant = app.calc_types['mutant']['qh']
            final_output.writeline(mut_str + ' MUTANT')
            final_output.writeline('ENTROPY RESULTS (QUASI-' +
                                   'HARMONIC APPROXIMATION) CALCULATED WITH PTRAJ:')
            final_output.add_section(qhmutant.print_summary())
        if INPUT['alarun'] and not INPUT['mutant_only']:
            final_output.add_section(('\nRESULT OF ALANINE SCANNING: (%s) DELTA ' +
                                      'DELTA S binding = %9.4f\n') % (mut_str,
                                                                      qhnorm.total_avg() - qhmutant.total_avg()))
    # end if INPUT['entropy']

    # Now print out the normal mode results
    if INPUT['nmoderun']:
        if not INPUT['mutant_only']:
            nm_norm = app.calc_types['nmode']['complex']
            final_output.write('ENTROPY RESULTS (HARMONIC APPROXIMATION) ' +
                               'CALCULATED WITH NMODE:\n')
            final_output.add_section('Complex:\n' + nm_norm.print_summary())
            # Now dump the energy vectors in CSV format
            if FILES.energyout:
                energyvectors.writerow(['NMODE entropy results'])
                nm_norm.print_vectors(energyvectors)
                energyvectors.writerow([])

        if INPUT['alarun']:
            nm_mut = app.calc_types['mutant']['nmode']['complex']
            final_output.write(mut_str + ' MUTANT\nENTROPY RESULTS (HARMONIC ' +
                               'APPROXIMATION) CALCULATED WITH NMODE:\n')
            final_output.add_section('Complex:\n' + nm_mut.print_summary())
            # Now dump the energy vectors in CSV format
            if FILES.energyout:
                energyvectors.writerow([mut_str + ' Mutant NMODE entropy results'])
                nm_mut.print_vectors(energyvectors)
                energyvectors.writerow([])

        # Now calculate the effect of alanine scanning
        if INPUT['alarun'] and not INPUT['mutant_only']:
            try:
                nm_diff_array = nm_norm.data['Total'] - nm_mut.data['Total']
                nm_davg = nm_diff_array.avg()
                nm_dstdev = nm_diff_array.stdev()
            except LengthError:
                nm_diff_array = []
                nm_davg = nm_norm.data['Total'].avg() - nm_mut.data['Total'].avg()
                nm_dstdev = sqrt(nm_norm.data['Total'].stdev() ** 2 +
                                 nm_mut.data['Total'].stdev() ** 2)
            final_output.add_section(('\nRESULT OF ALANINE SCANNING: (%s) DELTA ' +
                                      'DELTA S binding = %9.4f +/- %9.4f\n') % (mut_str,
                                                                                nm_davg, nm_dstdev))

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
            if INPUT['entropy'] == 1:
                final_output.add_section('Using Quasi-harmonic Entropy ' +
                                         'Approximation: FREE ENERGY (G) = %9.4f\n' %
                                         (com_norm.data['TOTAL'].avg() - qhnorm.total_avg()))
            if INPUT['nmoderun']:
                if len(com_norm.data['TOTAL']) != len(nm_norm.data['Total']):
                    davg = (com_norm.data['TOTAL'].avg() -
                            nm_norm.data['Total'].avg())
                    dstdev = sqrt(com_norm.data['TOTAL'].stdev() ** 2 +
                                  nm_norm.data['Total'].stdev() ** 2)
                else:
                    diff_array = com_norm.data['TOTAL'] - nm_norm.data['Total']
                    davg = diff_array.avg()
                    dstdev = diff_array.stdev()

                final_output.add_section('Using Normal Mode Entropy Approxima' +
                                         'tion: FREE ENERGY (G) =    %9.4f +/- %7.4f\n' % (davg, dstdev))

        if INPUT['alarun']:
            com_mut = app.calc_types['mutant'][key]['complex']
            final_output.write('%s MUTANT:%s' % (mut_str, headers[i]))
            final_output.add_section('Complex:\n' + com_mut.print_summary())
            # Dump energy vectors to a CSV
            if FILES.energyout:
                energyvectors.writerow([mut_str + ' Mutant ' + headers[i]])
                com_mut.print_vectors(energyvectors)
                energyvectors.writerow([])

            # Combine with the entropy(ies)
            if INPUT['entropy'] == 1:
                final_output.add_section('Using Quasi-harmonic Entropy ' +
                                         'Approximation: DELTA G binding = %9.4f\n' %
                                         (com_mut.data['TOTAL'].avg() - qhmutant.total_avg()))
            if INPUT['nmoderun']:
                if len(com_mut.data['TOTAL']) != len(nm_mut.data['Total']):
                    davg = (com_mut.data['TOTAL'].avg() -
                            nm_mut.data['Total'].avg())
                    dstdev = sqrt(com_mut.data['TOTAL'].stdev() ** 2 +
                                  nm_mut.data['Total'].stdev() ** 2)
                else:
                    diff_array = com_mut.data['TOTAL'] - nm_mut.data['Total']
                    davg = diff_array.avg()
                    dstdev = diff_array.stdev()

                final_output.add_section('Using Normal Mode Entropy Approxima' +
                                         'tion: FREE ENERGY (G) =    %9.4f +/- %7.4f\n' % (davg, dstdev))

        if INPUT['alarun'] and not INPUT['mutant_only']:
            diff_array = com_norm.data['TOTAL'] - com_mut.data['TOTAL']
            davg = diff_array.avg()
            dstdev = diff_array.stdev()
            final_output.write(('\nRESULT OF ALANINE SCANNING: (%s) DELTA ' +
                                'G = %9.4f  +/- %9.4f\n') % (mut_str, davg, dstdev))
            if INPUT['entropy'] == 1:
                final_output.write(('\n   (quasi-harmonic entropy) (%s) DELTA ' +
                                    'G = %9.4f\n') % (mut_str,
                                                      davg + qhnorm.total_avg() - qhmutant.total_avg()))
            if INPUT['nmoderun']:
                if len(nm_diff_array) == len(diff_array):
                    total_diffs = diff_array - nm_diff_array
                    davg1 = total_diffs.avg()
                    dstdev1 = total_diffs.stdev()
                else:
                    davg1 = davg - nm_davg
                    dstdev1 = sqrt(dstdev ** 2 - nm_dstdev ** 2)
                final_output.write(('\n      (normal mode entropy) (%s) DELTA ' +
                                    'G = %9.4f +/- %9.4f\n') % (mut_str, davg1, dstdev1))
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
    final_output.add_comment('Receptor mask:                  "%s"' %
                             INPUT['receptor_mask'])
    final_output.add_comment('Ligand mask:                    "%s"' %
                             INPUT['ligand_mask'])
    if prmtop_system.ligand_prmtop.ptr('nres') == 1:
        final_output.add_comment('Ligand residue name is "%s"' %
                                 prmtop_system.ligand_prmtop.parm_data['RESIDUE_LABEL'][0])
    final_output.add_comment('')
    final_output.add_comment('Calculations performed using %s complex frames.' %
                             app.numframes)
    if INPUT['nmoderun']:
        final_output.add_comment('NMODE calculations performed using %s frames.' %
                                 app.numframes_nmode)

    if INPUT['pbrun']:
        if INPUT['sander_apbs']:
            final_output.add_comment('Poisson Boltzmann calculations performed ' +
                                     'using iAPBS interface to sander (sander.APBS)')
        elif INPUT['use_sander'] or INPUT['decomprun']:
            final_output.add_comment('Poisson Boltzmann calculations performed ' +
                                     'using internal PBSA solver in sander.')
        else:
            final_output.add_comment('Poisson Boltzmann calculations performed ' +
                                     'using internal PBSA solver in mmpbsa_py_energy')
    final_output.add_comment('')

    if INPUT['gbrun']:
        if INPUT['molsurf']:
            final_output.add_comment('Generalized Born ESURF calculated using '
                                     '\'molsurf\' surface areas')
        else:
            final_output.add_comment('Generalized Born ESURF calculated using '
                                     '\'LCPO\' surface areas')
        final_output.add_comment('')

    final_output.add_comment('All units are reported in kcal/mole.')
    if INPUT['nmoderun'] or INPUT['entropy']:
        if INPUT['entropy'] == 2:
            final_output.add_comment('All entropy results have units kcal/mol ' +
                                     '(Temperature is %.2f K).' % INPUT['entropy_temp'])
        else:
            final_output.add_comment('All entropy results have units kcal/mol ' +
                                     '(Temperature is %.2f K).' % INPUT['temp'])
    if INPUT['ifqnt']:
        final_output.add_comment(('QM/MM: Residues %s are treated with the ' +
                                  'Quantum Hamiltonian %s') % (INPUT['qm_residues'], INPUT['qm_theory']))
    final_output.separate()

    # First do the entropies
    if INPUT['entropy'] == 1:
        if not INPUT['mutant_only']:
            qhnorm = app.calc_types['qh']
            final_output.writeline('ENTROPY RESULTS (QUASI-HARMONIC '
                                   'APPROXIMATION) CALCULATED WITH PTRAJ:')
            final_output.add_section(qhnorm.print_summary())
        if INPUT['alarun']:
            qhmutant = app.calc_types['mutant']['qh']
            final_output.writeline(mut_str + ' MUTANT')
            final_output.writeline('ENTROPY RESULTS (QUASI-' +
                                   'HARMONIC APPROXIMATION) CALCULATED WITH PTRAJ:')
            final_output.add_section(qhmutant.print_summary())
        if INPUT['alarun'] and not INPUT['mutant_only']:
            final_output.add_section(('\nRESULT OF ALANINE SCANNING: (%s) DELTA ' +
                                      'DELTA S binding = %9.4f\n') % (mut_str,
                                                                      qhnorm.total_avg() - qhmutant.total_avg()))
    # end if INPUT['entropy']

    # Now print out the normal mode results
    if INPUT['nmoderun']:
        if not INPUT['mutant_only']:
            nm_sys_norm = app.calc_types['nmode']['delta']
            final_output.write('ENTROPY RESULTS (HARMONIC APPROXIMATION) ' +
                               'CALCULATED WITH NMODE:\n\n')
            final_output.add_section(nm_sys_norm.print_summary())
            # Now dump the energy vectors in CSV format
            if FILES.energyout:
                energyvectors.writerow(['NMODE entropy results'])
                nm_sys_norm.print_vectors(energyvectors)
                energyvectors.writerow([])

        if INPUT['alarun']:
            nm_sys_mut = app.calc_types['mutant']['nmode']['delta']
            final_output.write(mut_str + ' MUTANT\nENTROPY RESULTS (HARMONIC ' +
                               'APPROXIMATION) CALCULATED WITH NMODE:\n')
            final_output.add_section(nm_sys_mut.print_summary())
            # Now dump the energy vectors in CSV format
            if FILES.energyout:
                energyvectors.writerow([mut_str + ' Mutant NMODE entropy results'])
                nm_sys_mut.print_vectors(energyvectors)
                energyvectors.writerow([])

        # Now calculate the effect of alanine scanning
        if INPUT['alarun'] and not INPUT['mutant_only']:
            davg, dstdev = nm_sys_norm.diff(nm_sys_mut, 'Total', 'Total')
            final_output.add_section(('\nRESULT OF ALANINE SCANNING: (%s) DELTA ' +
                                      'DELTA S binding = %9.4f +/- %9.4f\n') % (mut_str, davg,
                                                                                dstdev))

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
            if INPUT['entropy'] == 1:
                if isinstance(sys_norm.data['DELTA TOTAL'], EnergyVector):
                    final_output.add_section('Using Quasi-harmonic Entropy ' +
                                             'Approximation: DELTA G binding = %9.4f\n' %
                                             (sys_norm.data['DELTA TOTAL'].avg() - qhnorm.total_avg()))
                else:
                    final_output.add_section('Using Quasi-harmonic Entropy ' +
                                             'Approximation: DELTA G binding = %9.4f\n' %
                                             (sys_norm.data['DELTA TOTAL'][0] - qhnorm.total_avg()))
            elif INPUT['entropy'] == 2:
                b = 0
                if len(sys_norm.data['-TDS']) * app.INPUT['entropy_seg']/100 > 1:
                    b = ceil(len(sys_norm.data['-TDS']) * (1 - app.INPUT['entropy_seg']/100))
                norm_ts = sys_norm.data['-TDS'][b:].mean()
                if isinstance(sys_norm.data['DELTA TOTAL'], EnergyVector):
                    final_output.add_section('Interaction Entropy (-TDS): %9.4f\n' % norm_ts +
                                             'Using Interaction Entropy ' +
                                             'Approximation: DELTA G binding = %9.4f\n' %
                                             (sys_norm.data['DELTA TOTAL'].avg() + norm_ts))
                else:
                    final_output.add_section('Interaction Entropy (-TDS): %9.4f\n' % norm_ts +
                                             'Using Interaction Entropy ' +
                                             'Approximation: DELTA G binding = %9.4f\n' %
                                             (sys_norm.data['DELTA TOTAL'][0] + norm_ts))
            if INPUT['nmoderun']:
                davg, dstdev = sys_norm.diff(nm_sys_norm, 'DELTA TOTAL', 'Total')
                final_output.add_section('Using Normal Mode Entropy Approxima' +
                                         'tion: DELTA G binding =    %9.4f +/- %7.4f\n' % (davg, dstdev))

        if INPUT['alarun']:
            sys_mut = app.calc_types['mutant'][key]['delta']
            final_output.write('%s MUTANT:%s' % (mut_str, headers[i]))
            final_output.add_section(sys_mut.print_summary())
            # Dump energy vectors to a CSV
            if FILES.energyout:
                energyvectors.writerow([mut_str + ' Mutant ' + headers[i]])
                sys_mut.print_vectors(energyvectors)
                energyvectors.writerow([])

            # Combine with the entropy(ies)
            if INPUT['entropy'] == 1:
                if isinstance(sys_mut.data['DELTA TOTAL'], EnergyVector):
                    final_output.add_section('Using Quasi-harmonic Entropy ' +
                                             'Approximation: DELTA G binding = %9.4f\n' %
                                             (sys_mut.data['DELTA TOTAL'].avg() - qhmutant.total_avg()))
                else:
                    final_output.add_section('Using Quasi-harmonic Entropy ' +
                                             'Approximation: DELTA G binding = %9.4f\n' %
                                             (sys_mut.data['DELTA TOTAL'][0] - qhmutant.total_avg()))
            elif INPUT['entropy'] == 2:
                b = 0
                if len(sys_mut.data['-TDS']) * app.INPUT['entropy_seg']/100 > 1:
                    b = ceil(len(sys_mut.data['-TDS']) * app.INPUT['entropy_seg']/100 > 1)
                mut_ts = sys_mut.data['-TDS'][b:].mean()
                if isinstance(sys_mut.data['DELTA TOTAL'], EnergyVector):
                    final_output.add_section('Interaction Entropy (-TDS): %9.4f\n' % mut_ts +
                                             'Using Interaction Entropy ' +
                                             'Approximation: DELTA G binding = %9.4f\n' %
                                             (sys_mut.data['DELTA TOTAL'].avg() + mut_ts))
                else:
                    final_output.add_section('Interaction Entropy (-TDS): %9.4f\n' % mut_ts +
                                             'Using Interaction Entropy ' +
                                             'Approximation: DELTA G binding = %9.4f\n' %
                                             (sys_mut.data['DELTA TOTAL'][0] + mut_ts))
            if INPUT['nmoderun']:
                davg, dstdev = sys_mut.diff(nm_sys_mut, 'DELTA TOTAL', 'Total')
                final_output.add_section('Using Normal Mode Entropy Approxima' +
                                         'tion: DELTA G binding =    %9.4f +/- %7.4f\n' % (davg, dstdev))

        if INPUT['alarun'] and not INPUT['mutant_only']:
            davg, dstdev = sys_norm.diff(sys_mut, 'DELTA TOTAL', 'DELTA TOTAL')
            final_output.write(('\nRESULT OF ALANINE SCANNING: (%s) DELTA ' +
                                'DELTA G binding = %9.4f  +/- %9.4f\n') % (mut_str, davg, dstdev))
            if INPUT['entropy'] == 1:
                final_output.write(('\n   (quasi-harmonic entropy) (%s) DELTA ' +
                                    'DELTA G binding = %9.4f\n') % (mut_str,
                                                                    davg + qhnorm.total_avg() - qhmutant.total_avg()))
            elif INPUT['entropy'] == 2:
                final_output.write(('\n   (interaction entropy) (%s) DELTA ' +
                                    'DELTA G binding = %9.4f\n') % (mut_str, davg + norm_ts - mut_ts))
            if INPUT['nmoderun']:
                davg1, dstdev1 = nm_sys_norm.diff(nm_sys_mut, 'Total', 'Total')
                final_output.write(('\n      (normal mode entropy) (%s) DELTA ' +
                                    'DELTA G binding = %9.4f +/- %9.4f\n') % (mut_str, davg - davg1,
                                                                              sqrt(dstdev ** 2 + dstdev1 ** 2)))
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

    if INPUT['idecomp'] in [1, 2]:
        DecompClass = DecompOut
    else:
        DecompClass = PairDecompOut

    # Open up the files
    if INPUT['csv_format']:
        dec_out_file = OutputFile(FILES.decompout, 'wb', 0)
        decompout = writer(dec_out_file)
        decompout.writerow(['| Run on %s' % datetime.now().ctime()])
        if INPUT['gbrun']:
            decompout.writerow(['| GB non-polar solvation energies calculated '
                                'with gbsa=2'])
        decompout.writerow([idecompString[INPUT['idecomp']]])
    else:
        dec_out_file = OutputFile(FILES.decompout, 'w', 0)
        decompout = dec_out_file
        decompout.write_date()
        if INPUT['gbrun']:
            decompout.add_comment('| GB non-polar solvation energies calculated '
                                  'with gbsa=2')
        decompout.writeline(idecompString[INPUT['idecomp']])

    # Open up the energyvector CSV file
    if FILES.dec_energies: dec_energies = open(FILES.dec_energies, 'w')

    # GB first
    if INPUT['gbrun']:
        # See if we're dumping all energy terms
        csv_prefix = ''
        if FILES.dec_energies: csv_prefix = pre + 'gb'

        # See if we do normal system
        if not INPUT['mutant_only']:
            gb_com = DecompClass(pre + 'complex_gb.mdout',
                                 prmtop_system.complex_prmtop, INPUT['surften'], csv_prefix,
                                 size, INPUT['dec_verbose'])
            gb_com.fill_all_terms()
            if INPUT['csv_format']:
                decompout.writerow(['Energy Decomposition Analysis (All units ' +
                                    'kcal/mol): Generalized Born Solvent'])
                decompout.writerow([])
                gb_com.write_summary_csv(gb_com.numframes, decompout)
                decompout.writerow([])
            else:
                decompout.writeline('Energy Decomposition Analysis (All units ' +
                                    'kcal/mol): Generalized Born Solvent')
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
            gb_com = DecompClass(pre + 'mutant_complex_gb.mdout',
                                 mutant_system.complex_prmtop, INPUT['surften'], csv_prefix,
                                 size, INPUT['dec_verbose'])
            gb_com.fill_all_terms()
            if INPUT['csv_format']:
                decompout.writerow(['Energy Decomposition Analysis (All units ' +
                                    'kcal/mol): Generalized Born Solvent (%s)' % mutstr])
                decompout.writerow([])
                gb_com.write_summary_csv(gb_com.numframes, decompout)
                decompout.writerow([])
            else:
                decompout.writeline('Energy Decomposition Analysis (All units ' +
                                    'kcal/mol): Generalized Born Solvent (%s)' % mutstr)
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
        # See if we're dumping all energy terms
        csv_prefix = ''
        if FILES.dec_energies: csv_prefix = pre + 'pb'

        # See if we do normal system
        if not INPUT['mutant_only']:
            pb_com = DecompClass(pre + 'complex_pb.mdout',
                                 prmtop_system.complex_prmtop, INPUT['surften'], csv_prefix,
                                 size, INPUT['dec_verbose'])
            pb_com.fill_all_terms()
            if INPUT['csv_format']:
                decompout.writerow(['Energy Decomposition Analysis (All units ' +
                                    'kcal/mol): Generalized Born Solvent'])
                decompout.writerow([])
                pb_com.write_summary_csv(pb_com.numframes, decompout)
                decompout.writerow([])
            else:
                decompout.writeline('Energy Decomposition Analysis (All units ' +
                                    'kcal/mol): Poisson Boltzmann Solvent')
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
            pb_com = DecompClass(pre + 'mutant_complex_pb.mdout',
                                 mutant_system.complex_prmtop, INPUT['surften'], csv_prefix,
                                 size, INPUT['dec_verbose'])
            pb_com.fill_all_terms()
            if INPUT['csv_format']:
                decompout.writerow(['Energy Decomposition Analysis (All units ' +
                                    'kcal/mol): Poisson Boltzmann Solvent (%s)' % mutstr])
                decompout.writerow([])
                pb_com.write_summary_csv(pb_com.numframes, decompout)
                decompout.writerow([])
            else:
                decompout.writeline('Energy Decomposition Analysis (All units ' +
                                    'kcal/mol): Generalized Born Solvent (%s)' % mutstr)
                decompout.writeline('')
                pb_com.write_summary(pb_com.numframes, decompout)
                decompout.writeline('')
            # Dump the energy vectors if requested
            if FILES.dec_energies:
                del pb_com.csvwriter
                dec_energies.write('PB Decomposition Energies (%s mutant)' %
                                   mutstr + ls)
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
    else:  # Multiple trajectories
        # Per-residue
        if INPUT['idecomp'] in [1, 2]:
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
            decompout.writerow(['| GB non-polar solvation energies calculated '
                                'with gbsa=2'])
    else:
        dec_out_file = OutputFile(FILES.decompout, 'w', 0)
        decompout = dec_out_file
        decompout.write_date()
        if INPUT['gbrun']:
            decompout.add_comment('| GB non-polar solvation energies calculated '
                                  'with gbsa=2')

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
            if csv_prefix: csv_pre = csv_prefix % 'com'
            gb_com = SingleClass(pre + 'complex_gb.mdout',
                                 prmtop_system.complex_prmtop, INPUT['surften'],
                                 csv_pre, size, INPUT['dec_verbose'])
            if csv_prefix: csv_pre = csv_prefix % 'rec'
            gb_rec = SingleClass(pre + 'receptor_gb.mdout',
                                 prmtop_system.receptor_prmtop, INPUT['surften'],
                                 csv_pre, size, INPUT['dec_verbose'])
            if csv_prefix: csv_pre = csv_prefix % 'lig'
            gb_lig = SingleClass(pre + 'ligand_gb.mdout',
                                 prmtop_system.ligand_prmtop, INPUT['surften'],
                                 csv_pre, size, INPUT['dec_verbose'])
            if csv_prefix: csv_pre = pre + 'gb_bind'
            gb_bind = BindingClass(gb_com, gb_rec, gb_lig, prmtop_system,
                                   INPUT['idecomp'], INPUT['dec_verbose'], decompout,
                                   csv_pre, 'Energy Decomposition Analysis ' +
                                   '(All units kcal/mol): Generalized Born solvent')
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
            if csv_prefix: csv_pre = csv_prefix % 'com'
            gb_com = SingleClass(pre + 'mutant_complex_gb.mdout',
                                 mutant_system.complex_prmtop, INPUT['surften'],
                                 csv_pre, size, INPUT['dec_verbose'])
            if csv_prefix: csv_pre = csv_prefix % 'rec'
            gb_rec = SingleClass(pre + 'mutant_receptor_gb.mdout',
                                 mutant_system.receptor_prmtop, INPUT['surften'],
                                 csv_pre, size, INPUT['dec_verbose'])
            if csv_prefix: csv_pre = csv_prefix % 'lig'
            gb_lig = SingleClass(pre + 'mutant_ligand_gb.mdout',
                                 mutant_system.ligand_prmtop, INPUT['surften'],
                                 csv_pre, size, INPUT['dec_verbose'])
            if csv_prefix: csv_pre = pre + 'gb_bind'
            gb_bind = BindingClass(gb_com, gb_rec, gb_lig, mutant_system,
                                   INPUT['idecomp'], INPUT['dec_verbose'], decompout,
                                   csv_pre, 'Energy Decomposition Analysis ' +
                                   '(All units kcal/mol): Generalized Born solvent (%s)' %
                                   mutstr)
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
            if csv_prefix: csv_pre = csv_prefix % 'com'
            pb_com = SingleClass(pre + 'complex_pb.mdout',
                                 prmtop_system.complex_prmtop, INPUT['cavity_surften'],
                                 csv_pre, size, INPUT['dec_verbose'])
            if csv_prefix: csv_pre = csv_prefix % 'rec'
            pb_rec = SingleClass(pre + 'receptor_pb.mdout',
                                 prmtop_system.receptor_prmtop, INPUT['cavity_surften'],
                                 csv_pre, size, INPUT['dec_verbose'])
            if csv_prefix: csv_pre = csv_prefix % 'lig'
            pb_lig = SingleClass(pre + 'ligand_pb.mdout',
                                 prmtop_system.ligand_prmtop, INPUT['cavity_surften'],
                                 csv_pre, size, INPUT['dec_verbose'])
            if csv_prefix: csv_pre = pre + 'pb_bind'
            pb_bind = BindingClass(pb_com, pb_rec, pb_lig, prmtop_system,
                                   INPUT['idecomp'], INPUT['dec_verbose'], decompout,
                                   csv_pre, 'Energy Decomposition Analysis ' +
                                   '(All units kcal/mol): Poisson Boltzmann solvent')
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
            if csv_prefix: csv_pre = csv_prefix % 'com'
            pb_com = SingleClass(pre + 'mutant_complex_pb.mdout',
                                 mutant_system.complex_prmtop, INPUT['cavity_surften'],
                                 csv_pre, size, INPUT['dec_verbose'])
            if csv_prefix: csv_pre = csv_prefix % 'rec'
            pb_rec = SingleClass(pre + 'mutant_receptor_pb.mdout',
                                 mutant_system.receptor_prmtop, INPUT['cavity_surften'],
                                 csv_pre, size, INPUT['dec_verbose'])
            if csv_prefix: csv_pre = csv_prefix % 'lig'
            pb_lig = SingleClass(pre + 'mutant_ligand_pb.mdout',
                                 mutant_system.ligand_prmtop, INPUT['cavity_surften'],
                                 csv_pre, size, INPUT['dec_verbose'])
            if csv_prefix: csv_pre = pre + 'pb_bind'
            pb_bind = BindingClass(pb_com, pb_rec, pb_lig, mutant_system,
                                   INPUT['idecomp'], INPUT['dec_verbose'], decompout,
                                   csv_pre, 'Energy Decomposition Analysis ' +
                                   '(All units kcal/mol): Poisson Boltzmann solvent (%s)' %
                                   mutstr)
            # Write the data to the output file
            pb_bind.parse_all()
            if FILES.dec_energies:
                del pb_com.csvwriter, pb_rec.csvwriter, pb_lig.csvwriter
                del pb_bind.csvwriter
                dec_energies.write('PB Decomposition Energies (%s mutant)' %
                                   mutstr + ls)
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
        if self._binary:
            if isinstance(stuff, str):
                self._handle.write(stuff.encode())
            else:
                self._handle.write(stuff)
        else:
            if isinstance(stuff, str):
                self._handle.write(stuff)
            else:
                self._handle.write(stuff.decode())
        self._handle.flush()

    # ==================================================

    def separate(self):
        """ Delimiter between fields in output file """
        for i in range(2):
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
        import GMXMMPBSA
        stability = not FILES.receptor_prmtop

        self.writeline('|gmx_MMPBSA Version=%s based on MMPBSA.py v.%s' % (GMXMMPBSA.__version__,
                                                                              GMXMMPBSA.__mmpbsa_version__))

        # if FILES.solvated_prmtop:
        #    self.writeline('|Solvated complex topology file:  %s' %
        #               FILES.solvated_prmtop)
        self.writeline('|Complex topology file:           %s' %
                       FILES.complex_prmtop)

        if not stability:
            # if FILES.receptor_trajs:
            #    if FILES.solvated_receptor_prmtop:
            #       self.writeline('|Solvated receptor topology file: %s' %
            #                  FILES.solvated_receptor_prmtop)
            self.writeline('|Receptor topology file:          %s' %
                           FILES.receptor_prmtop)

            # if FILES.ligand_trajs:
            # if FILES.solvated_ligand_prmtop:
            #    self.writeline('|Solvated ligand topology file:   %s' %
            #               FILES.solvated_ligand_prmtop)
            self.writeline('|Ligand topology file:            %s' %
                           FILES.ligand_prmtop)

        if INPUT['alarun']:
            self.writeline('|Mutant complex topology file:    %s' %
                           FILES.mutant_complex_prmtop)
            if not stability:
                self.writeline('|Mutant receptor topology file:   %s' %
                               FILES.mutant_receptor_prmtop)
                self.writeline('|Mutant ligand topology file:     %s' %
                               FILES.mutant_ligand_prmtop)

        self.write('|Initial mdcrd(s):                ')
        for i in range(len(FILES.complex_trajs)):
            if i == 0:
                self.writeline(FILES.complex_trajs[i])
            else:
                self.writeline('|                                 %s' %
                               FILES.complex_trajs[i])

        if FILES.receptor_trajs:
            self.write('|Initial Receptor mdcrd(s):       ')
            for i in range(len(FILES.receptor_trajs)):
                if i == 0:
                    self.writeline(FILES.receptor_trajs[i])
                else:
                    self.writeline('|                                 %s' %
                                   FILES.receptor_trajs[i])

        if FILES.ligand_trajs:
            self.write('|Initial Ligand mdcrd(s):         ')
            for i in range(len(FILES.ligand_trajs)):
                if i == 0:
                    self.writeline(FILES.ligand_trajs[i])
                else:
                    self.writeline('|                                 %s' %
                                   FILES.ligand_trajs[i])

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
