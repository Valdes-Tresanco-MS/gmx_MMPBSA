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

from types import SimpleNamespace
from GMXMMPBSA import utils
from math import sqrt, ceil
from os import linesep as ls
import pickle


def data2pkl(app):
    info_file = SimpleNamespace(
        INPUT=app.INPUT,
        FILES=app.FILES,
        size=app.mpi_size,
        numframes=app.numframes,
        numframes_nmode=app.numframes_nmode,
        mutant_index=app.mutant_index,
        mut_str=app.mut_str,
        using_chamber=app.using_chamber,
        input_file=app.input_file_text,
        COM_PDB=''.join(open(app.FILES.complex_fixed).readlines()),
        output_file=''.join(open(app.FILES.output_file).readlines()),
        decomp_output_file=''.join(open(app.FILES.decompout).readlines()) if app.INPUT['decomprun']
        else None
    )

    with open('COMPACT_MMXSA_RESULTS.mmxsa', "wb") as f:
        pickle.dump(info_file, f)
        pickle.dump(app.calc_types, f)


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
        final_output.add_comment('Receptor mask:                  "%s"' % INPUT['receptor_mask'])
        final_output.add_comment('Ligand mask:                    "%s"' % INPUT['ligand_mask'])
        if prmtop_system.ligand_prmtop.ptr('nres') == 1:
            final_output.add_comment('Ligand residue name is:         "%s"' %
                                     prmtop_system.ligand_prmtop.parm_data['RESIDUE_LABEL'][0])
    final_output.add_comment('')
    final_output.add_comment('Calculations performed using %s complex frames' % app.numframes)
    if INPUT['nmoderun']:
        final_output.add_comment('NMODE calculations performed using %s frames' % app.numframes_nmode)
    if not stability:
        if INPUT['interaction_entropy']:
            final_output.add_comment('Interaction Entropy calculations performed using last %s frames' %
                                     ceil(app.numframes * (INPUT['ie_segment'] / 100)))
        if INPUT['c2_entropy']:
            final_output.add_comment('C2 Entropy Std. Dev. and Conf. Interv. (95%) have been obtained by '
                                     'bootstrapping with number of re-samplings = 2000')
    if INPUT['pbrun']:
        if INPUT['sander_apbs']:
            final_output.add_comment('Poisson Boltzmann calculations performed using iAPBS interface to sander '
                                     '(sander.APBS)')
        else:
            final_output.add_comment('Poisson Boltzmann calculations performed using internal PBSA solver in sander')
    final_output.add_comment('')

    if INPUT['gbrun']:
        if INPUT['molsurf']:
            final_output.add_comment('Generalized Born ESURF calculated using \'molsurf\' surface areas')
        else:
            final_output.add_comment('Generalized Born ESURF calculated using \'LCPO\' surface areas')
        final_output.add_comment('')

    final_output.add_comment('Using temperature = %.2f K)' % INPUT['temperature'])
    final_output.add_comment('All units are reported in kcal/mol')
    final_output.add_comment('')
    final_output.add_comment('SD - Sample standard deviation, SEM - Sample standard error of the mean')
    final_output.add_comment('SD(Prop.), SEM(Prop.) - SD and SEM obtained with propagation of uncertainty formula')
    final_output.add_comment('https://en.wikipedia.org/wiki/Propagation_of_uncertainty#Example_formulae')
    final_output.add_comment('')

    if INPUT['ifqnt']:
        final_output.add_comment(('QM/MM: Residues %s are treated with the ' +
                                  'Quantum Hamiltonian %s') % (INPUT['qm_residues'], INPUT['qm_theory']))
    final_output.separate()

    # First do the entropies
    if INPUT['qh_entropy']:
        if not INPUT['mutant_only']:
            if stability:
                qhnorm = app.calc_types.normal['qh']['complex']
            else:
                qhnorm = app.calc_types.normal['qh']['delta']
            final_output.writeline(f'Normal [ -TΔS ]')
            final_output.writeline('ENTROPY RESULTS (QH APPROXIMATION):')
            final_output.add_section(qhnorm.summary_output())
        if INPUT['alarun']:
            if stability:
                qhnorm = app.calc_types.mutant['qh']['complex']
            else:
                qhnorm = app.calc_types.mutant['qh']['delta']
            qhmutant = app.calc_types.mutant['qh']
            final_output.writeline(mut_str + ' Mutant [ -TΔS ]')
            final_output.writeline('ENTROPY RESULTS (QH APPROXIMATION):')
            final_output.add_section(qhmutant.summary_output())
        if INPUT['alarun'] and not INPUT['mutant_only']:
            qhmut_norm = app.calc_types.mut_norm['qh']
            final_output.writeline(f'Delta ( Mutant - Normal ) [ Δ(-TΔS) ]')
            final_output.writeline('ENTROPY RESULTS (QH APPROXIMATION):')
            final_output.add_section(qhmut_norm.summary_output())

    if not stability:
        # end if INPUT['entropy']
        if INPUT['interaction_entropy']:
            ie_inconsistent = False
            if not INPUT['mutant_only']:
                ienorm = app.calc_types.normal['ie']
                final_output.writeline('Normal [ -TΔS ]')
                final_output.writeline('ENTROPY RESULTS (INTERACTION ENTROPY):')
                final_output.add_section(ienorm.summary_output())
                # Now dump the energy vectors in CSV format
                if FILES.energyout:
                    energyvectors.writerow(['Interaction entropy results'])
                    ienorm._print_vectors(energyvectors)
                    energyvectors.writerow([])
                ie_inconsistent = ienorm['sigma'] > 3.6

            if INPUT['alarun']:
                iemutant = app.calc_types.mutant['ie']
                final_output.writeline(mut_str + ' Mutant [ -TΔS ]')
                final_output.writeline('ENTROPY RESULTS (INTERACTION ENTROPY):')
                final_output.add_section(iemutant.summary_output())
                ie_inconsistent = iemutant['sigma'] > 3.6

            if INPUT['alarun'] and not INPUT['mutant_only']:
                iemut_norm = app.calc_types.mut_norm['ie']
                final_output.writeline(f'Delta ( Mutant - Normal ) [ Δ(-TΔS) ]')
                final_output.writeline('ENTROPY RESULTS (INTERACTION ENTROPY):')
                final_output.add_section(iemut_norm.summary_output())
            if ie_inconsistent:
                final_output.writeline(
                    'WARNING: THE INTERACTION ENERGY STANDARD DEVIATION [ σ(Int. Energy) ]\n'
                    'IS GREATER THAN 3.6 kcal/mol (~15 kJ/mol). THUS, THE INTERACTION ENTROPY VALUES ARE\n'
                    'NOT RELIABLE. CHECK THIS PAPER FOR MORE INFO (https://doi.org/10.1021/acs.jctc.1c00374)\n\n')

        if INPUT['c2_entropy']:
            c2_inconsistent = False
            if not INPUT['mutant_only']:
                c2norm = app.calc_types.normal['c2']
                final_output.writeline('Normal [ -TΔS ]')
                final_output.writeline('ENTROPY RESULTS (C2 ENTROPY):')
                final_output.add_section(c2norm.summary_output())
                c2_inconsistent = c2norm['sigma'] > 3.6

            if INPUT['alarun']:
                c2mutant = app.calc_types.mutant['c2']
                final_output.writeline(mut_str + ' Mutant [ -TΔS ]')
                final_output.writeline('ENTROPY RESULTS (C2 ENTROPY):')
                final_output.add_section(c2mutant.summary_output())
                c2_inconsistent = c2mutant['sigma'] > 3.6

            if INPUT['alarun'] and not INPUT['mutant_only']:
                c2mut_norm = app.calc_types.mut_norm['c2']
                final_output.writeline(f'Delta ( Mutant - Normal ) [ Δ(-TΔS) ]')
                final_output.writeline('ENTROPY RESULTS (C2 ENTROPY):')
                final_output.add_section(c2mut_norm.summary_output())

            if c2_inconsistent:
                final_output.writeline(
                    'WARNING: THE INTERACTION ENERGY STANDARD DEVIATION [ σ(Int. Energy)]\n'
                    'IS GREATER THAN 3.6 kcal/mol (~15 kJ/mol). THUS, THE C2 ENTROPY VALUES ARE NOT\n'
                    'RELIABLE. CHECK THIS PAPER FOR MORE INFO (https://doi.org/10.1021/acs.jctc.1c00374)\n')

    # Now print out the normal mode results
    if INPUT['nmoderun']:
        if not INPUT['mutant_only']:
            if stability:
                nm_sys_norm = app.calc_types.normal['nmode']['complex']
            else:
                nm_sys_norm = app.calc_types.normal['nmode']['delta']
            final_output.writeline(f'Normal [ -TΔS ]')
            final_output.writeline('ENTROPY RESULTS (NMODE APPROXIMATION):')
            final_output.add_section(nm_sys_norm.summary_output())
            # Now dump the energy vectors in CSV format
            if FILES.energyout:
                energyvectors.writerow(['NMODE entropy results'])
                nm_sys_norm._print_vectors(energyvectors)
                energyvectors.writerow([])

        if INPUT['alarun']:
            if stability:
                nm_sys_mut = app.calc_types.mutant['nmode']['complex']
            else:
                nm_sys_mut = app.calc_types.mutant['nmode']['delta']
            final_output.writeline(mut_str + ' Mutant [ -TΔS ]')
            final_output.writeline('ENTROPY RESULTS (NMODE APPROXIMATION):')
            final_output.add_section(nm_sys_mut.summary_output())
            # Now dump the energy vectors in CSV format
            if FILES.energyout:
                energyvectors.writerow([mut_str + ' Mutant NMODE entropy results'])
                nm_sys_mut._print_vectors(energyvectors)
                energyvectors.writerow([])

        # Now calculate the effect of alanine scanning
        if INPUT['alarun'] and not INPUT['mutant_only']:
            nm_sys_mut_norm = app.calc_types.mut_norm['nmode']
            final_output.writeline(f'Delta ( Mutant - Normal ) [ Δ(-TΔS) ]')
            final_output.writeline('ENTROPY RESULTS (NMODE APPROXIMATION):')
            final_output.add_section(nm_sys_mut_norm.summary_output())

    # end if INPUT['nmoderun']

    triggers = ('gbrun', 'pbrun', 'rismrun_std', 'rismrun_gf', 'rismrun_pcplus')
    outkeys = ('gb', 'pb', 'rism std', 'rism gf', 'rism pcplus')
    headers = ('\nGENERALIZED BORN:\n\n', '\nPOISSON BOLTZMANN:\n\n',
               '\n3D-RISM:\n\n', '\n3D-RISM (Gauss. Fluct.):\n\n', '\n3D-RISM (PC+):\n\n')
    # Now print out the Free Energy results
    for i, key in enumerate(outkeys):
        if triggers[i] not in INPUT or not INPUT[triggers[i]]:
            continue
        if not INPUT['mutant_only']:
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
            if INPUT['qh_entropy']:
                dg_qh_davg, dg_qh_dstd = utils.calc_sum(sys_norm['TOTAL'], qhnorm['TOTAL'])
                final_output.add_section('Using Quasi-harmonic Entropy Approximation:\n'
                                         f'ΔG{"" if stability else " binding"} = {dg_qh_davg:9.2f} +'
                                         f'/- {dg_qh_dstd:7.2f}\n')
            if not stability:
                if INPUT['interaction_entropy']:
                    dg_ie_davg, dg_ie_dstd = utils.calc_sum(sys_norm['TOTAL'], ienorm['iedata'])
                    final_output.add_section(f"Using Interaction Entropy Approximation:\n"
                                             f"ΔG binding = {dg_ie_davg:9.2f} +/- {dg_ie_dstd:7.2f}\n")
                if INPUT['c2_entropy']:
                    dg_c2_davg, dh_dstd = utils.calc_sum(sys_norm['TOTAL'], c2norm['c2data'])
                    dg_c2_dstd = utils.get_std(dh_dstd, c2norm['c2_std'])
                    final_output.add_section(f"Using C2 Entropy Approximation:\n"
                                             f"ΔG binding = {dg_c2_davg:9.2f} +/- {dg_c2_dstd:7.2f}\n")
            if INPUT['nmoderun']:
                dg_nm_davg, dg_nm_dstd = utils.calc_sum(sys_norm['TOTAL'], nm_sys_norm['TOTAL'])
                final_output.add_section('Using Normal Mode Entropy Approximation:\n'
                                         f'ΔG{"" if stability else " binding"} = {dg_nm_davg:9.2f} +/-'
                                         f' {dg_nm_dstd:7.2f}\n')

        if INPUT['alarun']:
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

            if INPUT['qh_entropy']:
                mqh_davg, mqh_dstd = utils.calc_sum(sys_mut['TOTAL'], qhmutant['TOTAL'])
                final_output.add_section('Using Quasi-harmonic Entropy Approximation:\n'
                                         f'ΔG{"" if stability else " binding"} = {mqh_davg:9.2f} +/- {mqh_dstd:7.2f}\n')
            if not stability:
                if INPUT['interaction_entropy']:
                    mie_davg, mie_dstd = utils.calc_sum(sys_mut['TOTAL'], iemutant['iedata'])
                    final_output.add_section(f"Using Interaction Entropy Approximation:\n"
                                             f"ΔG binding = {mie_davg:9.2f} +/- {mie_dstd:7.2f}\n")
                if INPUT['c2_entropy']:
                    mc2_davg, mc2_dstd = utils.calc_sum(sys_mut['TOTAL'], c2mutant['c2data'])
                    final_output.add_section(f"Using C2 Entropy Approximation:\n"
                                             f"ΔG binding = {mc2_davg:9.2f} +/- {mc2_dstd:7.2f}\n")
            if INPUT['nmoderun']:
                mnm_davg, mnm_dstd = utils.calc_sum(sys_mut['TOTAL'], nm_sys_mut['TOTAL'])
                final_output.add_section('Using Normal Mode Entropy Approximation:\n'
                                         f'ΔG{"" if stability else " binding"} = {mnm_davg:9.2f} +/- {mnm_dstd:7.2f}\n')

        if INPUT['alarun'] and not INPUT['mutant_only']:
            mut_norm = app.calc_types.mut_norm[key]['delta']
            final_output.add_section(mut_norm.summary_output())
            ddh_davg = mut_norm['TOTAL'].mean()
            ddh_dstd = mut_norm['TOTAL'].std()

            final_output.write(f'\nRESULT OF ALANINE SCANNING ({mut_str}):\n' 
                               f"ΔΔH binding = {ddh_davg:9.2f} +/- {ddh_dstd:7.2f}\n")

            if INPUT['qh_entropy']:
                ddgqh_davg = ddh_davg + qhmut_norm.mean()
                # this std is the same of ΔΔH
                final_output.write('\n   (quasi-harmonic entropy)\n'
                                   f'ΔΔG{"" if stability else " binding"} = {ddgqh_davg:9.2f} +/- {ddh_dstd:7.2f}\n')
            if not stability:
                if INPUT['interaction_entropy']:
                    ddie_davg, ddie_dstd = iemut_norm['iedata'].mean(), iemut_norm['iedata'].std()
                    ddgie_davg = ddh_davg + ddie_davg
                    ddgie_dstd = utils.get_std(ddh_dstd, ddie_dstd)
                    final_output.write('\n   (interaction entropy)\n'
                                       f'ΔΔG binding = {ddgie_davg:9.2f} +/- {ddgie_dstd:7.2f}\n')
                if INPUT['c2_entropy']:
                    ddc2_davg, ddc2_dstd = c2mut_norm['c2data'], c2mut_norm['c2_std']
                    ddgc2_davg = ddh_davg + ddc2_davg
                    ddgc2_dstd = utils.get_std(ddh_dstd, ddc2_dstd)
                    final_output.write('\n   (C2 entropy)\n'
                                       f'ΔΔG binding = {ddgc2_davg:9.2f} +/- {ddgc2_dstd:7.2f}\n')
            if INPUT['nmoderun']:
                ddnm_davg, ddnm_dstd = nm_sys_mut_norm.mean(), nm_sys_mut_norm.std()
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
        if INPUT['gbrun']:
            decompout.writerow(['| GB non-polar solvation energies calculated with gbsa=2'])
    else:
        dec_out_file = OutputFile(FILES.decompout, 'w', 0)
        decompout = dec_out_file
        decompout.write_date()
        if INPUT['gbrun']:
            decompout.add_comment('| GB non-polar solvation energies calculated with gbsa=2')

    # Open up the CSV energy vector file
    if FILES.dec_energies:
        dec_energies_f = open(FILES.dec_energies, 'w')
        dec_energies = writer(dec_energies_f)

    for i, key in enumerate(outkeys):
        if triggers[i] not in INPUT or not INPUT[triggers[i]]:
            continue
        if not INPUT['mutant_only']:
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
        if INPUT['alarun']:
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
