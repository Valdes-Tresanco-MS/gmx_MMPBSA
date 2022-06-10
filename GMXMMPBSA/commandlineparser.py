# -*- coding: utf-8 -*-
"""
This module contains classes and such that are responsible for parsing
command-line arguments for gmx_MMPBSA.  All of the files specified for use
in gmx_MMPBSA will be assigned as attributes to the returned class.

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

import os
import sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, RawDescriptionHelpFormatter, RawTextHelpFormatter
from pathlib import Path
from GMXMMPBSA import __version__, __mmpbsa_version__, __ambertools_version__
from GMXMMPBSA.exceptions import GMXMMPBSA_ERROR
from types import SimpleNamespace

class OptionList(SimpleNamespace):
    """
    Just a container to hold the command-line options. Necessary when reading in
    a gmx_MMPBSA info file to have a container to load the results from the
    parser.
    """
    pass


if os.getenv('AMBERHOME'):
    # rism_xvv = os.path.join(os.getenv('AMBERHOME'), 'dat', 'mmpbsa', 'spc.xvv')
    path_to_file = os.path.join(os.getenv('AMBERHOME'), 'AmberTools', 'test', 'rism1d', 'tip3p-kh', 'tip3p.xvv.save')
    if os.path.exists(path_to_file):
        rism_xvv = path_to_file
    else:
        rism_xvv = os.path.join(Path(__file__).parent.joinpath('data'), 'xvv_files', 'tip3p.xvv')
else:
    rism_xvv = os.path.join(Path(__file__).parent.joinpath('data'), 'xvv_files', 'tip3p.xvv')


def check_arg(str_suffix, path=False):
    def decorator(f):
        def gmx_MMPBSA_file(*args):
            result = f(*args)
            if path and not Path(result).exists():
                GMXMMPBSA_ERROR(f'{result} do not exist or is inaccessible')
            if Path(result).suffix not in str_suffix:
                GMXMMPBSA_ERROR(f'{result} does not correspond to the required structure format {str_suffix}')
            return result
        return gmx_MMPBSA_file
    return decorator


@check_arg(['.tpr', '.pdb'], True)
def structure(arg):
    return arg


@check_arg(['.pdb'], True)
def pdb(arg):
    return arg


@check_arg(['.xtc', '.trr', '.pdb'], True)
def trajectory(arg):
    return arg


@check_arg(['.top'], True)
def topology(arg):
    return arg


@check_arg(['.mol2'], True)
def mol2(arg):
    return arg


@check_arg(['.ndx'], True)
def index(arg):
    return arg


class GMXMMPBSA_ArgParser(ArgumentParser):

    def exit(self, status=0, message=None):
        if message:
            GMXMMPBSA_ERROR(message)
        sys.exit(status)


description = ("gmx_MMPBSA is a new tool based on AMBER's MMPBSA.py aiming to perform end-state free energy calculations"
               " with GROMACS files. This program is an adaptation of Amber's MMPBSA.py and essentially works as such. "
               "gmx_MMPBSA works with any GROMACS version. This program will calculate binding free energies using "
               "end-state free energy methods on an ensemble of snapshots using a variety of implicit solvent models.")

complex_group_des = ("Complex files and info that are needed to perform the calculation. If the receptor and/or the "
                     "ligand info is not defined, we generate them from that of the complex.")

receptor_group_des = ("Receptor files and info that are needed to perform the calculation. If the receptor info is "
                      "not defined, we generate it from that of the complex.")

ligand_group_des = '''Ligand files and info that are needed to perform the calculation. If the ligand are not defined, 
                       we generate it from that of the complex.'''
# Set up the MM/PBSA parser here. It clutters up the MMPBSA_App to do it there

# noinspection PyTypeChecker
parser = GMXMMPBSA_ArgParser(epilog=f'''gmx_MMPBSA is an effort to implement the GB/PB and others calculations in GROMACS.
                                    \nBased on MMPBSA.py (version {__mmpbsa_version__}) and 
                                    AmberTools{__ambertools_version__}''',
                        description=(description + '''This is the core of gmx_MMPBSA and it will do all the 
                                    calculations'''),
                        formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument('-v', '--version', action='version',
                    version='''%%(prog)s %s based on MMPBSA version %s and AmberTools %s''' %
                            (__version__, __mmpbsa_version__, __ambertools_version__))
parser.add_argument('--input-file-help', dest='infilehelp', action='store_true',
                    help='Print all available options in the input file.',
                    default=False)
parser.add_argument('--create_input', dest='createinput', choices=['gb', 'pb', 'rism', 'ala', 'decomp', 'nmode', 'all'],
                    nargs='*', help='Create an new input file with selected calculation type.')
group = parser.add_argument_group('Miscellaneous Options')
group.add_argument('-O', '--overwrite', default=False, action='store_true',
                   help='Allow output files to be overwritten', dest='overwrite')
group.add_argument('-prefix', dest='prefix', default='_GMXMMPBSA_',
                   metavar='<file prefix>',
                   help='Prefix for intermediate files.')
group = parser.add_argument_group('Input and Output Files', '''These options specify the input files and optional 
output files.''')
group.add_argument('-i', dest='input_file', metavar='FILE', help='MM/PBSA input file.')
group.add_argument('-xvvfile', dest='xvvfile', help='XVV file for 3D-RISM.', default=rism_xvv)
group.add_argument('-o', dest='output_file', default='FINAL_RESULTS_MMPBSA.dat', metavar='FILE',
                   help='Output file with MM/PBSA statistics.')
group.add_argument('-do', dest='decompout', metavar='FILE', default='FINAL_DECOMP_MMPBSA.dat',
                   help='Output file for decomposition statistics summary.')
group.add_argument('-eo', dest='energyout', metavar='FILE',
                   help='''CSV-format output of all energy terms for every frame
                  in every calculation. File name forced to end in [.csv].
                  This file is only written when specified on the
                  command-line.''')
group.add_argument('-deo', dest='dec_energies', metavar='FILE',
                   help='''CSV-format output of all energy terms for each printed
                  residue in decomposition calculations. File name forced to end
                  in [.csv]. This file is only written when specified on the
                  command-line.''')
group.add_argument('-nogui', dest='gui', action='store_false', default=True,
                   help='No open gmx_MMPBSA_ana after all calculations finished')
group.add_argument('-s', '--stability', dest='stability', action='store_true', default=False,
                   help='''Perform stability calculation. Only the complex parameters are required. If
                         ligand is non-Protein (small molecule) type, then ligand *.mol2 file is 
                         required. In any other case receptor and ligand parameters will be ignored.
                         See description bellow''')

group = parser.add_argument_group('Complex', complex_group_des)
group.add_argument('-cs', dest='complex_tpr', metavar='<Structure File>', default=None, type=structure,
                   help='''Structure file of the complex. If it is Protein-Ligand 
                         (small molecule) complex and -cp is not defined, make 
                         sure that you define -lm option. See -lm description 
                         below. Allowed formats: *.tpr (recommended), *.pdb''')
group.add_argument('-ci', dest='complex_index', metavar='<Index File>', default=None, type=index,
                   help='Index file of the bound complex.')
group.add_argument('-cg', dest='complex_groups', metavar='index', nargs=2, default=None, type=int,
                   help='Groups of receptor and ligand in complex index file. The notation is as follows: "-cg '
                        '<Receptor group> <Ligand group>", ie. -cg 1 13')
group.add_argument('-ct', dest='complex_trajs', nargs='*', metavar='TRJ', type=trajectory,
                   help='''Complex trajectories. Make sure the trajectory is fitted and
                         pbc have been removed. Allowed formats: *.xtc (recommended), *.trr, *.pdb
                         (specify as many as you'd like).''')
group.add_argument('-cp', dest='complex_top', metavar='<Topology>', default=None, type=topology,
                   help='''The complex Topology file. When it is defined -lm option is not needed''')
group.add_argument('-cr', dest='reference_structure', metavar='<PDB File>', default=None, type=pdb,
                   help='''Complex Reference Structure file. This option is optional but recommended 
                         (Use the PDB file used to generate the topology in GROMACS). If not defined,
                         the chains ID assignment (if the structure used in -cs does not have chain
                         IDs) will be done automatically according to the structure (can generate
                         wrong mapping).''')

group = parser.add_argument_group('Receptor', receptor_group_des)
group.add_argument('-rs', dest='receptor_tpr', metavar='<Structure File>', default=None, type=structure,
                   help='''Structure file of the unbound receptor for multiple trajectory approach.
                         Allowed formats: *.tpr (recommended), *.pdb''')
group.add_argument('-ri', dest='receptor_index', metavar='<Index File>', default=None, type=index,
                   help='Index file of the unbound receptor.')
group.add_argument('-rg', dest='receptor_group', metavar='index', default=None, type=int,
                   help='''Receptor group in receptor index file. Notation: "-lg <Receptor group>", 
                         e.g. -rg 1''')
group.add_argument('-rt', dest='receptor_trajs', nargs='*', metavar='TRJ', type=trajectory,
                   help='''Input trajectories of the unbound receptor for multiple trajectory approach.
                         Allowed formats: *.xtc (recommended), *.trr, *.pdb, *.gro (specify as many
                         as you'd like).''')
group.add_argument('-rp', dest='receptor_top', metavar='<Topology>', default=None, type=topology,
                   help='''Topology file of the receptor.''')

group = parser.add_argument_group('Ligand', ligand_group_des)
group.add_argument('-lm', dest='ligand_mol2', metavar='<Structure File>', default=None, type=mol2,
                   help='''A *.mol2 file of the unbound ligand used to parametrize
                         ligand for GROMACS using Antechamber. Must be defined
                         if Protein-Ligand (small molecule) complex was define and -cp or -lp option are not defined.
                         No needed for Proteins, DNA, RNA, Ions, Glycans or any ligand parametrized in the Amber 
                         force fields. Must be the Antechamber output *.mol2.''')
group.add_argument('-ls', dest='ligand_tpr', metavar='<Structure File>', default=None, type=structure,
                   help='''Structure file of the unbound ligand. If ligand is a small molecule and -lp is not defined,
                   make sure that you define above -lm option. Allowed formats: *.tpr (recommended), *.pdb''')
group.add_argument('-li', dest='ligand_index', metavar='<Index File>', type=index,
                   default=None, help='Index file of the unbound ligand. Only if tpr file was define in -ls.')
group.add_argument('-lg', dest='ligand_group', metavar='index', default=None, type=int,
                   help='''Ligand group in ligand index file. Notation: "-lg <Ligand group>", 
                         e.g. -lg 13''')
group.add_argument('-lt', dest='ligand_trajs', nargs='*', metavar='TRJ', type=trajectory,
                   help='''Input trajectories of the unbound ligand for multiple trajectory approach. 
                         Allowed formats: *.xtc (recommended), *.trr, *.pdb, *.gro (specify as many
                         as you'd like).''')
group.add_argument('-lp', dest='ligand_top', metavar='<Topology>', default=None, type=topology,
                   help='''Topology file of the ligand.''')

group = parser.add_argument_group('Miscellaneous Actions')
group.add_argument('--rewrite-output', dest='rewrite_output', default=False,
                   action='store_true', help='''Do not re-run any calculations,
                  just parse the output files from the previous calculation and
                  rewrite the output files.''')
group.add_argument('--clean', dest='clean', action='store_true', default=False,
                   help='''Clean temporary files and quit.''')

# GUI parser
description = 'This program is part of gmx_MMPBSA and will show a workspace to analyze the gmx_MMPBSA results'
anaparser = ArgumentParser(epilog=f'gmx_MMPBSA is an effort to implement the GB/PB and others calculations in '
                                  f'GROMACS. \nBased on MMPBSA.py (version {__mmpbsa_version__}) and '
                                  f'AmberTools{__ambertools_version__}',
                           description=description,
                           formatter_class=ArgumentDefaultsHelpFormatter)
anaparser.add_argument('-v', '--version', action='version',
                       version='%%(prog)s %s based on MMPBSA version %s' % (__version__, __mmpbsa_version__))
group = anaparser.add_argument_group('Info file')
group.add_argument('-f', '--files', nargs='*', help='gmx_MMPBSA info files or container folder or list of them',
                   type=Path, default=None)
group.add_argument('-r', '--recursive', help='Search recursively in this folder at depth = 1', action='store_true',
                   default=False)

# tester parser
description = ('This program is part of gmx_MMPBSA and will allow you to run various gmx_MMPBSA examples easily.')
testparser = ArgumentParser(epilog=f'gmx_MMPBSA is an effort to implement the GB/PB and others calculations in '
                                  f'GROMACS. \nBased on MMPBSA.py (version {__mmpbsa_version__}) and '
                                  f'AmberTools{__ambertools_version__}',
                           description=description,
                           formatter_class=RawTextHelpFormatter)
testparser.add_argument('-v', '--version', action='version',
                       version='%%(prog)s %s based on MMPBSA version %s' % (__version__, __mmpbsa_version__))
group = testparser.add_argument_group('Test options')
group.add_argument('-t', dest='test', choices=list(range(19)) + [101], type=int, nargs='*', default=[2],
                   help='''\
The level the test is going to be run at. Multiple systems and analysis can be run at the same time.
      Nr. of Sys  
* 0      16     All -- Run all examples (Can take a long time!!!)
* 1      13     Minimal -- Does a minimal test with a set of systems and analyzes 
                that show that gmx_MMPBSA runs correctly. Only exclude 3drism, nmode
                protein-ligand MT because take a long time or are redundant
* 2       9     Fast -- Only the calculations that take a short time are run (Default)
[Systems]:
     Slow Frames
* 3    . | 10   Protein-Ligand (Single trajectory approximation)
* 4    . | 10   Protein-Protein
* 5    . | 10   Protein-DNA
* 6    x |  4   Protein-Membrane
* 7    . | 10   Protein-Glycan
* 8    x |  4   Metalloprotein-Peptide
* 9    . | 10   Protein-DNA-RNA-IONs-Ligand
* 10   x |  4   Protein-Ligand (CHARMM force field)
* 11   x |  4   Protein-ligand complex in membrane with CHARMMff 
[Analysis]:
     Slow Frames
* 12   . | 10   Alanine Scanning
* 13   . | 10   Stability calculation
* 14   . | 10   Decomposition Analysis
* 15   . | 16   Interaction Entropy approximation
* 16   . | 10   Protein-Ligand (Multiple trajectory approximation)
* 17   x |  4   Entropy calculation using Normal Mode approximation 
* 18   x |  4   Calculations using 3D-RISM approximation
''')
group.add_argument('-f', '--folder', help='Defines the folder to store all data', type=Path, default='.')
group.add_argument('-r', '--reuse', help='Defines the existing test forlder will be reuse', action='store_true')
group.add_argument('-ng', '--nogui', help='No open gmx_MMPBSA_ana after all calculations finished',
                   action='store_true',  default=False)
group.add_argument('-n', '--num_processors', type=int, default=4,
                   help='Defines the number of processor cores you want to use with MPI per calculation. If the number '
                        'of frames is less than the number of cpus defined, the calculation will be performed with '
                        'the number of processors = number of frames')
