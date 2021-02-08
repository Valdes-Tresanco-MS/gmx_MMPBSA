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
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from pathlib import Path
from GMXMMPBSA import __version__, __mmpbsa_version__, __ambertools_version__
from GMXMMPBSA.exceptions import GMXMMPBSA_ERROR


class OptionList(object):
    """
    Just a container to hold the command-line options. Necessary when reading in
    a gmx_MMPBSA info file to have a container to load the results from the
    parser.
    """
    pass


if os.getenv('AMBERHOME'):
    rism = os.path.join(os.getenv('AMBERHOME'), 'dat', 'mmpbsa', 'spc.xvv')
else:
    rism = None

def check_arg(str_suffix, path=False):
    def decorator(f):
        def gmx_MMPBSA_file(*args):
            result = f(*args)
            if path and not Path(result).exists():
                GMXMMPBSA_ERROR(f'{result} do not exist or is inaccessible')
            if not Path(result).suffix in str_suffix:
                GMXMMPBSA_ERROR(f'{result} does not correspond to the required structure format {str_suffix}')
            return result

        return gmx_MMPBSA_file

    return decorator


@check_arg(['.tpr', '.pdb', '.gro'], True)
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


description = '''gmx_MMPBSA is an effort to implement the GB/PB and others calculations in GROMACS. This program is 
                 an adaptation of Amber's MMPBSA.py and essentially works as such. As gmx_MMPBSA adapts MMPBSA.py, 
                 since it has all the resources of this script and work with any GROMACS version.'''

complex_group_des = '''Complex files and info that are needed to perform the calculation. If the receptor and / or the
                        ligand info is not defined, we generate them from that of the complex.'''

receptor_group_des = '''Receptor files and info that are needed to perform the calculation. If the receptor info is not 
                         defined, we generate it from that of the complex.'''

ligand_group_des = '''Ligand files and info that are needed to perform the calculation. If the ligand are not defined, 
                       we generate it from that of the complex.'''
# Set up the MM/PBSA parser here. It clutters up the MMPBSA_App to do it there

# noinspection PyTypeChecker
parser = ArgumentParser(epilog=f'''This program will calculate binding free energies using end-state free energy 
                                    methods on an ensemble of snapshots using a variety of implicit solvent 
                                    models.\nBased on MMPBSA.py (version {__mmpbsa_version__}) and 
                                    AmberTools{__ambertools_version__}''',
                        description=description, formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument('-v', '--version', action='version',
                    version='''%%(prog)s %s based on MMPBSA version %s and AmberTools %s''' %
                            (__version__, __mmpbsa_version__, __ambertools_version__))
parser.add_argument('--input-file-help', dest='infilehelp', action='store_true',
                    help='Print all available options in the input file.',
                    default=False)
group = parser.add_argument_group('Miscellaneous Options')
group.add_argument('-O', '--overwrite', default=False, action='store_true',
                   help='Allow output files to be overwritten', dest='overwrite')
group.add_argument('-prefix', dest='prefix', default='_GMXMMPBSA_',
                   metavar='<file prefix>',
                   help='Prefix for intermediate files.')
group = parser.add_argument_group('Input and Output Files', '''These options specify the input files and optional 
output files.''')
group.add_argument('-i', dest='input_file', metavar='FILE', help='MM/PBSA input file.')
group.add_argument('-xvvfile', dest='xvvfile', help='XVV file for 3D-RISM.', default=rism)
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
                   help='Open charts application (gmx_MMPBSA_ana) when all calculations finished')
group.add_argument('-s', '--stability', dest='stability', action='store_true', default=False,
                   help='''Perform stability calculation. Only the complex parameters are required. If
                         ligand is non-Protein (small molecule) type, then ligand *.mol2 file is 
                         required. In any other case receptor and ligand parameters will be ignored.
                         See description bellow''')

group = parser.add_argument_group('Complex', complex_group_des)
group.add_argument('-cs', dest='complex_tpr', metavar='<Structure File>', default=None, type=structure,
                   help='''Structure file of the complex. If it is Protein-Ligand (small molecule)
                         complex, make sure that you define -lm option. See -lm description below
                         Allowed formats: *.tpr (recommended), *.pdb, *.gro''')
group.add_argument('-ci', dest='complex_index', metavar='<Index File>', default=None, type=index,
                   help='Index file of the bound complex.')
group.add_argument('-cg', dest='complex_groups', metavar='index', nargs=2, default=None, type=int,
                   help='Groups of receptor and ligand in complex index file. The notation is as follows: "-cg '
                        '<Receptor group> <Ligand group>", ie. -cg 1 13')
group.add_argument('-ct', dest='complex_trajs', nargs='*', metavar='TRJ', type=trajectory,
                   help='''Input trajectories of the complex. Make sure the trajectory is fitted and
                         pbc have been removed. Allowed formats: *.xtc (recommended), *.trr, *.pdb
                         (specify as many as you'd like).''')
group.add_argument('-cp', dest='complex_top', metavar='<Topology>', default=None, type=topology,
                   help='''Topology file of the complex.''')
group.add_argument('-cr', dest='reference_structure', metavar='<PDB File>', default=None, type=pdb,
                   help='''Complex Reference Structure file. This option is optional but recommended 
                         (Use the PDB file used to generate the topology in GROMACS). If not defined,
                         the chains ID assignment (if the structure used in -cs does not have chain
                         IDs) will be done automatically according to the structure (can generate
                         inconsistencies).''')

group = parser.add_argument_group('Receptor', receptor_group_des)
group.add_argument('-rs', dest='receptor_tpr', metavar='<Structure File>', default=None, type=structure,
                   help='''Structure file of the unbound receptor for multiple trajectory approach.
                         Allowed formats: *.tpr (recommended), *.pdb, *.gro''')
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
                   help='''A *.mol2 file of the unbound ligand used to parametrize ligand for GROMACS
                         using Anetchamber. Must be defined if Protein-Ligand (small molecule) 
                         complex was define. No needed for Proteins, DNA, RNA, Ions and Glycans.
                         Antechamber output *.mol2 is recommended.''')
group.add_argument('-ls', dest='ligand_tpr', metavar='<Structure File>', default=None, type=structure,
                   help='''Structure file of the unbound ligand. If ligand is a small molecule, make 
                         sure that you define above -lm option. Allowed formats: *.tpr (recommended),
                         *.pdb, *.gro''')
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
group.add_argument('-make-mdins', dest='make_mdins', default=False,
                   action='store_true', help='''Create the input files for each
                  calculation and quit. This allows you to modify them and
                  re-run using -use-mdins''')
group.add_argument('-use-mdins', dest='use_mdins', default=False,
                   action='store_true', help='''Use existing input files for each
                  calculation. If they do not exist with the appropriate names,
                  %(prog)s will quit in error.''')
group.add_argument('-rewrite-output', dest='rewrite_output', default=False,
                   action='store_true', help='''Do not re-run any calculations,
                  just parse the output files from the previous calculation and
                  rewrite the output files.''')
group.add_argument('--clean', dest='clean', action='store_true', default=False,
                   help='''Clean temporary files and quit.''')

#### GUI parser
anaparser = ArgumentParser(epilog='''This program is part of gmx_MMPBSA and will show a workspace with 
                            charts to analyze the results''',
                           description=description,
                           formatter_class=ArgumentDefaultsHelpFormatter)
anaparser.add_argument('-v', '--version', action='version',
                       version='%%(prog)s %s based on MMPBSA version %s' % (__version__, __mmpbsa_version__))
group = anaparser.add_argument_group('Info file')
group.add_argument('-p', '--path', dest='path', help='Path to gmx_MMPBSA info file', required=True,
                   default=None)
