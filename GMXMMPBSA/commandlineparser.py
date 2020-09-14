# -*- coding: utf-8 -*-
"""
This module contains classes and such that are responsible for parsing
command-line arguments for MMPBSA.py.  All of the files specified for use
in MMPBSA.py will be assigned as attributes to the returned class.

                          GPL LICENSE INFO                             

  Copyright (C) 2009  Dwight McGee, Billy Miller III, and Jason Swails

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
   
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place - Suite 330,
  Boston, MA 02111-1307, USA.
"""

import os, sys
from GMXMMPBSA import __version__, __mmpbsa_version__


class OptionList(object):
    """
    Just a container to hold the command-line options. Necessary when reading in
    a MMPBSA.py info file to have a container to load the results from the
    parser.
    """
    pass


if os.getenv('AMBERHOME'):
    rism = os.path.join(os.getenv('AMBERHOME'), 'dat', 'mmpbsa', 'spc.xvv')
else:
    rism = None

description = """
                            🅶🅼🆇-🅼🅼🅿🅱🆂🅰 
                                🅥.{}
asdlml asdnllmasd asdjalsjdpasd asdljsadkpasd asdlkaspdkpoaskd askdaskdp1123123a
 """.format(__version__)
complex_group_des = """
Complex structure file (tpr), Gromacs index file (ndx), Receptor and Ligand 
groups names in the index file and forcefield to make Amber topology. 
If is necessary we split complex to generate Receptor and Ligand topology like 
ante-MMPBSA"""
receptor_group_des = """
Complex structure file (tpr), Gromacs index file (ndx), Receptor and Ligand 
groups names in the index file and forcefield to make Amber topology. 
If is necessary we split complex to generate Receptor and Ligand topology like 
ante-MMPBSA"""
ligand_group_des = """
Complex structure file (tpr), Gromacs index file (ndx), Receptor and Ligand 
groups names in the index file and forcefield to make Amber topology. 
If is necessary we split complex to generate Receptor and Ligand topology like 
ante-MMPBSA"""
# Set up the MM/PBSA parser here. It clutters up the MMPBSA_App to do it there
from argparse import ArgumentParser, RawDescriptionHelpFormatter

# noinspection PyTypeChecker
parser = ArgumentParser(epilog='''This program will calculate binding free
                        energies using end-state free energy methods on an
                        ensemble of snapshots using a variety of implicit
                        solvent models''',
                        description=description,
                        formatter_class=RawDescriptionHelpFormatter)
parser.add_argument('-v', '--version', action='version',
                    version='%%(prog)s %s based on MMPBSA version %s' % (__version__, __mmpbsa_version__))
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
group = parser.add_argument_group('Options',
                                  '''These options specify explicit calculation type and forcefield to prepare the 
Amber topologies''')
group.add_argument('-s', dest='stability', action='store_true', default=False,
                   help='Perform stability calculation. Only the complex parameter are required. Only, if ligand is '
                        'non-Protein (small molecule) type will required the ligand parameters. IN any other case '
                        'receptor and ligand parameters will be ignored')
group.add_argument('-pff', dest='protein_ff', default='amber14sb', metavar='<Forcefield>',
                   help='Used forcefield to make the protein MD in Gromacs. Allowed: amber14sb, amber99sb-ildn, '
                        'amber99sb, amber03, amber99, amber96, amber94. (Default: amber14sb)')
group.add_argument('-lff', dest='ligand_ff', metavar='<Forcefield>', default='gaff',
                   help='Used forcefield to make the ligand MD in Gromacs. Allowed: gaff, gaff2. (Default: gaff)')

group = parser.add_argument_group('Complex', complex_group_des)
group.add_argument('-cs', dest='complex_tpr', metavar='<Structure File>', default=None, required=True,
                   help='''Structure file of the bound complex''')
group.add_argument('-ci', dest='complex_index', metavar='<Index File>', default=None, required=True,
                   help='Index file of the bound complex.')
group.add_argument('-cg', dest='complex_groups', metavar='index', nargs=2, default=None, required=True,
                   help='Groups of the receptor and ligand in complex index file. The notation is as follows: "-cg '
                        '<Receptor group> <Ligand group>"')
group.add_argument('-ct', dest='complex_trajs', nargs='*', metavar='TRJ', required=True,
                   help='''Input trajectories of the complex. Allowed formats: *.xtc (recommended), *.trr, *.pdb
                  (specify as many as you'd like).''')

group = parser.add_argument_group('Receptor', receptor_group_des)
group.add_argument('-rs', dest='receptor_tpr', metavar='<Structure File>', default=None,
                   help='''Structure file of the unbound receptor. If omitted and stability is False, the structure 
                   from complex will be used.''')
group.add_argument('-ri', dest='receptor_index', metavar='<Index File>', default=None,
                   help='Index file of the unbound receptor.')
group.add_argument('-rg', dest='receptor_group', metavar='index', default=None,
                   help='Group of the receptor in receptor index file.')
group.add_argument('-rt', dest='receptor_trajs', nargs='*', metavar='TRJ',
                   help='''Input trajectories of the unbound receptor. Allowed formats: *.xtc (recommended), *.trr, 
                   *.pdb. (specify as many as you'd like).''')

group = parser.add_argument_group('Ligand', ligand_group_des)
group.add_argument('-ls', dest='ligand_tprOmol2', metavar='<Structure File>', default=None,
                   help='Structure file of the unbound ligand. Can be tpr file if ligand is protein-like or mol2 '
                        'file if is a small molecule. If omitted, is protein-like and stability is False the '
                        'structure from complex will be used. If ligand is a small molecule, this option always most '
                        'be defined and the rest (-li, -lg and -lt) are not needed')
group.add_argument('-li', dest='ligand_index', metavar='<Index File>',
                   default=None, help='Index file of the unbound ligand. Only if tpr file was define in -ls.')
group.add_argument('-lg', dest='ligand_group', metavar='index', default=None,
                   help='Group of the ligand in ligand index file. Only if tpr file was define in -ls.')
group.add_argument('-lt', dest='ligand_trajs', nargs='*', metavar='TRJ',
                   help='''Input trajectories of the unbound ligand. Allowed formats: *.xtc (recommended), *.trr, *.pdb
                  (specify as many as you'd like).''')


group = parser.add_argument_group('Mutant', ligand_group_des)
group.add_argument('-mc', dest='mutant_complex_tpr',
                   metavar='<Topology File>', help='''Complex topology file of
                  the mutant complex in which one residue has been mutated to
                  either a glycine or alanine to perform computational alanine
                  (or glycine) scanning.''')
group.add_argument('-mr', dest='mutant_receptor_tpr',
                   metavar='<Topology File>', help='''Receptor topology file of
                  the mutant receptor (see -mc above). If omitted, the mutation
                  is assumed to be in the ligand.''')
group.add_argument('-ml', dest='mutant_ligand_tpr',
                   metavar='<Topology File>', help='''Ligand topology file of the
                  mutant receptor (see -mc above). If omitted, the mutation is
                  assumed to be in the receptor.''')

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
# k = parser.parse_args(['sdfsd'])
# print k
