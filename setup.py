# ##############################################################################
#                           GPLv3 LICENSE INFO                                 #
#                                                                              #
#  Copyright (C) 2020  Mario S. Valdés-Tresanco and Mario E. Valdés-Tresanco   #
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

from setuptools import setup, find_packages
import versioneer
import sys

with open("README.md", "r") as f:
    LONG_DESCRIPTION = f.read()

if sys.version_info < (3, 11) or sys.version_info >= (3, 13):
    raise RuntimeError("gmx_MMPBSA requires python >=3.11,<3.13")

setup(
    name='gmx_MMPBSA',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    packages=find_packages(),
    package_data={"GMXMMPBSA": ["data/*", 'data/gmxMMPBSA/*', 'data/xvv_files/*', 'data/gmx_MMPBSA_test_manifest.json', 'analyzer/style/*', 'GMXMMPBSA.sh']},
    license='GPLv3',
    author='Mario S. Valdes-Tresanco and Mario E. Valdes-Tresanco ',
    author_email='mariosergiovaldes145@gmail.com',
    maintainer='Mario S. Valdes-Tresanco',
    maintainer_email='mariosergiovaldes145@gmail.com',
    url='https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA',
    description="gmx_MMPBSA is a new tool based on AMBER's MMPBSA.py aiming to perform end-state free energy  "
                "calculations with GROMACS files",
    long_description=LONG_DESCRIPTION,
    long_description_content_type='text/markdown',
    keywords=['gmx_MMPBSA', 'MMPBSA', 'MMGBSA', 'GROMACS', 'AmberTools'],
    # Dependency ranges allow the Python 3.12 stack required by newer GROMACS
    # packages while preserving Python 3.11 support.
    python_requires='>=3.11,<3.13',
    install_requires=['numpy>=1.26.4,<2',
                      'pandas>=2.2.0,<3',
                      'matplotlib>=3.8.0,<4',
                      'seaborn>=0.13.0,<0.14',
                      'scipy>=1.14.1,<2',
                      'mpi4py>=4.0.1,<5',
                      'parmed>=4.2.2,<5',
                      'tqdm',
                      'rich>=13,<15'],
    entry_points={
        "console_scripts": [
            "gmx_MMPBSA=GMXMMPBSA.app:gmxmmpbsa",
            "amber_MMPBSA=GMXMMPBSA.app:gmxmmpbsa_amber",
            "gmx_MMPBSA_ana=GMXMMPBSA.app:gmxmmpbsa_ana",
            "gmx_MMPBSA_test=GMXMMPBSA.app:gmxmmpbsa_test"]}
)
