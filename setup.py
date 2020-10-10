# ##############################################################################
#                           GPLv3 LICENSE INFO                                 #
#                                                                              #
#  Copyright (C) 2020  Mario S. Valdes-Tresanco and Mario E. Valdes-Tresanco   #
#  Copyright (C) 2014  Jason Swails, Bill Miller III, and Dwight McGee         #
#                                                                              #
#   Project: https://github.com/Valdes-Tresanco-MS/GMX-MMPBSA                  #
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

from setuptools import setup

with open("README.md", "r") as f:
    LONG_DESCRIPTION = f.read()

setup(
    name='GMX-MMPBSA',
    version='1.0.0',
    packages=['GMXMMPBSA'],
    license='GPLv3',
    author='Mario S. Valdes-Trasanco and Mario E. Valdes-Tresanco ',
    author_email='mariosergiovaldes145@gmail.com',
    maintainer='Mario S. Valdes-Trasanco',
    maintainer_email='mariosergiovaldes145@gmail.com',
    url='https://github.com/Valdes-Tresanco-MS/GMX-MMGBSA',
    description='Adaptation of MMPBSA.py (AMBER) to use Gromacs files',
    long_description=LONG_DESCRIPTION,
    long_description_content_type='text/markdown',
    keywords=['GMX-MMPBSA', 'MMPBSA', 'GROMACS', 'AmberTools'],
    entry_points={
        "console_scripts": [
            "gmx_mmpbsa=GMXMMPBSA.app:gmxmmpbsa",
            "gmx_mmpbsa_gui=GMXMMPBSA.app:gmxmmpbsa_gui",]}
)
