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
    name='gmx_MMPBSA',
    version='1.0.0',
    packages=['GMXMMPBSA'],
    license='GPLv3',
    author='Mario S. Valdes-Trasanco and Mario E. Valdes-Tresanco ',
    author_email='mariosergiovaldes145@gmail.com',
    maintainer='Mario S. Valdes-Trasanco',
    maintainer_email='mariosergiovaldes145@gmail.com',
    url='https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA',
    description='gmx_MMPBSA is a new tool aid to perform end-state free energy'
                'calculations with GROMACS files.',
    long_description=LONG_DESCRIPTION,
    long_description_content_type='text/markdown',
    project_urls={
        "Documentation": "https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA",
        "Issues": "https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/issues",
        "Publications": "Coming soon",
        "Source": "https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA",
    },
    classifiers=[
        "Programming Language :: Python :: 3.8",
        "License :: GPLv3",
        "Operating System :: OS Linux",
    ],
    keywords=['gmx_MMPBSA', 'MMPBSA', 'MMGBSA', 'GROMACS', 'AmberTools'],
    entry_points={
        "console_scripts": [
            "gmx_MMPBSA=GMXMMPBSA.app:gmxmmpbsa",
            "gmx_MMPBSA_gui=GMXMMPBSA.app:gmxmmpbsa_gui"]}
)
