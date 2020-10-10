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
    scripts=['GMX_MMPBSA.py']
)
