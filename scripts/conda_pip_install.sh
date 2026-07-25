#!/bin/bash
#
# reproducible script to install gmx_MMPBSA on Ubuntu 24.04.3 LTS
# might also work on some other Linux distributions

ENV=gmx_MMPBSA

# conda env and python
conda create -n $ENV "python=3.12" -y -q
conda activate $ENV

set -u
#set -x

# packages which are not available through pip
# (too bad, or we could get rid of conda-the-dumb-and-slow)
conda install -c conda-forge "ambertools>=24.8,<27" -y -q
conda install -c conda-forge "gromacs>=2022,<2027" pocl -y -q
conda install -c conda-forge pyqt6 -y -q

# pip for all other dependencies
pip install "matplotlib>=3.8,<4" \
            "mpi4py>=4.0.1,<5"  \
            "numpy>=1.26.4,<2"  \
            "pandas>=2.2,<3"    \
            "scipy>=1.14.1,<2"  \
            "seaborn>=0.13,<0.14" \
            gmx_MMPBSA

# very important: test that everything works
time gmx_MMPBSA_test -f tests -n 10
