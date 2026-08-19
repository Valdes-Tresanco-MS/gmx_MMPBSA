"""
 This package contains all of the functions for gmx_MMPBSA that it
 needs to run smoothly.
"""

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

__all__ = ['alamdcrd', 'amber_outputs', 'analyzer', 'API', 'app', 'calculation', 'commandlineparser', 'createinput',
           'exceptions', 'infofile', 'input_parser', 'main', 'make_top', 'make_trajs',
           'output_file', 'parm_setup', 'timer', 'utils', '__version__', '__mmpbsa_version__',
           '__ambertools_version__', '__gromacs_version__']

__author__ = "Mario S. Valdes Tresanco, Mario E. Valdes Tresanco, Pedro A. Valiente PhD and Ernesto Moreno PhD"
__license__ = "GPLv3"
__mmpbsa_author__ = "Jason Swails, Dwight McGee, and Bill Miller III"
__mmpbsa_version__ = "14.0"


def _installed_conda_version(package_name):
    """Return a package version from the active Conda environment, if present."""
    import json
    import os
    import sys
    from pathlib import Path

    prefixes = []
    for value in (os.environ.get('CONDA_PREFIX'), os.environ.get('AMBERHOME'), sys.prefix):
        if value and value not in prefixes:
            prefixes.append(value)

    for prefix in prefixes:
        metadata_dir = Path(prefix) / 'conda-meta'
        for metadata_file in sorted(metadata_dir.glob(f'{package_name}-*.json'), reverse=True):
            try:
                version = json.loads(metadata_file.read_text(encoding='utf-8')).get('version')
            except (OSError, ValueError):
                continue
            if version:
                return str(version)
    return None


def _installed_distribution_version(distribution_name):
    """Return a Python distribution version as a fallback for non-Conda installs."""
    try:
        from importlib.metadata import version
        return version(distribution_name)
    except Exception:
        return None


def _gromacs_command_version():
    """Return the version reported by the active GROMACS executable."""
    import re
    import shutil
    import subprocess

    executable = shutil.which('gmx')
    if not executable:
        return None
    try:
        result = subprocess.run([executable, '--version'], capture_output=True, text=True,
                                check=False, timeout=5)
    except (OSError, subprocess.TimeoutExpired):
        return None
    match = re.search(r'^GROMACS version:\s*(\S+)', result.stdout, re.MULTILINE)
    return match.group(1) if match else None


def _tool_version(package_name, distribution_name=None, command_version=None):
    return (
        _installed_conda_version(package_name)
        or (_installed_distribution_version(distribution_name) if distribution_name else None)
        or (command_version() if command_version else None)
        or 'unknown'
    )


__ambertools_version__ = _tool_version('ambertools', 'ambertools')
__gromacs_version__ = _tool_version('gromacs', 'gromacs', _gromacs_command_version)

from . import _version
__version__ = _version.get_versions()['version']
