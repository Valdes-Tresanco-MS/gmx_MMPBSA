#!/usr/bin/env python3
"""One-time helper to draft gmx_MMPBSA_test_manifest.json from example READMEs."""

from __future__ import annotations

import argparse
import json
import shlex
from pathlib import Path

TESTS = {
    3: ('Protein_ligand/ST', 'Protein_ligand/ST', 'Protein-Ligand (Single trajectory approximation)'),
    4: ('Protein_protein', 'Protein_protein', 'Protein-Protein'),
    5: ('Protein_DNA', 'Protein_DNA', 'Protein-DNA'),
    6: ('Protein_membrane', 'Protein_membrane', 'Protein-Membrane'),
    7: ('Protein_glycan', 'Protein_glycan', 'Protein-Glycan'),
    8: ('Metalloprotein_ligand', 'Metalloprotein_ligand', 'Metalloprotein-ligand'),
    9: ('Comp_receptor', 'Comp_receptor', 'Comp_receptor'),
    10: ('Protein_ligand_CHARMMff', 'Protein_ligand_CHARMMff', 'Protein-Ligand (CHARMM force field)'),
    11: ('Protein_membrane_CHARMMff', 'Protein_membrane_CHARMMff', 'Protein-ligand complex in membrane with CHARMMff'),
    12: ('Alanine_scanning', 'Alanine_scanning', 'Alanine Scanning'),
    13: ('Stability', 'Stability', 'Stability calculation'),
    14: ('Decomposition_analysis', 'Decomposition_analysis', 'Decomposition Analysis'),
    15: ('Entropy_calculations/Interaction_Entropy', 'Entropy_calculations/Interaction_Entropy', 'Interaction Entropy approximation'),
    16: ('Protein_ligand/MT', 'Protein_ligand/MT', 'Protein-Ligand (Multiple trajectory approximation)'),
    17: ('Entropy_calculations/nmode', 'Entropy_calculations/nmode', 'Entropy calculation using Normal Mode approximation'),
    18: ('3D-RISM', '3D-RISM', 'Calculations using 3D-RISM approximation'),
    19: ('Entropy_calculations/C2_Entropy', 'Entropy_calculations/C2_Entropy', 'C2 Entropy approximation'),
    20: ('Linear_PB_solver', 'Linear_PB_solver', 'LPB Calculation'),
    21: ('NonLinear_PB_solver', 'NonLinear_PB_solver', 'NLPB Calculation'),
    22: ('Protein_ligand_LPH_atoms_CHARMMff', 'Protein_ligand_LPH_atoms_CHARMMff', 'Protein-Ligand_LPH (CHARMM force field)'),
    23: ('QM_MMGBSA', 'QM_MMGBSA', 'QM/MMGBSA Calculation'),
    24: ('GBNSR6', 'GBNSR6', 'GBNSR6 Calculation'),
    25: ('AMBER', 'AMBER', 'AMBER input files'),
    26: ('Explicit_receptor_waters', 'Explicit_receptor_waters', 'ST MM/PB(GB)SA with explicit receptor waters'),
}

SLOW = {6, 8, 10, 11, 17, 18}
SUITE_MINIMAL = [3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15]
SUITE_FAST = [3, 4, 5, 7, 9, 12, 13, 14, 15]


def _suites_for(test_id: int) -> list[str]:
    suites = ['all']
    if test_id in SUITE_MINIMAL:
        suites.append('minimal')
    if test_id in SUITE_FAST:
        suites.append('fast')
    return suites


def _parse_serial_command(readme_path: Path) -> str | None:
    text = readme_path.read_text()
    in_serial = False
    for line in text.splitlines():
        if '=== "Serial"' in line:
            in_serial = True
            continue
        if in_serial and line.startswith('==='):
            break
        if in_serial and line.strip():
            command = line.strip()
            if command.startswith('gmx_MMPBSA ') or command.startswith('amber_MMPBSA '):
                return command
            if command.startswith('ggmx_MMPBSA '):
                return command.replace('ggmx_MMPBSA', 'gmx_MMPBSA', 1)
    for line in text.splitlines():
        command = line.strip()
        if (command.startswith('gmx_MMPBSA -O') or command.startswith('amber_MMPBSA -O')) and 'mpirun' not in command:
            return command
    return None


def _outputs_from_args(args: list[str]) -> list[str]:
    outputs = []
    index = 0
    while index < len(args):
        if args[index] in ('-o', '-eo', '-do', '-deo') and index + 1 < len(args):
            outputs.append(args[index + 1])
            index += 2
        else:
            index += 1
    return outputs


def build_manifest(examples_dir: Path) -> dict:
    tests = {}
    for test_id, (path, workdir, name) in TESTS.items():
        readme = examples_dir / path / 'README.md'
        command = _parse_serial_command(readme)
        if not command:
            raise SystemExit(f'No Serial command found for test {test_id} in {readme}')
        parts = shlex.split(command)
        executable = parts[0]
        command_args = parts[1:]
        entry = {
            'name': name,
            'path': path,
            'workdir': workdir,
            'input': 'mmpbsa.in',
            'executable': executable,
            'command_args': command_args,
            'suites': _suites_for(test_id),
            'slow': test_id in SLOW,
            'requires': [executable, 'mpirun'],
            'expected_outputs': _outputs_from_args(command_args),
        }
        if test_id == 18:
            entry['known_issues'] = ['rism_fortran']
        tests[str(test_id)] = entry

    return {
        'version': 1,
        'suites': {
            'all': {'id': 0, 'tests': list(range(3, 27))},
            'minimal': {'id': 1, 'tests': SUITE_MINIMAL},
            'fast': {'id': 2, 'tests': SUITE_FAST},
        },
        'aliases': {
            '101': 'all',
            'protein_ligand_st': 3,
            'protein_ligand_mt': 16,
            'explicit_receptor_waters': 26,
            'gbnsr6': 24,
            'amber': 25,
            'decomposition': 14,
            'rism': 18,
        },
        'tests': tests,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--examples-dir',
        type=Path,
        default=Path(__file__).resolve().parents[1] / 'examples',
        help='Path to the examples directory',
    )
    parser.add_argument(
        '--output',
        type=Path,
        default=Path(__file__).resolve().parents[1] / 'GMXMMPBSA' / 'data' / 'gmx_MMPBSA_test_manifest.json',
    )
    args = parser.parse_args()
    manifest = build_manifest(args.examples_dir)
    args.output.write_text(json.dumps(manifest, indent=2) + '\n')
    print(f'Wrote {args.output}')


if __name__ == '__main__':
    main()
