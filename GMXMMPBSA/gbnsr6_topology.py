"""Helpers for GBNSR6-specific topology compatibility."""

from pathlib import Path
import os
import shutil


DIHEDRAL_POINTER_INDEXES = (6, 7, 14, 17)
_LJ_FLAGS = ('LENNARD_JONES_ACOEF', 'LENNARD_JONES_BCOEF')


def prepare_gbnsr6_topology(prmtop, output_prmtop):
    """Create a GBNSR6-compatible prmtop copy.

    The legacy GBNSR6 reader has fixed-size atom-type and bonded-parameter
    tables. GROMACS topologies can contain many differently named atom types
    whose complete LJ interaction rows are nevertheless identical, and large
    protein topologies can exceed the reader's dihedral table. Merging exact
    LJ duplicates and clearing dihedral pointers keeps the GB-only input within
    those limits. Bonded MM terms are still calculated separately with the
    original topology.
    """
    prmtop = Path(prmtop)
    output_prmtop = Path(output_prmtop)

    tmp_prmtop = output_prmtop.with_name(f'{output_prmtop.name}.{os.getpid()}.tmp')
    shutil.copyfile(prmtop, tmp_prmtop)
    _compact_equivalent_lj_types(tmp_prmtop)
    _strip_dihedral_terms_in_place(tmp_prmtop)
    os.replace(tmp_prmtop, output_prmtop)
    return output_prmtop.as_posix()


def strip_dihedral_terms_for_gbnsr6(prmtop, output_prmtop):
    """Backward-compatible entry point for GBNSR6 topology preparation."""
    return prepare_gbnsr6_topology(prmtop, output_prmtop)


def _strip_dihedral_terms_in_place(prmtop):
    prmtop = Path(prmtop)
    lines = prmtop.read_text().splitlines(keepends=True)
    prmtop.write_text(''.join(_strip_dihedral_terms(lines)))


def _compact_equivalent_lj_types(prmtop):
    """Merge atom types with exactly identical LJ interaction matrix rows."""
    from parmed.amber import AmberFormat

    topology = AmberFormat(str(prmtop))
    data = topology.parm_data
    required = {'POINTERS', 'ATOM_TYPE_INDEX', 'NONBONDED_PARM_INDEX', *_LJ_FLAGS}
    if not required.issubset(data):
        return

    ntypes = data['POINTERS'][1]
    atom_type_indexes = data['ATOM_TYPE_INDEX']
    used_types = sorted(set(atom_type_indexes))
    nonbonded_indexes = data['NONBONDED_PARM_INDEX']
    acoef, bcoef = (data[flag] for flag in _LJ_FLAGS)

    def coefficient(atom_type_i, atom_type_j, values):
        matrix_index = ntypes * (atom_type_i - 1) + atom_type_j - 1
        return values[nonbonded_indexes[matrix_index] - 1]

    signatures = {
        atom_type_i: tuple(
            (coefficient(atom_type_i, atom_type_j, acoef),
             coefficient(atom_type_i, atom_type_j, bcoef))
            for atom_type_j in used_types
        )
        for atom_type_i in used_types
    }

    representatives = []
    compact_index = {}
    for atom_type_i in used_types:
        for index, representative in enumerate(representatives, start=1):
            if signatures[atom_type_i] == signatures[representative]:
                compact_index[atom_type_i] = index
                break
        else:
            representatives.append(atom_type_i)
            compact_index[atom_type_i] = len(representatives)

    compact_count = len(representatives)
    if compact_count == ntypes:
        return

    pair_indexes = {}
    compact_acoef = []
    compact_bcoef = []
    compact_nonbonded_indexes = []
    for type_i in range(1, compact_count + 1):
        for type_j in range(1, compact_count + 1):
            pair = tuple(sorted((type_i, type_j)))
            if pair not in pair_indexes:
                pair_indexes[pair] = len(compact_acoef) + 1
                original_i = representatives[type_i - 1]
                original_j = representatives[type_j - 1]
                compact_acoef.append(coefficient(original_i, original_j, acoef))
                compact_bcoef.append(coefficient(original_i, original_j, bcoef))
            compact_nonbonded_indexes.append(pair_indexes[pair])

    data['POINTERS'][1] = compact_count
    data['ATOM_TYPE_INDEX'] = [compact_index[index] for index in atom_type_indexes]
    data['NONBONDED_PARM_INDEX'] = compact_nonbonded_indexes
    data['LENNARD_JONES_ACOEF'] = compact_acoef
    data['LENNARD_JONES_BCOEF'] = compact_bcoef
    if 'SOLTY' in data:
        data['SOLTY'] = data['SOLTY'][:compact_count]

    compacted = Path(prmtop).with_name(f'{Path(prmtop).name}.compact')
    topology.write_parm(str(compacted))
    os.replace(compacted, prmtop)


def _strip_dihedral_terms(lines):
    output = []
    index = 0

    while index < len(lines):
        line = lines[index]
        if not line.startswith('%FLAG '):
            output.append(line)
            index += 1
            continue

        flag = line.split(None, 1)[1].strip()
        block, index = _read_flag_block(lines, index)

        if flag == 'POINTERS':
            output.extend(_strip_dihedral_pointers(block))
        else:
            output.extend(block)

    return output


def _read_flag_block(lines, start):
    index = start + 1
    while index < len(lines) and not lines[index].startswith('%FLAG '):
        index += 1
    return lines[start:index], index


def _strip_dihedral_pointers(block):
    header = _flag_header(block)
    values = []
    for line in block[len(header):]:
        values.extend(int(item) for item in line.split())

    for pointer_index in DIHEDRAL_POINTER_INDEXES:
        if pointer_index < len(values):
            values[pointer_index] = 0

    return header + _format_int_block(values)


def _flag_header(block):
    header = block[:1]
    if len(block) > 1 and block[1].startswith('%FORMAT'):
        header.append(block[1])
    return header


def _format_int_block(values):
    lines = []
    for index in range(0, len(values), 10):
        chunk = values[index:index + 10]
        lines.append(''.join(f'{value:8d}' for value in chunk) + '\n')
    return lines
