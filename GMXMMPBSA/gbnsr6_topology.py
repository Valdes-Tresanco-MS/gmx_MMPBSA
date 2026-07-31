"""Helpers for GBNSR6-specific topology compatibility."""

from pathlib import Path
import os


DIHEDRAL_POINTER_INDEXES = (6, 7, 14, 17)


def strip_dihedral_terms_for_gbnsr6(prmtop, output_prmtop):
    """Create a GBNSR6-only prmtop copy that reports no dihedral terms."""
    prmtop = Path(prmtop)
    output_prmtop = Path(output_prmtop)

    lines = prmtop.read_text().splitlines(keepends=True)
    stripped = _strip_dihedral_terms(lines)

    tmp_prmtop = output_prmtop.with_name(f'{output_prmtop.name}.{os.getpid()}.tmp')
    tmp_prmtop.write_text(''.join(stripped))
    os.replace(tmp_prmtop, output_prmtop)
    return output_prmtop.as_posix()


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
