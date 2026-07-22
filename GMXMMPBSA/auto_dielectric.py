"""
Helpers for automatic internal dielectric selection.
"""

from copy import deepcopy
from GMXMMPBSA.exceptions import InputError

AUTO_DIELECTRIC_VALUES = ('calculate', 'auto')

POSITIVE_AA = ['LYS', 'ARG', 'HIP']
NEGATIVE_AA = ['GLU', 'ASP']
NONPOLAR_AA = ['PHE', 'TRP', 'VAL', 'ILE', 'LEU', 'MET', 'PRO', 'CYX', 'ALA', 'GLY']
POLAR_AA = ['TYR', 'SER', 'THR', 'CYM', 'CYS', 'HIE', 'HID', 'GLN', 'ASN', 'ASH', 'GLH', 'LYN']

RESIDUE_CLASS_DIELECTRIC = {
    'nonpolar': 1.0,
    'polar': 3.0,
    'charged': 5.0,
}


def is_auto_dielectric(value):
    return isinstance(value, str) and value.strip().lower() in AUTO_DIELECTRIC_VALUES


def normalize_dielectric(value, name):
    if is_auto_dielectric(value):
        return value.strip().lower()
    try:
        return float(value)
    except (TypeError, ValueError) as exc:
        allowed = "', '".join(AUTO_DIELECTRIC_VALUES)
        raise InputError(f"{name} must be a number or one of '{allowed}'") from exc


def residue_class(resname):
    if resname in NONPOLAR_AA:
        return 'nonpolar'
    if resname in POLAR_AA:
        return 'polar'
    if resname in POSITIVE_AA or resname in NEGATIVE_AA:
        return 'charged'
    return None


def dielectric_for_residue(resname, INPUT):
    if resname in POLAR_AA:
        return INPUT['ala']['intdiel_polar'], 'polar'
    if resname in NONPOLAR_AA:
        return INPUT['ala']['intdiel_nonpolar'], 'nonpolar'
    if resname in POSITIVE_AA:
        return INPUT['ala']['intdiel_positive'], 'positive'
    if resname in NEGATIVE_AA:
        return INPUT['ala']['intdiel_negative'], 'negative'
    return None, None


def prepare_auto_dielectric(INPUT):
    requested = []
    fields = [('gb', 'intdiel'), ('pb', 'indi'), ('gbnsr6', 'epsin')]
    for nml, key in fields:
        value = normalize_dielectric(INPUT[nml][key], key)
        INPUT[nml][key] = value
        if is_auto_dielectric(value):
            requested.append((nml, key))

    if not requested:
        return

    INPUT['auto_dielectric'] = {
        'requested': requested,
        'done': False,
        'selected': None,
        'dominant': None,
        'totals': {},
        'print_res': None,
        'decomp': deepcopy(INPUT['decomp']),
    }

    INPUT['decomp']['decomprun'] = True
    INPUT['decomp']['idecomp'] = 2
    INPUT['decomp']['dec_verbose'] = 0
    if not INPUT['auto_dielectric']['decomp']['decomprun']:
        INPUT['decomp']['print_res'] = 'within 4'

    for nml, key in requested:
        INPUT[nml][key] = 5.0


def restore_auto_dielectric_decomp(INPUT):
    auto = INPUT.get('auto_dielectric')
    if not auto:
        return
    auto['print_res'] = INPUT['decomp']['print_res']
    INPUT['decomp'] = deepcopy(auto['decomp'])


def select_dielectric_from_decomp(decomp_delta):
    totals = dict.fromkeys(RESIDUE_CLASS_DIELECTRIC, 0.0)
    ignored = {}

    for res, terms in decomp_delta['TDC'].items():
        parts = res.split(':')
        resname = parts[2] if len(parts) >= 3 else ''
        cls = residue_class(resname)
        value = abs(float(terms['tot'].mean()))
        if cls is None:
            ignored[res] = value
            continue
        totals[cls] += value

    max_total = max(totals.values())
    if max_total == 0:
        raise InputError('Automatic dielectric selection found no classified residue contribution.')
    dominant_classes = [cls for cls, value in totals.items() if value == max_total]
    dominant = max(dominant_classes, key=lambda cls: RESIDUE_CLASS_DIELECTRIC[cls])
    selected = RESIDUE_CLASS_DIELECTRIC[dominant]

    return {
        'selected': selected,
        'dominant': dominant,
        'totals': totals,
        'ignored': ignored,
    }
