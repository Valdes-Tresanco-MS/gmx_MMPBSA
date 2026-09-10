"""Continuum-radius provenance and audit helpers.

The values used by Amber's continuum-solvation routines live in the final
prmtop arrays.  This module records those values without attempting to infer
or silently change a force-field/radius parameterization.
"""

import csv
import hashlib
import json
import logging
from pathlib import Path

import parmed


RADIUS_NAMES = {
    1: 'bondi',
    2: 'mbondi',
    3: 'mbondi2',
    4: 'mbondi3',
    5: 'mbondi_pb2',
    6: 'mbondi_pb3',
    7: 'charmm_radii',
}

GB_RECOMMENDED_RADII = {1: 'mbondi', 2: 'mbondi2', 5: 'mbondi2', 7: 'bondi', 8: 'mbondi3'}

_METAL_SYMBOLS = {
    'LI', 'BE', 'NA', 'MG', 'AL', 'K', 'CA', 'SC', 'TI', 'V', 'CR', 'MN', 'FE', 'CO', 'NI',
    'CU', 'ZN', 'GA', 'RB', 'SR', 'Y', 'ZR', 'NB', 'MO', 'TC', 'RU', 'RH', 'PD', 'AG', 'CD',
    'IN', 'SN', 'CS', 'BA', 'LA', 'CE', 'PR', 'ND', 'PM', 'SM', 'EU', 'GD', 'TB', 'DY', 'HO',
    'ER', 'TM', 'YB', 'LU', 'HF', 'TA', 'W', 'RE', 'OS', 'IR', 'PT', 'AU', 'HG', 'TL', 'PB',
    'BI', 'PO', 'AC', 'TH', 'PA', 'U', 'NP', 'PU', 'AM', 'CM', 'BK', 'CF', 'ES', 'FM', 'MD',
    'NO', 'LR',
}
_STANDARD_RESIDUES = {
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET',
    'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'ASH', 'CYM', 'CYX', 'HID', 'HIE', 'HIP',
    'LYN', 'GLH', 'HSD', 'HSE', 'HSP', 'ACE', 'NME', 'NHE', 'DA', 'DC', 'DG', 'DT', 'A', 'C',
    'G', 'U', 'RA', 'RC', 'RG', 'RU', 'WAT', 'HOH', 'SOL', 'TIP3P', 'TIP4P', 'OPC', 'SPC',
}
_ELEMENT_SYMBOLS = {
    1: 'H', 2: 'HE', 3: 'LI', 4: 'BE', 5: 'B', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 10: 'NE',
    11: 'NA', 12: 'MG', 13: 'AL', 14: 'SI', 15: 'P', 16: 'S', 17: 'CL', 18: 'AR',
    19: 'K', 20: 'CA', 26: 'FE', 27: 'CO', 28: 'NI', 29: 'CU', 30: 'ZN', 35: 'BR', 53: 'I',
}


def radius_name(radius_set):
    """Return the canonical radius name embedded in ``RADIUS_SET``."""
    text = str(radius_set or '').lower()
    for name in sorted(RADIUS_NAMES.values(), key=len, reverse=True):
        if f'({name})' in text or name in text:
            return name
    return 'unknown'


def radius_checksum(values):
    """Return a stable SHA-256 checksum for a numeric prmtop array."""
    payload = json.dumps([float(value) for value in values], separators=(',', ':'), allow_nan=False)
    return hashlib.sha256(payload.encode('ascii')).hexdigest()


def _input_value(input_data, key, default=None):
    general = input_data.get('general', {})
    if key in general:
        return general[key]
    # Accept the longer spelling in programmatic callers while keeping one
    # documented input keyword.
    return general.get('source_force_field_family', default) if key == 'source_force_field' else default


def source_force_field_family(input_data, parm=None):
    """Classify the source force-field family, honoring the explicit override."""
    override = str(_input_value(input_data, 'source_force_field', 'auto') or 'auto').strip().lower()
    if override != 'auto':
        return override

    if parm is not None and getattr(parm, 'chamber', False):
        return 'charmm'

    forcefields = _input_value(input_data, 'forcefields', []) or []
    if isinstance(forcefields, str):
        forcefields = [forcefields]
    text = ' '.join(str(value).lower() for value in forcefields)
    if 'gromos' in text:
        return 'gromos'
    if 'charmm' in text or 'cgenff' in text:
        return 'charmm'
    if 'opls' in text:
        return 'opls'
    if text:
        return 'amber'
    return 'unknown'


def atom_representation(parm):
    """Classify the topology conservatively as all-atom, united-atom, or unknown."""
    atoms = list(getattr(parm, 'atoms', []))
    if not atoms:
        return 'unknown'
    known = [getattr(atom, 'atomic_number', None) for atom in atoms]
    if all(number not in (None, 0) for number in known):
        return 'all_atom'
    if any(number in (None, 0) for number in known):
        return 'united_atom_or_unknown'
    return 'unknown'


def _element(atom):
    value = getattr(atom, 'element', None)
    if isinstance(value, str) and value.strip() and not value.strip().isdigit():
        symbol = value.strip().upper()
        return symbol if symbol in {s.upper() for s in parmed.periodic_table.Element[1:]} else ''
    number = value if value not in (None, '', 0) else getattr(atom, 'atomic_number', None)
    try:
        number = int(number)
        if 0 < number < len(parmed.periodic_table.Element):
            return parmed.periodic_table.Element[number].upper()
    except (TypeError, ValueError):
        pass
    return ''


def _lj_size(atom):
    for attribute in ('sigma', 'rmin', 'radius'):
        value = getattr(atom, attribute, None)
        if value is not None:
            try:
                return float(value)
            except (TypeError, ValueError):
                pass
    return None


def _atom_flags(atom):
    element = _element(atom)
    residue = getattr(getattr(atom, 'residue', None), 'name', '') or ''
    atom_name = str(getattr(atom, 'name', '') or '').strip().upper()
    atom_type = str(getattr(atom, 'type', '') or '').strip().upper()
    flags = []
    if not element:
        flags.append('unknown_element')
    if element in _METAL_SYMBOLS:
        flags.append('metal')
    if atom_name in {'EP', 'LP', 'DU', 'DUM', 'X', 'V'} or atom_type in {'EP', 'LP', 'DU', 'DUM'}:
        flags.append('dummy_or_extra_point')
    if str(residue).strip().upper() not in _STANDARD_RESIDUES:
        flags.append('nonstandard_residue')
    return flags


def _assignment_category(flags):
    if 'metal' in flags:
        return 'metal_unvalidated'
    if 'dummy_or_extra_point' in flags:
        return 'dummy_or_extra_point'
    if 'unknown_element' in flags:
        return 'unknown_element'
    if 'nonstandard_residue' in flags:
        return 'nonstandard_residue'
    return 'radius_set_rule'


def _models(input_data):
    models = {}
    if input_data.get('gb', {}).get('gbrun'):
        models['gb'] = f"igb={input_data['gb'].get('igb')}"
    if input_data.get('pb', {}).get('pbrun'):
        models['pb'] = f"ipb={input_data['pb'].get('ipb')}"
    if input_data.get('gbnsr6', {}).get('gbnsr6run'):
        models['gbnsr6'] = 'GBNSR6'
    return models


def _record(parm, component, requested, route, input_data, source_family):
    parm_data = getattr(parm, 'parm_data', {})
    radii = parm_data.get('RADII', [])
    screen = parm_data.get('SCREEN', [])
    radius_set = parm_data.get('RADIUS_SET', ['unknown'])
    radius_set = radius_set[0] if radius_set else 'unknown'
    models = _models(input_data)
    radius_known = radius_name(radius_set) != 'unknown'
    flags = [_atom_flags(atom) + ([] if radius_known else ['assignment_rule_unknown'])
             for atom in getattr(parm, 'atoms', [])]
    flat_flags = [flag for atom_flags in flags for flag in atom_flags]
    return {
        'component': component,
        'requested_radius_set': requested,
        'effective_radius_set': radius_name(radius_set),
        'RADIUS_SET': radius_set,
        'parmed_version': getattr(parmed, '__version__', 'unknown'),
        'assignment_route': route,
        'input_pbradii_applied': route == 'parmed_ChRad',
        'effective_radius_source': (
            'topology RADII/SCREEN preserved' if route == 'native_amber_topology_preserved'
            else 'Inherited radius family reapplied through ParmEd ChRad'
            if route == 'native_amber_mutant_inherited_ChRad'
            else 'GBNSR6-compatible topology copy derived from the prepared topology'
            if route == 'gbnsr6_prepared_copy'
            else 'PBRadii applied through ParmEd ChRad'
        ),
        'source_force_field_family': source_family,
        'atom_representation': atom_representation(parm),
        'models': models,
        'gb_model': models.get('gb'),
        'pb_model': models.get('pb'),
        'gbnsr6_model': models.get('gbnsr6'),
        'radii_sha256': radius_checksum(radii),
        'screen_sha256': radius_checksum(screen),
        'atom_count': len(getattr(parm, 'atoms', [])),
        'flag_counts': {flag: flat_flags.count(flag) for flag in sorted(set(flat_flags))},
    }


def _write_audit_csv(parm, component, output_path):
    fields = [
        'component', 'atom_index', 'residue_index', 'residue_name', 'atom_name', 'element',
        'source_atom_type', 'charge', 'lj_size', 'assigned_continuum_radius', 'screening_value',
        'assignment_category', 'flags',
    ]
    parm_data = getattr(parm, 'parm_data', {})
    radii = parm_data.get('RADII', [])
    screen = parm_data.get('SCREEN', [])
    radius_set = parm_data.get('RADIUS_SET', ['unknown'])
    radius_set = radius_set[0] if radius_set else 'unknown'
    radius_known = radius_name(radius_set) != 'unknown'
    with open(output_path, 'w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for index, atom in enumerate(getattr(parm, 'atoms', [])):
            residue = getattr(atom, 'residue', None)
            flags = _atom_flags(atom) + ([] if radius_known else ['assignment_rule_unknown'])
            writer.writerow({
                'component': component,
                'atom_index': index + 1,
                'residue_index': getattr(residue, 'number', getattr(residue, 'idx', '')),
                'residue_name': getattr(residue, 'name', ''),
                'atom_name': getattr(atom, 'name', ''),
                'element': _element(atom),
                'source_atom_type': getattr(atom, 'type', ''),
                'charge': getattr(atom, 'charge', ''),
                'lj_size': _lj_size(atom),
                'assigned_continuum_radius': radii[index] if index < len(radii) else '',
                'screening_value': screen[index] if index < len(screen) else '',
                'assignment_category': _assignment_category(flags),
                'flags': '|'.join(flags),
            })


def _load_prepared_topology_for_provenance(path, source_path=None):
    """Load a prepared topology while tolerating GBNSR6's reduced pointers.

    GBNSR6-compatible copies intentionally retain parameter arrays whose
    corresponding dihedral pointers are zeroed. ParmEd therefore rejects
    them as ordinary Amber topologies. The atom/radius data are unchanged by
    preparation, so use the source topology for atom metadata and replace its
    radius arrays with the raw arrays from the prepared copy.
    """
    path = str(path)
    source_path = str(source_path) if source_path is not None else None
    try:
        return parmed.load_file(path)
    except parmed.exceptions.AmberError:
        if source_path is None:
            raise
        source = parmed.load_file(source_path)
        raw = parmed.amber.AmberFormat(str(path))
        for key in ('RADIUS_SET', 'RADII', 'SCREEN'):
            if key in raw.parm_data:
                source.parm_data[key] = raw.parm_data[key]
        return source


def collect_radii_provenance(files, input_data, engine, additional_topologies=None):
    """Collect final topology provenance and write the requested artifacts."""
    requested = RADIUS_NAMES.get(input_data.get('general', {}).get('PBRadii'), 'unknown')
    route = 'native_amber_topology_preserved' if engine == 'amber' else 'parmed_ChRad'
    components = {
        'complex': getattr(files, 'complex_prmtop', None),
        'receptor': getattr(files, 'receptor_prmtop', None),
        'ligand': getattr(files, 'ligand_prmtop', None),
    }
    mutants = {
        'mutant_complex': getattr(files, 'mutant_complex_prmtop', None),
        'mutant_receptor': getattr(files, 'mutant_receptor_prmtop', None),
        'mutant_ligand': getattr(files, 'mutant_ligand_prmtop', None),
    }
    records = {}
    for component_group in (components, mutants):
        for component, path in component_group.items():
            if not path or not Path(path).exists():
                continue
            parm = parmed.load_file(path)
            family = source_force_field_family(input_data, parm)
            component_route = route
            if engine == 'amber' and component.startswith('mutant_'):
                normal_path = components.get(component.removeprefix('mutant_'))
                # Unchanged components may alias the normal topology after setup.
                unchanged = normal_path and Path(normal_path).exists() and Path(path).samefile(normal_path)
                if not unchanged:
                    component_route = 'native_amber_mutant_inherited_ChRad'
            records[component] = _record(parm, component, requested, component_route, input_data, family)
            _log_advisories(records[component], input_data)
            if input_data.get('general', {}).get('radii_audit', 0):
                audit_path = Path(f'GMXMMPBSA_radii_{component}.csv')
                _write_audit_csv(parm, component, audit_path)
                records[component]['audit_csv'] = str(audit_path)

    prepared_records = {}
    for component, topology in (additional_topologies or {}).items():
        if isinstance(topology, (tuple, list)):
            path, source_path = topology
        else:
            path, source_path = topology, None
        if not path or not Path(path).exists():
            continue
        parm = _load_prepared_topology_for_provenance(path, source_path)
        family = source_force_field_family(input_data, parm)
        prepared_records[component] = _record(
            parm, component, requested, 'gbnsr6_prepared_copy', input_data, family
        )
        prepared_records[component]['topology_path'] = str(path)
        if source_path is not None:
            prepared_records[component]['source_topology_path'] = str(source_path)
        _log_advisories(prepared_records[component], input_data)
        if input_data.get('general', {}).get('radii_audit', 0):
            audit_path = Path(f'GMXMMPBSA_radii_{component}.csv')
            _write_audit_csv(parm, component, audit_path)
            prepared_records[component]['audit_csv'] = str(audit_path)

    provenance = {
        'schema_version': 1,
        'requested_radius_set': requested,
        'source_force_field_override': _input_value(input_data, 'source_force_field', 'auto'),
        'input_pbradii_semantics': (
            'requested value is advisory; native topology RADII/SCREEN are effective for normal components; '
            'mutant topologies reapply the inherited radius family through ParmEd ChRad'
            if route == 'native_amber_topology_preserved'
            else 'requested value was applied while preparing the generated topology'
        ),
        'models': _models(input_data),
        'components': records,
        'prepared_topologies': prepared_records,
    }
    output_path = Path('GMXMMPBSA_radii.json')
    output_path.write_text(json.dumps(provenance, indent=2, sort_keys=True) + '\n', encoding='utf-8')
    logging.info('Continuum-radius provenance written to %s', output_path)
    for component, record in records.items():
        logging.info(
            'Radii %-15s requested=%s effective=%s RADIUS_SET=%r route=%s source=%s '
            'RADII=%s SCREEN=%s',
            component, record['requested_radius_set'], record['effective_radius_set'], record['RADIUS_SET'],
            record['assignment_route'], record['source_force_field_family'], record['radii_sha256'],
            record['screen_sha256'],
        )
    for component, record in prepared_records.items():
        logging.info(
            'Radii %-15s requested=%s effective=%s RADIUS_SET=%r route=%s source=%s '
            'RADII=%s SCREEN=%s',
            component, record['requested_radius_set'], record['effective_radius_set'], record['RADIUS_SET'],
            record['assignment_route'], record['source_force_field_family'], record['radii_sha256'],
            record['screen_sha256'],
        )
    return provenance


def _log_advisories(record, input_data):
    family = record['source_force_field_family']
    effective = record['effective_radius_set']
    models = record['models']
    if family == 'charmm' and 'gb' in models and effective.startswith(('bondi', 'mbondi')):
        logging.warning(
            'CHARMM topology with AMBER %s radii for GB is a cross-parameterization protocol; '
            'a native CHARMM implicit-solvent parameterization has not been established here.', effective
        )
    if family == 'charmm' and 'pb' in models and effective != 'charmm_radii':
        logging.info(
            'CHARMM PB topology uses %s radii; charmm_radii is available as the CHARMM-specific PB option '
            'but will not be selected automatically.', effective
        )
    if family == 'opls' and 'gb' in models and effective.startswith(('bondi', 'mbondi')):
        logging.warning(
            'OPLS topology with AMBER %s radii for GB is empirically unvalidated; interpret the calculation '
            'as a calibrated scoring protocol.', effective
        )
    if family == 'gromos' or record['atom_representation'] == 'united_atom_or_unknown':
        logging.warning(
            'GROMOS/united-atom or incompletely typed topology detected; standard all-atom continuum-radius '
            'rules have strong experimental-support limitations for this input.'
        )
    if 'gbnsr6' in models:
        logging.info(
            'GBNSR6 radius provenance is reported independently; pairwise-GB igb/radius compatibility rules '
            'are not applied to GBNSR6.'
        )
    igb = input_data.get('gb', {}).get('igb')
    recommended = GB_RECOMMENDED_RADII.get(igb)
    if 'gb' in models and recommended and effective != recommended:
        nuance = ''
        if igb == 8 and effective == 'mbondi2':
            nuance = ' Historical mbondi2 energy-only usage exists, but it is not the conventional igb=8 pairing.'
        logging.warning(
            'The %s topology uses %s radii, while igb=%s is conventionally used with %s radii.%s '
            'Stored topology values are preserved.',
            record['component'], effective, igb, recommended, nuance,
        )
