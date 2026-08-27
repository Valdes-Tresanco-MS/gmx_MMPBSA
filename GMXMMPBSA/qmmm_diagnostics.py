"""QM/MM diagnostics extracted from SANDER output."""

from dataclasses import dataclass
import re


@dataclass(frozen=True)
class QMMMDiagnostic:
    """A diagnostic extracted from a SANDER QM/MM output file."""

    code: str
    severity: str
    message: str
    remediation: str


def parse_qmmm_diagnostics(output_text, qm_theory=None):
    """Classify known QM/MM SANDER messages.

    The parser deliberately reports only messages that are actionable for a
    QM/MMGBSA user.  It does not infer universal method support from a single
    output file: parameter availability remains dependent on the SANDER build,
    its data files, and the elements in the selected QM region.
    """
    diagnostics = []

    scf_failures = re.findall(
        r'(?:No convergence in SCF after\s+\d+\s+steps?|'
        r'SCC-DFTB[^\n]*DID NOT CONVERGE[^\n]*)',
        output_text,
        re.IGNORECASE,
    )
    if scf_failures:
        diagnostics.append(QMMMDiagnostic(
            'scf_nonconvergence',
            'error',
            f'QM/MM SCF did not converge in {len(scf_failures)} step(s).',
            'Increase itrmax or adjust scfconv, then verify that all frames converge.',
        ))

    parameter_matches = list(re.finditer(
        r'There are no\s+(?P<method>.+?)\s+parameters for\s+'
        r'(?:(?:atomic number)\s+(?P<atomic_number>\d+)|this element)',
        output_text,
        re.IGNORECASE,
    ))
    for match in parameter_matches:
        method = match.group('method').strip()
        atomic_number = match.group('atomic_number')
        if atomic_number is None:
            context = output_text[max(0, match.start() - 500):match.start()]
            atomic_numbers = re.findall(r'has atomic number\s+(\d+)', context, re.IGNORECASE)
            atomic_number = atomic_numbers[-1] if atomic_numbers else None
        atom_text = f' for atomic number {atomic_number}' if atomic_number else ''
        diagnostics.append(QMMMDiagnostic(
            'qm_parameter_missing',
            'error',
            f'SANDER has no {method} QM parameters{atom_text}.',
            'Choose a method with parameters for every element in the QM region or change the QM region.',
        ))

    if not parameter_matches:
        unavailable = re.search(
            r'QM\s+(?P<method>[^\n]+?)\s+NOT AVAILABLE FOR THIS ATOM',
            output_text,
            re.IGNORECASE,
        )
        if unavailable:
            diagnostics.append(QMMMDiagnostic(
                'qm_parameter_missing',
                'error',
                f'SANDER reports QM method {unavailable.group("method").strip()} is not available for an atom in the QM region.',
                'Choose a method with parameters for every element in the QM region or change the QM region.',
            ))

    if re.search(
        r'Parameters for dispersion correction are not available for this atom',
        output_text,
        re.IGNORECASE,
    ):
        method_text = f' for {qm_theory}' if qm_theory else ''
        diagnostics.append(QMMMDiagnostic(
            'dispersion_parameter_missing',
            'error',
            f'Dispersion-correction parameters are missing{method_text} for an atom in the QM region.',
            'Use the correction only with a SANDER parameter set that covers the QM-region elements, or choose the base method.',
        ))

    dftb_file = re.search(
        r'Missing file:\s*(?:\n\s*)?(?P<path>[^\s]+\.skf)',
        output_text,
        re.IGNORECASE,
    )
    if dftb_file:
        diagnostics.append(QMMMDiagnostic(
            'dftb_parameter_file_missing',
            'error',
            f'DFTB Slater-Koster parameter file is missing: {dftb_file.group("path")}.',
            'Install the required DFTB parameter file for every QM-region element pair or choose another method.',
        ))

    if re.search(
        r'Analytical derivatives for d orbitals are not supported',
        output_text,
        re.IGNORECASE,
    ):
        diagnostics.append(QMMMDiagnostic(
            'numerical_qm_derivatives',
            'warning',
            'SANDER is using numerical derivatives because analytical d-orbital derivatives are unavailable.',
            'This is not a correctness failure, but it can substantially increase QM/MM runtime.',
        ))

    return diagnostics
