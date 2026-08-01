#!/usr/bin/env python3
"""Validate gmx_MMPBSA_test manifest against example READMEs and docs."""

from __future__ import annotations

import argparse
import shlex
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from GMXMMPBSA.test_manifest import load_manifest

EXAMPLES = REPO / 'examples'
DOCS = REPO / 'docs' / 'examples' / 'gmx_MMPBSA_test.md'


def _normalize_command(command: str) -> str:
    return ' '.join(shlex.split(command.strip()))


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
    return None


def _extract_gmx_test_tab(readme_path: Path) -> str:
    text = readme_path.read_text()
    in_tab = False
    lines = []
    for line in text.splitlines():
        if '=== "gmx_MMPBSA_test"' in line:
            in_tab = True
            continue
        if in_tab and line.startswith('==='):
            break
        if in_tab and line.strip():
            lines.append(line.strip())
    return '\n'.join(lines)


def _validate_readmes(manifest, examples_dir: Path, errors: list[str]) -> None:
    for test_id, test in manifest.tests.items():
        readme = examples_dir / test.workdir / 'README.md'
        if not readme.exists():
            errors.append(f'Test {test_id}: missing README {readme}')
            continue

        tab = _extract_gmx_test_tab(readme)
        expected_tab = f'gmx_MMPBSA_test -t {test_id}'
        if expected_tab not in tab:
            errors.append(f'Test {test_id}: README tab missing {expected_tab!r} in {readme}')

        serial = _parse_serial_command(readme)
        if not serial:
            errors.append(f'Test {test_id}: no Serial command in {readme}')
            continue

        manifest_cmd = test.executable + ' ' + ' '.join(shlex.quote(arg) for arg in test.command_args)
        if _normalize_command(serial) != _normalize_command(manifest_cmd):
            errors.append(
                f'Test {test_id}: Serial command mismatch\n'
                f'  README:   {_normalize_command(serial)}\n'
                f'  manifest: {_normalize_command(manifest_cmd)}'
            )


def _validate_index_footnotes(manifest, examples_dir: Path, errors: list[str]) -> None:
    index = examples_dir / 'README.md'
    if not index.exists():
        errors.append(f'Missing examples index: {index}')
        return

    suite_footnotes = {'all': '[^1]', 'minimal': '[^2]', 'fast': '[^3]'}
    lines = index.read_text().splitlines()
    for test_id, test in manifest.tests.items():
        expected = {suite_footnotes[suite] for suite in test.suites if suite in suite_footnotes}
        if not expected:
            continue

        link = f'({test.path}/README.md)'
        matches = [line for line in lines if link in line]
        if not matches:
            errors.append(f'Test {test_id}: missing examples index link {link!r} in {index}')
            continue

        if not any(expected.issubset(_footnotes_in_line(line)) for line in matches):
            errors.append(
                f'Test {test_id}: examples index link {link!r} missing suite footnotes '
                f'{", ".join(sorted(expected))}'
            )


def _footnotes_in_line(line: str) -> set[str]:
    return {f'[^{num}]' for num in ('1', '2', '3') if f'[^{num}]' in line}

def _validate_docs(errors: list[str]) -> None:
    if not DOCS.exists():
        errors.append(f'Missing docs file: {DOCS}')
        return

    text = DOCS.read_text()
    required_tokens = [str(i) for i in range(3, 27)] + ['101']
    usage_block = text.split('optional arguments:', 1)[0]
    option_block = text.split('Test options:', 1)[-1]

    for token in required_tokens:
        if token not in usage_block:
            errors.append(f'Docs usage block missing test selector {token!r}')
        if token not in option_block:
            errors.append(f'Docs option-help block missing test selector {token!r}')


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--examples-dir', type=Path, default=EXAMPLES)
    args = parser.parse_args()
    examples_dir = args.examples_dir.resolve()
    if not examples_dir.is_dir():
        print(f'Examples directory does not exist: {examples_dir}', file=sys.stderr)
        return 1

    manifest = load_manifest()
    errors: list[str] = []
    _validate_readmes(manifest, examples_dir, errors)
    _validate_index_footnotes(manifest, examples_dir, errors)
    _validate_docs(errors)

    if errors:
        print('gmx_MMPBSA_test documentation validation failed:', file=sys.stderr)
        for error in errors:
            print(f'- {error}', file=sys.stderr)
        return 1

    print('gmx_MMPBSA_test documentation validation passed.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
