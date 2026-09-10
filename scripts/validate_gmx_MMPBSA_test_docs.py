#!/usr/bin/env python3
"""Validate gmx_MMPBSA_test manifest against example READMEs and docs."""

from __future__ import annotations

import argparse
from collections.abc import Sequence
import re
import shlex
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from GMXMMPBSA.test_manifest import load_manifest

EXAMPLES = REPO / 'examples'
DOCS = REPO / 'docs' / 'examples' / 'gmx_MMPBSA_test.md'


def _command_tokens(command: str | Sequence[str]) -> list[str]:
    if isinstance(command, str):
        # A shell continuation is part of one command. Removing it before
        # shlex.split avoids treating a captured first line ending in '\\' as
        # an unterminated escape.
        command = re.sub(r'\\\s*\n', ' ', command.strip())
        return shlex.split(command)
    return list(command)


def _normalize_command(command: str | Sequence[str]) -> str:
    return shlex.join(_command_tokens(command))


def _command_signature(command: str | Sequence[str]) -> tuple[str, tuple[tuple[str, tuple[str, ...]], ...]]:
    """Return an order-independent executable/options representation."""
    tokens = _command_tokens(command)
    if not tokens:
        return '', ()

    executable = tokens[0]
    options: list[tuple[str, tuple[str, ...]]] = []
    index = 1
    while index < len(tokens):
        option = tokens[index]
        if not option.startswith('-'):
            raise ValueError(f'Unexpected positional argument {option!r}')
        index += 1
        values: list[str] = []
        while index < len(tokens) and not tokens[index].startswith('-'):
            values.append(tokens[index])
            index += 1
        options.append((option, tuple(values)))
    return executable, tuple(sorted(options))


def _section_text(text: str, heading: str) -> str:
    lines = text.splitlines()
    for index, line in enumerate(lines):
        if line.strip() != heading:
            continue
        section = []
        for candidate in lines[index + 1:]:
            stripped = candidate.strip()
            if (stripped.startswith('#') or stripped.startswith('=== ')) and stripped != heading:
                break
            section.append(candidate)
        return '\n'.join(section)
    return ''


def _command_from_lines(lines: list[str], prefixes: tuple[str, ...]) -> str | None:
    for index, line in enumerate(lines):
        candidate = line.strip()
        if not candidate.startswith(prefixes):
            continue

        command = [candidate]
        next_index = index + 1
        while command[-1].rstrip().endswith('\\') and next_index < len(lines):
            continuation = lines[next_index].strip()
            next_index += 1
            if continuation:
                command.append(continuation)
        return '\n'.join(command)
    return None


def _command_from_section(section: str, prefixes: tuple[str, ...]) -> str | None:
    fenced_blocks = re.findall(r'```[^\n]*\n(.*?)```', section, flags=re.DOTALL)
    for block in fenced_blocks:
        command = _command_from_lines(block.splitlines(), prefixes)
        if command:
            return command
    return _command_from_lines(section.splitlines(), prefixes)


def _parse_serial_command(readme_path: Path) -> str | None:
    text = readme_path.read_text()
    section = _section_text(text, '=== "Serial"')
    command = _command_from_section(section, ('gmx_MMPBSA ', 'amber_MMPBSA ', 'ggmx_MMPBSA '))
    if command and command.startswith('ggmx_MMPBSA '):
        return command.replace('ggmx_MMPBSA', 'gmx_MMPBSA', 1)
    return command


def _parse_bundled_test_command(readme_path: Path) -> str | None:
    text = readme_path.read_text()
    # Current guides use a dedicated heading and fenced command.
    section = _section_text(text, '### Run the bundled test')
    command = _command_from_section(section, ('gmx_MMPBSA_test ',))
    if command:
        return command

    # Retain compatibility with the two legacy guides that still use a tab.
    section = _section_text(text, '=== "gmx_MMPBSA_test"')
    return _command_from_section(section, ('gmx_MMPBSA_test ',))


def _test_selector(command: str) -> str | None:
    tokens = _command_tokens(command)
    if not tokens or tokens[0] != 'gmx_MMPBSA_test':
        return None
    for index, token in enumerate(tokens[:-1]):
        if token in {'-t', '--test'}:
            return tokens[index + 1]
    return None


def _validate_readmes(manifest, examples_dir: Path, errors: list[str]) -> None:
    for test_id, test in manifest.tests.items():
        readme = examples_dir / test.workdir / 'README.md'
        if not readme.exists():
            errors.append(f'Test {test_id}: missing README {readme}')
            continue

        bundled_test = _parse_bundled_test_command(readme)
        selector = _test_selector(bundled_test) if bundled_test else None
        if selector != str(test_id):
            errors.append(
                f'Test {test_id}: bundled-test command selects {selector!r}; '
                f'expected {test_id} in {readme}'
            )

        serial = _parse_serial_command(readme)
        if not serial:
            errors.append(f'Test {test_id}: no Serial command in {readme}')
            continue

        manifest_cmd = [test.executable, *test.command_args]
        if _command_signature(serial) != _command_signature(manifest_cmd):
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
