# ##############################################################################
#                           GPLv3 LICENSE INFO                                 #
# ##############################################################################

"""Load and validate the bundled gmx_MMPBSA_test manifest.

Selector resolution precedence:
1. Suite IDs 0, 1, 2 expand the suite test list.
2. Literal "101" resolves to suite "all" (same as 0).
3. Legacy test ID "11" resolves to the consolidated membrane test 6.
4. Numeric strings "3"..."26" resolve to individual tests.
5. Named aliases resolve to a test ID or suite name.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from importlib import resources
from typing import Any

__test__ = False


class ManifestError(ValueError):
    pass


@dataclass(frozen=True)
class TestEntry:
    id: int
    name: str
    path: str
    workdir: str
    input: str
    executable: str
    command_args: list[str]
    suites: list[str]
    slow: bool
    requires: list[str]
    expected_outputs: list[str]
    known_issues: tuple[str, ...] = ()


@dataclass(frozen=True)
class Manifest:
    version: int
    suites: dict[str, dict[str, Any]]
    aliases: dict[str, str | int]
    tests: dict[int, TestEntry]

    def get_test(self, test_id: int) -> TestEntry:
        try:
            return self.tests[test_id]
        except KeyError as exc:
            raise ManifestError(f'Unknown test id {test_id}') from exc

    def resolve_test_ids(self, selectors: list[str]) -> list[int]:
        if not selectors:
            raise ManifestError('No test was selected. Please define at least one test number')

        resolved: list[int] = []
        seen: set[int] = set()
        for selector in selectors:
            for test_id in self._expand_selector(selector):
                if test_id not in seen:
                    seen.add(test_id)
                    resolved.append(test_id)
        if not resolved:
            raise ManifestError('No test was selected. Please define at least one test number')
        return resolved

    def _expand_selector(self, selector: str) -> list[int]:
        token = str(selector).strip()
        if not token:
            raise ManifestError('Empty test selector')

        if token in {'0', '1', '2'}:
            return self._suite_tests_by_id(int(token))

        if token == '101':
            return self._suite_tests('all')

        if token in self.aliases:
            target = self.aliases[token]
            if isinstance(target, int):
                return [target]
            if isinstance(target, str):
                return self._suite_tests(target)
            raise ManifestError(f'Invalid alias target for {token!r}: {target!r}')

        if token.isdigit():
            test_id = int(token)
            if test_id in self.tests:
                return [test_id]
            raise ManifestError(f'Invalid test selector: {token}')

        raise ManifestError(f'Invalid test selector: {token}')

    def _suite_tests_by_id(self, suite_id: int) -> list[int]:
        for suite_name, suite in self.suites.items():
            if suite['id'] == suite_id:
                return list(suite['tests'])
        raise ManifestError(f'Unknown suite id {suite_id}')

    def _suite_tests(self, suite_name: str) -> list[int]:
        try:
            return list(self.suites[suite_name]['tests'])
        except KeyError as exc:
            raise ManifestError(f'Unknown suite {suite_name!r}') from exc

    def all_valid_choices(self) -> list[str]:
        choices = ['0', '1', '2', '101']
        choices.extend(str(test_id) for test_id in sorted(self.tests))
        choices.extend(name for name in self.aliases if name != '101')
        return list(dict.fromkeys(choices))


_MANIFEST: Manifest | None = None


def _manifest_resource_path() -> resources.Traversable:
    return resources.files('GMXMMPBSA').joinpath('data/gmx_MMPBSA_test_manifest.json')


def load_manifest() -> Manifest:
    global _MANIFEST
    if _MANIFEST is not None:
        return _MANIFEST

    with _manifest_resource_path().open('r', encoding='utf-8') as handle:
        raw = json.load(handle)

    version = raw.get('version')
    if version != 1:
        raise ManifestError(f'Unsupported manifest version: {version!r}')

    suites = raw.get('suites', {})
    aliases = raw.get('aliases', {})
    tests_raw = raw.get('tests', {})

    _validate_suites(suites)
    _validate_aliases(aliases, suites, tests_raw)

    tests: dict[int, TestEntry] = {}
    for key, entry in tests_raw.items():
        if not str(key).isdigit():
            raise ManifestError(f'Test keys must be numeric strings, got {key!r}')
        test_id = int(key)
        tests[test_id] = _parse_test_entry(test_id, entry)

    expected_ids = set(range(3, 27)) - {11}
    if set(tests) != expected_ids:
        missing = sorted(expected_ids - set(tests))
        extra = sorted(set(tests) - expected_ids)
        raise ManifestError(f'Manifest tests must be 3-26 except legacy alias 11. Missing={missing}, extra={extra}')

    for suite_name, suite in suites.items():
        for test_id in suite['tests']:
            if test_id not in tests:
                raise ManifestError(f'Suite {suite_name!r} references unknown test {test_id}')

    _MANIFEST = Manifest(version=version, suites=suites, aliases=aliases, tests=tests)
    return _MANIFEST


def _validate_suites(suites: dict[str, dict[str, Any]]) -> None:
    expected = {0, 1, 2}
    seen_ids: set[int] = set()
    for suite_name, suite in suites.items():
        suite_id = suite.get('id')
        tests = suite.get('tests')
        if not isinstance(suite_id, int):
            raise ManifestError(f'Suite {suite_name!r} has invalid id {suite_id!r}')
        if suite_id in seen_ids:
            raise ManifestError(f'Duplicate suite id {suite_id}')
        seen_ids.add(suite_id)
        if not isinstance(tests, list) or not tests:
            raise ManifestError(f'Suite {suite_name!r} must define a non-empty tests list')
    if seen_ids != expected:
        raise ManifestError(f'Suite ids must be exactly {sorted(expected)}')


def _validate_aliases(
    aliases: dict[str, str | int],
    suites: dict[str, dict[str, Any]],
    tests_raw: dict[str, Any],
) -> None:
    reserved = {'0', '1', '2', *tests_raw.keys()}
    suite_ids = {str(suite['id']) for suite in suites.values()}
    reserved |= suite_ids

    for alias, target in aliases.items():
        if alias in reserved:
            raise ManifestError(f'Alias {alias!r} collides with a suite or test id')
        if isinstance(target, int):
            if not (3 <= target <= 26):
                raise ManifestError(f'Alias {alias!r} target test id out of range: {target}')
        elif isinstance(target, str):
            if target not in suites:
                raise ManifestError(f'Alias {alias!r} references unknown suite {target!r}')
        else:
            raise ManifestError(f'Alias {alias!r} has invalid target {target!r}')


def _parse_test_entry(test_id: int, entry: dict[str, Any]) -> TestEntry:
    required = ('name', 'path', 'workdir', 'input', 'executable', 'command_args', 'expected_outputs')
    for key in required:
        if key not in entry or entry[key] in (None, '', []):
            raise ManifestError(f'Test {test_id} missing required field {key!r}')

    command_args = list(entry['command_args'])
    executable = entry['executable']
    if not isinstance(command_args, list) or not all(isinstance(arg, str) for arg in command_args):
        raise ManifestError(f'Test {test_id} command_args must be a list of strings')

    expected_outputs = list(entry['expected_outputs'])
    if not expected_outputs:
        raise ManifestError(f'Test {test_id} must define expected_outputs')

    requires = list(entry.get('requires') or [executable, 'mpirun'])
    known_issues = tuple(entry.get('known_issues') or ())

    return TestEntry(
        id=test_id,
        name=entry['name'],
        path=entry['path'],
        workdir=entry['workdir'],
        input=entry['input'],
        executable=executable,
        command_args=command_args,
        suites=list(entry.get('suites') or []),
        slow=bool(entry.get('slow', False)),
        requires=requires,
        expected_outputs=expected_outputs,
        known_issues=known_issues,
    )


def build_help_text() -> str:
    manifest = load_manifest()
    lines = [
        'The level the test is going to be run at. Multiple systems and analysis can be run at the same time.',
        '      Nr. of Sys  ',
        '* 0      22     All -- Run all examples (Can take a long time!!!)',
        '* 1      11     Minimal -- Does a minimal test with a set of systems and analyzes',
        '                that show that gmx_MMPBSA runs correctly. Comp_receptor is excluded',
        '                until its required GROMACS topology is shipped; slow or redundant',
        '                cases are omitted from this suite',
        '* 2       8     Fast -- Only the calculations that take a short time are run (Default)',
        '[Systems]:',
        '     Slow Frames',
    ]

    system_ids = [3, 4, 5, 6, 7, 8, 9, 10]
    for test_id in system_ids:
        lines.append(_format_help_line(manifest.get_test(test_id)))
    lines.append('* 11       |     Legacy alias for test 6 (consolidated membrane example)')

    lines.extend([
        '[Analysis]:',
        '     Slow Frames',
    ])
    for test_id in range(12, 27):
        lines.append(_format_help_line(manifest.get_test(test_id)))

    return '\n'.join(lines)


def _format_help_line(test: TestEntry) -> str:
    slow_mark = 'x' if test.slow else '.'
    frame_count = _estimate_frame_count_label(test)
    return f'* {test.id:<4} {slow_mark} | {frame_count:<3} {test.name}'


def _estimate_frame_count_label(test: TestEntry) -> str:
    # Help text historically shows approximate frame counts; keep simple labels.
    defaults = {
        6: 4, 8: 10, 10: 4, 11: 4, 15: 10, 17: 10, 18: 4, 25: 5,
    }
    return str(defaults.get(test.id, 10))


def snapshot_outputs(work_dir, expected_outputs: list[str]) -> dict[str, tuple[int, int, int] | None]:
    """Capture metadata used to distinguish fresh outputs from stale files."""
    from pathlib import Path

    base = Path(work_dir)
    snapshot = {}
    for name in expected_outputs:
        path = base / name
        try:
            stat = path.stat()
        except OSError:
            snapshot[name] = None
        else:
            snapshot[name] = (stat.st_ino, stat.st_size, stat.st_mtime_ns)
    return snapshot


def verify_outputs(
    work_dir,
    expected_outputs: list[str],
    previous: dict[str, tuple[int, int, int] | None] | None = None,
) -> list[str]:
    """Return expected outputs that are missing or unchanged since ``previous``."""
    from pathlib import Path

    base = Path(work_dir)
    missing = []
    for name in expected_outputs:
        path = base / name
        try:
            stat = path.stat()
        except OSError:
            missing.append(name)
            continue
        if not path.is_file():
            missing.append(name)
            continue
        if previous is not None and previous.get(name) is not None:
            current = (stat.st_ino, stat.st_size, stat.st_mtime_ns)
            if current == previous[name]:
                missing.append(name)
    return missing
