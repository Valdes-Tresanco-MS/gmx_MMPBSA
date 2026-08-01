#!/usr/bin/env python3
"""Validate docs/examples README copies match transformed examples/ READMEs."""

from __future__ import annotations

import argparse
import difflib
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
SCRIPTS = REPO / 'scripts'
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from example_readme_sync import (  # noqa: E402
    docs_readme_path,
    is_skipped_docs_target,
    iter_example_readmes,
    transform_readme_content,
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--examples-dir', type=Path, default=REPO / 'examples')
    parser.add_argument('--docs-dir', type=Path, default=REPO / 'docs' / 'examples')
    args = parser.parse_args()

    examples_root = args.examples_dir.resolve()
    docs_root = args.docs_dir.resolve()
    errors: list[str] = []

    for examples_readme in iter_example_readmes(examples_root):
        docs_readme = docs_readme_path(examples_root, docs_root, examples_readme)
        rel = docs_readme.relative_to(docs_root)
        if is_skipped_docs_target(docs_readme, docs_root):
            continue

        expected = transform_readme_content(examples_readme.read_text(encoding='utf-8'))
        if not docs_readme.exists():
            errors.append(f'{rel}: missing docs copy')
            continue

        actual = docs_readme.read_text(encoding='utf-8')
        if actual != expected:
            diff = difflib.unified_diff(
                actual.splitlines(),
                expected.splitlines(),
                fromfile=f'docs/{rel.as_posix()}',
                tofile=f'expected/{rel.as_posix()}',
                lineterm='',
            )
            errors.append(f'{rel}: content mismatch\n' + '\n'.join(diff))

    if errors:
        print('Example README parity validation failed:', file=sys.stderr)
        for error in errors:
            print(f'- {error}', file=sys.stderr)
        return 1

    print('Example README parity validation passed.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
