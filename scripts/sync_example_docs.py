#!/usr/bin/env python3
"""Sync example READMEs from examples/ to docs/examples/ for MkDocs."""

from __future__ import annotations

import argparse
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


def sync_readmes(examples_root: Path, docs_root: Path, check_only: bool) -> list[str]:
    errors: list[str] = []
    changed = 0

    for examples_readme in iter_example_readmes(examples_root):
        docs_readme = docs_readme_path(examples_root, docs_root, examples_readme)
        if is_skipped_docs_target(docs_readme, docs_root):
            continue

        source = transform_readme_content(examples_readme.read_text(encoding='utf-8'))
        current = docs_readme.read_text(encoding='utf-8') if docs_readme.exists() else None

        if current == source:
            continue

        changed += 1
        rel = docs_readme.relative_to(docs_root)
        if check_only:
            errors.append(f'{rel}: docs copy is out of date (run scripts/sync_example_docs.py)')
        else:
            docs_readme.parent.mkdir(parents=True, exist_ok=True)
            docs_readme.write_text(source, encoding='utf-8')
            print(f'Updated {docs_readme.relative_to(REPO)}')

    if check_only and not errors:
        print(f'Example README sync check passed ({len(iter_example_readmes(examples_root))} files, {changed} would change).')
    elif not check_only:
        print(f'Sync complete ({changed} file(s) updated).')

    return errors


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--examples-dir', type=Path, default=REPO / 'examples')
    parser.add_argument('--docs-dir', type=Path, default=REPO / 'docs' / 'examples')
    parser.add_argument('--check', action='store_true', help='Exit 1 if docs copies are stale')
    args = parser.parse_args()

    examples_root = args.examples_dir.resolve()
    docs_root = args.docs_dir.resolve()

    if not examples_root.is_dir():
        print(f'Examples directory does not exist: {examples_root}', file=sys.stderr)
        return 1

    errors = sync_readmes(examples_root, docs_root, args.check)
    if errors:
        print('Example README sync check failed:', file=sys.stderr)
        for error in errors:
            print(f'- {error}', file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
