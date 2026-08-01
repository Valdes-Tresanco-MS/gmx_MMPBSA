#!/usr/bin/env python3
"""Run the Phase 2 symlink spike and print a link audit summary."""

from __future__ import annotations

import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]

AUDIT_PAGES = [
    (
        'examples/AMBER/index.html',
        '../../docs/amber_MMPBSA.md',
        'broken under symlink layout',
    ),
    (
        'examples/Protein_ligand/ST/index.html',
        '../../../gmx_MMPBSA_command-line/',
        'MkDocs-oriented footnote depth',
    ),
    (
        'examples/Entropy_calculations/nmode/index.html',
        '../../../gmx_MMPBSA_command-line/',
        'MkDocs-oriented footnote depth',
    ),
    (
        'examples/psf_dcd/protein_protein/index.html',
        '../../../gmx_MMPBSA_command-line/',
        'MkDocs-oriented footnote depth',
    ),
]


def extract_hrefs(html: str) -> set[str]:
    return {match.strip('"\'') for match in re.findall(r'href=([^\s>]+)', html)}


def restore_docs_examples(backup: Path, docs_examples: Path) -> None:
    if docs_examples.is_symlink():
        docs_examples.unlink()
    elif docs_examples.exists():
        shutil.rmtree(docs_examples)
    shutil.copytree(backup, docs_examples)


def main() -> int:
    docs = REPO / 'docs'
    examples = REPO / 'examples'
    docs_examples = docs / 'examples'
    site_dir = REPO / 'site_symlink_spike'
    gmx_test_tmp: Path | None = None
    build_code = 1
    link_issues: list[str] = []

    if not docs_examples.is_dir() or docs_examples.is_symlink():
        print('docs/examples must be a regular directory before running the spike', file=sys.stderr)
        return 1

    with tempfile.TemporaryDirectory(prefix='gmx_mmpbsa_examples_spike_') as tmp:
        backup = Path(tmp) / 'examples_backup'
        gmx_test_backup = backup / 'gmx_MMPBSA_test.md'

        shutil.copytree(docs_examples, backup)
        shutil.rmtree(docs_examples)
        docs_examples.symlink_to('../examples', target_is_directory=True)

        if gmx_test_backup.exists():
            gmx_test_tmp = examples / 'gmx_MMPBSA_test.md'
            shutil.copy2(gmx_test_backup, gmx_test_tmp)

        if site_dir.exists():
            shutil.rmtree(site_dir)

        build = subprocess.run(
            ['mkdocs', 'build', '--strict', '-d', str(site_dir)],
            cwd=REPO,
            capture_output=True,
            text=True,
        )
        build_code = build.returncode
        build_output = '\n'.join(part for part in (build.stdout, build.stderr) if part)

        print('mkdocs build --strict exit code:', build_code)
        broken_link_warnings = [
            line.strip()
            for line in build_output.splitlines()
            if 'WARNING -  Doc file' in line and 'is not found among documentation files' in line
        ]
        if broken_link_warnings:
            print('\nBroken doc links reported by mkdocs:')
            for line in broken_link_warnings:
                print(f'  {line}')

        if build_code == 0 or site_dir.exists():
            for page, expected_href, note in AUDIT_PAGES:
                html_path = site_dir / page
                if not html_path.exists():
                    link_issues.append(f'{page}: HTML page missing')
                    continue
                hrefs = extract_hrefs(html_path.read_text(encoding='utf-8', errors='replace'))
                if expected_href in hrefs:
                    status = 'FOUND' if 'broken' in note else 'OK'
                    print(f'{status}: {page} contains {expected_href!r} ({note})')
                    if 'broken' in note:
                        link_issues.append(f'{page}: footnote href {expected_href!r} is broken')
                else:
                    print(f'MISSING: {page} does not contain {expected_href!r} ({note})')
                    if 'broken' not in note:
                        link_issues.append(f'{page}: expected footnote href {expected_href!r}')

        restore_docs_examples(backup, docs_examples)

    if gmx_test_tmp and gmx_test_tmp.exists():
        gmx_test_tmp.unlink()
    if site_dir.exists():
        shutil.rmtree(site_dir)

    print('\nSpike summary:')
    if build_code != 0:
        print('- mkdocs build --strict FAILED with symlink layout')
    else:
        print('- mkdocs build --strict succeeded with symlink layout')

    if link_issues:
        print('- link audit issues:')
        for item in link_issues:
            print(f'  - {item}')
        print('- recommendation: keep Phase 1 sync; do not merge symlink PR')
        return 1

    print('- sampled footnote links look correct in built HTML')
    print('- recommendation: manual HTML review still advised before any structural PR')
    return 0 if build_code == 0 else 1


if __name__ == '__main__':
    sys.exit(main())
