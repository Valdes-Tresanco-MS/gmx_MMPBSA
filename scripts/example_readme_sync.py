"""Shared helpers for syncing example READMEs to docs/examples/."""

from __future__ import annotations

from pathlib import Path

EXCLUDE_DIR_NAMES = {
    'API (copy)',
    '__pycache__',
    '.gmx_mmpbsa_temp',
}

EXCLUDE_FILE_SUFFIXES = (
    '.log',
    '.kate-swp',
)

SKIP_DOCS_RELATIVE = {
    Path('gmx_MMPBSA_test.md'),
}


def should_skip_examples_dir(path: Path) -> bool:
    return any(part in EXCLUDE_DIR_NAMES for part in path.parts)


def iter_example_readmes(examples_root: Path) -> list[Path]:
    readmes: list[Path] = []
    for readme in sorted(examples_root.rglob('README.md')):
        rel = readme.relative_to(examples_root)
        if should_skip_examples_dir(rel):
            continue
        readmes.append(readme)
    return readmes


def transform_readme_content(content: str) -> str:
    # Canonical READMEs use repository-root `docs/` targets so their links work
    # when viewed directly on GitHub.  The mirrored `docs/examples/` tree is
    # one directory deeper, so remove the docs prefix for each supported depth.
    transformed = content.replace('../../../docs/examples/gmx_MMPBSA_test.md', '../../gmx_MMPBSA_test.md')
    transformed = transformed.replace('../../docs/examples/gmx_MMPBSA_test.md', '../gmx_MMPBSA_test.md')
    transformed = transformed.replace('../../../docs/', '../../../')
    transformed = transformed.replace('../../docs/', '../../')
    transformed = transformed.replace('../docs/', '../')
    transformed = transformed.replace('\u2013', '-').replace('\u2014', '-')
    return transformed


def docs_readme_path(examples_root: Path, docs_root: Path, examples_readme: Path) -> Path:
    rel = examples_readme.relative_to(examples_root)
    return docs_root / rel


def is_skipped_docs_target(docs_readme: Path, docs_root: Path) -> bool:
    rel = docs_readme.relative_to(docs_root)
    return rel in SKIP_DOCS_RELATIVE
