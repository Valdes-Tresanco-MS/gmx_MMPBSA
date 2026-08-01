# Example docs symlink spike (Phase 2)

Date: 2026-07-31

## Procedure

On a local throwaway layout (restored automatically by `scripts/run_example_docs_symlink_spike.py`):

1. Back up `docs/examples/` to a temp directory outside `docs/`.
2. Replace `docs/examples/` with `ln -s ../examples docs/examples`.
3. Copy `gmx_MMPBSA_test.md` into `examples/` for the spike only (MkDocs nav entry).
4. Run `mkdocs build --strict -d site_symlink_spike`.
5. Inspect built HTML footnote hrefs for four nested READMEs.
6. Restore the synced `docs/examples/` tree from backup.

Re-run:

```bash
python scripts/run_example_docs_symlink_spike.py
```

## mkdocs build --strict

**Result: FAILED** (exit code 1)

MkDocs reported a broken documentation link:

```
WARNING - Doc file 'examples/AMBER/README.md' contains a link '../../docs/amber_MMPBSA.md',
but the target 'docs/amber_MMPBSA.md' is not found among documentation files.
```

Under the symlink layout, MkDocs resolves relative links from the logical path under `docs/`. From `docs/examples/AMBER/README.md`, `../../docs/amber_MMPBSA.md` becomes `docs/docs/amber_MMPBSA.md`, which does not exist.

This confirms the external review finding: a single `../../docs/...` footnote convention cannot serve both GitHub (`examples/`) and MkDocs (`docs/examples/`) when those trees are unified by symlink.

## Footnote HTML audit (4 nested READMEs)

| Page | Source footnote | Built HTML href | Verdict |
|------|-----------------|-----------------|---------|
| `examples/AMBER/README.md` | `[1]: ../../docs/amber_MMPBSA.md` | `../../docs/amber_MMPBSA.md` | **Broken** |
| `examples/Protein_ligand/ST/README.md` | `[1]: ../../../gmx_MMPBSA_command-line.md#...` | `../../../gmx_MMPBSA_command-line/` | OK |
| `examples/Entropy_calculations/nmode/README.md` | `[1]: ../../../gmx_MMPBSA_command-line.md#...` | `../../../gmx_MMPBSA_command-line/` | OK |
| `examples/psf_dcd/protein_protein/README.md` | (analyzer footnote) | `../../../gmx_MMPBSA_command-line/` nav links OK | OK |

Most example READMEs already use MkDocs-relative footnotes (`../../../...` from nested paths). Only **AMBER** uses a GitHub-oriented `../../docs/...` link in the canonical `examples/` tree.

## Recommendation

**Do not merge a symlink-based structural dedup (PR2).** Stay on **Phase 1**:

- Canonical edits in `examples/**/README.md`
- `python scripts/sync_example_docs.py` rewrites `../../docs/` → `../../` for the published `docs/examples/` copy
- CI `--check` + parity validator prevent drift

If structural dedup is revisited later, safer options are:

1. MkDocs build-time link rewrite on a symlinked tree
2. `pymdownx.snippets` thin wrappers in `docs/examples/`
3. Absolute site URLs (worse GitHub UX)

Relocating `docs/examples/gmx_MMPBSA_test.md` into `examples/` remains optional and blocked until a winning structural strategy exists.
