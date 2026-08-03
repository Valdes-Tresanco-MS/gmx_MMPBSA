from contextlib import redirect_stderr
from io import StringIO
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

from scripts import run_example_docs_symlink_spike as spike


class ExampleDocsSymlinkSpikeTest(unittest.TestCase):
    def test_missing_mkdocs_restores_docs_examples_and_preserves_site(self):
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            docs_examples = root / 'docs' / 'examples'
            examples = root / 'examples'
            docs_examples.mkdir(parents=True)
            examples.mkdir()
            (docs_examples / 'README.md').write_text('original docs')
            (docs_examples / 'gmx_MMPBSA_test.md').write_text('test docs')
            (examples / 'gmx_MMPBSA_test.md').write_text('existing example file')
            site = root / 'site_symlink_spike'
            site.mkdir()
            (site / 'sentinel').write_text('keep me')

            with patch.object(spike, 'REPO', root), patch.object(
                spike.subprocess, 'run', side_effect=FileNotFoundError('mkdocs')
            ):
                diagnostics = StringIO()
                with redirect_stderr(diagnostics):
                    status = spike.main()

            self.assertEqual(status, 1)
            self.assertIn('Could not run mkdocs', diagnostics.getvalue())
            self.assertFalse(docs_examples.is_symlink())
            self.assertEqual((docs_examples / 'README.md').read_text(), 'original docs')
            self.assertEqual((examples / 'gmx_MMPBSA_test.md').read_text(), 'existing example file')
            self.assertEqual((site / 'sentinel').read_text(), 'keep me')


if __name__ == '__main__':
    unittest.main()
