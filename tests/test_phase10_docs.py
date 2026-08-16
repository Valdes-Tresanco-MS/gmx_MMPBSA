import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class LoggingDocumentationTest(unittest.TestCase):
    def test_logging_guide_documents_supported_modes_and_cluster_monitoring(self):
        guide = (ROOT / 'docs/logging.md').read_text()
        for text in ('gmx_MMPBSA.log', 'tail -f gmx_MMPBSA.log', '`auto`', '`rich`', '`classic`', '`plain`', '`none`'):
            with self.subTest(text=text):
                self.assertIn(text, guide)

    def test_logging_guide_is_in_mkdocs_navigation(self):
        nav = (ROOT / 'mkdocs.yml').read_text()
        self.assertIn('Logging and progress: logging.md', nav)

    def test_changelog_mentions_record_based_logging_changes(self):
        changelog = (ROOT / 'docs/changelog.md').read_text()
        self.assertIn('record-based warning/error totals', changelog)
        self.assertIn('not a stable machine-readable interface', changelog)


if __name__ == '__main__':
    unittest.main()
