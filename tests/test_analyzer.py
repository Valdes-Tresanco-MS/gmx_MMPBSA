import os
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')

from PyQt6.QtWidgets import QApplication

from GMXMMPBSA.analyzer.dialogs import InitDialog


class AnalyzerInfoLoadingTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.qt_app = QApplication.instance() or QApplication([])

    def test_mutant_only_info_file_does_not_create_normal_entry(self):
        with TemporaryDirectory() as tmpdir:
            info_file = Path(tmpdir) / '_GMXMMPBSA_info'
            info_file.write_text(
                "INPUT['general']['exp_ki'] = [15.0]\n"
                "INPUT['general']['sys_name'] = 'H15'\n"
                "INPUT['ala']['mutant_only'] = 1\n"
                "FILES.stability = False\n"
                "mut_str = 'A/15 - HIDxALA'\n"
            )

            dialog = InitDialog(None)
            dialog.get_files_info([info_file])

            system = dialog.f_item.child(0)
            self.assertEqual(system.childCount(), 1)
            self.assertEqual(system.child(0).text(2), 'A/15 - HIDxALA')

            dialog.close()
