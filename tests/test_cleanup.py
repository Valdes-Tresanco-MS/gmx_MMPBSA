import os
import tempfile
import unittest
from pathlib import Path

from GMXMMPBSA.utils import remove


class CleanupTest(unittest.TestCase):
    def test_clean_removes_default_csv_outputs_with_summary_outputs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            for filename in (
                'FINAL_RESULTS_MMPBSA.dat', 'FINAL_RESULTS_MMPBSA.csv',
                'FINAL_DECOMP_MMPBSA.dat', 'FINAL_DECOMP_MMPBSA.csv',
            ):
                (root / filename).write_text('results\n')

            old_cwd = os.getcwd()
            os.chdir(tmpdir)
            try:
                remove(-1)
            finally:
                os.chdir(old_cwd)

            for filename in (
                'FINAL_RESULTS_MMPBSA.dat', 'FINAL_RESULTS_MMPBSA.csv',
                'FINAL_DECOMP_MMPBSA.dat', 'FINAL_DECOMP_MMPBSA.csv',
            ):
                self.assertFalse((root / filename).exists())

    def test_minimal_cleanup_preserves_info_metadata(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            (root / '_GMXMMPBSA_info').write_text("INPUT['pb']['inp'] = 1\n")
            (root / '_GMXMMPBSA_pb.mdin').write_text('inp=1\n')
            (root / 'GMXMMPBSA_membrane_parameters.csv').write_text('diagnostic\n')
            (root / 'GMXMMPBSA_membrane_parameters.png').write_bytes(b'diagnostic')
            (root / 'COM.prmtop').write_text('temporary\n')
            (root / 'COMPACT_MMXSA_RESULTS.mmxsa').write_text('results\n')

            old_cwd = os.getcwd()
            os.chdir(tmpdir)
            try:
                remove(0)
            finally:
                os.chdir(old_cwd)

            self.assertTrue((root / '_GMXMMPBSA_info').exists())
            self.assertTrue((root / 'GMXMMPBSA_membrane_parameters.csv').exists())
            self.assertTrue((root / 'GMXMMPBSA_membrane_parameters.png').exists())
            self.assertTrue((root / 'COMPACT_MMXSA_RESULTS.mmxsa').exists())
            self.assertFalse((root / '_GMXMMPBSA_pb.mdin').exists())
            self.assertFalse((root / 'COM.prmtop').exists())


if __name__ == '__main__':
    unittest.main()
