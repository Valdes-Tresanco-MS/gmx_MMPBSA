import tempfile
import unittest
from pathlib import Path

try:
    from GMXMMPBSA.createinput import SanderGBInput, SanderPBSAInput
except ModuleNotFoundError as exc:
    if exc.name == 'parmed':
        SanderGBInput = None
        SanderPBSAInput = None
    else:
        raise


class PBInputDefaultTest(unittest.TestCase):
    @unittest.skipIf(SanderPBSAInput is None, 'ParmEd is required to generate PB input files')
    def test_pb_mdin_fallback_uses_inp_one(self):
        pb = SanderPBSAInput({
            'general': {'netcdf': False},
            'pb': {'istrng': 0.0},
            'decomp': {'idecomp': 0, 'dec_verbose': 0},
        })
        pb.make_mdin()

        with tempfile.TemporaryDirectory() as tmpdir:
            output = Path(tmpdir) / 'pb.mdin'
            pb.write_input(output)
            text = output.read_text()

        self.assertIn('inp=1', text.replace(' ', ''))


class QMMMInputTest(unittest.TestCase):
    @unittest.skipIf(SanderGBInput is None, 'ParmEd is required to generate GB input files')
    def test_itrmax_is_written_to_qmmm_namelist(self):
        gb = SanderGBInput({
            'general': {'netcdf': False},
            'gb': {'ifqnt': 1, 'qmmask': "':1'", 'qm_theory': "'PM6'", 'qmcharge': 0, 'itrmax': 5000},
            'decomp': {'idecomp': 0, 'dec_verbose': 0},
        })
        gb.make_mdin()

        with tempfile.TemporaryDirectory() as tmpdir:
            output = Path(tmpdir) / 'qmmm.mdin'
            gb.write_input(output)
            text = output.read_text().replace(' ', '')

        self.assertIn('&qmmm', text)
        self.assertIn('itrmax=5000', text)


if __name__ == '__main__':
    unittest.main()
