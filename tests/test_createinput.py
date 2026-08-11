import tempfile
import unittest
from pathlib import Path

try:
    from GMXMMPBSA.createinput import SanderPBSAInput
except ModuleNotFoundError as exc:
    if exc.name == 'parmed':
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


if __name__ == '__main__':
    unittest.main()
