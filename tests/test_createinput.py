import tempfile
import unittest
from pathlib import Path

try:
    from GMXMMPBSA.createinput import SanderAPBSInput, SanderGBInput, SanderMMInput, SanderPBSAInput
except ModuleNotFoundError as exc:
    if exc.name == 'parmed':
        SanderAPBSInput = SanderGBInput = SanderMMInput = None
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


class DecompositionMdinTest(unittest.TestCase):
    @unittest.skipIf(SanderMMInput is None, 'ParmEd is required to generate MM input files')
    def _write_mdin(self, decomprun, idecomp, dec_verbose):
        mm = SanderMMInput({
            'general': {'netcdf': False},
            'decomp': {
                'decomprun': decomprun,
                'idecomp': idecomp,
                'dec_verbose': dec_verbose,
            },
        })
        mm.make_mdin()

        with tempfile.TemporaryDirectory() as tmpdir:
            output = Path(tmpdir) / 'mm.mdin'
            mm.write_input(output)
            return output.read_text().replace(' ', '')

    def test_normal_mdin_disables_decomposition(self):
        text = self._write_mdin(decomprun=False, idecomp=2, dec_verbose=1)

        self.assertNotIn('idecomp=2', text)
        self.assertIn('dec_verbose=0', text)

    def test_decomposition_mdin_preserves_nonzero_settings(self):
        text = self._write_mdin(decomprun=True, idecomp=2, dec_verbose=1)

        self.assertIn('idecomp=2', text)
        self.assertIn('dec_verbose=1', text)


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

    def test_ndiis_attempts_is_optional_and_written_when_set(self):
        gb = SanderGBInput({
            'general': {'netcdf': False},
            'gb': {
                'ifqnt': 1, 'qmmask': "':1'", 'qm_theory': "'AM1'", 'qmcharge': -2,
                'ndiis_attempts': 700,
            },
            'decomp': {'idecomp': 0, 'dec_verbose': 0},
        })
        gb.make_mdin()

        with tempfile.TemporaryDirectory() as tmpdir:
            output = Path(tmpdir) / 'qmmm.mdin'
            gb.write_input(output)
            text = output.read_text().replace(' ', '')

        self.assertIn('ndiis_attempts=700', text)

    def test_ndiis_attempts_is_not_written_when_omitted(self):
        gb = SanderGBInput({
            'general': {'netcdf': False},
            'gb': {'ifqnt': 1, 'qmmask': "':1'", 'qm_theory': "'AM1'", 'qmcharge': -2},
            'decomp': {'idecomp': 0, 'dec_verbose': 0},
        })
        gb.make_mdin()

        with tempfile.TemporaryDirectory() as tmpdir:
            output = Path(tmpdir) / 'qmmm.mdin'
            gb.write_input(output)
            text = output.read_text().replace(' ', '')

        self.assertNotIn('ndiis_attempts=', text)


class APBSInputTest(unittest.TestCase):
    @unittest.skipIf(SanderAPBSInput is None, 'ParmEd is required to generate APBS input files')
    def test_apbs_mdin_applies_custom_pb_settings(self):
        pb = SanderAPBSInput({
            'general': {'netcdf': False},
            'pb': {
                'scale': 2.0,
                'istrng': 0.15,
                'cavity_surften': 0.0378,
                # APBS/ParmEd defaults suppress matching values; use non-defaults.
                'inp': 1,           # APBS cntrl default is 2
                'radiopt': 1,       # APBS default is 0
                'prbrad': 1.6,      # APBS default is 1.4
                'indi': 3.5,        # APBS pdie default is 2.0
                'exdi': 70.0,       # APBS sdie default is 78.4
                'pb_maxcyc': 7,
            },
            'decomp': {'decomprun': False, 'idecomp': 2, 'dec_verbose': 1},
        })
        pb.make_mdin()

        with tempfile.TemporaryDirectory() as tmpdir:
            output = Path(tmpdir) / 'apbs.mdin'
            pb.write_input(output)
            text = output.read_text().replace(' ', '')

        self.assertEqual(pb.namelist, 'pb')
        self.assertIn('&apbs', text)
        self.assertNotIn('&pb', text)
        self.assertIn("grid='2.0,2.0,2.0'", text)
        self.assertIn("ionc='0.15,0.15'", text)
        self.assertIn('gamma=0.1581552', text)
        self.assertIn('pdie=3.5', text)
        self.assertIn('sdie=70.0', text)
        self.assertIn('radiopt=1', text)
        self.assertIn('srad=1.6', text)
        self.assertIn('inp=1', text)
        self.assertIn('maxcyc=7', text)
        # decomprun is off, so nonzero template idecomp must not leak
        self.assertNotIn('idecomp=2', text)
        self.assertEqual(pb.mdin.cntrl_nml['idecomp'], 0)


if __name__ == '__main__':
    unittest.main()
