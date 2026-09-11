import copy
from io import StringIO
import re
import tempfile
import unittest
from pathlib import Path

from GMXMMPBSA.exceptions import InputError
from GMXMMPBSA.input_parser import input_file


class InputParserTest(unittest.TestCase):
    def _parse(self, contents):
        with tempfile.TemporaryDirectory() as tmpdir:
            infile = Path(tmpdir) / 'mmpbsa.in'
            infile.write_text(contents)
            parser = copy.deepcopy(input_file)
            for namelist in parser.namelists.values():
                namelist.open = False
            return parser.Parse(infile.as_posix())

    def test_generated_none_value_is_accepted_for_optional_ndiis_attempts(self):
        parsed = self._parse('&gb\n  ndiis_attempts = None\n/\n')

        self.assertIsNone(parsed['gb']['ndiis_attempts'])

    def test_invalid_numeric_value_reports_input_variable(self):
        with self.assertRaisesRegex(InputError, r"Invalid value 'not-an-int' for itrmax"):
            self._parse('&gb\n  itrmax = not-an-int\n/\n')

    def test_indented_inline_comments_are_accepted(self):
        parsed = self._parse('&general\n  startframe = 7 # first production frame\n/\n')

        self.assertEqual(parsed['general']['startframe'], 7)

    def test_implicit_solvent_defaults_match_170_contract(self):
        parsed = self._parse('&general\n/\n&gb\n/\n&pb\n/\n')

        self.assertEqual(parsed['general']['PBRadii'], 4)
        self.assertEqual(parsed['gb']['igb'], 8)
        self.assertEqual(parsed['pb']['exdi'], 78.5)

    def test_pbradii_accepts_named_sets_case_insensitively(self):
        expected = {
            'bondi': 1,
            'mbondi': 2,
            'mbondi2': 3,
            'mbondi3': 4,
            'mbondi_pb2': 5,
            'mbondi_pb3': 6,
            'charmm_radii': 7,
        }

        for name, number in expected.items():
            parsed = self._parse(f'&general\n  PBRadii = "{name.upper()}"\n/\n')
            self.assertEqual(parsed['general']['PBRadii'], number)

    def test_pbradii_invalid_named_set_reports_allowed_names(self):
        with self.assertRaisesRegex(InputError, r'expected int or one of .*mbondi3'):
            self._parse('&general\n  PBRadii = not_a_radius_set\n/\n')

    def test_unterminated_namelist_is_rejected(self):
        with self.assertRaisesRegex(InputError, r'Unterminated namelist general'):
            self._parse('&general\n  startframe = 7\n')

    def test_membrane_parameters_accept_automatic_and_atom_names(self):
        parsed = self._parse('&pb\n  memopt=1, mthick=automatic, mctrdz=automatic, membrane_atoms="P;N"\n/\n')

        self.assertEqual(parsed['pb']['mthick'], 'automatic')
        self.assertEqual(parsed['pb']['mctrdz'], 'automatic')
        self.assertEqual(parsed['pb']['membrane_atoms'], 'P;N')

    def test_membrane_parameters_default_to_automatic_phosphorus_detection(self):
        parsed = self._parse('&pb\n/\n')

        self.assertEqual(parsed['pb']['mthick'], 'automatic')
        self.assertEqual(parsed['pb']['mctrdz'], 'automatic')
        self.assertEqual(parsed['pb']['membrane_atoms'], 'P')

    def test_create_input_comments_unset_optional_values(self):
        parser = copy.deepcopy(input_file)
        for namelist in parser.namelists.values():
            namelist.open = False
        output = Path(tempfile.mkdtemp()) / 'mmpbsa.in'
        try:
            parser.print_contents(output, ('general', 'gb'))
            text = output.read_text()

            self.assertRegex(text, r'(?m)^\s*#\s*ndiis_attempts\s*=\s*None\s')
            self.assertIsNone(re.search(r'(?m)^\s+ndiis_attempts\s*=', text))
            parsed = parser.Parse(output)
            self.assertIsNone(parsed['gb']['ndiis_attempts'])
            self.assertFalse(any(n.open for n in parser.namelists.values()))
        finally:
            output.unlink(missing_ok=True)
            output.parent.rmdir()

    def test_create_input_pb_membrane_template_uses_membrane_defaults(self):
        parser = copy.deepcopy(input_file)
        output = StringIO()

        parser.print_contents(output, ('general', 'pb_mem'))
        text = output.getvalue()

        expected_values = {
            'memopt': '1',
            'emem': '7.0',
            'indi': '4.0',
            'mctrdz': '"automatic"',
            'mthick': '"automatic"',
            'poretype': '1',
            'radiopt': '0',
            'istrng': '0.15',
            'fillratio': '1.25',
            'inp': '2',
            'sasopt': '0',
            'solvopt': '2',
            'ipb': '1',
            'bcopt': '10',
            'nfocus': '1',
            'linit': '1000',
            'eneopt': '1',
            'cutfd': '7.0',
            'cutnb': '99.0',
            'maxarcdot': '15000',
            'npbverb': '1',
            'membrane_atoms': '"P"',
        }
        for name, value in expected_values.items():
            self.assertRegex(text, rf'(?m)^\s+{name}\s*=\s*{re.escape(value)}\s+#')

        self.assertIn('&pb\n', text)
        self.assertNotIn('&pb_mem', text)
        self.assertEqual(input_file.namelists['pb'].variables['memopt'].value, 0)


    def test_failed_parse_releases_namelist_open_state(self):
        parser = copy.deepcopy(input_file)
        for namelist in parser.namelists.values():
            namelist.open = False

        with tempfile.TemporaryDirectory() as tmpdir:
            infile = Path(tmpdir) / 'mmpbsa.in'
            infile.write_text('&general\norphan_value,\n/\n')

            with self.assertRaises(InputError):
                parser.Parse(infile.as_posix())

        self.assertFalse(any(n.open for n in parser.namelists.values()))


if __name__ == '__main__':
    unittest.main()
