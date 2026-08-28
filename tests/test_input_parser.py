import copy
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
