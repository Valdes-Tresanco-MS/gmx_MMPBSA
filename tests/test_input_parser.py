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

    def test_create_input_comments_unset_optional_values(self):
        parser = copy.deepcopy(input_file)
        output = Path(tempfile.mkdtemp()) / 'mmpbsa.in'
        try:
            parser.print_contents(output, ('general', 'gb'))
            text = output.read_text()

            self.assertRegex(text, r'(?m)^\s*#\s*ndiis_attempts\s*=\s*None\s')
            self.assertIsNone(re.search(r'(?m)^\s+ndiis_attempts\s*=', text))
            parsed = parser.Parse(output)
            self.assertIsNone(parsed['gb']['ndiis_attempts'])
        finally:
            output.unlink(missing_ok=True)
            output.parent.rmdir()


if __name__ == '__main__':
    unittest.main()
