import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from scripts.validate_gmx_MMPBSA_test_docs import (
    _command_signature,
    _parse_bundled_test_command,
    _parse_serial_command,
    _test_selector,
)


class GmxMmpbsaTestDocumentationValidatorTest(unittest.TestCase):
    def test_parses_current_bundled_test_section(self):
        with TemporaryDirectory() as tmpdir:
            readme = Path(tmpdir) / 'README.md'
            readme.write_text(
                '### Run the bundled test\n\n'
                '```bash\n'
                'gmx_MMPBSA_test -t 4\n'
                '```\n\n'
                '### Run it manually\n'
            )
            command = _parse_bundled_test_command(readme)

        self.assertEqual(command, 'gmx_MMPBSA_test -t 4')
        self.assertEqual(_test_selector(command), '4')

    def test_parses_multiline_serial_command(self):
        with TemporaryDirectory() as tmpdir:
            readme = Path(tmpdir) / 'README.md'
            readme.write_text(
                '=== "Serial"\n\n'
                '    ```bash\n'
                '    gmx_MMPBSA -O \\\n'
                '      -i mmpbsa.in \\\n'
                '      -cs com.tpr\n'
                '    ```\n\n'
                '=== "With MPI"\n'
            )
            command = _parse_serial_command(readme)

        self.assertEqual(
            _command_signature(command),
            _command_signature(['gmx_MMPBSA', '-cs', 'com.tpr', '-i', 'mmpbsa.in', '-O']),
        )


if __name__ == '__main__':
    unittest.main()
