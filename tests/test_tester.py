import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from GMXMMPBSA.tester import run_process


class TestOutputValidation(unittest.TestCase):
    def _run(self, root, command, expected=('result.dat',), skip=False):
        return run_process(
            root,
            'Synthetic test',
            'synthetic',
            command,
            root / 'synthetic.log',
            list(expected),
            skip,
        )

    def test_unchanged_existing_output_is_reported_stale(self):
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            (root / 'result.dat').write_text('old result')

            result = self._run(root, ['/bin/true'])

            self.assertEqual(result, ('synthetic', True, ['result.dat']))

    def test_new_output_is_accepted(self):
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)

            result = self._run(root, ['/bin/sh', '-c', 'printf new > result.dat'])

            self.assertEqual(result, ('synthetic', False, []))

    def test_changed_output_is_accepted(self):
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            (root / 'result.dat').write_text('old result')

            result = self._run(root, ['/bin/sh', '-c', 'printf new > result.dat'])

            self.assertEqual(result, ('synthetic', False, []))

    def test_skip_output_check_preserves_bypass(self):
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            (root / 'result.dat').write_text('old result')

            result = self._run(root, ['/bin/true'], skip=True)

            self.assertEqual(result, ('synthetic', False, []))


if __name__ == '__main__':
    unittest.main()
