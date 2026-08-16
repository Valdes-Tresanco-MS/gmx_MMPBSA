import io
import logging
import unittest
from types import SimpleNamespace
from unittest.mock import patch

from GMXMMPBSA.exceptions import GMXMMPBSA_ERROR, InputError, MMPBSA_Error


class Phase9ErrorReportingTest(unittest.TestCase):
    def test_fatal_helper_logs_one_concise_record_and_exception(self):
        with self.assertLogs(level='ERROR') as logs:
            with self.assertRaises(InputError) as exc:
                GMXMMPBSA_ERROR('bad input value', InputError)

        self.assertEqual(len(logs.output), 1)
        self.assertIn('InputError: bad input value', logs.output[0])
        self.assertNotIn('Check the gmx_MMPBSA.log file', logs.output[0])
        self.assertEqual(str(exc.exception), 'bad input value')
        self.assertTrue(exc.exception._gmxmmpbsa_logged)

    def test_domain_error_is_not_logged_twice_at_application_boundary(self):
        from GMXMMPBSA import app

        error = InputError('invalid input')
        with self.assertLogs(level='ERROR') as logs:
            app._log_uncaught_exception(error)
            app._log_uncaught_exception(error)

        self.assertEqual(len(logs.output), 1)
        self.assertIn('InputError: invalid input', logs.output[0])

    def test_unexpected_error_includes_traceback_context(self):
        from GMXMMPBSA import app

        try:
            raise RuntimeError('unexpected failure')
        except RuntimeError as error:
            with self.assertLogs(level='ERROR') as logs:
                app._log_uncaught_exception(error)

        text = '\n'.join(logs.output)
        self.assertIn('Unexpected internal error: unexpected failure', text)
        self.assertIn('Traceback (most recent call last)', text)

    def test_bundle_path_is_recorded_without_creating_an_error_record(self):
        from GMXMMPBSA import app

        fake_app = SimpleNamespace(
            master=True,
            FILES=SimpleNamespace(no_error_bundle=False),
        )
        with patch.object(app, 'create_error_bundle', return_value='/tmp/error-bundle.zip'):
            with patch('sys.stderr', new_callable=io.StringIO):
                with self.assertLogs(level='INFO') as logs:
                    app._maybe_create_error_bundle(fake_app, MMPBSA_Error('failure'))

        self.assertEqual(sum('Diagnostic error bundle created' in message for message in logs.output), 1)
        self.assertIn('/tmp/error-bundle.zip', '\n'.join(logs.output))


if __name__ == '__main__':
    unittest.main()
