import unittest

from GMXMMPBSA.utils import create_input_args


class CreateInputArgsTest(unittest.TestCase):
    def test_gbnsr6_with_decomp_is_kept(self):
        self.assertEqual(
            create_input_args(['gbnsr6', 'decomp']),
            ['general', 'gbnsr6', 'decomp'],
        )

    def test_decomp_without_compatible_method_is_ignored(self):
        with self.assertLogs(level='WARNING') as logs:
            result = create_input_args(['nmode', 'decomp'])
        self.assertEqual(result, ['general', 'nmode'])
        self.assertTrue(any('gbnsr6' in message for message in logs.output))


if __name__ == '__main__':
    unittest.main()
