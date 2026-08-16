import unittest

from GMXMMPBSA.utils import reconcile_qm_charges


class SharedTopologyMessageTest(unittest.TestCase):
    def test_automatic_qm_charge_assignment_is_shared_and_info_level(self):
        gb_input = {'qmcharge_com': 0, 'qmcharge_rec': 0, 'qmcharge_lig': 0}

        with self.assertLogs(level='INFO') as logs:
            reconcile_qm_charges(gb_input, -1, 2)

        self.assertEqual(gb_input, {
            'qmcharge_com': 1,
            'qmcharge_rec': -1,
            'qmcharge_lig': 2,
        })
        self.assertEqual(sum(message.startswith('INFO:') for message in logs.output), 3)
        self.assertFalse(any(message.startswith('WARNING:') for message in logs.output))

    def test_user_qm_charge_mismatch_remains_visible_once_per_value(self):
        gb_input = {'qmcharge_com': 0, 'qmcharge_rec': 0, 'qmcharge_lig': 2}

        with self.assertLogs(level='INFO') as logs:
            reconcile_qm_charges(gb_input, -1, 2)

        self.assertEqual(gb_input, {
            'qmcharge_com': 0,
            'qmcharge_rec': 0,
            'qmcharge_lig': 2,
        })
        self.assertEqual(sum(message.startswith('INFO:') for message in logs.output), 3)
        self.assertEqual(sum(message.startswith('WARNING:') for message in logs.output), 2)


if __name__ == '__main__':
    unittest.main()
