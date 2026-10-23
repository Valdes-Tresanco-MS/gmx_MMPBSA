import unittest
from unittest.mock import Mock

from GMXMMPBSA.make_top import CheckMakeTop
from GMXMMPBSA.make_top_amber import CheckAmberTop
from GMXMMPBSA.utils import Residue, reconcile_qm_charges


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


class QMMMTopologyMessageTest(unittest.TestCase):
    @staticmethod
    def _checker(builder, gbrun=True, ifqnt=0, mask=':1'):
        checker = object.__new__(builder)
        checker.INPUT = {
            'decomp': {'decomprun': False},
            'gb': {'gbrun': gbrun, 'ifqnt': ifqnt, 'com_qmmask': mask, 'qm_residues': 'all',
                   'qmcharge_com': 0, 'qmcharge_rec': 0, 'qmcharge_lig': 0},
        }
        checker.explicit_waters = 0
        checker.explicit_waters_mask = ''
        checker.gmx2pdb = Mock()
        checker.amber2pdb = Mock()
        checker.gmxtop2prmtop = Mock(return_value=('COM.prmtop', 'REC.prmtop', 'LIG.prmtop'))
        checker.ambertop2prmtop = checker.gmxtop2prmtop
        checker._resolve_explicit_waters_mask = Mock()
        residues = [Residue(1, 1, 'A', 'R', 1, 'ALA')]
        checker.get_selected_residues = Mock(return_value=(residues, (-1, 2)))
        checker.cleanup_trajs = Mock()
        return checker

    def test_disabled_qmmm_does_not_report_charge_overrides(self):
        for builder in (CheckMakeTop, CheckAmberTop):
            with self.subTest(builder=builder.__name__):
                checker = self._checker(builder)
                with self.assertNoLogs(level='INFO'):
                    checker.buildTopology()

    def test_disabled_gb_does_not_report_qmmm_setup(self):
        for builder in (CheckMakeTop, CheckAmberTop):
            for mask in ('', ':1'):
                with self.subTest(builder=builder.__name__, mask=mask):
                    checker = self._checker(builder, gbrun=False, ifqnt=1, mask=mask)
                    with self.assertNoLogs(level='INFO'):
                        checker.buildTopology()

    def test_active_qmmm_with_user_masks_still_reports_charge_override(self):
        for builder in (CheckMakeTop, CheckAmberTop):
            with self.subTest(builder=builder.__name__):
                checker = self._checker(builder, ifqnt=1)
                with self.assertLogs(level='WARNING') as logs:
                    checker.buildTopology()
                self.assertIn('qmcharge_com', logs.output[0])
                self.assertEqual(len(logs.output), 1)

    def test_active_qmmm_without_masks_still_assigns_residues_and_charges(self):
        for builder in (CheckMakeTop, CheckAmberTop):
            with self.subTest(builder=builder.__name__):
                checker = self._checker(builder, ifqnt=1, mask='')
                with self.assertLogs(level='INFO'):
                    checker.buildTopology()
                self.assertEqual(checker.INPUT['gb']['qm_residues'], '1')
                self.assertEqual(checker.INPUT['gb']['qmcharge_com'], 1)
                self.assertEqual(checker.INPUT['gb']['qmcharge_rec'], -1)
                self.assertEqual(checker.INPUT['gb']['qmcharge_lig'], 2)


if __name__ == '__main__':
    unittest.main()
