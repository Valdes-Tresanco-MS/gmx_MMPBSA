import logging
import tempfile
import unittest
from pathlib import Path
from unittest.mock import mock_open, patch

from GMXMMPBSA.fake_mpi import MPI as FakeMPI
from GMXMMPBSA.logging_utils import (
    enable_file_logging,
    format_command_line,
    setup_logging as _setup_logging,
)
from GMXMMPBSA import utils


class LoggingOwnershipTest(unittest.TestCase):
    def setUp(self):
        self.root = logging.getLogger()
        self.handlers = self.root.handlers[:]
        self.level = self.root.level
        self.root.handlers.clear()

    def tearDown(self):
        for handler in self.root.handlers[:]:
            self.root.removeHandler(handler)
            handler.close()
        self.root.handlers[:] = self.handlers
        self.root.setLevel(self.level)

    def _file_handlers(self):
        return [handler for handler in self.root.handlers if isinstance(handler, logging.FileHandler)]

    def test_master_rank_owns_the_calculation_log(self):
        with tempfile.TemporaryDirectory() as directory:
            log_file = Path(directory) / 'gmx_MMPBSA.log'

            _setup_logging(log_file, master=True, rank=0, force=True)
            logging.info('master record')

            self.assertEqual(len(self._file_handlers()), 1)
            self.assertTrue(log_file.exists())
            self.assertIn('master record', log_file.read_text())

    def test_terminal_only_setup_does_not_create_calculation_log(self):
        with tempfile.TemporaryDirectory() as directory:
            log_file = Path(directory) / 'gmx_MMPBSA.log'

            _setup_logging(log_file, master=True, rank=0, force=True, file_enabled=False)
            logging.info('terminal-only record')

            self.assertFalse(log_file.exists())
            self.assertEqual(self._file_handlers(), [])

    def test_file_records_use_compact_format(self):
        with tempfile.TemporaryDirectory() as directory:
            log_file = Path(directory) / 'gmx_MMPBSA.log'

            _setup_logging(log_file, master=True, rank=0, force=True, file_enabled=False)
            enable_file_logging(log_file, rank=0)
            logging.info('metadata record')

            text = log_file.read_text()
            self.assertEqual(text, '[INFO   ] metadata record\n')

    def test_warning_counter_does_not_count_final_summary(self):
        content = (
            '[WARNING] warning\n'
            '[ERROR  ] error\n'
            '[INFO   ] '
            'Finalizing gmx_MMPBSA: [ERROR  ] = 0; [WARNING] = 1\n'
        )
        with patch('builtins.open', mock_open(read_data=content)):
            self.assertEqual(utils.get_warnings(), {'warning': 1, 'error': 1})

    def test_command_line_quoting_and_reconstructed_mpi_label(self):
        command = format_command_line(
            ['-i', 'path with spaces/mmpbsa.in', '-cg', ':A & !@', 'mpi'],
            mpi_size=2,
            mpi_requested=True,
        )

        self.assertIn("'path with spaces/mmpbsa.in'", command)
        self.assertIn("':A & !@'", command)
        self.assertNotIn(' mpi ', command)
        self.assertTrue(command.endswith(' (reconstructed)'))

    def test_non_master_rank_does_not_open_or_modify_master_log(self):
        with tempfile.TemporaryDirectory() as directory:
            log_file = Path(directory) / 'gmx_MMPBSA.log'
            log_file.write_text('previous master output\n')

            _setup_logging(log_file, master=False, rank=1, force=True)
            logging.info('rank one record')

            self.assertEqual(self._file_handlers(), [])
            self.assertEqual(log_file.read_text(), 'previous master output\n')
            self.assertFalse(Path(directory, 'gmx_MMPBSA.rank-1.log').exists())

            logging.error('rank one failure')
            rank_log = Path(directory, 'gmx_MMPBSA.rank-1.log')
            self.assertTrue(rank_log.exists())
            self.assertIn('rank one failure', rank_log.read_text())
            self.assertEqual(log_file.read_text(), 'previous master output\n')

    def test_fake_mpi_serial_rank_owns_the_log(self):
        with tempfile.TemporaryDirectory() as directory:
            log_file = Path(directory) / 'gmx_MMPBSA.log'
            communicator = FakeMPI.COMM_WORLD

            _setup_logging(
                log_file,
                master=communicator.Get_rank() == 0,
                rank=communicator.Get_rank(),
                force=True,
            )

            self.assertEqual(communicator.Get_size(), 1)
            self.assertEqual(len(self._file_handlers()), 1)
            self.assertTrue(log_file.exists())

    def test_reinitialization_does_not_duplicate_handlers(self):
        with tempfile.TemporaryDirectory() as directory:
            log_file = Path(directory) / 'gmx_MMPBSA.log'

            _setup_logging(log_file, master=True, rank=0, force=True)
            _setup_logging(log_file, master=True, rank=0, force=True)

            self.assertEqual(len(self._file_handlers()), 1)
            self.assertEqual(
                len([handler for handler in self.root.handlers if isinstance(handler, logging.StreamHandler)]),
                2,
            )


if __name__ == '__main__':
    unittest.main()
