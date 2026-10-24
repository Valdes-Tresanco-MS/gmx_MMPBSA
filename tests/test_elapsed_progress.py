import io
import logging
import os
import threading
import unittest
from contextlib import contextmanager
from types import SimpleNamespace
from unittest.mock import patch

from GMXMMPBSA.logging_utils import elapsed_progress, setup_logging
from GMXMMPBSA.make_top import CheckMakeTop
from GMXMMPBSA.make_top_amber import CheckAmberTop


class ObservedStream(io.StringIO):
    def __init__(self, interactive):
        super().__init__()
        self.interactive = interactive
        self.updated = threading.Event()
        self.flushed = threading.Event()

    def isatty(self):
        return self.interactive

    def write(self, value):
        result = super().write(value)
        if '[|]' in value:
            self.updated.set()
        return result

    def flush(self):
        super().flush()
        self.flushed.set()


class ElapsedPreparationProgressTest(unittest.TestCase):
    def test_terminal_and_redirected_streams_receive_live_updates_before_exit(self):
        for interactive in (True, False):
            with self.subTest(interactive=interactive):
                stream = ObservedStream(interactive)
                with self.assertLogs(level='INFO') as records:
                    with elapsed_progress('Preparing test structure', stream=stream, interval=.01) as timer:
                        self.assertTrue(stream.updated.wait(2), 'no live preparation update')
                        self.assertTrue(stream.flushed.wait(2))
                self.assertFalse(timer._thread.is_alive())
                output = stream.getvalue()
                self.assertIn('This can take a while; please wait... [|]', output)
                self.assertNotIn('\x1b', output)
                if interactive:
                    self.assertTrue(output.startswith('\r'))
                    self.assertTrue(output.endswith('\r'))
                else:
                    self.assertNotIn('\r', output)
                    self.assertTrue(output.endswith('\n'))
                self.assertEqual(len(records.output), 1)
                self.assertIn('Preparing test structure completed in', records.output[0])

    def test_mpi_forwarding_redraws_one_line_and_clears_before_completion(self):
        class LineBufferedPipe(ObservedStream):
            def __init__(self):
                super().__init__(False)
                self.pending = ''
                self.forwarded = []
                self.second_update = threading.Event()

            def write(self, value):
                self.pending += value
                while '\n' in self.pending:
                    line, self.pending = self.pending.split('\n', 1)
                    self.forwarded.append(line)
                    super().write(line + '\n')
                    if '[/]' in line:
                        self.second_update.set()
                return len(value)

        for rank_variable in ('OMPI_COMM_WORLD_RANK', 'PMI_RANK', 'PMIX_RANK', 'MPI_LOCALRANKID'):
            with self.subTest(launcher=rank_variable), patch.dict(os.environ, {rank_variable: '0'}, clear=True):
                stream = LineBufferedPipe()
                console = logging.StreamHandler(stream)
                with self.assertLogs(level='INFO') as records:
                    logger = logging.getLogger()
                    logger.addHandler(console)
                    try:
                        with elapsed_progress('MPI preparation', stream=stream, interval=.01) as timer:
                            self.assertTrue(stream.second_update.wait(2), 'MPI launcher buffered live redraws')
                            self.assertTrue(timer._thread.is_alive())
                            self.assertIn('\x1b[1A\r', stream.getvalue())
                    finally:
                        logger.removeHandler(console)
                self.assertFalse(timer._thread.is_alive())
                self.assertEqual(stream.pending, '')
                self.assertTrue(stream.forwarded[-1].startswith('\x1b[1A\r'))
                self.assertIn('MPI preparation completed in', stream.forwarded[-1])
                self.assertEqual(len(records.output), 1)

    def test_mpi_fast_step_does_not_move_cursor(self):
        with patch.dict(os.environ, {'OMPI_COMM_WORLD_RANK': '0'}, clear=True):
            stream = ObservedStream(False)
            with self.assertLogs(level='INFO'):
                with elapsed_progress('Fast MPI step', stream=stream, interval=60):
                    pass
            self.assertEqual(stream.getvalue(), '')

    def test_exception_cleans_up_and_logs_stopped_duration_without_swallowing_error(self):
        for interactive in (True, False):
            with self.subTest(interactive=interactive):
                stream = ObservedStream(interactive)
                with self.assertLogs(level='INFO') as records:
                    with self.assertRaisesRegex(ValueError, 'preparation failed'):
                        with elapsed_progress('Loading topology', stream=stream, interval=.01) as timer:
                            self.assertTrue(stream.updated.wait(2))
                            raise ValueError('preparation failed')
                self.assertFalse(timer._thread.is_alive())
                self.assertIn('Loading topology stopped after', records.output[0])
                self.assertNotIn('\x1b', stream.getvalue())

    def test_fast_steps_do_not_write_empty_redraws(self):
        for interactive in (True, False):
            with self.subTest(interactive=interactive), self.assertLogs(level='INFO'):
                stream = ObservedStream(interactive)
                with elapsed_progress('Fast operation', stream=stream, interval=60):
                    pass
                self.assertEqual(stream.getvalue(), '')

    def test_custom_message_and_stream_without_isatty(self):
        class PlainStream:
            def __init__(self):
                self.output = io.StringIO()
                self.updated = threading.Event()

            def write(self, text):
                self.output.write(text)
                self.updated.set()

            def flush(self):
                pass

        stream = PlainStream()
        with self.assertLogs(level='INFO'):
            with elapsed_progress('Preparation', message='Working', stream=stream, interval=.01):
                self.assertTrue(stream.updated.wait(2))
        self.assertIn('Working [|]', stream.output.getvalue())
        self.assertNotIn('\r', stream.output.getvalue())

    def test_animation_is_not_written_to_the_calculation_log(self):
        import tempfile
        from pathlib import Path

        logger = logging.getLogger()
        old_handlers, old_level = logger.handlers[:], logger.level
        try:
            with tempfile.TemporaryDirectory() as directory:
                path = Path(directory) / 'calculation.log'
                setup_logging(path, force=True)
                stream = ObservedStream(False)
                with elapsed_progress('Structure preparation', stream=stream, interval=.01):
                    self.assertTrue(stream.updated.wait(2))
                self.assertIn('[|]', stream.getvalue())
                text = path.read_text()
                self.assertIn('Structure preparation completed in', text)
                self.assertNotIn('[|]', text)
                self.assertNotIn('\x1b', text)
        finally:
            for handler in logger.handlers:
                if handler not in old_handlers:
                    handler.close()
            logger.handlers[:] = old_handlers
            logger.setLevel(old_level)


class SharedPreparationFeedbackTest(unittest.TestCase):
    def test_both_engines_time_the_same_conversion_for_each_component(self):
        @contextmanager
        def record(label):
            labels.append(label)
            yield

        class ParmFactory:
            @staticmethod
            def from_structure(structure):
                return structure

        for builder in (CheckMakeTop, CheckAmberTop):
            for component in ('complex', 'receptor', 'ligand'):
                with self.subTest(builder=builder.__name__, component=component):
                    labels = []
                    structure = object()
                    with patch('GMXMMPBSA.make_top.elapsed_progress', record):
                        actual = builder._convert_amber_topology(structure, ParmFactory, component)
                    self.assertIs(actual, structure)
                    self.assertEqual(labels, [f'Converting selected {component} topology to AMBER ParmFactory'])

    def test_both_engines_report_validation_for_all_three_components(self):
        @contextmanager
        def record(label):
            labels.append(label)
            yield

        for builder in (CheckMakeTop, CheckAmberTop):
            with self.subTest(builder=builder.__name__):
                labels = []
                checker = object.__new__(builder)
                checker.FILES = SimpleNamespace(reference_structure=None, prefix='unused')
                checker.INPUT = {'general': {'assign_chainID': 0}}
                structure = SimpleNamespace(residues=[SimpleNamespace(chain='A')], save=lambda *a, **k: None)
                module = builder.__module__
                with patch(module + '.elapsed_progress', record), patch(module + '.check_str'):
                    checker.check_structures(structure, structure, structure)
                self.assertEqual(labels, ['Complex structure validation', 'Receptor structure validation',
                                          'Ligand structure validation'])


if __name__ == '__main__':
    unittest.main()
