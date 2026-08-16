import io
import logging
import os
import re
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from GMXMMPBSA.progress import MAX_RICH_WIDTH, FrameCounter, _StallNotifier, monitor_progress, resolve_progress_style


class _Stream(io.StringIO):
    def __init__(self, is_terminal):
        super().__init__()
        self.is_terminal = is_terminal

    def isatty(self):
        return self.is_terminal


class ProgressStyleTest(unittest.TestCase):
    def test_auto_uses_classic_when_mpi_hides_terminal(self):
        self.assertEqual(resolve_progress_style('auto', _Stream(False)), 'classic')

    def test_auto_uses_rich_for_terminal(self):
        with patch.dict(os.environ, {'TERM': 'xterm-256color'}):
            self.assertEqual(resolve_progress_style('auto', _Stream(True)), 'rich')

    def test_explicit_styles_are_preserved(self):
        for style in ('classic', 'plain', 'none'):
            self.assertEqual(resolve_progress_style(style, _Stream(False)), style)


class FrameCounterTest(unittest.TestCase):
    marker = '                    FINAL RESULTS\n'

    def test_counts_only_new_frames_across_rank_files(self):
        with tempfile.TemporaryDirectory() as directory:
            template = str(Path(directory, 'complex_gb.mdout.%d'))
            Path(template % 0).write_text(self.marker)
            Path(template % 1).write_text(self.marker * 2)
            counter = FrameCounter(template, mpi_size=2)
            self.assertEqual(counter.count(), 3)

            with Path(template % 0).open('a') as output:
                output.write(self.marker * 2)
            self.assertEqual(counter.count(), 5)
            self.assertEqual(counter.count(), 5)

    def test_resets_a_truncated_output(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory, 'complex_gb.mdout.0')
            output.write_text(self.marker * 2)
            counter = FrameCounter(str(output), mpi_size=1)
            self.assertEqual(counter.count(), 2)

            output.write_text(self.marker)
            self.assertEqual(counter.count(), 1)

    def test_nmode_marker(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory, 'complex_nm.out.0')
            output.write_text('Total: 1.0\nTotal: 2.0\n')
            self.assertEqual(FrameCounter(str(output), nmode=True).count(), 2)


class PlainProgressTest(unittest.TestCase):
    def test_stall_notices_are_rate_limited(self):
        notifier = _StallNotifier(total=10, timeout=5, started=0)
        self.assertIsNone(notifier.check(0, now=4))
        self.assertEqual(notifier.check(0, now=5), 5)
        self.assertIsNone(notifier.check(0, now=9))
        self.assertEqual(notifier.check(0, now=10), 10)
        self.assertIsNone(notifier.check(1, now=11))

    def test_plain_monitor_reports_completion(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory, 'complex_gb.mdout.0')
            output.write_text(FrameCounterTest.marker * 2)
            with self.assertLogs(level=logging.INFO) as messages:
                monitor_progress(str(output), nframes=2, style='plain', label='Complex', poll_interval=0)
            text = '\n'.join(messages.output)
            self.assertIn('Complex completed: 2 frames in', text)
            self.assertEqual(text.count('Complex completed:'), 1)

    def test_rich_monitor_logs_completion_summary(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory, 'receptor_gb.mdout.0')
            output.write_text(FrameCounterTest.marker * 3)
            with self.assertLogs(level=logging.INFO) as messages:
                monitor_progress(
                    str(output), nframes=3, style='rich', label='Receptor',
                    stream=_Stream(True), poll_interval=0,
                )
            self.assertIn('Receptor completed: 3 frames in', '\n'.join(messages.output))

    def test_rich_monitor_logs_cluster_progress_checkpoints(self):
        stream = _Stream(True)
        with patch('GMXMMPBSA.progress.FrameCounter.count', side_effect=[5, 10]):
            with self.assertLogs(level=logging.DEBUG) as messages:
                monitor_progress(
                    'complex_gb.mdout.%d', nframes=10, mpi_size=4, style='rich',
                    label='Complex', stream=stream, poll_interval=0,
                )
        text = '\n'.join(messages.output)
        self.assertIn('Complex progress: 5/10 frames (50%)', text)
        self.assertIn('[4 MPI ranks]', text)
        self.assertIn('Complex completed: 10 frames in', text)

    def test_none_does_not_read_outputs(self):
        monitor_progress('/does/not/exist.%d', nframes=2, style='none', poll_interval=0)

    def test_rich_monitor_renders_to_terminal_stream(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory, 'complex_gb.mdout.0')
            output.write_text(FrameCounterTest.marker)
            stream = _Stream(True)
            monitor_progress(
                str(output), nframes=1, style='rich', label='Complex', stream=stream,
                poll_interval=0,
            )
            self.assertIn('Complex', stream.getvalue())

    def test_forced_rich_renders_when_mpi_hides_terminal(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory, 'complex_gb.mdout.0')
            output.write_text(FrameCounterTest.marker)
            stream = _Stream(False)
            monitor_progress(
                str(output), nframes=1, style='rich', label='Complex', stream=stream,
                poll_interval=0,
            )
            self.assertIn('Complex', stream.getvalue())
            self.assertNotIn('frame/s\n', stream.getvalue())

    def test_rich_output_is_capped_on_wide_terminals(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory, 'complex_gb.mdout.0')
            output.write_text(FrameCounterTest.marker)
            stream = _Stream(True)
            with patch('GMXMMPBSA.progress.shutil.get_terminal_size', return_value=os.terminal_size((240, 40))):
                monitor_progress(
                    str(output), nframes=1, style='rich', label='Complex', stream=stream,
                    poll_interval=0,
                )
            visible_lines = [
                re.sub(r'\x1b\[[0-?]*[ -/]*[@-~]', '', line)
                for line in stream.getvalue().splitlines() if line
            ]
            self.assertTrue(visible_lines)
            self.assertTrue(all(len(line) <= MAX_RICH_WIDTH for line in visible_lines))


if __name__ == '__main__':
    unittest.main()
