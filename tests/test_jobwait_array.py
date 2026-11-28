"""Slurm can show an array under individual or compressed task identifiers."""
import pytest

from smj.junction.jobwait import wait_many


@pytest.mark.parametrize("listed", ["712_3", "712_[0-11%2]"])
def test_active_array_survives_two_queue_polls_until_sentinel(listed):
    class Queue:
        polls = 0

        def run(self, command, timeout):
            self.polls += 1
            if self.polls < 3:
                return 0, f"@@QRC 0\n@@ID {listed}\n@@END", ""
            return 0, "@@SENT 0 0 today\n@@QRC 0\n@@END", ""

    matched, unrelated = wait_many(
        Queue(), [{"remote_dir": "/a", "job_id": "712"},
                  {"remote_dir": "/b", "job_id": "71"}],
        poll_s=0, sleep=lambda _: None)
    assert matched["finished"] and not matched["failed"] and matched["polls"] == 3
    assert unrelated["failed"] and unrelated["polls"] == 2
