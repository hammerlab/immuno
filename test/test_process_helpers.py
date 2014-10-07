import tempfile
from subprocess import CalledProcessError
from os import remove

from immuno.process_helpers import (
    AsyncProcess,
    run_command,
    run_multiple_commands,
    run_multiple_commands_redirect_stdout
)
from nose.tools import raises
def test_async_ls():
    process = AsyncProcess(["ls"])
    process.wait()

def test_run_ls():
    run_command(["ls", "-als"])

@raises(OSError)
def test_run_bad_command():
    run_command(["__DSH#*#*&SHJ"])

@raises(CalledProcessError)
def test_run_invalid_args():
    run_command(["ls", "-Z_z"], suppress_stderr = True)

def test_run_mulitple_ls():
    run_multiple_commands([["ls"], ["ls"]])

def test_run_multiple_redirect():
    """
    Create two temporary files, run ls and redirect its output
    to the temporary files, check that they're not empty
    """

    t1 = tempfile.NamedTemporaryFile(mode='w', delete=False)
    t2 = tempfile.NamedTemporaryFile(mode='w', delete=False)
    run_multiple_commands_redirect_stdout({
        t1: ['ls'],
        t2: ['ls', '-als'],
    })
    t1.close()
    t2.close()
    with open(t1.name, 'r') as t1:
        assert len(t1.read()) > 0
    remove(t1.name)

    with open (t2.name, 'r') as t2:
        assert len(t2.read()) > 0
    remove(t2.name)
