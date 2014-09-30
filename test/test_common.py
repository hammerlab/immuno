from os.path import split 
from subprocess import CalledProcessError
from nose.tools import eq_, raises

from immuno.common import (
    splitext_permissive, find_paths,
    AsyncProcess, run_command, run_multiple_commands
)

def test_splitext_permissive():
    base, ext = splitext_permissive("", [".txt", ".gz"])
    eq_(base, "")
    eq_(ext, "")
    base, ext = splitext_permissive("apple", [".txt", ".gz"])
    eq_(base, "apple")
    eq_(ext, "")
    base, ext = splitext_permissive("apple.banana.tar.gz", [".txt", "png"])
    eq_(base, "apple.banana.tar")
    eq_(ext, ".gz")
    base, ext = splitext_permissive("apple.banana.tar.gz", [".txt", ".gz"])
    eq_(base, "apple.banana")
    eq_(ext, ".tar")
    base, ext = splitext_permissive("apple.banana.tar.gz.txt", 
            [".txt", ".tar"])
    eq_(base, "apple.banana.tar")
    eq_(ext, ".gz")
    base, ext = splitext_permissive("apple.banana.tar.gz.txt", 
            [".txt", ".tar", ".gz"])
    eq_(base, "apple")
    eq_(ext, ".banana")
    base, ext = splitext_permissive("apple.banana.tar.gz.txt", 
            ["apple", ".banana", ".tar", ".gz", ".txt"])
    eq_(base, "apple")
    eq_(ext, "")

@raises(ValueError)    
def test_splitext_permissive_error1():
    base, ext = splitext_permissive("", "")

@raises(ValueError)
def test_splitext_permissive_error2():
    base, ext = splitext_permissive("test", [".tar", ""])

def test_find_paths():
    curr_dir = split(__file__)[0]
    test_files = find_paths(directory_string = curr_dir, extensions = [".py"])
    assert len(test_files) > 0

@raises(OSError)
def test_find_paths_wrong_dir():
    curr_dir = split(__file__)[0]
    wrong_dir = curr_dir + "_NONSENSE_!!!!!"
    test_files = find_paths(directory_string = wrong_dir)

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
