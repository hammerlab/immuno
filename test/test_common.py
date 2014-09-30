from os import remove
from os.path import split
from subprocess import CalledProcessError

from nose.tools import eq_, raises

from immuno.common import (
    normalize_chromosome_name,
    is_valid_peptide,
    splitext_permissive,
    find_paths,
)

def test_normalize_chromosome_name():
    assert normalize_chromosome_name("1") == "1"
    assert normalize_chromosome_name(1) == "1"
    assert normalize_chromosome_name("chr1") == "1"
    assert normalize_chromosome_name("chrX") == "X"
    assert normalize_chromosome_name("M") == "M"
    assert normalize_chromosome_name("MT") == "M"

def test_is_valid_peptide():
    assert is_valid_peptide("SYFPTHEI")
    assert not is_valid_peptide("X")
    assert not is_valid_peptide("XSYFPTEHI")
    assert not is_valid_peptide("Z")
    assert not is_valid_peptide("")
    assert not is_valid_peptide(1)
    assert not is_valid_peptide(None)
    assert not is_valid_peptide("_")
    assert not is_valid_peptide(".")

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
