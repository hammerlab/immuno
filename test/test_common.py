from immuno.common import splitext_permissive
from nose.tools import eq_, raises

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
