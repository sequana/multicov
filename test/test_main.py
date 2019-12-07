import easydev
import os
import tempfile
import subprocess
import sys


sequana_path = easydev.get_package_location('sequana_coverage')
sharedir = os.sep.join([sequana_path , "sequana_pipelines/coverage/data"])
annotation = sharedir + os.sep + "JB409847.gbk"
reference = sharedir + os.sep + "JB409847.fa"

def test_standalone_subprocess():
    directory = tempfile.TemporaryDirectory()
    cmd = "sequana_pipelines_coverage --input-directory {} "
    cmd += "--working-directory {} --force --annotation {} "
    cmd += " --reference {} -o "
    cmd = cmd.format(sharedir, directory.name, annotation, reference)
    subprocess.call(cmd.split())

def test_standalone_script():
    directory = tempfile.TemporaryDirectory()
    import sequana_pipelines.coverage.main as m
    sys.argv = ["test", "--input-directory", directory.name, 
         "--force"]
    m.main()

def test_wrong_reference():
    directory = tempfile.TemporaryDirectory()
    import sequana_pipelines.coverage.main as m
    sys.argv = ["test", "--input-directory", directory.name, 
        "--force", "--reference", "wrong"]
    try:
        m.main()
        assert False
    except IOError:
        assert True

def test_wrong_genbank():
    directory = tempfile.TemporaryDirectory()
    import sequana_pipelines.coverage.main as m
    sys.argv = ["test", "--input-directory", directory.name, 
        "--force", "--genbank", "wrong"]
    try:
        m.main()
        assert False
    except IOError:
        assert True


# For travis, you may want to add this with snakemake:
#if "TRAVIS_PYTHON_VERSION" in os.environ:
#    cmd += ["--snakemake-jobs", "1"]

def test_check_output():

    with tempfile.TemporaryDirectory() as wk:

        # create the wokring directory and script
        cmd = "sequana_pipelines_coverage --input-directory {} "
        cmd += "--working-directory {} --force --annotation {} "
        cmd += " --reference {} -o "
        cmd = cmd.format(sharedir, wk, annotation, reference)
        subprocess.call(cmd.split())

        subprocess.call("sh coverage.sh".split(), cwd=wk)


