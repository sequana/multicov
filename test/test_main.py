import os
import subprocess
import sys

from . import test_dir

sharedir = f"{test_dir}/data"

annotation = f"{test_dir}/data/JB409847.gbk"
reference = f"{test_dir}/data/JB409847.fa"

def test_standalone_subprocess(tmpdir):
    directory = tmpdir.mkdir("wkdir")

    cmd = "sequana_multicov --input-directory {} "
    cmd += "--working-directory {} --force --annotation {} "
    cmd += " --reference {} -o "
    cmd = cmd.format(sharedir, directory, annotation, reference)
    subprocess.call(cmd.split())

def test_standalone_script(tmpdir):
    directory = tmpdir.mkdir("wkdir")
    import sequana_pipelines.multicov.main as m
    sys.argv = ["test", "--input-directory", sharedir, 
         "--force"]
    m.main()

def test_wrong_reference(tmpdir):
    import sequana_pipelines.multicov.main as m
    directory = tmpdir.mkdir("wkdir")
    sys.argv = ["test", "--input-directory", str(directory), 
        "--force", "--reference", "wrong"]
    try:
        m.main()
        assert False
    except IOError:
        assert True

def test_wrong_genbank(tmpdir):
    directory = tmpdir.mkdir("wkdir")
    import sequana_pipelines.multicov.main as m
    sys.argv = ["test", "--input-directory", str(directory), 
        "--force", "--genbank", "wrong"]
    try:
        m.main()
        assert False
    except IOError:
        assert True

def test_check_output(tmpdir):
    wkdir = tmpdir.mkdir("wkdir")

    # create the wokring directory and script
    cmd = f"sequana_multicov --input-directory {test_dir}/data "
    cmd += f"--working-directory {wkdir} --force --annotation {annotation} "
    cmd += f" --reference {reference} -o "
    subprocess.call(cmd.split())

    subprocess.call("sh multicov.sh".split(), cwd=wkdir)


