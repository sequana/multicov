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
    cmd += "--working-directory {} --force --annotation-file {} "
    cmd += " --reference-file {} --circular "
    cmd = cmd.format(sharedir, directory, annotation, reference)
    subprocess.call(cmd.split())


def test_standalone_script(tmpdir):
    directory = tmpdir.mkdir("wkdir")
    import sequana_pipelines.multicov.main as m
    sys.argv = ["test", "--input-directory", sharedir, "--force"]
    m.main(standalone_mode=False)


def test_wrong_reference(tmpdir):
    import sequana_pipelines.multicov.main as m
    directory = tmpdir.mkdir("wkdir")
    sys.argv = ["test", "--input-directory", str(directory),
        "--force", "--reference-file", "wrong"]
    try:
        m.main(standalone_mode=False)
        assert False
    except IOError:
        assert True


def test_wrong_annotation(tmpdir):
    directory = tmpdir.mkdir("wkdir")
    import sequana_pipelines.multicov.main as m
    sys.argv = ["test", "--input-directory", str(directory),
        "--force", "--annotation-file", "wrong"]
    try:
        m.main(standalone_mode=False)
        assert False
    except IOError:
        assert True


def test_check_output(tmpdir):
    wkdir = tmpdir.mkdir("wkdir")

    # create the working directory and script
    cmd = f"sequana_multicov --input-directory {test_dir}/data "
    cmd += f"--working-directory {wkdir} --force --annotation-file {annotation} "
    cmd += f" --reference-file {reference} --circular "
    subprocess.call(cmd.split())

    subprocess.call("bash multicov.sh".split(), cwd=wkdir)
