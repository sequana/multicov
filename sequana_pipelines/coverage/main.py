import shutil
import sys
import os
import argparse

from sequana.pipelines_common import *
from sequana.snaketools import Module
from sequana import logger
logger.level = "INFO"

col = Colors()

NAME = "coverage"
m = Module(NAME)
m.is_executable()


class Options(argparse.ArgumentParser):
    def __init__(self, prog=NAME):
        usage = col.purple(
            """This script prepares the sequana pipeline coverage layout to
            include the Snakemake pipeline and its configuration file ready to
            use.

            In practice, it copies the config file and the pipeline into a
            directory (coverage) together with an executable script

            For a local run, use :

                sequana_pipelines_coverage --input-directory PATH_TO_DATA --run-mode local

            For a run on a SLURM cluster:

                sequana_pipelines_coverage --input-directory PATH_TO_DATA --run-mode slurm

        """
        )
        super(Options, self).__init__(usage=usage, prog=prog, description="")

        # add a new group of options to the parser
        so = SlurmOptions()
        so.add_options(self)

        # add a snakemake group of options to the parser
        so = SnakemakeOptions(working_directory=NAME)
        so.add_options(self)

        so = InputOptions(input_pattern="*.bed")
        so.add_options(self)

        so = GeneralOptions()
        so.add_options(self)

        pipeline_group = self.add_argument_group("pipeline")

        pipeline_group.add_argument("-o", "--circular", action="store_true")
        pipeline_group.add_argument("--double-threshold", default=0.5)
        pipeline_group.add_argument("--genbank", default=None,
            help="the genbank to annotate the events found")
        pipeline_group.add_argument("--reference", default=None,
            help="the genome reference used to plot GC content")
        pipeline_group.add_argument("--high-threshold", default=4)
        pipeline_group.add_argument("--low-threshold", default=-4)
        pipeline_group.add_argument("--mixture-models", default=2, type=int,
            help="""Number of models to use in the mixture model. (default 2).
                 No need to change this value. Possibly, you may want to set 
                 to 1 or 3 in some rate occasions. """)
        pipeline_group.add_argument("--window", default=20000, type=int, 
            help="""Length of the running median window. Keep to 20000 as much as
            possible. This allows the detection of CNV up to 10kb. If longer
            event are present, increase this window size.""")
        pipeline_group.add_argument("--chunksize", default=5000000, type=int)
        pipeline_group.add_argument("--binning", default=-1, type=int)
        pipeline_group.add_argument("--cnv-clustering", default=-1)


def main(args=None):

    if args is None:
        args = sys.argv

    options = Options(NAME).parse_args(args[1:])

    manager = PipelineManager(options, NAME)

    # create the beginning of the command and the working directory
    manager.setup()

    # fill the config file with input parameters
    cfg = manager.config.config

    cfg.input_directory = os.path.abspath(options.input_directory)
    cfg.input_pattern = options.input_pattern


    cfg.coverage.circular = options.circular
    cfg.coverage.double_threshold = options.double_threshold

    if options.genbank:
        genbank = os.path.abspath(options.genbank)
        cfg.coverage.genbank_file = genbank
        if os.path.exists(genbank):
            shutil.copy(genbank, manager.workdir)
        else:
            raise IOError("{} not found".format(options.genbank))

    if options.reference:
        reference = os.path.abspath(options.reference)
        cfg.coverage.reference_file = reference
        if os.path.exists(reference):
            shutil.copy(reference, manager.workdir)
        else:
            raise IOError("{} not found".format(options.reference))

    cfg.coverage.high_threshold = options.high_threshold
    cfg.coverage.low_threshold = options.low_threshold
    cfg.coverage.mixture_models = options.mixture_models
    cfg.coverage.window = options.window
    cfg.coverage.chunksize = options.chunksize
    cfg.coverage.binning = options.binning
    cfg.coverage.cnv_clustering = options.cnv_clustering




    # finalise the command and save it; copy the snakemake. update the config
    # file and save it.
    manager.teardown()


if __name__ == "__main__":
    main()
