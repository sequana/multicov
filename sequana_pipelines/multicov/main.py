#
#  This file is part of Sequana software
#
#  Copyright (c) 2016-2021 - Sequana Development Team
#
#  File author(s):
#      Thomas Cokelaer <thomas.cokelaer@pasteur.fr>
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
import shutil
import sys
import os
import argparse
import subprocess

from sequana_pipetools.options import *
from sequana_pipetools.options import before_pipeline
from sequana_pipetools.misc import Colors
from sequana_pipetools.info import sequana_epilog, sequana_prolog
from sequana_pipetools import SequanaManager

col = Colors()

NAME = "multicov"


class Options(argparse.ArgumentParser):
    def __init__(self, prog=NAME, epilog=None):
        usage = col.purple(sequana_prolog.format(**{"name": NAME}))
        super(Options, self).__init__(
            usage=usage,
            prog=prog,
            description="",
            epilog=epilog,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )

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

        self.add_argument("--run", default=False, action="store_true",
            help="execute the pipeline directly")

    def parse_args(self, *args):
        args_list = list(*args)
        if "--from-project" in args_list:
            if len(args_list) > 2:
                msg = (
                    "WARNING [sequana]: With --from-project option, "
                    + "pipeline and data-related options will be ignored."
                )
                print(col.error(msg))
            for action in self._actions:
                if action.required is True:
                    action.required = False
        options = super(Options, self).parse_args(*args)
        return options


def main(args=None):

    if args is None:
        args = sys.argv

    # whatever needs to be called by all pipeline before the options parsing
    before_pipeline(NAME)

    # option parsing including common epilog
    options = Options(NAME, epilog=sequana_epilog).parse_args(args[1:])

    # the real stuff is here
    manager = SequanaManager(options, NAME)

    # create the beginning of the command and the working directory
    manager.setup()
    from sequana import logger

    logger.setLevel(options.level)
    logger.name = "sequana_rnaseq"
    logger.info(f"#Welcome to sequana_multicov pipeline.")

    # fill the config file with input parameters
    if options.from_project is None:
        cfg = manager.config.config

        cfg.input_directory = os.path.abspath(options.input_directory)
        cfg.input_pattern = options.input_pattern


        cfg.sequana_coverage.circular = options.circular
        cfg.sequana_coverage.double_threshold = options.double_threshold

        if options.genbank:
            genbank = os.path.abspath(options.genbank)
            cfg.sequana_coverage.genbank_file = genbank
            if os.path.exists(genbank):
                shutil.copy(genbank, manager.workdir)
            else:
                raise IOError("{} not found".format(options.genbank))

        if options.reference:
            reference = os.path.abspath(options.reference)
            cfg.sequana_coverage.reference_file = reference
            if os.path.exists(reference):
                shutil.copy(reference, manager.workdir)
            else:
                raise IOError("{} not found".format(options.reference))

        cfg.sequana_coverage.high_threshold = options.high_threshold
        cfg.sequana_coverage.low_threshold = options.low_threshold
        cfg.sequana_coverage.mixture_models = options.mixture_models
        cfg.sequana_coverage.window_size = options.window
        cfg.sequana_coverage.chunksize = options.chunksize
        cfg.sequana_coverage.binning = options.binning
        cfg.sequana_coverage.cnv_clustering = options.cnv_clustering




    # finalise the command and save it; copy the snakemake. update the config
    # file and save it.
    manager.teardown()

    if options.run:
        subprocess.Popen(["sh", "{}.sh".format(NAME)], cwd=options.workdir)


if __name__ == "__main__":
    main()
