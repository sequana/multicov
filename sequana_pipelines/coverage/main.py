import shutil
import sys
import os
import argparse

from sequana_pipetools.options import *
from sequana_pipetools.misc import Colors
from sequana_pipetools.info import sequana_epilog, sequana_prolog

col = Colors()

NAME = "coverage"




class Options(argparse.ArgumentParser):
    def __init__(self, prog=NAME, epilog=None):
        usage = col.purple(sequana_prolog.format(**{"name": NAME}))
        super(Options, self).__init__(usage=usage, prog=prog, description="",
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

    def parse_args(self, *args):
        args_list = list(*args)
        if "--from-project" in args_list:
            if len(args_list)>2:
                msg = "WARNING [sequana]: With --from-project option, " + \
                        "pipeline and data-related options will be ignored."
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
    from sequana_pipetools.options import before_pipeline
    before_pipeline(NAME)

    # option parsing including common epilog
    options = Options(NAME, epilog=sequana_epilog).parse_args(args[1:])


    from sequana.pipelines_common import SequanaManager

    # the real stuff is here
    manager = SequanaManager(options, NAME)

    # create the beginning of the command and the working directory
    manager.setup()

    # fill the config file with input parameters
    if options.from_project is None:
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
        cfg.coverage.window_size = options.window
        cfg.coverage.chunksize = options.chunksize
        cfg.coverage.binning = options.binning
        cfg.coverage.cnv_clustering = options.cnv_clustering




    # finalise the command and save it; copy the snakemake. update the config
    # file and save it.
    manager.teardown()


if __name__ == "__main__":
    main()
