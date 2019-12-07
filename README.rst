This is is the **coverage** pipeline from the `Sequana <https://sequana.readthedocs.org>`_ projet

:Overview: Parallelised version of sequana_coverage for large eukaryotes genome.
:Input: A set of BAM or BED files. BED file must have 3 or 4 columns. First column is
    the chromosome/contig name, second column stored positions and third the
    coverage. Fourth optional columns contains a filtered coverage (not used in
    the analysis but shown in the HTML reports)
:Output: a set of HTML reports for each chromosomes and a multiqc report
:Status: production
:Citation: 
    - Sequana coverage (https://doi.org/10.1101/092478)
    - Cokelaer et al, (2017), ‘Sequana’: a Set of Snakemake NGS pipelines, Journal of Open Source Software, 2(16), 352, JOSS DOI doi:10.21105/joss.00352


Installation
~~~~~~~~~~~~

You must install Sequana first::

    pip install sequana

Then, just install this package::

    pip install sequana_coverage

This gives an executable called sequana_pipelines_coverage. Note that is should
not be confused with the original sequana_coverage standalone from Sequana
library. Indeed, this pipeline calls sequana_coverage behund the scene. 

Usage
~~~~~

::

    sequana_pipelines_coverage --help
    sequana_pipelines_coverage --input-directory DATAPATH 

This creates a directory with the pipeline and configuration file. You will then need 
to execute the pipeline::

    cd coverage
    sh coverage.sh  # for a local run

This launch a snakemake pipeline. If you are familiar with snakemake, you can 
retrieve the pipeline itself and its configuration files and then execute the pipeline yourself with specific parameters::

    snakemake -s coverage.rules -c config.yaml --cores 4 --stats stats.txt

Or use `sequanix <https://sequana.readthedocs.io/en/master/sequanix.html>`_ interface as follows::

    sequanix -w analysis -i . -p coverage

Go to the second panel, in Input data and then in Input directory. There, you
must modify the pattern (empty field by default meaning search for fastq files)
and set the field to either::

    *.bed

or::

    *.bam


You are ready to go. Save the project and press Run. Once done, open the HTML report.




Requirements
~~~~~~~~~~~~

This pipelines requires the following executable(s):

- sequana_coverage from **Sequana**, which should be installed automatically.

.. image:: https://raw.githubusercontent.com/sequana/sequana_coverage/master/sequana_pipelines/coverage/dag.png


Details
~~~~~~~~~

This pipeline runs **coverage** in parallel on the input BAM files (or BED file). 


The coverage tool takes as input a BAM or a BED file. The BED file must have 3
or 4 columns as explained in the standalone application (sequana_coverage) 
`documentation <http://sequana.readthedocs.io/en/master/applications.html?highlight=coverage#sequana-coverage>`_. 
In short, the first column is the chromosome name, the second column is the
position (sorted) and the third column is the coverage (an optional fourth
column would contain a coverage signal, which could be high quality coverage for
instance).

If you have only BAM files, you can convert them using **bioconvert** tool or
the command::

    samtools depth -aa input.bam > output.bed

If you have a CRAM file::

    samtools view -@ 4 -T reference.fa -b -o out.bam  in.cram

For very large BAM/BED files, we recommend to split the BED file by
chromosomes. For instance for the chromosome  chr1, type::

    # samtools index in.bam
    samtools depth -aa input.bam -r chr1 in.bam > chr1.bed

The standalone or Snakemake application can also take as input your BAM file and
will convert it automatically into a BED file.


Rules and configuration details
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here is the `latest documented configuration file <https://raw.githubusercontent.com/sequana/sequana_coverage/master/sequana_pipelines/coverage/config.yaml>`_
to be used with the pipeline. Each rule used in the pipeline may have a section in the configuration file. 

