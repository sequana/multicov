This is is the **coverage** pipeline from the `Sequana <https://sequana.readthedocs.org>`_ project


.. image:: https://badge.fury.io/py/sequana-multicov.svg
     :target: https://pypi.python.org/pypi/sequana_multicov

.. image:: http://joss.theoj.org/papers/10.21105/joss.00352/status.svg
    :target: http://joss.theoj.org/papers/10.21105/joss.00352
    :alt: JOSS (journal of open source software) DOI

.. image:: https://github.com/sequana/multicov/actions/workflows/main.yml/badge.svg
   :target: https://github.com/sequana/multicov/actions/workflows    


:Overview: Parallelised version of sequana_coverage for large eukaryotes genome.
:Input: A set of BAM or BED files. BED file must have 3 or 4 columns. First column is
    the chromosome/contig name, second column stored positions and third the
    coverage. Fourth optional columns contains a filtered coverage (not used in
    the analysis but shown in the HTML reports)
:Output: a set of HTML reports for each chromosomes and a multiqc report
:Status: production
:Citation: 
    Dimitri Desvillechabrol, Christiane Bouchier, Sean Kennedy, Thomas Cokelaer
    *Sequana coverage: detection and characterization of genomic variations 
    using running median and mixture models*
    GigaScience, Volume 7, Issue 12, December 2018, giy110, 
    https://doi.org/10.1093/gigascience/giy110

    and 

    Cokelaer et al, (2017), ‘Sequana’: a Set of Snakemake NGS pipelines, Journal of Open Source Software, 2(16), 352, JOSS DOI https://doi:10.21105/joss.00352


Installation
~~~~~~~~~~~~


sequana_multicov is based on Python3, just install the package as follows::

    pip install sequana_multicov --upgrade


Usage
~~~~~

::

    sequana_multicov --help
    sequana_multicov --input-directory DATAPATH 

By default, this looks for BED file. WARNING. This are BED3 meaning a 3-columns
tabulated file like this one::

    chr1 1 10
    chr1 2 11
    ...
    chr1 N1 10
    chr2 1 20
    chr2 2 21
    ...
    chr2 N2 20

where the first column stored the chromosome name, the second is the position
and the third is the coverage itself. See sequana_coverage documentation for
details. If you have BAM files as input, we will do the conversion for you. In
such case, use this option::

    --input-pattern "*.bam"

The sequana_coverage script creates a directory with the pipeline and 
its configuration file. You will then need 
to execute the pipeline::

    cd coverage
    sh coverage.sh  # for a local run

This launch a snakemake pipeline. If you are familiar with snakemake, you can 
retrieve the pipeline itself and its configuration files and then execute the pipeline yourself with specific parameters::

    snakemake -s multicov.rules -c config.yaml --cores 4 --stats stats.txt

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
- multiqc

.. .. image:: https://raw.githubusercontent.com/sequana/multicov/master/sequana_pipelines/multicov/dag.png


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

Here is the `latest documented configuration file <https://raw.githubusercontent.com/sequana/multicov/main/sequana_pipelines/multicov/config.yaml>`_
to be used with the pipeline. Each rule used in the pipeline may have a section in the configuration file. 


Changelog
~~~~~~~~~

========= ====================================================================
Version   Description
========= ====================================================================
1.1.0     * set apptainer containers and use wrappers
1.0.0     * renamed into multicov.
          * update to use latest sequana_pipetools (v0.9.2)
0.9.1     * rename genbank field into annotation, window into window_size
0.9.0     * first version
========= ====================================================================

