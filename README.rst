
.. image:: https://badge.fury.io/py/sequana-multicov.svg
     :target: https://pypi.python.org/pypi/sequana_multicov

.. image:: https://github.com/sequana/multicov/actions/workflows/main.yml/badge.svg
   :target: https://github.com/sequana/multicov/actions/workflows/main.yml

.. image:: https://img.shields.io/badge/python-3.11%20%7C%203.12-blue.svg
    :target: https://pypi.python.org/pypi/sequana_multicov
    :alt: Python 3.11 | 3.12

.. image:: http://joss.theoj.org/papers/10.21105/joss.00352/status.svg
    :target: http://joss.theoj.org/papers/10.21105/joss.00352
    :alt: JOSS (journal of open source software) DOI

This is the **multicov** pipeline from the `Sequana <https://sequana.readthedocs.org>`_ project

:Overview: Parallelised version of sequana_coverage for multi-sample genomic coverage analysis and CNV detection
:Input: A set of BED files (3 or 4 columns: chromosome, position, coverage, optional filtered coverage)
:Output: Per-sample HTML coverage reports, a MultiQC report, and a summary.html with links to all reports
:Status: Production
:Documentation: This README file and https://sequana.readthedocs.io
:Citation:
    Dimitri Desvillechabrol, Christiane Bouchier, Sean Kennedy, Thomas Cokelaer
    *Sequana coverage: detection and characterization of genomic variations
    using running median and mixture models*
    GigaScience, Volume 7, Issue 12, December 2018, giy110,
    https://doi.org/10.1093/gigascience/giy110

    and

    Cokelaer et al, (2017), 'Sequana': a Set of Snakemake NGS pipelines, Journal of Open Source Software, 2(16), 352, JOSS DOI https://doi:10.21105/joss.00352


Installation
~~~~~~~~~~~~

If you already have all requirements, install the package using pip::

    pip install sequana_multicov --upgrade


Usage
~~~~~

Scan BED files in a directory and set up the pipeline (replace ``DATAPATH`` with your input directory)::

    sequana_multicov --input-directory DATAPATH

To provide a reference FASTA file for GC content plots::

    sequana_multicov --input-directory DATAPATH --reference-file genome.fa

To provide a GenBank annotation file for event annotation::

    sequana_multicov --input-directory DATAPATH --annotation-file genome.gbk

This creates a ``multicov/`` directory with the pipeline and configuration file. Execute the pipeline locally::

    cd multicov
    sh multicov.sh

If you are familiar with Snakemake, you can also run the pipeline directly::

    snakemake -s multicov.rules --cores 4 --stats stats.txt

See ``.sequana/profile/config.yaml`` to tune Snakemake behaviour (cores, cluster settings, etc.).

Usage with apptainer
~~~~~~~~~~~~~~~~~~~~~

With apptainer, initiate the working directory as follows::

    sequana_multicov --input-directory DATAPATH  --apptainer-prefix ~/.sequana/apptainers

Images are downloaded in the shared location. Then::

    cd multicov
    sh multicov.sh


Input format
~~~~~~~~~~~~~

BED files must have 3 or 4 tab-separated columns::

    chr1    1    10
    chr1    2    11
    ...
    chr2    1    20
    chr2    2    21
    ...

where the first column is the chromosome/contig name, the second is the position (1-based, sorted), and the third is the coverage depth. An optional fourth column may contain a filtered coverage signal (shown in reports but not used in the analysis).

If you only have BAM files, convert them with::

    samtools depth -aa input.bam > output.bed

For a specific chromosome only::

    samtools depth -aa -r chr1 input.bam > chr1.bed

For CRAM files, convert to BAM first::

    samtools view -@ 4 -T reference.fa -b -o out.bam in.cram


Requirements
~~~~~~~~~~~~

This pipeline requires the following executables:

- **sequana_coverage** — from the `Sequana <https://sequana.readthedocs.io>`_ package (installed automatically)
- **multiqc** — aggregated HTML report across samples

Install all dependencies at once::

    mamba env create -f environment.yml


Details
~~~~~~~~~

This pipeline runs **sequana_coverage** in parallel across all input BED files. For each sample it produces a standalone HTML report with:

- coverage plots and running-median normalisation
- ROI (region of interest) detection using z-score thresholds
- CNV clustering
- GC content overlay (when a reference FASTA is provided)
- Event annotation (when a GenBank file is provided)

On success, a ``summary.html`` is generated listing all samples with direct links to their individual reports, plus a MultiQC report aggregating key statistics across samples.

For very large genomes the ``--binning`` and ``--chunksize`` options can be used to reduce memory usage.


Rules and configuration details
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here is the `latest documented configuration file <https://raw.githubusercontent.com/sequana/multicov/main/sequana_pipelines/multicov/config.yaml>`_
to be used with the pipeline. Each rule used in the pipeline may have a section in the configuration file.


Changelog
~~~~~~~~~

========= ====================================================================
Version   Description
========= ====================================================================
1.2.0     * convert packaging from setup.py to pyproject.toml (Poetry)
          * add apptainer container for sequana_coverage rule
          * add summary.html report with sample count and per-sample links
1.1.0     * set apptainer containers and use wrappers
1.0.0     * renamed into multicov
          * update to use latest sequana_pipetools (v0.9.2)
0.9.1     * rename genbank field into annotation, window into window_size
0.9.0     * first version
========= ====================================================================


Contribute & Code of Conduct
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To contribute to this project, please take a look at the
`Contributing Guidelines <https://github.com/sequana/sequana/blob/main/CONTRIBUTING.rst>`_ first. Please note that this project is released with a
`Code of Conduct <https://github.com/sequana/sequana/blob/main/CONDUCT.md>`_. By contributing to this project, you agree to abide by its terms.
