Sequana coverage pipeline documentation
#####################################################

|version|, |today|, status: production

We strongly recommend to first test the standalone **sequana_coverage**,
already in the **Sequana** suite. If this is not fast enough, then this 
pipeline may be useful.

The **coverage** pipeline is a `Sequana <https://github.com/sequana/sequana>`_ pipeline. You can find the source code
on  `https://github.com/sequana/sequana_coverage <https://github.com/sequana/sequana_coverage/>`_. Would you have issues
about the code, usage or lack of information, please fill a report
on `Sequana itself <https://github.com/sequana/sequana/issues>`_ indicating the pipeline name (We centralized all
pipelines issues on **Sequana** repository only so as to be more responsive).

If you use **Sequana**, please do not forget to cite us:

    Cokelaer et al, (2017), 'Sequana': a Set of Snakemake NGS pipelines, Journal of
    Open Source Software, 2(16), 352, `JOSS DOI doi:10.21105/joss.00352 <http://www.doi2bib.org/bib/10.21105/joss.00352>`_

This pipeline has also its own citation:

    Dimitri Desvillechabrol, Christiane Bouchier, Sean Kennedy, Thomas Cokelaer
    *Sequana coverage: detection and characterization of genomic variations 
    using running median and mixture models*
    GigaScience, Volume 7, Issue 12, December 2018, giy110, 
    https://doi.org/10.1093/gigascience/giy110



.. contents::
   :depth: 2

What is Sequana ?
=====================

**Sequana** is a versatile tool that provides

#. A Python library dedicated to NGS analysis (e.g., tools to visualise standard NGS formats).
#. A set of Pipelines dedicated to NGS in the form of Snakefiles
#. Standalone applications
    #. sequana_coverage ease the
       extraction of genomic regions of interest and genome coverage information
    #. sequana_taxonomy performs a quick
       taxonomy of your FastQ. This requires dedicated databases to be downloaded.
    #. Sequanix, a GUI for Snakemake workflows (hence Sequana pipelines as well)

To join the project, please let us know on `github <https://github.com/sequana/sequana/issues/306>`_.

For more information, please see `github <https://sequana.readthedocs.io>`_.

The Sequana coverage pipeline
==============================================

.. include:: ../README.rst


Rules
=====

.. snakemakerule:: sequana_coverage

