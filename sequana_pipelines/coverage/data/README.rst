The genbank, fasta and fastq can be found online as desribed in the sequana
coverage paper from giga science, and the notebook published oneline.

The BAM file can be generated using bwa (a wrapper available in sequana_mapping
is also possible):

    sequana_mapping --file1 JB409847_R1_clean.fastq.gz --reference JB409847.fa
    # the previous command creates a BAM file
    samtools depth -aa PREV_OUTPUT > JB409847.bed

