params.bam = "."

bam_set = Channel.fromPath("${params.bam}/*.bam")
kmer_set = Channel.from(2..12)
bam_set.combine(kmer_set).into { bam_kmer }

process kmer_count {

    cpus 8
    memory '16 GB'

    tag { "${bam} - ${kmer} "}

    input:
        set file(bam), val(kmer) from bam_kmer
    output:
        file("out.tsv.gz") into out

    script:
        strain_name = bam.getName().replace(".bam", "")
    """
        samtools fastq ${bam} | fastq-kmers -k ${kmer} | \\
        awk ' { print \$0 "\\t${kmer}\\t${strain_name}\\tfastq-kmer" }' | gzip > out.tsv.gz
    """
}

process kmer_join {

    publishDir 'results'

    input:
        file("out*.tsv.gz") from out.toSortedList()

    output:
        file("kmers.tsv.gz")

    """
        zcat *.tsv.gz | gzip > kmers.tsv.gz
    """
}
