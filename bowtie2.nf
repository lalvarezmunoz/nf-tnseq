process bowtie2 {
    // directives
    container 'quay.io/biocontainers/bowtie2:2.5.3--py310ha0a81b8_0'
    
    input:
        tuple val(prefix), path(trimmed_reads) 
        path(reference)
        val(organism)

    output:
        tuple val(prefix), path("${prefix}_BowtieMAP.txt"), emit: bam
        tuple val(prefix), path("${prefix}_BowtieMAP.log"), emit: log

    script:
    """
    (bowtie2 -p 8 -x ${reference}/${organism} ${trimmed_reads} "${prefix}_BowtieMAP.txt") 2> "${prefix}_BowtieMAP.log"
    
    """
}