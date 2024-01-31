process cutadapt {
    // directives
    container 'quay.io/biocontainers/cutadapt:4.6--py39hf95cd2a_1'

    input:
        tuple val(prefix), path(reads), val(paramaters), val(step) // muestraA, /mnt/c/muestraA, -d XXX | -a XXX | -m XXX

    output:
        tuple val(prefix), path("${prefix}_step${step}_cutadapt.fastq"), emit: sequence
        tuple val(prefix), path("${prefix}_step${step}_cutadapt.log"), emit: log

    script:
    """
    cutadapt ${paramaters} -o "${prefix}_step${step}_cutadapt.fastq" ${reads} > "${prefix}_step${step}_cutadapt.log"
    
    """
}