process tamapper {
    // directives

     output:
        tuple val(prefix), path("${prefix}_BowtieMAP.txt"), emit: bam
        tuple val(prefix), path("${prefix}_BowtieMAP.log"), emit: log

    input:
        tuple val(prefix), path(aligned_reads) 
        //path(reference)
        //val(organism)

    output:
    // 7 files:
        // counter_file=prefix+"_TAmap.txt"
        // percent_file=prefix+"_TApercent.txt"
        // counts_file=prefix+"_TAHitFraction.txt"
        // artemis1_file=prefix+"_artemis_chrm1.txt"
        // artemis2_file=prefix+"_artemis_chrm2.txt"
        // artemis3_file=prefix+"_artemis_pAt.txt"
        // artemis4_file=prefix+"_artemis_pTi.txt"


        tuple val(prefix), path("${prefix}_BowtieMAP.txt"), emit: bam
        tuple val(prefix), path("${prefix}_BowtieMAP.log"), emit: log

    script:
    """
    (python TAmapperVC2019.py xxxxBowtieMAP.txt) 2> "${prefix}_BowtieMAP.log"
    
    """
}