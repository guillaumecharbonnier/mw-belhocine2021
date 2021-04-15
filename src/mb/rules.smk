rule createWorkingTree:
	input:
	WDIR:"Mohamed/projects/Broad"

	shell:
	"
		mkdir -p  {input.WDIR}/PROGS  ;
		mkdir -p  {input.WDIR}/LOGS  ;
		mkdir -p  {input.WDIR}/ANNOTATIONS/TXT  ;
		mkdir -p  {input.WDIR}/ANNOTATIONS/GTF  ;
		mkdir -p  {input.WDIR}/DATA  ;
		mkdir -p  {input.WDIR}/DATA/RNA-SEQ  ;
		mkdir -p  {input.WDIR}/DATA/CHIP-SEQ  ;
		mkdir -p  {input.WDIR}/DATA/MICROARRAYS  ;
		mkdir -p  {input.WDIR}/RESULTS  ;
		mkdir -p  {input.WDIR}/DOCS ;
	"
	
rule createRNASeqDirectories4Sample:
	input:
	WDIR:"Mohamed/projects/Broad"
	shell:
	"
		mkdir -p  {input.WDIR}/DATA/RNA-SEQ/{sample} ;
		mkdir -p  {input.WDIR}/DATA/RNA-SEQ/{sample}/BAM ;
		mkdir -p  {input.WDIR}/DATA/RNA-SEQ/{sample}/BWIG ;
		mkdir -p  {input.WDIR}/DATA/RNA-SEQ/{sample}/FASTQ ;
		mkdir -p  {input.WDIR}/DATA/RNA-SEQ/{sample}/TRIMMED ;
		mkdir -p  {input.WDIR}/DATA/RNA-SEQ/{sample}/TOPHAT ;
		mkdir -p  {input.WDIR}/DATA/RNA-SEQ/{sample}/CUFFLINKS ;
		mkdir -p  {input.WDIR}/DATA/RNA-SEQ/{sample}/DESEQ2 ;
		mkdir -p  {input.WDIR}/DATA/RNA-SEQ/{sample}/PROGS ;
		mkdir -p  {input.WDIR}/DATA/RNA-SEQ/{sample}/LOGS ;
		mkdir -p  {input.WDIR}/DATA/RNA-SEQ/{sample}/QUAL_CHECK ;
	"

rule createChIPSeqDirectories4Sample:
	input:
	WDIR:"Mohamed/projects/Broad"
	shell:
	"
		mkdir -p  {input.WDIR}/DATA/CHIP-SEQ/{sample} ;
		mkdir -p  {input.WDIR}/DATA/CHIP-SEQ/{sample}/BAM ;
		mkdir -p  {input.WDIR}/DATA/CHIP-SEQ/{sample}/BWIG ;
		mkdir -p  {input.WDIR}/DATA/CHIP-SEQ/{sample}/FASTQ ;
		mkdir -p  {input.WDIR}/DATA/CHIP-SEQ/{sample}/TRIMMED ;
		mkdir -p  {input.WDIR}/DATA/CHIP-SEQ/{sample}/MACS2 ;
		mkdir -p  {input.WDIR}/DATA/CHIP-SEQ/{sample}/PROGS ;
		mkdir -p  {input.WDIR}/DATA/CHIP-SEQ/{sample}/LOGS ;
		mkdir -p  {input.WDIR}/DATA/CHIP-SEQ/{sample}/QUAL_CHECK ;
	"


######################################################################################################3


rule chipseq_fastqc:
    input:
        "DATA/CHIP-SEQ/{sample}/FASTQ/{sample}.{ext}"
    output:
	    html = "DATA/CHIP-SEQ/{sample}/QUAL_CHECK/{sample}_fastqc.html",
        zip  = "DATA/CHIP-SEQ/{sample}/QUAL_CHECK/{sample}_fastqc.zip"
    log:
               "DATA/CHIP-SEQ/{sample}/LOGS/fastqc_{sample}.log"

    wildcard_constraints:
        ext="fastq|fastq.gz|bam|sam"
    threads:
        MAX_THREADS
    shell:
        "fastqc --threads {threads} --outdir `dirname {output.html}` {input} &> {log}"
		

rule sickle_se:
    input:
        fastq = "DATA/CHIP-SEQ/{sample}/FASTQ/{sample}.{ext}"
    output:
        fastq = "DATA/CHIP-SEQ/{sample}/TRIMMED/{sample}_TRIM.{ext}"
    params:
        "-g -q 20 -t sanger"
    log:
        "DATA/CHIP-SEQ/{sample}/LOGS/sickle_{sample}.log"
    wildcard_constraints:
        ext="fastq|fastq.gz"
    threads:
        1
    shell:
        "sickle se -f {input.fastq} {params} -o {output.fastq} &> {log}"

rule sickle_pe:
    input:
        m1="DATA/CHIP-SEQ/{sample}/FASTQ/{sample}_1.{ext}",
        m2="DATA/CHIP-SEQ/{sample}/FASTQ/{sample}_2.{ext}"
    output:
        m1=    "DATA/CHIP-SEQ/{sample}/TRIMMED/{sample}_TRIM_1.{ext}",
        m2=    "DATA/CHIP-SEQ/{sample}/TRIMMED/{sample}_TRIM_2.{ext}",
        single="DATA/CHIP-SEQ/{sample}/TRIMMED/{sample}_TRIM_single.{ext}"
    params:
        gz = "-g"
    log:
        "DATA/CHIP-SEQ/{sample}/LOGS/sickle_{sample}.log"
    wildcard_constraints:
        ext="fastq|fastq.gz"
    shell:
        "sickle pe -f {input.m1} -r {input.m2} -q 20 -t sanger"
        "-o {output.m1} -p {output.m2} "
        "-s {output.single} {params} &> {log}"


rule star_se:
    input:
        fwd="DATA/CHIP-SEQ/{sample}/TRIMMED/{sample}_TRIM.fastq.gz",
        index="ANNOTATIONS/index/star/build_index/{index}/Genome",
        gtf="ANNOTATIONS/GTF/annotation_ensembl/{index}.gtf"
    output:
        bam="DATA/CHIP-SEQ/{sample}/BAM/{sample}.bam",
		log="DATA/CHIP-SEQ/{sample}/LOGS/STAR_{sample}.log"
    params:
        genomedir="ANNOTATIONS/index/star/build_index/{index}/",
        outdir="DATA/CHIP-SEQ/{sample}/STAR/"
    threads:
        12
    shell:
        """
        STAR \
            --genomeDir {params.genomedir} \
            --readFilesCommand zcat \
            -c --readFilesIn {input.fwd}\
            --runThreadN {threads} \
            --sjdbGTFfile {input.gtf} \
            --outFilterMismatchNoverLmax 0.05 \
            --outSAMtype BAM SortedByCoordinate \
            --outFilterMultimapNmax 1 \
            --genomeLoad NoSharedMemory
        mv Aligned.sortedByCoord.out.bam {output.bam}
        """

rule star_pe:
    input:
        fwd="DATA/CHIP-SEQ/{sample}/TRIMMED/{sample}_TRIM_1.fastq.gz",
        rev="DATA/CHIP-SEQ/{sample}/TRIMMED/{sample}_TRIM_2.fastq.gz",
        index="ANNOTATIONS/index/star/build_index/{index}/Genome",
        gtf="ANNOTATIONS/GTF/annotation_ensembl/{index}.gtf"
    output:
        bam="DATA/CHIP-SEQ/{sample}/BAM/{sample}.bam",
		log="DATA/CHIP-SEQ/{sample}/LOGS/STAR_{sample}.log"
    params:
        genomedir="ANNOTATIONS/index/star/build_index/{index}/",
        outdir="DATA/CHIP-SEQ/{sample}/STAR/"
    threads:
        12
    shell:
        """
        STAR \
            --genomeDir {params.genomedir} \
            --readFilesIn {input.fwd} {input.rev} \
            --readFilesCommand zcat \
            --runThreadN {threads} \
            --sjdbGTFfile {input.gtf} \
            --outFilterMismatchNoverLmax 0.05 \
            --outSAMtype BAM SortedByCoordinate \
            --outFilterMultimapNmax 1 \
            --genomeLoad NoSharedMemory
        mv Aligned.sortedByCoord.out.bam {output.bam}
        """

rule macs2_callpeak_broad_without_control:
    input:
        bam="DATA/CHIP-SEQ/{sample}/BAM/{sample}.sort.bam"
        bai="DATA/CHIP-SEQ/{sample}/BAM/{sample}.sort.bam.bai"
    output:
        bed =  "DATA/CHIP-SEQ/{sample}/MACS2/{sample}_peaks.broadPeak",
        xls =  "DATA/CHIP-SEQ/{sample}/MACS2/{sample}_peaks.xls"
    params:
        outdir="DATA/CHIP-SEQ/{sample}/MACS2/"
	shell:"""
    macs2 callpeak\
        --broad\
        -t {input.bam}\
        -f BAM\
        -g {gsize}\
		--broad-cutoff 0.01\
		--nomodel\
        --name {file}\
        --outdir {params.outdir}
    """

rule macs2_callpeak_broad_with_control:
    input:
        bam="DATA/CHIP-SEQ/{sample}/BAM/{sample}.sort.bam"
        bai="DATA/CHIP-SEQ/{sample}/BAM/{sample}.sort.bam.bai"
        bamc="DATA/CHIP-SEQ/{sample}/BAM/{control}.sort.bam"
        baic="DATA/CHIP-SEQ/{sample}/BAM/{control}.sort.bam.bai"
    output:
        bed =  "DATA/CHIP-SEQ/{sample}/MACS2/{sample}_peaks.broadPeak"
        xls =  "DATA/CHIP-SEQ/{sample}/MACS2/{sample}_peaks.xls"
    params:
        outdir="DATA/CHIP-SEQ/{sample}/MACS2/"
	shell:"""
    macs2 callpeak\
        --broad\
        -t {input.bam}\
		-c {input.bamc}\
        -f BAM\
        -g {gsize}\
		--broad-cutoff 0.01\
		--nomodel\
        --name {file}\
        --outdir {params.outdir}
    """

	
rule bedtools_sort:
	input:
		bed =  "DATA/CHIP-SEQ/{sample}/MACS2/{sample}_peaks.broadPeak"
	output:
		bed= "DATA/CHIP-SEQ/{sample}/MACS2/{sample}_peaks.broadPeak.sort"
	shell:"""
	bedtools sort -i {input.bed} > {output.bed}
	"""

rule closest_features:
	input:
		bed= "DATA/CHIP-SEQ/{sample}/MACS2/{sample}_peaks.broadPeak.sort"
	output:
		txt= "DATA/CHIP-SEQ/{sample}/MACS2/{sample}_peaks.broadPeak.sort.txt"
	params:
		ref="ANNOTATIONS/TXT/BED/{index}.bed"
	shell:"""
	closest-features --dist --delim '\t' {input.bed} {params.ref} > {output.txt}
	"""	

#####################################################################################################
		
rule rnaseq_fastqc:
    input:
        "DATA/RNA-SEQ/{sample}/FASTQ/{sample}.{ext}"
    output:
	    html = "DATA/RNA-SEQ/{sample}/QUAL_CHECK/{sample}_fastqc.html",
        zip  = "DATA/RNA-SEQ/{sample}/QUAL_CHECK/{sample}_fastqc.zip"
    log:
               "DATA/RNA-SEQ/{sample}/LOGS/fastqc_{sample}.log"

    wildcard_constraints:
        ext="fastq|fastq.gz|bam|sam"
    threads:
        MAX_THREADS
    shell:
        "fastqc --threads {threads} --outdir `dirname {output.html}` {input} &> {log}"

rule bowtie2_paired_end:
        input:
        read1="DATA/RNA-SEQ/{sample}/TRIMMED/{sample}_TRIM_1.fastq.gz",
        read2="DATA/RNA-SEQ/{sample}/TRIMMED/{sample}_TRIM_2.fastq.gz",
        
    output:
        sam                = "DATA/RNA-SEQ/{sample}/BAM/{sample}.sam",
        unmapped_single    = "DATA/RNA-SEQ/{sample}/BAM/unmapped/single.fastq.gz",
        unmapped_pair1     = "DATA/RNA-SEQ/{sample}/BAM/unmapped/pair.fastq.1.gz",
        unmapped_pair2     = "DATA/RNA-SEQ/{sample}/BAM/unmapped/pair.fastq.2.gz"
    log:
        "DATA/RNA-SEQ/{sample}/LOGS/BOWTIE_{sample}.log"
    params:
        unmapped_pair = "DATA/RNA-SEQ/{sample}/BAM/unmapped/pair.fastq.gz",
        index="ANNOTATIONS/index/Bowtie2Index/{index}/Genome",
    threads:
        MAX_THREADS
    shell:
        """
        bowtie2 -p {threads} -x {params.index}\
            -1 {input.read1} -2 {input.read2}\
            --un-gz {output.unmapped_single}\
            --un-conc-gz {params.unmapped_pair}\
            -S {output.sam}\
            2> {log}
        """
	
rule cufflinks:
    input:
	    sam="DATA/RNA-SEQ/{sample}/BAM/{sample}.sam",
        bam="DATA/RNA-SEQ/{sample}/BAM/{sample}.bam",
        gtfGuide=gtf="ANNOTATIONS/GTF/annotation_ensembl/{index}.gtf"
    output:
        gtf="DATA/RNA-SEQ/{sample}/CUFFLINKS/{sample}_transcripts.gtf"
    log:
        "DATA/RNA-SEQ/{sample}/LOGS/CUFFLINKS_{sample}.log"
    params:
        outdir="DATA/RNA-SEQ/{sample}/CUFFLINKS/"
    wildcard_constraints:
        libraryType="fr-unstranded|fr-firststrand|ff-firststrand|ff-secondstrand|fr-secondstrand|ff-unstranded|transfrags"
    threads:
        MAX_THREADS
    shell:
        """
		cufflinks \
            --library-type {wildcards.libraryType} \
            --GTF-guide {input.gtfGuide} \
            -p {threads}  \
            -o {params.outdir} \
            {input.bam} 2> {log}
        """

ALLSAMPLES=$(shell cd BAMs  ; ls --color=none *bam | perl -npe 's/accepted_hits.bam//g'| perl -npe 's/\n/,/g' )
ALLBAMS=$(shell ls --color=none BAMs/*bam )

run_cuffdiff:
    input:
		gtf="ANNOTATIONS/GTF/gencodeV19Lnc_mergeWithRefSeq.gtf"	
    output:
		dir="RESULTS"

	log:
		"RESULTS/CUFFDIFF.log"
    threads:
        MAX_THREADS

    shell:
	"""
	cuffdiff -o {output}\
	-L $(ALLSAMPLES) \
	--FDR 0.05 -p {threads} --no-diff --library-type fr-firststrand --library-norm-method classic-fpkm \
	{input.gtf} \
	$(ALLBAMS) \
	2> {log}
	"""
	
rule samtools_view:
    input:
        sam="DATA/RNA-SEQ/{sample}/BAM/{sample}.sam"
    output:
        bam="DATA/RNA-SEQ/{sample}/BAM/{sample}.bam"
	log:
        "DATA/RNA-SEQ/{sample}/LOGS/SAMTOOLS_View_{sample}.log"
    threads:
        MAX_THREADS
    shell:
        "samtools view -@ {threads} -S -b {input.sam} > {input.bam} &> {log}"


rule samtools_sort:
    input:
        bam="DATA/RNA-SEQ/{sample}/BAM/{sample}.bam"
    output:
        bam="DATA/RNA-SEQ/{sample}/BAM/{sample}.sort.bam"
    log:
        "DATA/RNA-SEQ/{sample}/LOGS/SAMTOOLS_Sort_{sample}.log"

    threads:
        MAX_THREADS
    shell:
        "samtools sort -@ {threads} {input.bam} -o {output.bam} &> {log}"


rule samtools_index:
    input:
        bam="DATA/RNA-SEQ/{sample}/BAM/{sample}.sort.bam"
    output:
        bam="DATA/RNA-SEQ/{sample}/BAM/{sample}.sort.bam"
        bai="DATA/RNA-SEQ/{sample}/BAM/{sample}.sort.bam.bai"
    log:
            "DATA/RNA-SEQ/{sample}/LOGS/SAMTOOLS_index_{sample}.log"
    threads:
        1
    shell:
        """
        samtools index -@ {threads} {input.bam}) &> {log}
        """

rule doBigWigScaled:
	input:
		bam="DATA/RNA-SEQ/{sample}/BAM/{sample}.sort.bam"
		chromsize="ANNOTATIONS/genomes/Homo_sapiens/UCSC/hg19/Annotation/Genes/ChromInfo.txt"
	output:
		bg="DATA/RNA-SEQ/{sample}/BAM/{sample}.sort.bg"
		bw="DATA/RNA-SEQ/{sample}/BWIG/{sample}.sort.bw"
    shell:
    """
	genomeCoverageBed -split -bg  \
	-ibam {input.bam} \
	-g {input.chromsize} > {output.bg}
	
	bedGraphToBigWig {output.bg}  \
	{input.chromsize} {output.bw}
	"""

rule deseq2:
	output:
		txt= "DATA/RNA-SEQ/RESULTS/DESeq2_results.txt"
	shell:"""
	Rscript deseq2_script.R
	"""
	
	
rule inflexion_point:
	output:
		txt= "DATA/CHIP-SEQ/{sample}/MACS2/Inflexion.txt"
	shell:"""
	Rscript inflexion_point_script.R
	"""

rule affy_processing:
	output:
		txt= "DATA/MICROARRAY/{project}/Expression.txt"
	shell:"""
	Rscript Affy_script.R
	"""	



EXTRA CODE CHECK gtftk tool for more information :
##################################################


rule gtftk_coverage:
    input:
        gtf = "",
        chrominfo = "",
        bw = ""
    output:
        txt="project/output/gtftk_coverage/{file}/{file}.txt"
    log:
        "project/output/gtftk_coverage/{file}/{file}.log"
    threads:
        MAX_THREADS
    shell:
        """
        gtftk coverage --inputfile {input.gtf} --outputfile {output.txt} --chrom-info {input.chrominfo} --nb-proc {threads} {input.bw} &> {log}
        """
##################################################
		
rule gtftk_mk_matrix:	
		input:
        gtf = "",
        bw = ""
    output:
        txt="project/output/gtftk_coverage/{file}/{file}.txt"
    log:
        "project/output/gtftk_coverage/{file}/{file}.log"
    threads:
        MAX_THREADS
    shell:
        """
        gtftk mk_matrix \
		--inputfile {input.gtf}\
		--outputfile {output.txt}\
		--bin-around-frac 0.5\
		-t transcript\
		-d 5000\
		-u 5000\
		-w 200\
		-c {index}\
		-l {file}\
		-y {input.bw}\
		--nb-proc {threads}\
		{input.bw} &> {log}
        """
##################################################

rule gtftk_profile :	
		input:
        covergae = "",
    output:
        txt="project/output/gtftk_coverage/{file}/{file}.txt"
    log:
        "project/output/gtftk_coverage/{file}/{file}.log"
    threads:
        MAX_THREADS
    shell:
        """
        gtftk profile \
		-D \
		-i {input.coverage}\
		-c 'red' \
		-d {file} \
		-o profile \
		-pf png -if  \
		{file}_profile.png
        """
##################################################

rule tophat:
	input:
		gtf=""
		m1=""
		m2=""
		index=""
	output:
		bam=""
		dir="DATA/RNA-SEQ/{sample}/TOPHAT"
	shell:
        """
		tophat2 -p 4 -C -x 1 -g 1 --bowtie2 --keep-tmp --library-type fr-firststrand \
		-G {input.gtf} \
		-o {input.dir}\
		{input.index} \
		{input.m1} \
		{input.m2} ;
		mv {output.dir}/accepted_hits.bam {output.bam} ;
		samtools index {output.bam}
		"""
		
		