rule r_compute_pausing_index:
    input:
        'out/deepTools/bamCoverage_--binSize_20_--minMappingQuality_0_--normalizeUsing_RPKM/samtools/index/samtools/sort/samtools/view_sam_to_bam_-q_30/bowtie2/se_hg19/sickle/se_-t_sanger_-q_20/ln/alias/sst/all_samples/fastq/Jurkat_SRR2157609_Pol2.bw',
        'out/r/extract_control_sharp_subsets_from_Atak2013/NM_Broad.bed',
        'out/r/extract_control_sharp_subsets_from_Atak2013/NM_Sharp1.bed',
        'out/r/extract_control_sharp_subsets_from_Atak2013/NM_Sharp2.bed',
        'out/r/extract_control_sharp_subsets_from_Atak2013/NM_Sharp3.bed'
    output:
        'out/r/compute_pausing_index/pi_boxplot.pdf',
        'out/r/compute_pausing_index/pi_boxplot.png',
    conda:
        "../../../../mw-lib/src/snakemake/envs/r_brgenomics.yaml"
    script:
        "../../r/script/compute_pausing_index.R"

rule r_compute_pausing_index_from_our_rnaseq:
    input:
        'out/deepTools/bamCoverage_--binSize_20_--minMappingQuality_0_--normalizeUsing_RPKM/samtools/index/samtools/sort/samtools/view_sam_to_bam_-q_30/bowtie2/se_hg19/sickle/se_-t_sanger_-q_20/ln/alias/sst/all_samples/fastq/Jurkat_SRR2157609_Pol2.bw',
        'out/r/extract_control_sharp_subsets_from_our_rnaseq/NM_genes_associated_with_Broad.bed',
        'out/r/extract_control_sharp_subsets_from_our_rnaseq/NM_genes_associated_with_Sharp1.bed',
        'out/r/extract_control_sharp_subsets_from_our_rnaseq/NM_genes_associated_with_Sharp2.bed',
        'out/r/extract_control_sharp_subsets_from_our_rnaseq/NM_genes_associated_with_Sharp3.bed'
    output:
        'out/r/compute_pausing_index_in_broad_and_sharp_subsets_from_our_rnaseq/pi_boxplot.pdf',
        'out/r/compute_pausing_index_in_broad_and_sharp_subsets_from_our_rnaseq/pi_boxplot.png',
        'out/r/compute_pausing_index_in_broad_and_sharp_subsets_from_our_rnaseq/pi_genebodies_trim_500b_boxplot.png',
    conda:
        "../../../../mw-lib/src/snakemake/envs/r_brgenomics.yaml"
    script:
        "../../r/script/compute_pausing_index_in_broad_and_sharp_subsets_from_our_rnaseq.R"


