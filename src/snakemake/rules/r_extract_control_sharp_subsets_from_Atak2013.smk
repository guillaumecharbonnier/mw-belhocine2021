JURKAT_BROAD_AND_SHARP_SUBSETS=['Broad', 'Sharp1','Sharp2','Sharp3']

rule r_extract_control_sharp_subsets_from_Atak2013:
    input:
        'inp/KalenderAtak2013/Global_Data_RNA-seq_TALL_Isoforms_ATAK.csv.bz2',
        'out/gtftk/get_5p_3p_coords/gtftk/convert_ensembl/gunzip/to-stdout/wget/https/ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_genomic.bed',
        'out/sort/_-u_-k1,1_-k2,2n/ln/updir/mw-tall/inp/Peaks_hg19_Mohammed/Sharp/NM_Sharp_Jurkat.bed',
        'out/sort/_-u_-k1,1_-k2,2n/ln/updir/mw-tall/inp/Peaks_hg19_Mohammed/Broad/NM_Broad_Jurkat.bed'
    output:
        bed=expand('out/r/extract_control_sharp_subsets_from_Atak2013/NM_{subset}{variants}.bed',
        subset=JURKAT_BROAD_AND_SHARP_SUBSETS,
        variants=[
        "",
        "_length_gt_1500b",
        "_length_gt_1500b_trim_500b_at_TSS_TES"])
    conda:
        "../../../../mw-lib/src/snakemake/envs/r_genomicranges.yaml"
    script:
        "../../r/script/extract_control_sharp_subsets_from_Atak2013.R"

rule r_extract_control_sharp_subsets_from_our_rnaseq:
    input:
        'out/subread/featureCounts_-O_gtf-GRCh37-ensembl-r87-protein-coding_bam-hg19-rna-jurkat-loucy.tsv',
        'out/sort/_-u_-k1,1_-k2,2n/ln/updir/mw-tall/inp/Peaks_hg19_Mohammed/Sharp/NM_Sharp_Jurkat.bed',
        'out/sort/_-u_-k1,1_-k2,2n/ln/updir/mw-tall/inp/Peaks_hg19_Mohammed/Broad/NM_Broad_Jurkat.bed'
    output:
        bed=expand('out/r/extract_control_sharp_subsets_from_our_rnaseq/NM_{feature}{subset}{variants}.bed',
        feature=["", "genes_associated_with_"],
        subset=JURKAT_BROAD_AND_SHARP_SUBSETS + ["Sharp"],
        variants=[
        "",
        "_trim_500b_at_TSS_TES"])
    conda:
        "../../../../mw-lib/src/snakemake/envs/r_genomicranges.yaml"
    script:
        "../../r/script/extract_control_sharp_subsets_from_our_rnaseq.R"


