rule meta_broad_tall_rev3_q4:
    """
    Aim:
        Gather desired heatmaps and profiles into bookown folder so they can be published with the report and linked from.
    """
    input:
        expand("out/deepTools/plotHeatmap_--colorList_blueCyanYellowOrangeRed_--heatmapWidth_6_--whatToShow_cph{averageTypeSummaryPlot}{legendLocation}/deepTools/computeMatrix_{mode}_--maxThreshold_1000_{bed}_{bw}.pdf",
            averageTypeSummaryPlot=[
                "",
                # With --maxThreshold_1000, mean profiles look smoother with no artefact
                #"_--averageTypeSummaryPlot_median",
                ],
            legendLocation = list(pandas.read_csv('../mw-belhocine2021/src/snakemake/tables/meta_broad_tall_rev3_q4_legendLocation.tsv', sep='\t',na_filter=False)['id']),
            mode = list(pandas.read_csv('../mw-belhocine2021/src/snakemake/tables/meta_broad_tall_rev3_q4_mode.tsv', sep='\t',na_filter=False)['id']),
            bed = list(pandas.read_csv('../mw-belhocine2021/src/snakemake/tables/meta_broad_tall_rev3_q4_bed.tsv', sep='\t',na_filter=False)['id']),
            bw = list(pandas.read_csv("../mw-belhocine2021/src/snakemake/tables/meta_broad_tall_rev3_q4_bw.tsv", sep='\t',na_filter=False)['id']))

rule meta_broad_tall_rev3_q4_single:
    input:
        expand("out/deepTools/plotHeatmap_--colorList_blueCyanYellowOrangeRed_--heatmapWidth_6_--whatToShow_cph{legendLocation}/deepTools/computeMatrix_{mode}_--maxThreshold_1000_{bed}/ln/alias/sst/all_samples/hg19/bw/{sample}.pdf",
            legendLocation = list(pandas.read_csv('../mw-belhocine2021/src/snakemake/tables/meta_broad_tall_rev3_q4_legendLocation.tsv', sep='\t',na_filter=False)['id']),
            mode = list(pandas.read_csv('../mw-belhocine2021/src/snakemake/tables/meta_broad_tall_rev3_q4_mode.tsv', sep='\t',na_filter=False)['id']),
            bed = list(pandas.read_csv('../mw-belhocine2021/src/snakemake/tables/meta_broad_tall_rev3_q4_bed.tsv', sep='\t',na_filter=False)['id']),
            sample = [x.strip() for x in open("../mw-belhocine2021/src/snakemake/lists/samples_broad_tall_rev3_q4.txt","r")])

rule meta_broad_tall_rev3_q4_spikein_single:
    input:
        expand("out/deepTools/plotHeatmap_--colorList_blueCyanYellowOrangeRed_--heatmapWidth_6_--whatToShow_cph{legendLocation}/deepTools/computeMatrix_{mode}_--maxThreshold_1000_{bed}/deepTools/spikeIn_bamCoverage_BDGP6/samtools/sam_to_bam_bai_-q_30/bowtie2/pe_hg19/sickle/pe_-t_sanger_-q_20/ln/alias/sst/all_samples/fastq/{sample}.pdf",
            legendLocation = list(pandas.read_csv('../mw-belhocine2021/src/snakemake/tables/meta_broad_tall_rev3_q4_legendLocation.tsv', sep='\t',na_filter=False)['id']),
            bed = list(pandas.read_csv('../mw-belhocine2021/src/snakemake/tables/meta_broad_tall_rev3_q4_bed.tsv', sep='\t',na_filter=False)['id']),
            mode = list(pandas.read_csv('../mw-belhocine2021/src/snakemake/tables/meta_broad_tall_rev3_q4_mode.tsv', sep='\t',na_filter=False)['id']),
            sample=[x.strip() for x in open("../mw-belhocine2021/src/snakemake/lists/meta_broad_tall_rev3_q4_single_spikein_samples.txt","r")])

rule meta_broad_tall_thymocytes_single:
    input:
        expand("out/deepTools/plotHeatmap_--colorList_blueCyanYellowOrangeRed_--heatmapWidth_6_--whatToShow_cph{legendLocation}/deepTools/computeMatrix_{mode}_--maxThreshold_1000_{bed}/ln/alias/sst/all_samples/hg19/bw/{sample}.pdf",
            legendLocation = list(pandas.read_csv('../mw-belhocine2021/src/snakemake/tables/meta_broad_tall_rev3_q4_legendLocation.tsv', sep='\t',na_filter=False)['id']),
            mode = list(pandas.read_csv('../mw-belhocine2021/src/snakemake/tables/meta_broad_tall_rev3_q4_mode.tsv', sep='\t',na_filter=False)['id']),
            bed = list(pandas.read_csv('../mw-belhocine2021/src/snakemake/tables/meta_broad_tall_rev3_q4_bed.tsv', sep='\t',na_filter=False)['id']),
            sample = [x.strip() for x in open("../mw-belhocine2021/src/snakemake/lists/samples_broad_tall_thymocytes.txt","r")])

rule meta_qc_insert_size:
    input:
        expand("out/picard/CollectInsertSizeMetrics/samtools/sam_to_bam_bai_-q_30/bowtie2/pe_{genome}/sickle/pe_-t_sanger_-q_20/ln/alias/sst/all_samples/fastq/{sample}.pdf",
            sample=[x.strip() for x in open("../mw-belhocine2021/src/snakemake/lists/meta_broad_tall_rev3_q4_single_spikein_samples.txt","r")],
            genome=['BDGP6','hg19'])
