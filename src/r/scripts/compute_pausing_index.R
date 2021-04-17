# This function installs and loads a package
loadLibrary <- function(package){
        if (!require(basename(package), character.only=TRUE)) BiocManager::install(package, update=FALSE)
    library(basename(package), character.only=TRUE)
}

packages <- c('BRGenomics',
              'ggpubr',
              'GenomicRanges'
              )

invisible(lapply(packages, loadLibrary))

smi <- 'out/deepTools/bamCoverage_--binSize_20_--minMappingQuality_0_--normalizeUsing_RPKM/samtools/index/samtools/sort/samtools/view_sam_to_bam_-q_30/bowtie2/se_hg19/sickle/se_-t_sanger_-q_20/ln/alias/sst/all_samples/fastq/Jurkat_SRR2157609_Pol2.bw'
Pol2 <- import(smi)

smi <- 'out/r/extract_control_sharp_subsets_from_Atak2013/NM_Broad.bed'
NM_Broad <- import(smi)

smi <- 'out/r/extract_control_sharp_subsets_from_Atak2013/NM_Sharp1.bed'
NM_Sharp1 <- import(smi)

smi <- 'out/r/extract_control_sharp_subsets_from_Atak2013/NM_Sharp2.bed'
NM_Sharp2 <- import(smi)

smi <- 'out/r/extract_control_sharp_subsets_from_Atak2013/NM_Sharp3.bed'
NM_Sharp3 <- import(smi)


#narrow(NM_Broad, start=500, end=-500)
#narrow(NM_Broad, start=250, end=-250)
#Error in .Call2("solve_user_SEW", refwidths, start, end, width, translate.negative.coord,  : 
#                solving row 163: the supplied start/end lead to a negative width

PI_Broad <- getPausingIndices(dataset.gr   = Pol2,
                              promoters.gr = promoters(NM_Broad, upstream=250, downstream=250),
                              genebodies.gr= narrow(NM_Broad, start=200, end=-200))

PI_Sharp1 <- getPausingIndices(dataset.gr   = Pol2,
                              promoters.gr = promoters(NM_Sharp1, upstream=250, downstream=250),
                              genebodies.gr= narrow(NM_Sharp1, start=200, end=-200))

PI_Sharp2 <- getPausingIndices(dataset.gr   = Pol2,
                              promoters.gr = promoters(NM_Sharp2, upstream=250, downstream=250),
                              genebodies.gr= narrow(NM_Sharp2, start=200, end=-200))

PI_Sharp3 <- getPausingIndices(dataset.gr   = Pol2,
                              promoters.gr = promoters(NM_Sharp3, upstream=250, downstream=250),
                              genebodies.gr= narrow(NM_Sharp3, start=200, end=-200))

d <- rbind(data.frame(PI=PI_Broad , n_gene=1:length(PI_Broad), peak_type='Broad', subset='Broad'),
           data.frame(PI=PI_Sharp1, n_gene=1:length(PI_Broad), peak_type='Sharp', subset='Sharp1'),
           data.frame(PI=PI_Sharp2, n_gene=1:length(PI_Broad), peak_type='Sharp', subset='Sharp2'),
           data.frame(PI=PI_Sharp3, n_gene=1:length(PI_Broad), peak_type='Sharp', subset='Sharp3'))


outdir="out/r/compute_pausing_index"
dir.create(outdir)

comparisons_list <- list(c('Broad','Sharp1'),
                         c('Broad','Sharp2'),
                         c('Broad','Sharp3'))

p <- ggplot(d, aes(x=subset, y=PI, color=peak_type))
p <- p + geom_violin()
p <- p + geom_boxplot(width=0.3, outlier.alpha=0.2)
p <- p + stat_summary(fun=mean, colour="darkred", geom="point", shape=18, size=3,show.legend = FALSE)
p <- p + stat_compare_means(label =  "p.signif",
                            label.x = 1.5,
                            method = 'wilcox',
                            method.args = list(alternative = "less"),
                            comparisons = comparisons_list)
p <- p + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
p <- p + theme_bw()

ggsave(filename=file.path(outdir,'pi_boxplot.pdf'), p)
ggsave(filename=file.path(outdir,'pi_boxplot.png'), p)
