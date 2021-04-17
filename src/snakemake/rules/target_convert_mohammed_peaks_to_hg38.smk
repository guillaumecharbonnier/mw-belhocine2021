rule target_convert_mohammed_peaks_to_hg38:
    """
    Note:
        This target was changed on 2020-06-07 10:31:40
        to add the "sort/_-u" because I identified later
        that Mohammed peaks were duplicated.
    """
    input:
        expand("out/bedtools/merge_-d_5/sort/_-k1,1_-k2,2n/crossmap/chain-hg19-to-hg38/sort/_-u/ln/updir/mw-tall/inp/Peaks_hg19_Mohammed/{files}", files=[x.strip() for x in open("../mw-belhocine2021/src/snakemake/lists/mohammed_peaks.txt","r")])

rule tmp_debug_crossmap:
    "out/crossmap/chain-hg19-to-hg38/ln/updir/mw-tall/inp/Peaks_hg19_Mohammed/Broad/NM_Broad_CD4_TH91_H3K4me3_S010R5H1.bed"
    input:
        expand("out/sort/_-u/ln/updir/mw-tall/inp/Peaks_hg19_Mohammed/{files}", files=[x.strip() for x in open("../mw-belhocine2021/src/snakemake/lists/mohammed_peaks.txt","r")])
