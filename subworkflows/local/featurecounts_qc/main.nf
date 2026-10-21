//
// Gene-level (gene_id-grouped) featureCounts + custom QC report.
//
// Independent of the biotype-QC featureCounts call in bam_qc_rnaseq (which
// groups by gene_biotype/gene_type purely to feed the MultiQC biotype plot).
// This subworkflow exists to feed a bespoke downstream QC report and is not
// gated by --skip_qc/--use_rustqc/--skip_biotype_qc.
//

include { SUBREAD_FEATURECOUNTS as SUBREAD_FEATURECOUNTS_GENE_ID } from '../../../modules/nf-core/subread/featurecounts/main'
include { FEATURECOUNTS_REPORT                                   } from '../../../modules/local/featurecounts_report/main'

workflow FEATURECOUNTS_QC {

    take:
    ch_genome_bam   // channel: [ val(meta), path(bam) ]
    ch_gtf          // channel: path(gtf)

    main:
    SUBREAD_FEATURECOUNTS_GENE_ID (
        ch_genome_bam.combine(ch_gtf)
    )

    FEATURECOUNTS_REPORT (
        SUBREAD_FEATURECOUNTS_GENE_ID.out.counts
            .map { _meta, tsv -> tsv }
            .collect()
    )

    emit:
    counts         = SUBREAD_FEATURECOUNTS_GENE_ID.out.counts  // channel: [ val(meta), path(tsv) ]
    summary        = SUBREAD_FEATURECOUNTS_GENE_ID.out.summary // channel: [ val(meta), path(tsv) ]
    merged_counts  = FEATURECOUNTS_REPORT.out.merged_counts    // channel: path(tsv)
    qc_metrics     = FEATURECOUNTS_REPORT.out.qc_metrics       // channel: path(tsv)
    pca_scores     = FEATURECOUNTS_REPORT.out.pca_scores       // channel: path(tsv)
    report_pdf     = FEATURECOUNTS_REPORT.out.pdf              // channel: path(pdf)
    report_png     = FEATURECOUNTS_REPORT.out.png              // channel: path(png)
}
