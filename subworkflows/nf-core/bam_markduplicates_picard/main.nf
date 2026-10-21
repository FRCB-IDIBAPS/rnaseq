//
// Picard MarkDuplicates, index BAM file and run samtools stats, flagstat, idxstats and collectrnaseqmetrics
//

include { PICARD_MARKDUPLICATES       } from '../../../modules/nf-core/picard/markduplicates/main'
include { PICARD_COLLECTRNASEQMETRICS } from '../../../modules/nf-core/picard/collectrnaseqmetrics/main'
include { SAMTOOLS_INDEX              } from '../../../modules/nf-core/samtools/index/main'
include { BAM_STATS_SAMTOOLS          } from '../bam_stats_samtools/main'

workflow BAM_MARKDUPLICATES_PICARD {
    take:
    ch_reads          // channel: [ val(meta), path(reads) ]
    ch_fasta_fai // channel: [ val(meta), path(fasta), path(fai)]
    ch_flatfile       // channel [ path(flatfile) ]
    ch_rrna_intervals // channel: [ path(rrna_intervals) ]
    main:
    PICARD_MARKDUPLICATES(ch_reads, ch_fasta_fai)

    ch_markdup = PICARD_MARKDUPLICATES.out.bam.mix(PICARD_MARKDUPLICATES.out.cram)

    SAMTOOLS_INDEX(ch_markdup)

    ch_reads_index = ch_markdup.join(SAMTOOLS_INDEX.out.index, by: [0])

    BAM_STATS_SAMTOOLS(ch_reads_index, ch_fasta_fai)

    ch_per_sample_mqc_bundle = BAM_STATS_SAMTOOLS.out.stats
        .join(BAM_STATS_SAMTOOLS.out.flagstat,   remainder: true)
        .join(BAM_STATS_SAMTOOLS.out.idxstats,   remainder: true)
        .join(PICARD_MARKDUPLICATES.out.metrics, remainder: true)
        .map { row -> [row[0], row.drop(1).findAll { f -> f != null }.collectMany { e -> (e instanceof List) ? e : [e] }] }

    PICARD_COLLECTRNASEQMETRICS ( ch_markdup, ch_fasta_fai.map { meta, fasta, _fai -> [meta, fasta] }, ch_flatfile, ch_rrna_intervals)
    emit:
    bam                   = PICARD_MARKDUPLICATES.out.bam // channel: [ val(meta), path(bam) ]
    cram                  = PICARD_MARKDUPLICATES.out.cram // channel: [ val(meta), path(cram) ]
    metrics               = PICARD_MARKDUPLICATES.out.metrics // channel: [ val(meta), path(metrics) ]
    rnaseqmetrics         = PICARD_COLLECTRNASEQMETRICS.out.metrics // channel: [ val(meta), path(metrics) ]
    index                 = SAMTOOLS_INDEX.out.index // channel: [ val(meta), path(index) ]
    stats                 = BAM_STATS_SAMTOOLS.out.stats // channel: [ val(meta), path(stats) ]
    flagstat              = BAM_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), path(flagstat) ]
    idxstats              = BAM_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), path(idxstats) ]
    per_sample_mqc_bundle = ch_per_sample_mqc_bundle // channel: [ val(meta), list(files) ]
}
