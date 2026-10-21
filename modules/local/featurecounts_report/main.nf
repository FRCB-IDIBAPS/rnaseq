process FEATURECOUNTS_REPORT {
    label "process_medium"

    conda "${moduleDir}/environment.yml"
    container "featurecounts_report"

    input:
    path counts, stageAs: 'input/*'

    output:
    path "merged_counts.tsv" , emit: merged_counts
    path "qc_metrics.tsv"    , emit: qc_metrics
    path "pca_scores.tsv"    , emit: pca_scores
    path "qc_report.pdf"     , emit: pdf
    path "*.png"             , emit: png
    tuple val("${task.process}"), val('r-base'), eval("Rscript -e 'cat(as.character(getRversion()))'"), emit: versions_r_base, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    featurecounts_report.R input . $args
    """

    stub:
    """
    touch merged_counts.tsv
    touch qc_metrics.tsv
    touch pca_scores.tsv
    touch qc_report.pdf
    touch library_size.png
    """
}
