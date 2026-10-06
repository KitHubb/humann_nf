nextflow.enable.dsl = 2

/*
 * Minimal HUMAnN 3 workflow for paired, host-filtered shotgun reads.
 * R1/R2 are concatenated only inside each Nextflow task work directory.
 */

process HUMANN {
    tag "${meta.id}"
    label 'humann'
    maxForks params.max_parallel_samples as int

    input:
    tuple val(meta), path(read1), path(read2), path(tax_profile)

    output:
    tuple val(meta), path("${meta.id}_genefamilies.${params.output_format}"), emit: genefamilies
    tuple val(meta), path("${meta.id}_pathabundance.${params.output_format}"), emit: pathabundance
    tuple val(meta), path("${meta.id}_pathcoverage.${params.output_format}"), emit: pathcoverage
    tuple val(meta), path("${meta.id}.humann.log"), emit: logs

    script:
    def removeTemp = params.remove_temp_output ? '--remove-temp-output' : ''
    def removeStratified = params.remove_stratified_output ? '--remove-stratified-output' : ''
    def removeDescription = params.remove_column_description_output ? '--remove-column-description-output' : ''
    def bypassNucleotideIndex = params.bypass_nucleotide_index ? '--bypass-nucleotide-index' : ''
    def bypassNucleotide = params.bypass_nucleotide_search ? '--bypass-nucleotide-search' : ''
    def bypassTranslated = params.bypass_translated_search ? '--bypass-translated-search' : ''
    def translatedIdentity = params.translated_identity_threshold == null ? '' : "--translated-identity-threshold ${params.translated_identity_threshold}"
    def pathwaysDatabase = params.pathways_database ? "--pathways-database ${params.pathways_database}" : ''
    def idMapping = params.id_mapping ? "--id-mapping ${params.id_mapping}" : ''

    """
    cat ${read1} ${read2} > ${meta.id}.fastq.gz

    humann \
        --input ${meta.id}.fastq.gz \
        --input-format fastq.gz \
        --output . \
        --output-basename ${meta.id} \
        --taxonomic-profile ${tax_profile} \
        --nucleotide-database ${params.chocophlan} \
        --protein-database ${params.uniref} \
        --threads ${task.cpus} \
        --memory-use ${params.memory_use} \
        --search-mode ${params.search_mode} \
        --prescreen-threshold ${params.prescreen_threshold} \
        --bowtie-options='${params.bowtie_options}' \
        --nucleotide-identity-threshold ${params.nucleotide_identity_threshold} \
        --nucleotide-query-coverage-threshold ${params.nucleotide_query_coverage_threshold} \
        --nucleotide-subject-coverage-threshold ${params.nucleotide_subject_coverage_threshold} \
        --translated-alignment ${params.translated_alignment} \
        --diamond-options='${params.diamond_options}' \
        --evalue ${params.evalue} \
        ${translatedIdentity} \
        --translated-query-coverage-threshold ${params.translated_query_coverage_threshold} \
        --translated-subject-coverage-threshold ${params.translated_subject_coverage_threshold} \
        --gap-fill ${params.gap_fill} \
        --minpath ${params.minpath} \
        --pathways ${params.pathways} \
        --xipe ${params.xipe} \
        --log-level ${params.log_level} \
        --o-log ${meta.id}.humann.log \
        --output-format ${params.output_format} \
        --output-max-decimals ${params.output_max_decimals} \
        ${pathwaysDatabase} \
        ${idMapping} \
        ${bypassNucleotideIndex} \
        ${bypassNucleotide} \
        ${bypassTranslated} \
        ${removeTemp} \
        ${removeStratified} \
        ${removeDescription} \
        ${params.humann_extra_args}
    """
}

process JOIN_TABLES {
    tag "${table_type}"
    label 'utility'

    input:
    tuple val(table_type), path(tables)

    output:
    tuple val(table_type), path("humann_${table_type}.tsv"), emit: joined

    script:
    """
    mkdir input_tables
    cp ${tables} input_tables/
    humann_join_tables \
        --input input_tables \
        --file_name ${table_type} \
        --output humann_${table_type}.tsv
    """
}

process RENORM_TABLE {
    tag "${table_type}:${params.renorm_units}"
    label 'utility'

    input:
    tuple val(table_type), path(table)

    output:
    tuple val(table_type), path("humann_${table_type}_${params.renorm_units}.tsv"), emit: normalized

    when:
    params.run_renorm

    script:
    """
    humann_renorm_table \
        --input ${table} \
        --units ${params.renorm_units} \
        --mode ${params.renorm_mode} \
        --special ${params.renorm_special} \
        --update-snames \
        --output humann_${table_type}_${params.renorm_units}.tsv
    """
}

process SPLIT_STRATIFIED {
    tag "${table_type}"
    label 'utility'

    input:
    tuple val(table_type), path(table)

    output:
    path "split_${table_type}/*", emit: split

    when:
    params.run_split_stratified && !params.remove_stratified_output

    script:
    """
    mkdir split_${table_type}
    humann_split_stratified_table \
        --input ${table} \
        --output split_${table_type}
    """
}

workflow {
    if (!params.input) {
        error "--input samplesheet.csv is required"
    }

    reads = Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            if (!row.sample || !row.fastq_1 || !row.fastq_2 || !row.taxonomic_profile) {
                error "Samplesheet requires: sample,fastq_1,fastq_2,taxonomic_profile"
            }
            def meta = [id: row.sample.toString()]
            def read1 = file(row.fastq_1, checkIfExists: true)
            def read2 = file(row.fastq_2, checkIfExists: true)
            def profile = file(row.taxonomic_profile, checkIfExists: true)
            tuple(meta, read1, read2, profile)
        }

    HUMANN(reads)

    if (params.run_join) {
        genes = HUMANN.out.genefamilies.map { meta, table -> tuple('genefamilies', table) }.groupTuple()
        abundance = HUMANN.out.pathabundance.map { meta, table -> tuple('pathabundance', table) }.groupTuple()
        coverage = HUMANN.out.pathcoverage.map { meta, table -> tuple('pathcoverage', table) }.groupTuple()
        join_inputs = genes.mix(abundance, coverage)
        JOIN_TABLES(join_inputs)
        // Pathway coverage is a 0-1 completeness score and must not be CPM/relative-abundance normalized.
        abundance_tables = JOIN_TABLES.out.joined.filter { table_type, table -> table_type != 'pathcoverage' }
        coverage_tables = JOIN_TABLES.out.joined.filter { table_type, table -> table_type == 'pathcoverage' }
        RENORM_TABLE(abundance_tables)
        split_source = params.run_renorm ? RENORM_TABLE.out.normalized.mix(coverage_tables) : JOIN_TABLES.out.joined
        SPLIT_STRATIFIED(split_source)
    }
}
