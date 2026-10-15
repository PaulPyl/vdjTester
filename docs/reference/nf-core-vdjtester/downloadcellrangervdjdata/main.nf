// TODO nf-core: If in doubt look at other nf-core/modules to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/modules/nf-core/
//               You can also ask for help via your pull request or on the #modules channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A module file SHOULD only define input and output files as command-line parameters.
//               All other parameters MUST be provided using the "task.ext" directive, see here:
//               https://www.nextflow.io/docs/latest/process.html#ext
//               where "task.ext" is a string.
//               Any parameters that need to be evaluated in the context of a particular sample
//               e.g. single-end/paired-end data MUST also be defined and evaluated appropriately.
// TODO nf-core: Software that can be piped together SHOULD be added to separate module files
//               unless there is a run-time, storage advantage in implementing this way
//               e.g. it's ok to have a single module for bwa to output BAM instead of SAM:
//                 bwa mem | samtools view -B -T ref.fasta
// TODO nf-core: Optional inputs are not currently supported by Nextflow. However, using an empty
//               list (`[]`) instead of a file can be used to work around this issue.

process DOWNLOADCELLRANGERVDJDATA {
    tag '$meta.project_name'
    label 'process_low'

    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bash':
        'biocontainers/biocontainers:latest' }"

    input:
    // TODO nf-core: Where applicable all sample-specific information e.g. "id", "single_end", "read_group"
    //               MUST be provided as an input via a Groovy Map called "meta".
    //               This information may not be required in some instances e.g. indexing reference genome files:
    //               https://github.com/nf-core/modules/blob/master/modules/nf-core/bwa/index/main.nf
    // TODO nf-core: Where applicable please provide/convert compressed files as input/output
    //               e.g. "*.fastq.gz" and NOT "*.fastq", "*.bam" and NOT "*.sam" etc.
    val meta

    output:
    // TODO nf-core: Named file extensions MUST be emitted for ALL output channels
    path "${meta.project_name}/sc5p_v2_hs_B_1k_multi_5gex_b_fastqs/*.fastq.gz", emit: fastqs
    path "${meta.project_name}/refdata-cellranger-vdj-GRCh38-alts-ensembl-${meta.ref_version}-5.0.0/*.{fa,fasta,json}", emit: vdj_reference
    // TODO nf-core: List additional required output channels/values here
    // TODO nf-core: Update the command here to obtain the version number of the software used in this module
    // TODO nf-core: If multiple software packages are used in this module, all MUST be added here
    //               by copying the line below and replacing the current tool with the extra tool(s)
    tuple val("${task.process}"), val('downloadcellrangervdjdata'), eval("echo 'Cell Ranger V(D)J Reference Data'"), topic: versions, emit: versions_downloadcellrangervdjdata

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    set -e  # Exit on error
    
    echo "=== nf-core-vdjtester Example Data Download ==="
    
    cd ${meta.project_name}
    downloadcellrangervdjdata "${meta.ref_version}"
    
    echo ""
    echo "Download complete!"
    ls -laR .
    """

    stub:
    """
    echo "=== nf-core-vdjtester Example Data Download (Stub Mode) ==="
    mkdir -p ${prefix}/dataset-vdj-practice
    cd ${prefix}/dataset-vdj-practice
    curl --fail -L \
        https://cf.10xgenomics.com/samples/cell-vdj/6.0.0/sc5p_v2_hs_B_1k_multi_5gex_b_Multiplex/sc5p_v2_hs_B_1k_multi_5gex_b_Multiplex_fastqs.tar \
        -o sc5p_v2_hs_B_1k_multi_5gex_b_Multiplex_fastqs.tar
    tar -xf sc5p_v2_hs_B_1k_multi_5gex_b_Multiplex_fastqs.tar
    ls -laR .
    """
}