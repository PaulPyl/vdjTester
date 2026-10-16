#!/usr/bin/env nextflow

// vdjTester v0.1: one 10x sample through Cell Ranger multi, Immcantation and MiXCR,
// a cell-by-cell comparison of their VDJ calls, and a check of those calls against
// B/T identity from the matched gene expression. See docs/PLAN.md.

include { FETCH_TAR        } from './modules/local/fetch_tar'
include { CELLRANGER_MULTI } from './modules/local/cellranger_multi'
include { IMMCANTATION     } from './modules/local/immcantation'
include { MIXCR            } from './modules/local/mixcr'
include { COMPARE_VDJ      } from './modules/local/compare'
include { EXPRESSION_CHECK } from './modules/local/compare'

// Directory name a 10x reference tarball unpacks to.
def tarName(String url) {
    return url.tokenize('/').last().replaceAll(/\.tar(\.gz)?$/, '')
}

workflow {
    log.info """
    vdjTester v0.1
      sample   : ${params.sample_id}
      libraries: ${params.libraries.collect { lib -> lib.feature_types }.join(', ')}
      outdir   : ${params.outdir}
    """.stripIndent()

    // --- Inputs: use a given path, otherwise download into store_dir --------------
    ch_inputs = channel.of(
        ['fastqs',  params.fastq_dir,     "${params.sample_id}_fastqs", params.fastq_url],
        ['gex_ref', params.gex_reference, tarName(params.gex_ref_url),   params.gex_ref_url],
        ['vdj_ref', params.vdj_reference, tarName(params.vdj_ref_url),   params.vdj_ref_url],
    ).branch { row ->
        given: row[1]
        fetch: true
    }

    FETCH_TAR(ch_inputs.fetch.map { row -> tuple(row[0], row[2], row[3]) })

    ch_cr_input = ch_inputs.given
        .map { row -> tuple(row[0], file(row[1], checkIfExists: true)) }
        .mix(FETCH_TAR.out)
        .collect(flat: false)
        .map { pairs ->
            def by_key = pairs.collectEntries { pair -> [(pair[0]): pair[1]] }
            tuple(params.sample_id, by_key.fastqs, by_key.gex_ref, by_key.vdj_ref)
        }

    // --- Callers -------------------------------------------------------------------
    CELLRANGER_MULTI(ch_cr_input, params.libraries)

    ch_cr_vdj = CELLRANGER_MULTI.out.vdj_b.map { sample, dir -> tuple(sample, 'vdj_b', dir) }
        .mix(CELLRANGER_MULTI.out.vdj_t.map { sample, dir -> tuple(sample, 'vdj_t', dir) })
    IMMCANTATION(ch_cr_vdj)

    ch_mixcr_calls = channel.empty()
    if (!params.skip_mixcr) {
        def vdj_libraries = params.libraries.findAll { lib -> lib.feature_types.startsWith('VDJ') }
        ch_mixcr_input = ch_cr_input.flatMap { sample, fastqs, _gex_ref, _vdj_ref ->
            vdj_libraries.collect { lib -> tuple(sample, lib.fastq_id, fastqs, lib.subdir) }
        }
        MIXCR(ch_mixcr_input)
        ch_mixcr_calls = MIXCR.out.calls
    }

    // --- Comparison and expression check ---------------------------------------------
    ch_calls = CELLRANGER_MULTI.out.calls.flatten()
        .mix(IMMCANTATION.out.calls, ch_mixcr_calls)
        .collect()
    COMPARE_VDJ(params.sample_id, ch_calls)

    EXPRESSION_CHECK(
        CELLRANGER_MULTI.out.filtered_h5,
        COMPARE_VDJ.out.cell_classes,
        file(params.markers, checkIfExists: true),
    )
}
