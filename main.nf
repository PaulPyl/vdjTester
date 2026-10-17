#!/usr/bin/env nextflow

// vdjTester: 10x 5' samples through Cell Ranger multi, Immcantation and MiXCR, a
// cell-by-cell comparison of their VDJ calls, a check of those calls against B/T
// identity from the matched gene expression, and clone tracking across timepoints
// of the same patient. See docs/PLAN.md.
//
// Input, one of:
//   --samplesheet sheet.csv   sample_id,patient,timepoint,feature_types,fastq_id,runs
//                             (runs: ENA/SRA run accessions, one per lane, ';'-separated)
//   (default)                 the 10x 1k B-cell example, params.fastq_url / fastq_dir

include { FETCH_TAR        } from './modules/local/fetch_tar'
include { FETCH_ENA        } from './modules/local/fetch_ena'
include { CELLRANGER_MULTI } from './modules/local/cellranger_multi'
include { IMMCANTATION     } from './modules/local/immcantation'
include { MIXCR            } from './modules/local/mixcr'
include { MIXCR_IMGT       } from './modules/local/mixcr'
include { TRUST4           } from './modules/local/trust4'
include { FETCH_FILES      } from './modules/local/fetch_files'
include { COMPARE_VDJ      } from './modules/local/compare'
include { EXPRESSION_CHECK } from './modules/local/compare'
include { CLONE_TRACKING   } from './modules/local/compare'
include { COMPARE_CLONES   } from './modules/local/compare'
include { REPORT           } from './modules/local/compare'
include { FETCH_IMGT; MKVDJREF_IMGT } from './modules/local/imgt_reference'
include { FETCH_FILES as FETCH_MIXCR_LIBRARY } from './modules/local/fetch_files'

// Directory name a 10x tarball unpacks to.
def tarName(String url) {
    return url.tokenize('/').last().replaceAll(/\.tar(\.gz)?$/, '')
}

workflow {
    log.info """
    vdjTester
      input    : ${params.samplesheet ?: "10x example ${params.sample_id}"}
      reference: ${params.reference_set}
      outdir   : ${params.outdir}
    """.stripIndent()
    if (!(params.reference_set in ['default', 'imgt'])) {
        error "--reference_set must be 'default' or 'imgt', not '${params.reference_set}'"
    }
    def imgt = params.reference_set == 'imgt'

    // --- References and the example FASTQs: a given path, or download into store_dir
    ch_tars = channel.of(
        ['gex_ref', params.gex_reference, tarName(params.gex_ref_url), params.gex_ref_url],
        ['vdj_ref', params.vdj_reference, tarName(params.vdj_ref_url), params.vdj_ref_url],
    )
    if (!params.samplesheet) {
        ch_tars = ch_tars.mix(channel.of(['fastqs', params.fastq_dir, "${params.sample_id}_fastqs", params.fastq_url]))
    }
    ch_tars = ch_tars.branch { row ->
        given: row[1]
        fetch: true
    }
    FETCH_TAR(ch_tars.fetch.map { row -> tuple(row[0], row[2], row[3]) })
    ch_resolved = ch_tars.given
        .map { row -> tuple(row[0], file(row[1], checkIfExists: true)) }
        .mix(FETCH_TAR.out)

    ch_refs = ch_resolved
        .filter { key, _path -> key in ['gex_ref', 'vdj_ref'] }
        .collect(flat: false)
        .map { pairs ->
            def by_key = pairs.collectEntries { pair -> [(pair[0]): pair[1]] }
            tuple(by_key.gex_ref, by_key.vdj_ref)
        }

    // TRUST4's germline FASTAs: TRUST4's references, and the IMGT V genes used for
    // paralog-aware V agreement, so they are fetched even when TRUST4 is skipped.
    FETCH_FILES(channel.of(tuple(
        "trust4_ref_${params.trust4_ref_version}",
        ['hg38_bcrtcr.fa', 'human_IMGT+C.fa'].collect { f ->
            [f, "https://raw.githubusercontent.com/liulab-dfci/TRUST4/${params.trust4_ref_version}/${f.replace('+', '%2B')}"]
        },
    )))

    // --- Germline references: each tool's own, or IMGT through each tool's own route --
    // (modules/local/imgt_reference.nf). Immcantation and TRUST4 use IMGT germlines
    // either way; with 'imgt', Cell Ranger gets a fetch-imgt reference and MiXCR
    // MiLaboratories' IMGT library.
    ch_vdj_refs  = ch_refs
    ch_mixcr_lib = channel.empty()
    if (imgt) {
        FETCH_IMGT()
        MKVDJREF_IMGT(FETCH_IMGT.out)
        ch_vdj_refs = ch_refs.combine(MKVDJREF_IMGT.out).map { gex, _vdj, vdj_imgt -> tuple(gex, vdj_imgt) }
        if (!params.skip_mixcr) {
            def lib_file = params.mixcr_imgt_library.tokenize('/').last()
            FETCH_MIXCR_LIBRARY(channel.of(tuple("mixcr_library_${lib_file.replaceAll(/\.json(\.gz)?$/, '')}", [[lib_file, params.mixcr_imgt_library]])))
            ch_mixcr_lib = FETCH_MIXCR_LIBRARY.out.map { dir -> dir.resolve(lib_file) }.first()
        }
    }

    // --- Libraries: one item per library, (sample_id, library map, FASTQ dir) -------
    if (params.samplesheet) {
        ch_rows = channel.fromPath(params.samplesheet, checkIfExists: true)
            .splitCsv(header: true)
        ch_meta = ch_rows
            .map { row -> tuple(row.sample_id, [patient: row.patient, timepoint: row.timepoint as int]) }
            .unique()
        FETCH_ENA(ch_rows.map { row -> tuple(row.sample_id, row.fastq_id, row.runs.tokenize(';')) })
        // FETCH_ENA emits (sample_id, fastq_id, dir); attach feature_types by fastq_id.
        ch_types = ch_rows.map { row -> tuple(row.fastq_id, row.feature_types) }
        ch_lib_items = FETCH_ENA.out
            .map { sid, fid, dir -> tuple(fid, sid, dir) }
            .join(ch_types)
            .map { fid, sid, dir, ft -> tuple(sid, [fastq_id: fid, feature_types: ft, dir: dir.name], dir) }
    } else {
        ch_meta = channel.of(tuple(params.sample_id, [patient: params.sample_id, timepoint: 1]))
        ch_lib_items = ch_resolved
            .filter { key, _path -> key == 'fastqs' }
            .flatMap { _key, root ->
                params.libraries.collect { lib ->
                    tuple(params.sample_id, [fastq_id: lib.fastq_id, feature_types: lib.feature_types, dir: lib.subdir], root.resolve(lib.subdir))
                }
            }
    }

    // Libraries per sample, known up front, so each sample starts as soon as its own
    // libraries are in rather than after every download (groupTuple without a size
    // waits for the whole channel).
    def libs_per_sample = params.samplesheet
        ? file(params.samplesheet).splitCsv(header: true).countBy { row -> row.sample_id }
        : [(params.sample_id): params.libraries.size()]
    def vdj_per_sample = params.samplesheet
        ? file(params.samplesheet).splitCsv(header: true).findAll { row -> row.feature_types.startsWith('VDJ') }.countBy { row -> row.sample_id }
        : [(params.sample_id): params.libraries.count { lib -> lib.feature_types.startsWith('VDJ') }]

    // --- Cell Ranger multi, one task per sample ------------------------------------
    ch_samples = ch_lib_items
        .map { sid, lib, dir -> tuple(groupKey(sid, libs_per_sample[sid]), lib, dir) }
        .groupTuple()
        .map { key, libs, dirs ->
            // Fixed library order (VDJ first, then GEX; by name), so the Cell Ranger
            // config, and so -resume, do not depend on which download finished first.
            def sorted = [libs, dirs].transpose().sort { pair ->
                [pair[0].feature_types.startsWith('VDJ') ? 0 : 1, pair[0].fastq_id]
            }
            // ArrayBag, as groupTuple returns: -resume hashes bags order-independently,
            // so tasks from runs before this sorting still match.
            tuple(key.toString(), new nextflow.util.ArrayBag(sorted.collect { pair -> pair[0] }), new nextflow.util.ArrayBag(sorted.collect { pair -> pair[1] }))
        }
        .combine(ch_vdj_refs)
    CELLRANGER_MULTI(ch_samples)

    ch_cr_vdj = CELLRANGER_MULTI.out.vdj_b.map { sid, dir -> tuple(sid, 'vdj_b', dir) }
        .mix(CELLRANGER_MULTI.out.vdj_t.map { sid, dir -> tuple(sid, 'vdj_t', dir) })
    IMMCANTATION(ch_cr_vdj)

    // MiXCR needs only the reads, so it starts as soon as a VDJ library is there.
    ch_mixcr_calls = channel.empty()
    if (!params.skip_mixcr) {
        ch_mixcr_in = ch_lib_items
            .filter { _sid, lib, _dir -> lib.feature_types.startsWith('VDJ') }
            .map { sid, lib, dir -> tuple(sid, lib.fastq_id, dir) }
        if (imgt) {
            MIXCR_IMGT(ch_mixcr_in, ch_mixcr_lib)
            ch_mixcr_calls = MIXCR_IMGT.out.calls
        } else {
            MIXCR(ch_mixcr_in)
            ch_mixcr_calls = MIXCR.out.calls
        }
    }

    // TRUST4: a second assembler from raw reads.
    ch_trust4_calls = channel.empty()
    if (!params.skip_trust4) {
        TRUST4(
            ch_lib_items
                .filter { _sid, lib, _dir ->
                    lib.feature_types.startsWith('VDJ') || (params.trust4_gex && lib.feature_types == 'Gene Expression')
                }
                .map { sid, lib, dir -> tuple(sid, lib.fastq_id, dir) }
                .combine(FETCH_FILES.out)
        )
        ch_trust4_calls = TRUST4.out.calls
    }

    // --- Per-sample comparison and expression check --------------------------------
    // Per VDJ library: one call file each from Cell Ranger, Immcantation, and MiXCR
    // and TRUST4 unless skipped.
    def calls_per_vdj = 2 + (params.skip_mixcr ? 0 : 1) + (params.skip_trust4 ? 0 : 1)
    def gex_per_sample = params.samplesheet
        ? file(params.samplesheet).splitCsv(header: true).findAll { row -> row.feature_types == 'Gene Expression' }.countBy { row -> row.sample_id }
        : [(params.sample_id): params.libraries.count { lib -> lib.feature_types == 'Gene Expression' }]
    def gex_calls = !params.skip_trust4 && params.trust4_gex
    def gex_ids = params.samplesheet
        ? file(params.samplesheet).splitCsv(header: true).findAll { row -> row.feature_types == 'Gene Expression' }.collect { row -> row.fastq_id }
        : params.libraries.findAll { lib -> lib.feature_types == 'Gene Expression' }.collect { lib -> lib.fastq_id }
    ch_calls = CELLRANGER_MULTI.out.calls
        .transpose()
        .mix(IMMCANTATION.out.calls, ch_mixcr_calls, ch_trust4_calls)
        .map { sid, f -> tuple(groupKey(sid, vdj_per_sample[sid] * calls_per_vdj + (gex_calls ? (gex_per_sample[sid] ?: 0) : 0)), f) }
        .groupTuple()
        .map { key, files -> tuple(key.toString(), files) }
    // The GEX matrix gives the cell barcodes TRUST4 calls are restricted to; the VDJ
    // reference gives germline V sequences for paralog-aware agreement.
    COMPARE_VDJ(
        ch_calls.join(CELLRANGER_MULTI.out.filtered_h5).map { sid, calls, h5 -> tuple(sid, calls, h5, gex_ids) },
        ch_vdj_refs.map { _gex, vdj -> vdj }.first(),
        FETCH_FILES.out.first(),   // value channel: reused by every sample
    )

    EXPRESSION_CHECK(
        CELLRANGER_MULTI.out.filtered_h5.join(COMPARE_VDJ.out.cell_classes),
        file(params.markers, checkIfExists: true),
    )

    // --- Clone-level comparison of the methods, per sample --------------------------
    COMPARE_CLONES(COMPARE_VDJ.out.calls_long)

    // --- Clone tracking across a patient's timepoints ------------------------------
    def samples_per_patient = params.samplesheet
        ? file(params.samplesheet).splitCsv(header: true).collect { row -> [row.patient, row.sample_id] }.unique().countBy { pair -> pair[0] }
        : [(params.sample_id): 1]
    ch_tracking = COMPARE_VDJ.out.calls_long
        .join(ch_meta)
        .map { sid, calls, meta -> tuple(groupKey(meta.patient, samples_per_patient[meta.patient]), sid, meta.timepoint, calls) }
        .groupTuple()
        .map { key, sids, tps, calls -> tuple(key.toString(), sids, tps, calls) }
    CLONE_TRACKING(ch_tracking)

    // --- HTML report ------------------------------------------------------------------
    ch_samples_tsv = ch_meta
        .map { sid, meta -> "${sid}\t${meta.patient}\t${meta.timepoint}" }
        .collectFile(name: 'samples.tsv', newLine: true, sort: true, seed: "sample_id\tpatient\ttimepoint")
    // Each method's germline reference, for the run overview.
    ch_ref_info = imgt
        ? MKVDJREF_IMGT.out.map { ref ->
            def date = ref.resolve('imgt_downloaded.txt')
            "IMGT via fetch-imgt + mkvdjref, downloaded ${date.exists() ? date.text.trim() : 'on an unrecorded date'}"
        }
        : channel.value("${tarName(params.vdj_ref_url)} (Ensembl-based)")
    ch_run_info = CELLRANGER_MULTI.out.versions.first()
        .combine(ch_ref_info)
        .map { f, ref_info ->
            [
                ['run_name', workflow.runName], ['date', new Date().format('yyyy-MM-dd')],
                ['nextflow', nextflow.version.toString()], ['profile', workflow.profile],
                ['input', params.samplesheet ?: "10x example ${params.sample_id}"],
                ['reference set', params.reference_set],
                ['germlines: Cell Ranger', ref_info],
                ['germlines: MiXCR', params.skip_mixcr ? 'skipped' : (imgt ? params.mixcr_imgt_library.tokenize('/').last() : 'built-in library')],
                ['germlines: Immcantation', "IMGT germlines in ${params.immcantation_container}"],
                ['germlines: TRUST4', "IMGT file of TRUST4 ${params.trust4_ref_version}"],
                ['cellranger', f.text.trim()],
                ['immcantation', params.immcantation_container],
                ['mixcr', params.skip_mixcr ? 'skipped' : params.mixcr_container],
                ['trust4', params.skip_trust4 ? 'skipped' : "${params.trust4_container} (germlines ${params.trust4_ref_version}; on GEX: ${params.trust4_gex})"],
                ['trust4 filter', "min UMIs ${params.trust4_min_umis} (GEX ${params.trust4_gex_min_umis}), min fraction ${params.trust4_min_frac}"],
                ['expression scores', "${params.score_method}, min markers ${params.min_markers}"],
                ['python container', params.python_container],
            ].collect { kv -> "${kv[0]}\t${kv[1]}" }.join('\n')
        }
        .collectFile(name: 'run_info.tsv', newLine: true, seed: "key\tvalue")
    ch_metrics = CELLRANGER_MULTI.out.metrics
        .collectFile { sid, f -> ["${sid}.metrics_summary.csv", f.text] }
    ch_report_data = COMPARE_VDJ.out.report.flatten()
        .mix(EXPRESSION_CHECK.out.report.flatten(), CLONE_TRACKING.out.report.flatten(), COMPARE_CLONES.out.report.flatten(), ch_metrics, ch_samples_tsv, ch_run_info)
        .collect()
    REPORT(ch_report_data, file("${projectDir}/assets/report/vdjtester_report.qmd", checkIfExists: true))
}
