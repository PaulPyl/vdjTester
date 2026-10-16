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
include { COMPARE_VDJ      } from './modules/local/compare'
include { EXPRESSION_CHECK } from './modules/local/compare'
include { CLONE_TRACKING   } from './modules/local/compare'

// Directory name a 10x tarball unpacks to.
def tarName(String url) {
    return url.tokenize('/').last().replaceAll(/\.tar(\.gz)?$/, '')
}

workflow {
    log.info """
    vdjTester
      input : ${params.samplesheet ?: "10x example ${params.sample_id}"}
      outdir: ${params.outdir}
    """.stripIndent()

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
        .combine(ch_refs)
    CELLRANGER_MULTI(ch_samples)

    ch_cr_vdj = CELLRANGER_MULTI.out.vdj_b.map { sid, dir -> tuple(sid, 'vdj_b', dir) }
        .mix(CELLRANGER_MULTI.out.vdj_t.map { sid, dir -> tuple(sid, 'vdj_t', dir) })
    IMMCANTATION(ch_cr_vdj)

    // MiXCR needs only the reads, so it starts as soon as a VDJ library is there.
    ch_mixcr_calls = channel.empty()
    if (!params.skip_mixcr) {
        MIXCR(
            ch_lib_items
                .filter { _sid, lib, _dir -> lib.feature_types.startsWith('VDJ') }
                .map { sid, lib, dir -> tuple(sid, lib.fastq_id, dir) }
        )
        ch_mixcr_calls = MIXCR.out.calls
    }

    // --- Per-sample comparison and expression check --------------------------------
    // Per VDJ library: one Cell Ranger file, one Immcantation file, one MiXCR file.
    def calls_per_vdj = params.skip_mixcr ? 2 : 3
    ch_calls = CELLRANGER_MULTI.out.calls
        .transpose()
        .mix(IMMCANTATION.out.calls, ch_mixcr_calls)
        .map { sid, f -> tuple(groupKey(sid, vdj_per_sample[sid] * calls_per_vdj), f) }
        .groupTuple()
        .map { key, files -> tuple(key.toString(), files) }
    COMPARE_VDJ(ch_calls)

    EXPRESSION_CHECK(
        CELLRANGER_MULTI.out.filtered_h5.join(COMPARE_VDJ.out.cell_classes),
        file(params.markers, checkIfExists: true),
    )

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
}
