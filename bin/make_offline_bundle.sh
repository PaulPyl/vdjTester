#!/usr/bin/env bash
# Collect everything vdjTester needs on an air-gapped cluster (e.g. UPPMAX Bianca) into
# one directory with a checksum manifest. Run on a machine with internet, from the
# repository. See docs/OFFLINE.md.
#
#   bin/make_offline_bundle.sh --dry-run            # list everything, with sizes
#   bin/make_offline_bundle.sh /path/to/bundle      # collect it
#
# The bundle:
#   pipeline/vdjTester-<commit>.tar.gz   the pipeline code at --commit (git archive)
#   nextflow/nextflow-<version>-dist     Nextflow, self-contained (needs Java 17+)
#   nextflow/plugins/                    plugins given with --plugins (the pipeline needs none)
#   containers/*.img                     Apptainer images, under the names Nextflow's cache uses
#   store/<name>.tar[.gz]                one archive per store_dir entry: 10x references,
#                                        TRUST4 germlines, IMGT reference, MiXCR IMGT library
#   tools/                               Cell Ranger, only if given with --cellranger-tar
#   unpack.sh, README.txt, bundle_info.tsv, MANIFEST.sha256
#
# References and images are copied from an existing store_dir and image cache where
# they are (so the bundle has exactly what was used here, e.g. the IMGT release that
# fetch-imgt downloaded), and downloaded or pulled otherwise.
#
# Written for bash 3.2 (macOS) as well as Linux.
set -euo pipefail

repo=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)

usage() {
    cat <<EOF
Usage: $(basename "$0") [options] [OUTDIR]

  -n, --dry-run              list what would be collected, with sizes; write nothing
  --commit REF               pipeline commit to bundle (default: HEAD)
  --store DIR                existing store_dir to take references from
                             (default: $repo/data/store)
  --container-cache DIR      existing Nextflow image cache to take images from
                             (default: $repo/container_cache/nextflow)
  --reference-set SET        default, imgt or both (default: both)
  --skip-mixcr               leave out MiXCR's image and IMGT library
  --nextflow-version V       Nextflow to bundle (default: the local \`nextflow -v\`)
  --plugins 'a@1.0 b@2.0'    Nextflow plugins to bundle (default: none; the pipeline uses none)
  --cellranger-tar FILE      Cell Ranger tarball you downloaded from 10x (not redistributable;
                             see README.txt in the bundle)
  --arch ARCH                image architecture to pull (default: amd64)
  -h, --help
EOF
}

dry_run=0; commit=HEAD; store="$repo/data/store"; cache="$repo/container_cache/nextflow"
reference_set=both; skip_mixcr=0; nf_version=''; plugins=''; cellranger_tar=''; arch=amd64; outdir=''
while [ $# -gt 0 ]; do
    case $1 in
        -n|--dry-run)        dry_run=1 ;;
        --commit)            commit=$2; shift ;;
        --store)             store=$2; shift ;;
        --container-cache)   cache=$2; shift ;;
        --reference-set)     reference_set=$2; shift ;;
        --skip-mixcr)        skip_mixcr=1 ;;
        --nextflow-version)  nf_version=$2; shift ;;
        --plugins)           plugins=$2; shift ;;
        --cellranger-tar)    cellranger_tar=$2; shift ;;
        --arch)              arch=$2; shift ;;
        -h|--help)           usage; exit 0 ;;
        -*)                  echo "unknown option: $1" >&2; usage >&2; exit 2 ;;
        *)                   outdir=$1 ;;
    esac
    shift
done
case $reference_set in default|imgt|both) ;; *) echo "--reference-set: default, imgt or both" >&2; exit 2 ;; esac
if [ $dry_run -eq 0 ] && [ -z "$outdir" ]; then
    echo "OUTDIR is required (or --dry-run)" >&2; usage >&2; exit 2
fi
if [ -n "$cellranger_tar" ] && [ ! -f "$cellranger_tar" ]; then
    echo "--cellranger-tar: no such file: $cellranger_tar" >&2; exit 2
fi

command -v nextflow >/dev/null || { echo "nextflow is needed on this machine (to read the pipeline's settings)" >&2; exit 1; }
sha=$(git -C "$repo" rev-parse --verify "$commit^{commit}")
short=${sha:0:10}
# `nextflow -v` adds a build number (26.04.6.12646); releases are tagged without it.
[ -n "$nf_version" ] || nf_version=$(nextflow -v | awk '{print $3}' | cut -d. -f1-3)

tmp=$(mktemp -d "${TMPDIR:-/tmp}/vdjbundle.XXXXXX")
trap 'rm -rf "$tmp"' EXIT

# --- The pipeline's settings at that commit -------------------------------------------
git -C "$repo" archive --format=tar --prefix=vdjTester/ "$sha" | tar -x -C "$tmp"
flat=$(NXF_OFFLINE=true nextflow config "$tmp/vdjTester" -flat 2>/dev/null)
param() {
    printf '%s\n' "$flat" | awk -v k="params.$1" '$1 == k { sub(/^[^=]*= /, ""); gsub(/^'"'"'|'"'"'$/, ""); print; exit }'
}
tar_name() { basename "$1" | sed -E 's/\.tar(\.gz)?$//'; }
# The file name Nextflow gives an image in its Apptainer/Singularity cache (as in main.nf).
image_name() { printf '%s.img\n' "$(printf '%s' "${1#*://}" | tr ':/' '--')"; }

gex_url=$(param gex_ref_url); vdj_url=$(param vdj_ref_url)
trust4_version=$(param trust4_ref_version); mixcr_lib_url=$(param mixcr_imgt_library)
images="$(param immcantation_container) $(param trust4_container) $(param python_container)"
[ $skip_mixcr -eq 1 ] || images="$images $(param mixcr_container)"
for v in gex_url vdj_url trust4_version mixcr_lib_url; do
    [ -n "${!v}" ] || { echo "could not read $v from the pipeline config at $short" >&2; exit 1; }
done

# --- Sizes ------------------------------------------------------------------------------
url_size() {  # bytes; -1 if the final response is not a success
    curl -sIL --max-time 30 "$1" 2>/dev/null | tr -d '\r' | awk '
        /^HTTP\// { ok = ($2 ~ /^2/); n = 0 }
        tolower($1) == "content-length:" { n = $2 }
        END { print (ok ? n + 0 : -1) }'
}
url_row() {  # kind, bundle path, url
    local s; s=$(url_size "$3")
    if [ "$s" -lt 0 ]; then missing "$1" "$2" "not reachable: $3"; else row "$1" "$2" "$3" "$s"; fi
}
local_size() { du -skL "$1" 2>/dev/null | awk '{ print $1 * 1024 }'; }
human() { awk -v b="$1" 'BEGIN { if (b <= 0) { print "?"; exit } split("B KB MB GB TB", u); i = 1; while (b >= 1024 && i < 5) { b /= 1024; i++ } printf (i == 1 ? "%d %s" : "%.1f %s"), b, u[i] }'; }

total=0; missing=0
row() {  # kind, bundle path, source, size in bytes (0 = unknown)
    printf '%-10s %-62s %10s  %s\n' "$1" "$2" "$(human "$4")" "$3"
    total=$((total + $4))
}
missing() {  # kind, bundle path, why
    printf '%-10s %-62s %10s  %s\n' "$1" "$2" "MISSING" "$3"
    missing=$((missing + 1))
}

# --- What to collect: kind|bundle path|source ---------------------------------------
# Sources: git:<sha>, url:<url>, dir:<local dir to archive>, file:<local file>,
# image:<uri> (pulled), files:<name>=<url>;... (downloaded into one directory, archived)
plan="$tmp/plan.txt"
: > "$plan"
add() { printf '%s|%s|%s\n' "$1" "$2" "$3" >> "$plan"; }

add pipeline "pipeline/vdjTester-$short.tar.gz" "git:$sha"
add nextflow "nextflow/nextflow-$nf_version-dist" "url:https://github.com/nextflow-io/nextflow/releases/download/v$nf_version/nextflow-$nf_version-dist"
for p in $plugins; do
    add plugin "nextflow/plugins/$p" "plugin:$p"
done
for uri in $images; do
    name=$(image_name "$uri")
    if [ -f "$cache/$name" ]; then
        add container "containers/$name" "file:$cache/$name"
    else
        add container "containers/$name" "image:$uri"
    fi
done

# 10x references: the store copy (an unpacked directory) or 10x's tarball.
for url in "$gex_url" "$vdj_url"; do
    name=$(tar_name "$url")
    if [ -d "$store/$name" ]; then
        add reference "store/$name.tar" "dir:$store/$name"
    else
        add reference "store/$name.tar.gz" "url:$url"
    fi
done
trust4_name="trust4_ref_$trust4_version"
if [ -d "$store/$trust4_name" ]; then
    add reference "store/$trust4_name.tar" "dir:$store/$trust4_name"
else
    base="https://raw.githubusercontent.com/liulab-dfci/TRUST4/$trust4_version"
    add reference "store/$trust4_name.tar" "files:hg38_bcrtcr.fa=$base/hg38_bcrtcr.fa;human_IMGT+C.fa=$base/human_IMGT%2BC.fa"
fi
if [ "$reference_set" != default ]; then
    # IMGT serves only its current release, so the reference is the one built here
    # (fetch-imgt + mkvdjref by a --reference_set imgt run), never a fresh download.
    add reference "store/imgt.tar" "dir:$store/imgt"
    if [ $skip_mixcr -eq 0 ]; then
        lib_file=$(basename "$mixcr_lib_url")
        lib_name="mixcr_library_$(printf '%s' "$lib_file" | sed -E 's/\.json(\.gz)?$//')"
        if [ -d "$store/$lib_name" ]; then
            add reference "store/$lib_name.tar" "dir:$store/$lib_name"
        else
            add reference "store/$lib_name.tar" "files:$lib_file=$mixcr_lib_url"
        fi
    fi
fi
if [ -n "$cellranger_tar" ]; then
    add tool "tools/$(basename "$cellranger_tar")" "file:$cellranger_tar"
fi

# --- Dry run: the list, with sizes ---------------------------------------------------
engine=$(command -v apptainer || command -v singularity || true)
printf 'vdjTester offline bundle: commit %s, Nextflow %s, reference set %s%s\n\n' \
    "$short" "$nf_version" "$reference_set" "$([ $skip_mixcr -eq 1 ] && echo ', without MiXCR')"
printf '%-10s %-62s %10s  %s\n' KIND PATH SIZE SOURCE
while IFS='|' read -r kind path src; do
    case $src in
        git:*)    row "$kind" "$path" "git archive $short" "$(git -C "$repo" archive --format=tar.gz "$sha" | wc -c | tr -d ' ')" ;;
        url:*)    url_row "$kind" "$path" "${src#url:}" ;;
        file:*)   row "$kind" "$path" "${src#file:}" "$(local_size "${src#file:}")" ;;
        plugin:*) row "$kind" "$path" "nextflow plugin install ${src#plugin:}" 0 ;;
        dir:*)
            if [ -d "${src#dir:}" ]; then
                row "$kind" "$path" "${src#dir:}/" "$(local_size "${src#dir:}")"
            else
                missing "$kind" "$path" "${src#dir:} not found: run the pipeline once with --reference_set imgt here, or pass --store"
            fi ;;
        image:*)
            if [ -n "$engine" ]; then
                row "$kind" "$path" "$(basename "$engine") pull --arch $arch ${src#image:}" 0
            else
                missing "$kind" "$path" "not in $cache and no apptainer/singularity here to pull ${src#image:}"
            fi ;;
        files:*)
            s=0; list=${src#files:}
            IFS=';' read -r -a pairs <<< "$list"
            for pair in "${pairs[@]}"; do
                n=$(url_size "${pair#*=}")
                if [ "$n" -lt 0 ]; then s=-1; break; fi
                s=$((s + n))
            done
            if [ "$s" -lt 0 ]; then missing "$kind" "$path" "not reachable: ${pair#*=}"; else row "$kind" "$path" "${list//;/ }" "$s"; fi ;;
    esac
done < "$plan"
printf '\nTotal (known sizes): %s\n' "$(human "$total")"
if [ -z "$cellranger_tar" ]; then
    echo
    echo "Not included: Cell Ranger. 10x does not allow redistributing it: download it yourself"
    echo "(https://www.10xgenomics.com/support/software/cell-ranger/downloads, accepting the licence)"
    echo "and pass --cellranger-tar, or use a Cell Ranger installed on the target system."
    for f in "$repo"/bin/cellranger-*.tar.gz "$repo"/tools/cellranger-*.tar.gz; do
        [ -f "$f" ] && echo "  found locally: $f"
    done
fi
if [ $skip_mixcr -eq 0 ]; then
    echo
    echo "MiXCR: its licence is checked online at every run, which fails without internet."
    echo "Ask MiLaboratories for a licence that works offline; otherwise run with --skip_mixcr."
fi
if [ $missing -gt 0 ]; then
    echo; echo "$missing item(s) MISSING (see above)." >&2
    [ $dry_run -eq 1 ] && exit 0 || exit 1
fi
[ $dry_run -eq 1 ] && exit 0

# --- Collect ---------------------------------------------------------------------------
mkdir -p "$outdir"
outdir=$(cd "$outdir" && pwd)
fetch() { curl -fL --retry 5 --retry-all-errors -C - -o "$2" "$1"; }
while IFS='|' read -r kind path src; do
    dest="$outdir/$path"
    if [ -e "$dest" ] && [ "$kind" != plugin ]; then
        echo "have    $path"; continue
    fi
    echo "collect $path"
    mkdir -p "$(dirname "$dest")"
    case $src in
        git:*)  git -C "$repo" archive --format=tar.gz --prefix=vdjTester/ -o "$dest" "$sha" ;;
        url:*)  fetch "${src#url:}" "$dest.part" && mv "$dest.part" "$dest" ;;
        file:*) cp "${src#file:}" "$dest.part" && mv "$dest.part" "$dest" ;;
        dir:*)  d=${src#dir:}; tar -chf "$dest.part" -C "$(dirname "$d")" "$(basename "$d")" && mv "$dest.part" "$dest" ;;
        image:*)
            uri=${src#image:}; case $uri in *://*) ;; *) uri="docker://$uri" ;; esac
            "$engine" pull --arch "$arch" "$dest.part" "$uri" && mv "$dest.part" "$dest" ;;
        plugin:*)
            NXF_HOME="$tmp/nxf_home" nextflow plugin install "${src#plugin:}"
            mkdir -p "$outdir/nextflow/plugins"
            cp -R "$tmp/nxf_home/plugins/." "$outdir/nextflow/plugins/" ;;
        files:*)
            name=$(basename "$path" .tar); rm -rf "$tmp/$name"; mkdir -p "$tmp/$name"
            IFS=';' read -r -a pairs <<< "${src#files:}"
            for pair in "${pairs[@]}"; do fetch "${pair#*=}" "$tmp/$name/${pair%%=*}"; done
            tar -cf "$dest.part" -C "$tmp" "$name" && mv "$dest.part" "$dest" ;;
    esac
done < "$plan"
chmod +x "$outdir/nextflow/nextflow-$nf_version-dist"
cp "$repo/bin/unpack_offline_bundle.sh" "$outdir/unpack.sh"
chmod +x "$outdir/unpack.sh"

{
    printf 'key\tvalue\n'
    printf 'commit\t%s\n' "$sha"
    printf 'created\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    printf 'nextflow\t%s\n' "$nf_version"
    printf 'reference_set\t%s\n' "$reference_set"
    printf 'mixcr\t%s\n' "$([ $skip_mixcr -eq 1 ] && echo skipped || echo included)"
    printf 'images\t%s\n' "$images"
    printf 'image_arch\t%s\n' "$arch"
    while IFS='|' read -r kind path src; do printf 'source:%s\t%s\n' "$path" "$src"; done < "$plan"
} > "$outdir/bundle_info.tsv"

cat > "$outdir/README.txt" <<EOF
vdjTester offline bundle (commit $short, Nextflow $nf_version)

Move this whole directory into the air-gapped system (on Bianca through the wharf),
then there:

    sha256sum -c MANIFEST.sha256        # every file intact?
    ./unpack.sh /proj/<project>/vdjtester

and follow vdjTester/docs/OFFLINE.md in the unpacked pipeline.

Cell Ranger: not redistributable (10x licence), so only here if whoever made the
bundle added their own download (tools/). Otherwise download it from 10x yourself,
accepting its licence, and bring it in the same way, or use a Cell Ranger installed on
the system; pass its directory as --cellranger_dir.

MiXCR: checks its licence online at every run, which fails without internet. Ask
MiLaboratories for a licence that works offline and run
with --mixcr_offline_licence; otherwise run with --skip_mixcr.
EOF

(
    cd "$outdir"
    sum=$(command -v sha256sum >/dev/null && echo sha256sum || echo 'shasum -a 256')
    find . -type f ! -name MANIFEST.sha256 ! -name '*.part' | LC_ALL=C sort | sed 's#^\./##' \
        | while IFS= read -r f; do $sum "$f"; done > MANIFEST.sha256
)
echo
echo "Bundle in $outdir ($(du -sh "$outdir" | cut -f1)); verify with: (cd $outdir && sha256sum -c MANIFEST.sha256)"
