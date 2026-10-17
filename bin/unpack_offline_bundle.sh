#!/usr/bin/env bash
# Unpack a vdjTester offline bundle (bin/make_offline_bundle.sh), on the air-gapped
# system. Shipped in the bundle as unpack.sh; see docs/OFFLINE.md.
#
#   ./unpack.sh DEST
#
# Checks MANIFEST.sha256 first, then lays out DEST as -profile offline/bianca expect
# with --bundle_dir DEST:
#   DEST/vdjTester/            the pipeline
#   DEST/nextflow/nextflow     Nextflow (and DEST/nextflow/plugins, if bundled)
#   DEST/containers/*.img      images (container_cache)
#   DEST/store/<name>/         references (store_dir)
#   DEST/tools/cellranger      Cell Ranger, if the bundle has it (cellranger_dir)
set -euo pipefail

bundle=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
[ $# -eq 1 ] || { echo "Usage: $0 DEST" >&2; exit 2; }
mkdir -p "$1"
dest=$(cd "$1" && pwd)

echo "Checking checksums in $bundle ..."
(cd "$bundle" && sha256sum --quiet -c MANIFEST.sha256)

mkdir -p "$dest/store" "$dest/containers" "$dest/tools" "$dest/nextflow"

tar -xzf "$bundle"/pipeline/vdjTester-*.tar.gz -C "$dest"
cp "$bundle"/nextflow/nextflow-*-dist "$dest/nextflow/nextflow"
chmod +x "$dest/nextflow/nextflow"
if [ -d "$bundle/nextflow/plugins" ]; then
    cp -R "$bundle/nextflow/plugins" "$dest/nextflow/"
fi
for img in "$bundle"/containers/*.img; do
    [ -e "$img" ] || continue
    echo "image     $(basename "$img")"
    cp "$img" "$dest/containers/"
done

# One archive per store_dir entry, unpacked under its own name (top directory stripped,
# as the pipeline's own download step does).
for archive in "$bundle"/store/*.tar "$bundle"/store/*.tar.gz; do
    [ -e "$archive" ] || continue
    name=$(basename "$archive" | sed -E 's/\.tar(\.gz)?$//')
    echo "reference $name"
    rm -rf "$dest/store/$name.partial"
    mkdir -p "$dest/store/$name.partial"
    tar -xf "$archive" -C "$dest/store/$name.partial" --strip-components=1
    rm -rf "$dest/store/$name"
    mv "$dest/store/$name.partial" "$dest/store/$name"
done

for crtar in "$bundle"/tools/cellranger-*.tar.gz; do
    [ -e "$crtar" ] || continue
    echo "tool      $(basename "$crtar")"
    tar -xzf "$crtar" -C "$dest/tools"
    crdir=$(basename "$crtar" .tar.gz)
    ln -sfn "$crdir" "$dest/tools/cellranger"
done

cat <<EOF

Unpacked into $dest. To run (see $dest/vdjTester/docs/OFFLINE.md):

    module load <java 17 or later>      # Nextflow needs Java 17+
    export NXF_OFFLINE=true NXF_HOME=$dest/nextflow/home
    [ -d $dest/nextflow/plugins ] && export NXF_PLUGINS_DIR=$dest/nextflow/plugins
    $dest/nextflow/nextflow run $dest/vdjTester -profile bianca \\
        --bundle_dir $dest --slurm_account <project> \\
        --samplesheet <samplesheet.csv> --outdir <results> --skip_mixcr
EOF
[ -e "$dest/tools/cellranger" ] || echo "
No Cell Ranger in the bundle: pass --cellranger_dir <unpacked cellranger-x.y.z>."
