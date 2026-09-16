process DOWNLOAD_GENE_BED {
    tag { genome }
    storeDir "${params.outdir}/reference"

    // Needs wget, so it cannot run bare under a container-only profile.
    conda "${moduleDir}/environment.yml"
    container params.echidna_container

    input:
    val genome

    output:
    path "${genome}_refGene.bed", emit: bed

    script:
    """
    set -euo pipefail

    # Compute nodes on many clusters have no outbound network. If this fails,
    # fetch the BED once on a login node and pass --gene_bed instead.
    wget -qO- "https://hgdownload.soe.ucsc.edu/goldenPath/${genome}/database/refGene.txt.gz" \\
        | gunzip -c \\
        | awk 'BEGIN{OFS="\\t"} {print \$3, \$5, \$6, \$13}' \\
        | sort -k1,1 -k2,2n \\
        > ${genome}_refGene.bed \\
        || { echo "ERROR: could not download the refGene table for '${genome}' from UCSC." >&2
             echo "       Compute nodes often have no outbound network. Download it once on a" >&2
             echo "       login node and re-run with: --gene_bed /path/to/${genome}_refGene.bed" >&2
             exit 1; }

    [[ -s "${genome}_refGene.bed" ]] || {
        echo "ERROR: the gene BED for '${genome}' came back empty." >&2
        echo "       Check that '${genome}' is a valid UCSC genome name, or supply an existing" >&2
        echo "       BED with: --gene_bed /path/to/${genome}_refGene.bed" >&2
        exit 1; }
    """

    stub:
    """
    echo -e 'chr1\\t1\\t2\\tGENE1' > ${genome}_refGene.bed
    """
}
