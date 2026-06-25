# Neighbors -> UMAP -> Leiden -> marker genes -> (optional) annotation.
rule cluster_annotate:
    input:  f"{RESULTS}/integrate/integrated.h5ad"
    output:
        h5ad=f"{RESULTS}/cluster/annotated.h5ad",
        markers=f"{RESULTS}/cluster/markers.csv",
        umap=f"{RESULTS}/figures/umap_clusters.png",
    conda:  "../envs/scanpy.yaml"
    script: "../scripts/cluster_annotate.py"
