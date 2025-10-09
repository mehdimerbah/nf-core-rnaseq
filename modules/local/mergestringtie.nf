process MERGESTRINGTIE {
    debug true
   
conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/pip_pandas:3f82ac0fa1619d28':
        'community.wave.seqera.io/library/pip_pandas:99d2c619734f6c46' }"

    input:
    path(infiles)
    // val metamaps

    output:
    path "*.tsv", emit: tsv
    // path tsv, emit: tsv

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    #!/bin/env python

    import pandas as pd
    import os

    merged_df = None

    files = ${infiles.collect { '"' + it.name + '"' }}  # local staged copies

    for file in files:
        prefix = os.path.basename(file).replace(".gene.abundance.txt", "")
        df = pd.read_csv(file, sep="\t")
        
        # Keep only relevant columns
        df = df[["Gene ID", "Gene Name", "TPM"]].copy()

        df.rename(columns={"Gene ID": "gene_id", "Gene Name": "gene_name", "TPM": prefix}, inplace=True)
    
        if merged_df is None:
            merged_df = df
        else:
            merged_df = pd.merge(merged_df, df, on=["gene_id", "gene_name"], how="outer")

    # Sort by Gene ID
    merged_df.sort_values("gene_id", inplace=True)

    # Save result
    merged_df.to_csv("${params.genetable_outfile}.tsv", sep="\t", index=False)
    """
}
