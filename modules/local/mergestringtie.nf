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
    
    print(f"Processing {len(files)} files: {files}")

    for file in files:
        # Handle both regular and aggregated abundance files
        if file.endswith(".gene.abundance.aggregated.txt"):
            prefix = os.path.basename(file).replace(".gene.abundance.aggregated.txt", "")
        else:
            prefix = os.path.basename(file).replace(".gene.abundance.txt", "")
            
        print(f"Processing file: {file} -> sample: {prefix}")
        df = pd.read_csv(file, sep="\t")
        print(f"  Loaded {len(df)} genes")
        
        # Keep only relevant columns
        df = df[["Gene ID", "Gene Name", "TPM"]].copy()

        df.rename(columns={"Gene ID": "gene_id", "Gene Name": "gene_name", "TPM": prefix}, inplace=True)
    
        if merged_df is None:
            merged_df = df
            print(f"  Initial dataframe: {len(merged_df)} genes")
        else:
            before_merge = len(merged_df)
            merged_df = pd.merge(merged_df, df, on=["gene_id"], how="left")
            print(f"  After merge: {before_merge} -> {len(merged_df)} genes")

    # Sort by Gene ID
    merged_df.sort_values("gene_id", inplace=True)

    # Sort columns by sample
    merged_df = merged_df[['gene_id', 'gene_name'] + sorted(merged_df.columns.difference(['gene_id', 'gene_name']))]

    print(f"Final merged dataframe: {len(merged_df)} genes x {len(merged_df.columns)} columns")
    print(f"Sample columns: {sorted(merged_df.columns.difference(['gene_id', 'gene_name']))}")
    
    # Save result
    output_file = "${params.genetable_outfile}.tsv"
    merged_df.to_csv(output_file, sep="\t", index=False)
    print(f"Saved TPM matrix to: {output_file}")
    """
}
