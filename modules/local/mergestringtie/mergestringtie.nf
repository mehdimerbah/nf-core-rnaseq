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
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
python3 - <<EOF

import pandas as pd
import os

merged_df = None

files = ${infiles.collect { '"' + it.name + '"' }}  # local staged copies


for file in files:
    # Handle both regular and aggregated abundance files
    if file.endswith(".gene.abundance.aggregated.txt"):
        prefix = os.path.basename(file).replace(".gene.abundance.aggregated.txt", "")
    else:
        prefix = os.path.basename(file).replace(".gene.abundance.txt", "")
        
    df = pd.read_csv(file, sep="\t")
    
    # Keep only relevant columns
    df = df[["Gene ID", "Gene Name", "TPM"]].copy()

    df.rename(columns={"Gene ID": "gene_id", "Gene Name": "gene_name", "TPM": prefix}, inplace=True)

    if merged_df is None:
        merged_df = df
    else:
        merged_df = pd.merge(merged_df, df[["gene_id", prefix]], on=["gene_id"], how="left")


# Sort by Gene ID
merged_df.sort_values("gene_id", inplace=True)

# Sort columns by sample (gene_name should still be present from the first dataframe)
sample_columns = sorted([col for col in merged_df.columns if col not in ['gene_id', 'gene_name']])
merged_df = merged_df[['gene_id', 'gene_name'] + sample_columns]

print(f"Final merged dataframe: {len(merged_df)} genes x {len(merged_df.columns)} columns")
print(f"Sample columns: {sample_columns}")

# Save result
output_file = "${params.genetable_outfile}.tsv"
merged_df.to_csv(output_file, sep="\t", index=False)
EOF

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    python: \$(python --version | sed 's/Python //g')
    pandas: \$(python -c "import pandas; print(pandas.__version__)")
END_VERSIONS
    """
}
