process AGGREGATESTRINGTIE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/pip_pandas:3f82ac0fa1619d28':
        'community.wave.seqera.io/library/pip_pandas:99d2c619734f6c46' }"

    input:
    tuple val(meta), path(abundance)

    output:
    tuple val(meta), path("*.aggregated.txt"), emit: abundance
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
python3 - <<EOF
import pandas as pd
import sys

def aggregate_stringtie_genes(input_file, output_file):
    # Read the StringTie gene abundance file
    df = pd.read_csv(input_file, sep='\\t')
    
    
    # Check for duplicates
    duplicates = df['Gene ID'].duplicated(keep=False)
    if duplicates.any():
        duplicate_count = duplicates.sum()
        unique_duplicate_genes = len(df[duplicates]['Gene ID'].unique())
        print(f"Found {duplicate_count} duplicate entries for {unique_duplicate_genes} genes")
        
        # Show which genes have duplicates
        duplicate_genes = df[duplicates]['Gene ID'].unique()
        for gene in duplicate_genes[:5]:  # Show first 5 only
            gene_entries = df[df['Gene ID'] == gene]
            print(f"Gene {gene}: {len(gene_entries)} entries")
    
    # Aggregate duplicates using first non-zero strategy
    def first_non_zero_or_first(series):
        \"\"\"Return first non-zero value, or first value if all are zero\"\"\"
        non_zero = series[series > 0]
        if len(non_zero) > 0:
            return non_zero.iloc[0]
        else:
            return series.iloc[0]
    
    aggregated = df.groupby('Gene ID').agg({
        'Gene Name': 'first',  # Take first occurrence
        'Reference': 'first',  # Take first occurrence  
        'Strand': 'first',     # Take first occurrence
        'Start': 'min',        # Take minimum start position
        'End': 'max',          # Take maximum end position
        'Coverage': first_non_zero_or_first,  # First non-zero coverage
        'FPKM': first_non_zero_or_first,      # First non-zero FPKM
        'TPM': first_non_zero_or_first        # First non-zero TPM
    }).reset_index()
    
    # Reorder columns to match original
    aggregated = aggregated[['Gene ID', 'Gene Name', 'Reference', 'Strand', 'Start', 'End', 'Coverage', 'FPKM', 'TPM']]
    
    print(f"After aggregation: {len(aggregated)} rows")
    
    # Save aggregated results
    aggregated.to_csv(output_file, sep='\\t', index=False)
    
    return len(df) - len(aggregated)  # Return number of duplicates removed

# Process the file
duplicates_removed = aggregate_stringtie_genes("${abundance}", "${prefix}.gene.abundance.aggregated.txt")
print(f"Successfully aggregated ${abundance}, removed {duplicates_removed} duplicate entries")
EOF

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    python: \$(python --version | sed 's/Python //g')
    pandas: \$(python -c "import pandas; print(pandas.__version__)")
END_VERSIONS
    """
}