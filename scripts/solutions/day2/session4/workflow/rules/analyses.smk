'''
This Snakefile contains rules used to gather files containing read counts into
a single table and perform Differential Expression Analyses.
'''


# Input function used in rule count_table
def get_gene_counts(wildcards):
    '''
    This function lists all the gene count tables of samples in the config file
    '''
    return [f"results/{sample}/{sample}_genes_read_quantification.tsv"
            for sample in config['samples']]


rule count_table:
    '''
    This rule merges all the gene count tables of an assembly into one table.
    '''
    input:
        get_gene_counts
    output:
        table = 'results/total_count_table.tsv'
    log:
        'logs/total_count_table.log'
    benchmark:
        'benchmarks/total_count_table.txt'
    conda:
        '../envs/py.yaml'
    resources:
        mem_mb = 500
    threads: 1
    script:
        '../scripts/count_table.py'


rule differential_expression:
    '''
    This rule detects DEGs and plots associated visual control graphs (PCA,
    heatmaps...).
    '''
    input:
        table = rules.count_table.output.table
    output:
        deg = 'results/deg_list.tsv',
        pdf = 'results/deg_plots.pdf'
    log:
        'logs/differential_expression.log'
    benchmark:
        'benchmarks/differential_expression.txt'
    container:
        'docker://geertvangeest/deseq2:v1'
    resources:
        mem_gb = 1
    threads: 2
    script:
        '../scripts/DESeq2.R'
