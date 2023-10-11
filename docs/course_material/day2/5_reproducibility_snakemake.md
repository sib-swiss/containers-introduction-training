## Learning outcomes

**After having completed this chapter you will be able to:**

* Create and use an input function
* Deploy a conda environment and run scripts within it
* Deploy a Docker/Singularity container and run scripts within it

## Exercises

In this series of exercises, we will create the last two rules of the workflow. Each rule will execute a script, one in Python and one in R, and both rules will have dedicated environment that you will need to take into account in the snakefiles.

!!! note "Development and back-up"
    During this session, we will modify our snakefiles quite heavily, so it may be a good idea to start by making a back-up: `cp -r worklow/ worklow_backup`. As a general rule, if you have a doubt on the code you are developing, do not hesitate to make a back-up.

!!! hint
    This is not a programming course, so you won't need to write the scripts: they were already prepared for you!

### Creating a rule to gather read count files







### Creating a rule to detect Differentially Expressed Genes









<!-- ### Optimising a workflow by multi-threading

When working with real datasets, most processes are very long and computationally expensive. Fortunately, they can be parallelised very efficiently to decrease the computation time by using several [threads](https://en.wikipedia.org/wiki/Thread_(computing)) for a single job.

**Exercise:** Parallelise as much processes as possible using the `threads` directive and test its effect:

1. Identify which software can make use of parallelisation
1. Identify in each software the parameter that controls multi-threading
1. Implement the multi-threading

!!! hint
    * Check the software documentation and parameters with the `-h/--help` flags
    * Remember that multi-threading only applies to software that can make use of a threads parameters, Snakemake itself cannot parallelize a software automatically
    * Remember that you need to add threads to the Snakemake rule but also to the commands! Just increasing the number of threads in Snakemake will not magically run a command with multiple threads
    * Remember that you have 4 threads in total, so even if you ask for more in a rule, Snakemake will cap this value at 4. And if you use 4 threads in a rule, that means that no other job can run parallel!

??? done "Answer"
    It turns out that all the software except `samtools index` can handle multi-threading:

    * `atropos trim`, `hisat2`, `samtools view`, and `samtools sort` use the `--threads` option
    * `featureCounts` uses the `-T` option

    Let's use 4 threads for the mapping step and 2 for the other steps. Your Snakefile should look like this:
    ```
    rule fastq_trim:
        '''
        This rule trims paired-end reads to improve their quality. Specifically, it removes:
        - Low quality bases
        - A stretches longer than 20 bases
        - N stretches
        '''
        input:
            reads1 = 'data/{sample}_1.fastq',
            reads2 = 'data/{sample}_2.fastq',
        output:
            trim1 = 'results/{sample}/{sample}_atropos_trimmed_1.fastq',
            trim2 = 'results/{sample}/{sample}_atropos_trimmed_2.fastq'
        log:
            'logs/{sample}/{sample}_atropos_trimming.log'
        benchmark:
            'benchmarks/{sample}/{sample}_atropos_trimming.txt'
        resources:
            mem_mb = 500
        threads: 2
        shell:
            '''
            echo "Trimming reads in <{input.reads1}> and <{input.reads2}>" > {log}
            atropos trim -q 20,20 --minimum-length 25 --trim-n --preserve-order --max-n 10 --no-cache-adapters -a "A{{20}}" -A "A{{20}}" --threads {threads} -pe1 {input.reads1} -pe2 {input.reads2} -o {output.trim1} -p {output.trim2} &>> {log}
            echo "Trimmed files saved in <{output.trim1}> and <{output.trim2}> respectively" >> {log}
            echo "Trimming report saved in <{log}>" >> {log}
            '''

!!! note "Explicit is better than implicit"
    Even if a software cannot multi-thread, it is useful to add `threads: 1` in the rule to keep the rule consistency and clearly state that the software works with a single thread.


Your DAG should resemble this:
<figure align="center">
  <img src="../../assets/images/all_samples_dag.png" width="100%"/>
</figure>

And your filegraph, this:

<figure align="center">
  <img src="../../assets/images/all_samples_filegraph.png" height=500/>
</figure>
 -->
