## Learning outcomes

**After having completed this chapter you will be able to:**

* Understand the structure of a Snakemake workflow
* Write rules and Snakefiles to produce the desired outputs
* Chain rules together
* Run a Snakemake workflow

## Exercises

!!! note "Command &lt;cmd_name&gt; not found"  <!-- AT. Check if special characters worked -->
    If you try to rum a command and get an error such as _Command 'snakemake' not found_, you are probably not in the right environment. To list them, use `mamba env list`. Then activate the right environment with `mamba activate <env_name>`. You can deactivate an environment with `mamba deactivate`.

### Workflow structure

It is strongly advised to implement your answers in a directory called `workflow` (the reason for this will be explained later). You are free to chose the names and location of files for the different steps of your workflow, but we recommend that you at least group all outputs from the workflow in a `results` directory within the `workflow` directory.

### Creating a basic rule

Rules are the basic blocks of a Snakemake workflow. A **rule** is like a recipe indicating how to produce a specific **output**; the actual application of a rule to create an output is called a **job**. A rule is defined in a Snakefile with the _keyword_ `rule`, and contains _directives_ which indicate the rule's properties. We will learn about other directives later in the course.

To create the simplest rule possible, we need at least two directives: 
- `output`: path of the output file for this rule
- `shell`: shell commands to execute in order to generate the output

**Exercise:** The following example shows the minimal syntax to implement a rule. What do you think it does? Does it create a file? If so, how is it called?

```python
rule first_step:
    output:
        'results/first_step.txt'
    shell:
        'echo “snakemake” > results/first_step.txt'
```

??? done "Answer"
    This rule uses the `echo` shell command to print the line "snakemake" in an output file called `first_step.txt`, located in the `results` folder

Rules are defined and written in a file called **Snakefile** (note the capital `S` and the absence of extension in the filename). This file should be located at the root of the workflow directory (here, `workflow/Snakefile`).

!!! note "Paths in Snakemake"
    All the paths in Snakefile are relative to the directory containing the Snakefile.

**Exercise:** Create a Snakefile and copy the rule in it. Because the Snakemake language is basically Python, _do not forget to keep the indentation as is and use space characters in the indents instead of tabs._

### Executing a workflow with a precise output

It is now time to execute your first worklow!

**Exercise:** Execute the workflow with `snakemake --cores 1 <output_name>`. What value should you use as _output_name_? Do you see the file?

??? done "Answer"
    Execute the workflow: `snakemake --cores 1 results/first_step.txt`
    Visualise your directory content: `ls -alh results/`
    Check the output content with `cat results/first_step.txt`

Note that during the execution of the workflow, Snakemake automatically created the **missing folder** (`results/`) in the output path. If several folders are missing (for example, here, `test1/test2/test3/first_step.txt`), Snakemake will create **all of them**.

**Exercise:** Rerun the exact same command. What happens?

??? done "Answer"
<!-- AT. Check how this looks -->
    Well... Nothing! You get a message saying that Snakemake did not do anything:

    ```
    Building DAG of jobs...
    Nothing to be done (all requested files are present and up to date).
    ```

    By default, existing outputs are only generated again if the input of the rule that generates them is newer than them, based on file modification dates. You can change this behaviour and force the rerun by using the `-f` option: `snakemake --cores 1 -f results/first_step.txt` or the `-F` option to force recreate ALL the outputs of the workflow.

Note that in the previous example, values for these two directives are **strings**. For the `shell` directive (we will see other types of directive values later in the course), long string can be written on multiple lines for clarity, simply using a set of quotes for each line:
<!-- AT. Check how next lines look -->
```python
shell:
    'echo "I want to print a very very very very very very '
    'very very very very long string in my output" > results/first_step.txt'
```

### Using the input directive

The next directive used by most rules is `input`. Like `output`, `input` indicates the path to a file that is required by the rule to generate the output. In the following example, we modified the previous rule to use the file previously created `results/first_step.tsv` as an input and copy this file to `results/second_step.txt`:

```python
rule second_step:
    input:
        'results/first_step.txt'
    output:
        'results/second_step.txt'
    shell:
        'cp results/first_step.txt results/second_step.txt'
```

Note that with this rule definition, Snakemake **will not run** if `data/first_step.tsv` does not exist!

**Exercise:** Modify your first rule to add an input and execute the workflow. Check that the output was created and that the files are identical.

??? done "Answer"
    Execute the workflow: `snakemake --cores 1 results/second_step.txt`
    Visualise your directory content: `ls -alh results/`
    Check that the files are identical `diff results/first_step.txt results/second_step.txt`

### Using several rules in a workflow

Creating one Snakefile per rule doesn't seem like a good solution, so let's try to improve this.

**Exercise:** Delete the `results/` folder, gather the two previous rules in the same Snakefile (place the `first_step` rule first) and try to run the workflow **without specifying an output**. What happens?

??? done "Answer"
    Execute the workflow without outputs: `snakemake --cores 1`.
    When executed, Snakemake tries to generate a specific output called **target**, and resolves all dependencies based on this target. A target can be any output that can be generated by any rule in the workflow. When you do not specify a target, the default one is the output of the first rule in the Snakefile, here `results/first_step.txt`. If you had placed the `second_step` rule in first position, Snakemake would have crashed because the input for this rule does not exist. If you are in advance, feel free to try this!

**Exercise:** With this in mind, use a space-separated list of targets (instead of one filename) in your command to generate multiple targets. Use the `-F` to force the rerun of the whole workflow or delete your `results/` folder beforehand.

??? done "Answer"
    Execute the workflow with multiple targets: `snakemake --cores 1 -F results/first_step.txt results/second_step.txt`
    You should now see Snakemake execute the 2 rules and produce both targets/outputs.

### Chaining rules

Once again, writing all the outputs in the `snakemake` command is very time-consuming, error-prone and annoying! Imagine what happens when your workflow generate tens of outputs? Fortunately, there is a way to simplify this, which relies on rules dependency. The core principle of Snakemake's execution is to compute a Directed Acyclic Graph (DAG) that summarizes dependencies between all inputs and outputs required to generate the final desired output. For each job, starting from the jobs generating the final output, Snakemake checks if the required inputs exist. If they do not, the software looks for a rule that generates the input. This process is repeated until all dependencies are resolved. This is why Snakemake is said to have a 'bottom-up' approach: it starts from the last outputs and go back to the first inputs.

!!! hint
    Your Snakefile should look like this:

    ```python
    rule first_step:
        output:
            'results/first_step.txt'
        shell:
            'echo “snakemake” > results/first_step.txt'

    rule second_step:
        input:
            'results/first_step.txt'
        output:
            'results/second_step.txt'
        shell:
            'cp results/first_step.txt results/second_step.txt'

    ```

**Exercise:** Delete the `results/` folder, identify your final output and execute the workflow **specifying only this output** in the command.

??? done "Answer"
    Execute the workflow: `snakemake --cores 1 results/second_step.txt`
    Visualise your directory content: `ls -alh results/`
    You should now see Snakemake execute the 2 rules and produce both targets/outputs. To generate the output `results/second_step.txt`, Snakemake requires the input `results/first_step.txt`. Before the workflow is executed, this file does not exist, therefore, Snakemake looks for a rule that generates `results/first_step.txt`, in this case the first defined rule `first_step`. The process is then repeated for `first_step`. In this case, the rule does not require an input, so all dependencies are resolved, and Snakemake can generate the DAG.

### Important notes on chaining rules

#### Rules produce unique outputs

Because of the rules dependency process, by default, an output can only be generated by a single rule. Otherwise, Snakemake cannot decide which rule to use to generate this output, and the rules are considered **ambiguous**. In practice, there are ways to deal with ambiguous rules, but we will not cover them in this course (see [the relevant section in the official documentation](https://snakemake.readthedocs.io/en/v7.32.3/snakefiles/rules.html#handling-ambiguous-rules)).

#### Rules dependency can be written more easily

It is possible to refer to the output of a rule directly in another rule with the syntax `rules.<rule_name>.output`. Note that you don't need quotes around this statement, because it is a Snakemake object. The following example implements this syntax for the two rule defined above:

```python
rule first_step:
    output:
        'results/first_step.txt'
    shell:
        'echo “snakemake” > results/first_step.txt'

rule second_step:
    input:
        rules.first_step.output
    output:
        'results/second_step.txt'
    shell:
        'cp results/first_step.txt results/second_step.txt'
```

This method has several advantages, among which:
* It limits the risk of error because you do not have to write the same filename at several locations
* A change in output name will be automatically propagated to rules that depend on it, *i.e.* the name only has to be changed once
* This makes the code much clearer and easier to understand: with this syntax, you instantly know the object type (`rule`), how it is created (`first_step`) and what it is (`output`)
