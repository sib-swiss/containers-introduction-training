## Learning outcomes

**After having completed this chapter you will be able to:**

* Login to a remote machine with `ssh`
* Use `apptainer pull` to convert an image from dockerhub to the 'apptainer image format' (`.sif`)
* Execute a apptainer container
* Explain the difference in default mounting behaviour between `docker` and `apptainer`
* Use `apptainer shell` to generate an interactive shell inside a `.sif` image
* Search and use images with both `docker` and `apptainer` from [bioconda](https://bioconda.github.io/index.html)

## Material

[:fontawesome-solid-file-pdf: Download the presentation](../../assets/pdf/apptainer.pdf){: .md-button }

<iframe width="560" height="315" src="https://www.youtube.com/embed/d3kxtzUutjk" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

* [Apptainer documentation](https://apptainer.org/docs/user/latest/)
* [An article on Docker vs Apptainer](https://pythonspeed.com/articles/containers-filesystem-data-processing/)
* [Using conda and containers with snakemake](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#containerization-of-conda-based-workflows)

## Exercises

!!! info "Important"
    This last part should be done on the remote server that has apptainer installed. If you haven't set up a VScode SSH connection to the server yet, you can find instructions [here](../../precourse.md#ssh-connections).

### Pulling an image

Apptainer can take several image formats (e.g. a `docker` image), and convert them into it's own `.sif` format. Unlike `docker` this image doesn't live in a local image cache, but it's stored as an actual file.

**Exercise:** On the remote server, pull the docker image that has the adjusted default `CMD` that we have pushed to dockerhub [in this exercise](dockerfiles.md#using-cmd) (`ubuntu-figlet-df:v3`) with `apptainer pull`. The syntax is:

```sh
apptainer pull docker://[USER NAME]/[IMAGE NAME]:[TAG]
```

??? done "Answer"
    ```sh
    apptainer pull docker://[USER NAME]/ubuntu-figlet:v3
    ```
    This will result in a file called `ubuntu-figlet_v3.sif`

!!! note
    If you weren't able to push the image in the previous exercises to your docker hub, you can use `geertvangeest` as username to pull the image.

### Executing an image

These `.sif` files can be run as standalone executables:

```sh
./ubuntu-figlet_v3.sif
```

!!! note
    This is shorthand for:

    ```sh
    apptainer run ubuntu-figlet_v3.sif
    ```

And you can overwrite the default command like this:

```sh
apptainer run [IMAGE NAME].sif [COMMAND]
```

!!! note
    In this case, you can also use

    ```sh
    ./[IMAGE NAME].sif [COMMAND]
    ```

    However, most applications require `apptainer run`. Especially if you want to provide options like `--bind` (for mounting directories). 


**Exercise:** Run the `.sif` file without a command, and with a command that runs `figlet`. Do you get expected output? Do the same for the python or R image you've created in the previous chapter.

!!! note "Entrypoint and apptainer"
    The `daterange` image has an entrypoint set, and `apptainer run` does not overwrite it. In order to ignore both the entrypoint and cmd use `apptainer exec`.  

??? done "Answer"
    Running it without a command (`./ubuntu-figlet_v3.sif`) should give:

    ```
    __  __         _                                                 _        _
    |  \/  |_   _  (_)_ __ ___   __ _  __ _  ___  __      _____  _ __| | _____| |
    | |\/| | | | | | | '_ ` _ \ / _` |/ _` |/ _ \ \ \ /\ / / _ \| '__| |/ / __| |
    | |  | | |_| | | | | | | | | (_| | (_| |  __/  \ V  V / (_) | |  |   <\__ \_|
    |_|  |_|\__, | |_|_| |_| |_|\__,_|\__, |\___|   \_/\_/ \___/|_|  |_|\_\___(_)
           |___/                     |___/
    ```
    Which is the default command that we changed in the `Dockerfile`.

    Running with a another `figlet` command:

    ```sh
    ./ubuntu-figlet_v3.sif figlet 'Something else'
    ```

    Should give:

    ```
    ____                       _   _     _                    _
    / ___|  ___  _ __ ___   ___| |_| |__ (_)_ __   __ _    ___| |___  ___
    \___ \ / _ \| '_ ` _ \ / _ \ __| '_ \| | '_ \ / _` |  / _ \ / __|/ _ \
    ___) | (_) | | | | | |  __/ |_| | | | | | | | (_| | |  __/ \__ \  __/
    |____/ \___/|_| |_| |_|\___|\__|_| |_|_|_| |_|\__, |  \___|_|___/\___|
                                                 |___/

    ```
    === "R"
        Pulling the `deseq2` image:

        ```sh
        apptainer pull docker://[USER NAME]/deseq2:v1
        ```

        Running it without command:

        ```sh
        ./deseq2_v1.sif
        ```

        Running with a command:

        ```sh
        ./deseq2_v1.sif --rows 200
        ```

        To overwrite both entrypoint and the command:

        ```sh
        apptainer exec deseq2_v1.sif test_deseq2.R --rows 200
        ```
        
    === "Python"
        Pulling the `daterange.py` image:

        ```sh
        apptainer pull docker://[USER NAME]/daterange:v1
        ```

        Running it without command:

        ```sh
        ./daterange_v1.sif
        ```

        Running with a command:

        ```sh
        ./daterange_v1.sif --date 20221005
        ```

        To overwrite both entrypoint and the command:

        ```sh
        apptainer exec daterange_v1.sif daterange.py --date 20221005
        ```

### Mounting with Apptainer

Apptainer is also different from Docker in the way it handles mounting. By default, Apptainer binds your home directory and a number of paths in the root directory to the container. This results in behaviour that is almost like if you are working on the directory structure of the host.  

!!! note "If your directory is not mounted by default"
    It depends on the apptainer settings whether most directories are mounted by default to the container. If your directory is not mounted, you can do that with the `--bind` option of `apptainer exec`:

    ```sh
    apptainer exec --bind /my/dir/to/mount/ [IMAGE NAME].sif [COMMAND]
    ```


Running the command `pwd` (full name of current working directory) will therefore result in a path on the host machine:

```sh
./ubuntu-figlet_v3.sif pwd
```

**Exercise:** Run the above command. What is the output? How would the output look like if you would run a similar command with Docker?

!!! hint
    A similar Docker command would look like (run this on your local computer):

    ```sh
    docker run --rm ubuntu-figlet:v3 pwd
    ```

??? done "Answer"
    The output of `./ubuntu-figlet_v3.sif pwd` is the current directory on the host: i.e. `/home/username` if you have it in your home directory. The output of `docker run --rm ubuntu-figlet:v3 pwd` (on the local host) would be `/`, which is the default workdir (root directory) of the container. As we did not mount any host directory, this directory exists only within the container (i.e. separated from the host).

### Interactive shell

If you want to debug or inspect an image, it can be helpful to have a shell inside the container. You can do that with `apptainer shell`:

```sh
apptainer shell ubuntu-figlet_v3.sif
```

!!! note
    To exit the shell type `exit`.

**Exercise:** Can you run `figlet` inside this shell?

??? done "Answer"
    Yes:
    ```
    Apptainer> figlet test
     _            _
    | |_ ___  ___| |_
    | __/ _ \/ __| __|
    | ||  __/\__ \ |_
     \__\___||___/\__|

    ```

During the lecture you have learned that apptainer takes over the user privileges of the user on the host. You can get user information with command like `whoami`, `id`, `groups` etc.

**Exercise:** Run the `figlet` container interactively. Do you have the same user privileges as if you were on the host? How is that with `docker`?

??? done "Answer"
    A command like `whoami` will result in your username printed at stdout:

    ```
    Apptainer> whoami
    myusername
    Apptainer> id
    uid=1030(myusername) gid=1031(myusername) groups=1031(myusername),1001(condausers)
    Apptainer> groups
    myusername condausers
    ```

    With apptainer, you have the same privileges inside the apptainer container as on the host. If you do this in the docker container (based on the same image), you'll get output like this:

    ```
    root@a3d6e59dc19d:/# whoami
    root
    root@a3d6e59dc19d:/# groups
    root
    root@a3d6e59dc19d:/# id
    uid=0(root) gid=0(root) groups=0(root)
    ```

### A bioinformatics example (extra)

All bioconda packages also have a pre-built container. Have a look at the [bioconda website](https://bioconda.github.io/index.html), and search for `fastqc`. In the search results, click on the appropriate record (i.e. package 'fastqc'). Now, scroll down and find the docker command to obtain the fastqc image. 

**Exercise:** Check out the container image at quay.io, by following [quay.io/biocontainers/fastqc](https://quay.io/biocontainers/fastqc). Choose a tag, and pull it with `apptainer pull`.

??? done "Answer"
    It's up to you which tag you choose. The tag with the latest version is `0.12.1--hdfd78af_0`.

    ```sh
    apptainer pull docker://quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0
    ```

!!! tip "Apptainer images at galaxy.org"
    Most (if not all) biocontainer images are available as apptainer (singularity) image at [https://depot.galaxyproject.org/singularity/](https://depot.galaxyproject.org/singularity/). You can simply download them with `wget` or `curl`, e.g.:

    ```sh
    wget https://depot.galaxyproject.org/singularity/fastqc%3A0.12.1--hdfd78af_0
    ```

Let's test the image. Download some sample reads first:

```sh
mkdir reads
cd reads
wget https://introduction-containers.s3.eu-central-1.amazonaws.com/ecoli_reads.tar.gz
tar -xzvf ecoli_reads.tar.gz
rm ecoli_reads.tar.gz
```

Now you can simply run the image as an executable preceding the commands you would like to run within the container. E.g. running `fastqc` would look like:

```sh
cd
./fastqc_0.12.1--hdfd78af_0.sif fastqc ./reads/ecoli_*.fastq.gz
```

This will result in `html` files in the directory `./reads`. These are quality reports for the sequence reads. If you'd like to view them, you can download them with `scp` or e.g. [FileZilla](https://filezilla-project.org/), and view them with your local browser.
