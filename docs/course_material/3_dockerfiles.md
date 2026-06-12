# Working with Dockerfiles

## Learning outcomes

**After having completed this chapter you will be able to:**

* Build an image based on a Dockerfile
* Use the basic Dockerfile syntax
* Change the default command of an image and validate the change
* Map ports to a container to display interactive content through a browser

## Material

* [Official `Dockerfile` reference](https://docs.docker.com/engine/reference/builder/)
* [Ten simple rules for writing Dockerfiles](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1008316)

## Exercises

!!! failure "Important"
    This part should be done locally **on your laptop**.

To make your images shareable and adjustable, it's good practice to work with a `Dockerfile`. Think of a `Dockerfile` as a script with a set of instructions to build your image from an existing image.

### Basic `Dockerfile`

You can generate an image from a `Dockerfile` using the command `docker build`. A `Dockerfile` has its own syntax for giving instructions. Luckily, they are rather simple. The script always contains a line starting with `FROM` that takes the image name from which the new image will be built. After that you usually want to run some commands to e.g. configure and/or install software. The instruction to run these commands during building starts with `RUN`.  In our `figlet` example that would be:

```dockerfile
FROM ubuntu:jammy-20250415.1
RUN apt update
RUN apt install figlet
```

!!! note "On writing reproducible `Dockerfiles`"
    At the `FROM` statement in the the above `Dockerfile` you see that we have added a specific tag to the image (i.e. `jammy-20250415.1`). We could also have written:

    ```dockerfile
    FROM ubuntu
    RUN apt update
    RUN apt install figlet
    ```

    This will automatically pull the image with the tag `latest`. However, if the maintainer of the `ubuntu` images decides to tag another `ubuntu` version as `latest`, rebuilding with the above `Dockerfile` will not give you the same result. Therefore it's always good practice to add the (stable) tag to the image in a `Dockerfile`. More rules on making your `Dockerfiles` more reproducible [here](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1008316).

**Exercise:** Create a file on your computer called `Dockerfile`, and paste the above instruction lines in that file. Make the directory containing the `Dockerfile` your current directory. Build a new image based on that `Dockerfile` with:

=== "x86_64 / AMD64"
    ```sh
    docker build .
    ```
=== "ARM64 (MacOS M1 chip)"
    ```sh
    export DOCKER_DEFAULT_PLATFORM=linux/amd64
    docker build .
    ```

!!! warning "If using an Apple M chip (newer Macs)"
    If you are using a computer with an Apple M chip, you have the less common ARM system architecture, which can limit transferability of images to (more common) `x86_64/AMD64` machines. When building images on a Mac with an M chip (especially if you have sharing in mind), it's best to set the `DOCKER_DEFAULT_PLATFORM` to `linux/amd64` with `export DOCKER_DEFAULT_PLATFORM=linux/amd64`. 

!!! note "The argument of `docker build`"
    The command `docker build` takes a directory as input (providing `.` means the current directory). This directory should contain the `Dockerfile`, but it can also contain more of the build context, e.g. (Python, R, shell) scripts that are required to build the image.

What has happened? What is the name of the build image?

??? success "Answer"
    A new image was created based on the `Dockerfile`. You can check it with: `docker image ls`, which gives something like:

    ```
    REPOSITORY                        TAG       IMAGE ID       CREATED             SIZE
    <none>                            <none>    92c980b09aad   7 seconds ago       101MB
    ubuntu-figlet                     latest    e08b999c7978   About an hour ago   101MB
    ubuntu                            latest    f63181f19b2f   30 hours ago        72.9MB
    ```

    It has created an image without a name or tag. That's a bit inconvenient.

**Exercise:** Build a new image with a specific name. You can do that with adding the option `-t` to `docker build`. Before that, remove the nameless image.

!!! info
    An image without a name is usually a "dangling image". You can remove those with `docker image prune`.

??? success "Answer"
    Remove the nameless image with `docker image prune`.

    After that, rebuild an image with a name:

    === "x86_64 / AMD64"
        ```sh
        docker build -t ubuntu-figlet:v2 .
        ```
    === "ARM (MacOS M1 chip)"
        ```sh
        export DOCKER_DEFAULT_PLATFORM=linux/amd64
        docker build -t ubuntu-figlet:v2 .
        ```

### Using `CMD`

As you might remember the second positional argument of `docker run` is a command (i.e. `docker run IMAGE [CMD]`). If you leave it empty, it uses the default command. You can change the default command in the `Dockerfile` with an instruction starting with `CMD`. For example:

```dockerfile
FROM ubuntu:jammy-20250415.1
RUN apt update
RUN apt install figlet
CMD figlet My image works!
```

**Exercise:** Build a new image based on the above `Dockerfile`. Can you validate the change using `docker image inspect`? Can you overwrite this default with `docker run`?

??? success "Answer"
    Copy the new line to your `Dockerfile`, and build the new image like this:

    === "x86_64 / AMD64"
        ```sh
        docker build -t ubuntu-figlet:v3 .
        ```
    === "ARM64 (MacOS M1 chip)"
        ```sh
        export DOCKER_DEFAULT_PLATFORM=linux/amd64
        docker build -t ubuntu-figlet:v3 .
        ```

    The command `docker inspect ubuntu-figlet:v3` will give:

    ```
    "Cmd": [
        "/bin/sh",
        "-c",
        "figlet My image works!"
    ]
    ```

    So the default command (`/bin/bash`) has changed to `figlet My image works!`

    Running the image (with clean-up (`--rm`)):

    ```sh
    docker run --rm ubuntu-figlet:v3
    ```

    Will result in:

    ```
    __  __         _                                                 _        _
    |  \/  |_   _  (_)_ __ ___   __ _  __ _  ___  __      _____  _ __| | _____| |
    | |\/| | | | | | | '_ ` _ \ / _` |/ _` |/ _ \ \ \ /\ / / _ \| '__| |/ / __| |
    | |  | | |_| | | | | | | | | (_| | (_| |  __/  \ V  V / (_) | |  |   <\__ \_|
    |_|  |_|\__, | |_|_| |_| |_|\__,_|\__, |\___|   \_/\_/ \___/|_|  |_|\_\___(_)
           |___/                     |___/
    ```

    And of course you can overwrite the default command:

    ```sh
    docker run --rm ubuntu-figlet:v3 figlet another text
    ```

    Resulting in:

    ```
    _   _                 _            _
    __ _ _ __   ___ | |_| |__   ___ _ __  | |_ _____  _| |_
    / _` | '_ \ / _ \| __| '_ \ / _ \ '__| | __/ _ \ \/ / __|
    | (_| | | | | (_) | |_| | | |  __/ |    | ||  __/>  <| |_
    \__,_|_| |_|\___/ \__|_| |_|\___|_|     \__\___/_/\_\\__|

    ```

!!! note "Two flavours of `CMD`"
    You have seen in the output of `docker inspect` that docker translates the command (i.e. `figlet "my image works!"`) into this: `["/bin/sh", "-c", "figlet 'My image works!'"]`. The notation we used in the `Dockerfile` is the *shell notation* while the notation with the square brackets (`[]`) is the *exec-notation*. You can use both notations in your `Dockerfile`. Altough the *shell notation* is more readable, the *exec notation* is directly used by the image, and therefore less ambiguous.

    A `Dockerfile` with shell notation:

    ```dockerfile
    FROM ubuntu:jammy-20250415.1
    RUN apt update
    RUN apt install figlet
    CMD figlet My image works!
    ```

    A `Dockerfile` with exec notation:

    ```dockerfile
    FROM ubuntu:jammy-20250415.1
    RUN apt update
    RUN apt install figlet
    CMD ["/bin/sh", "-c", "figlet My image works!"]
    ```

**Exercise:** Now push our created image (with a version tag) to Docker Hub. We will use it later for the [`apptainer` exercises](4_apptainer.md).

??? success "Answer"
    ```sh
    docker tag ubuntu-figlet:v3 [USER NAME]/ubuntu-figlet:v3
    docker push [USER NAME]/ubuntu-figlet:v3
    ```

### Build an image for your own script

Often containers are built for a specific purpose. For example, you can use a container to ship all dependencies together with your developed set of scripts/programs. For that you will need to add your scripts to the container. That is quite easily done with the instruction `COPY`. However, in order to make your container more user-friendly, there are several additional instructions that can come in useful. We will treat the most frequently used ones below. Depending on your preference, either choose **R** or **Python** below. 

=== "R"
    In the exercises, we will use a simple script called `test_deseq2.R`. You can download it [here](https://raw.githubusercontent.com/sib-swiss/containers-introduction-training/main/docker/exercise_r_script/test_deseq2.R), or copy-paste it:

    ```R title="test_deseq2.R"
    #!/usr/bin/env Rscript

    # Load packages required for this script
    write("Loading packages required for this script", stderr())
    suppressPackageStartupMessages({
        library(DESeq2)
        library(optparse)
    })

    # Workaround for issue 112: https://github.com/thelovelab/DESeq2/issues/112
    # This can probably be removed in the future
    setOldClass("ExpData")

    # Load dependency packages for testing installations
    write("Loading dependency packages for testing installations", stderr())
    suppressPackageStartupMessages({
        library(apeglm)
        library(IHW)
        library(limma)
        library(data.table)
        library(ggplot2)
        library(ggrepel)
        library(pheatmap)
        library(RColorBrewer)
        library(scales)
        library(stringr)
    })

    # Create parsing options list
    option_list <- list(
        make_option(c("--rows"),
            type = "integer",
            help = "Number of rows in dummy matrix [default = %default]",
            default = 100
        )
    )

    # Implement parser with optparse
    opt_parser <- OptionParser(
        option_list = option_list,
        description = "Runs DESeq2 on dummy data"
    )

    # Parse options with optparse
    opt <- parse_args(opt_parser)

    # Create a random dummy count matrix
    cnts <- matrix(rnbinom(n = opt$row * 10, mu = 100, size = 1 / 0.5), ncol = 10)
    cond <- factor(rep(1:2, each = 5))

    # Object construction
    dds <- DESeqDataSetFromMatrix(cnts, DataFrame(cond), ~cond)

    # Standard analysis
    dds <- DESeq(dds)
    res <- results(dds)

    # Print results to stdout
    print(res)
    ```
    
    After you have downloaded it, make sure to set the permissions to executable:

    ```sh
    chmod +x test_deseq2.R
    ```
    It is a relatively simple script that runs DESeq2 on a dummy dataset. To execute it with default parameters, you can use:

    ```sh
    ./test_deseq2.R
    ```

    With this command, the dummy dataset will contain 100 rows, the default value as written in the `# Create parsing options list` code block of the script. If you want to run it on a different number of rows, you can use the `--rows` optional argument that specifies the number of rows generated in the input count matrix:

    ```sh
    ./test_deseq2.R --rows 75
    ```

    When running the script, it will return a bunch of messages and at the end an overview of differential gene expression analysis results:

    ```
        baseMean log2FoldChange     lfcSE         stat    pvalue      padj
        <numeric>      <numeric> <numeric>    <numeric> <numeric> <numeric>
    1     66.1249       0.281757  0.727668     0.387206  0.698604  0.989804
    2     76.9682       0.305763  0.619209     0.493796  0.621451  0.989804
    3     64.7843      -0.694525  0.479445    -1.448603  0.147448  0.931561
    4    123.0252       0.631247  0.688564     0.916758  0.359269  0.931561
    5     93.2002      -0.453430  0.686043    -0.660936  0.508653  0.941951
    ...       ...            ...       ...          ...       ...       ...
    96    64.0177    0.757585137  0.682683  1.109718054  0.267121  0.931561
    97   114.3689   -0.580010850  0.640313 -0.905823841  0.365029  0.931561
    98    79.9620    0.000100617  0.612442  0.000164288  0.999869  0.999869
    99    92.6614    0.563514308  0.716109  0.786910869  0.431334  0.939106
    100   96.4410   -0.155268696  0.534400 -0.290547708  0.771397  0.989804
    ```

    From the script you can see it has `DESeq2` and `optparse` as dependencies. If we want to run the script inside a container, we would have to install them. We do this in the `Dockerfile` below, with the following instructions:

    - Use the [r2u base image](https://hub.docker.com/r/rocker/r2u) version jammy
    - Install the package `DESeq2`, `optparse` and some additional packages we will need later on. We perform the installations with `install2.r`, which is a helper command that is present inside most rocker images. More info [here](https://rocker-project.org/use/extending.html#install2.r).
    - Copy the script `test_deseq2.R` to `/opt` inside the container:

    ```dockerfile
    FROM rocker/r2u:jammy

    RUN install2.r \
        DESeq2 \
        optparse \
        apeglm \
        IHW \
        limma \
        data.table \
        ggrepel \
        pheatmap \
        stringr

    COPY test_deseq2.R /opt 
    ```

    !!! note
        In order to use `COPY`, the file that needs to be copied should be in the same directory as the `Dockerfile` or one of its subdirectories.

    !!! note "R image stack"
        The most used R image stack is from the [rocker project](https://rocker-project.org/). It contains many different base images (e.g. with Shiny, Rstudio, tidyverse etc.). The possible installation methods (for example, with `apt` or `install2.r`) depend on the type of image you are using. To understand more about how to install R packages in different containers, check this [cheat sheet](https://raw.githubusercontent.com/sib-swiss/containers-introduction-training/main/r-docker-cheatsheet/r-docker-cheatsheet.pdf), or visit [rocker-project.org](https://rocker-project.org/).

    **Exercise:** Download the `test_deseq2.R` and build the image with `docker build`. Name the image `deseq2`. After that, start an interactive session and execute the script inside the container. 

    !!! info
        Make an interactive session with the options `-i` and `-t` and use `/bin/bash` as the command. 

    ??? success "Answer"
        Build the container:

        === "x86_64 / AMD64"
            ```sh
            docker build -t deseq2 .
            ```
        === "ARM64 (MacOS M1 chip)"
            ```sh
            export DOCKER_DEFAULT_PLATFORM=linux/amd64
            docker build -t deseq2 .
            ```

        Run the container:

        ```sh
        docker run -it --rm deseq2 /bin/bash
        ```

        Inside the container we look up the script:

        ```sh
        cd /opt
        ls
        ```

        This should return `test_deseq2.R`. 

        Now you can execute it from inside the container:

        ```sh
        ./test_deseq2.R --rows 100
        ```

    That's kind of nice. We can ship our R script inside our container. However, we don't want to run it interactively every time. So let's make some changes to make it easy to run it as an executable. For example, we can add `/opt` to the global `$PATH` variable with `ENV`. 

    !!! note "The `$PATH` variable"
        The path variable is a special variable that consists of a list of path seperated by colons (`:`). These paths are searched if you are trying to run an executable. More info this topic at e.g. [wikipedia](https://en.wikipedia.org/wiki/PATH_(variable)). 

    ```dockerfile
    FROM rocker/r2u:jammy

    RUN install2.r \
        DESeq2 \
        optparse \
        apeglm \
        IHW \
        limma \
        data.table \
        ggrepel \
        pheatmap \
        stringr

    COPY test_deseq2.R /opt 

    ENV PATH=/opt:$PATH
    ```

    !!! note
        The `ENV` instruction can be used to set any variable. 

    **Exercise**: Rebuild the image and start an interactive bash session inside the new image. Is the path variable updated? (i.e. can we execute `test_deseq2.R` from anywhere?)

    ??? success "Answer"
        After re-building we start an interactive session:

        ```sh
        docker run -it --rm deseq2 /bin/bash
        ```

        The path is upated, `/opt` is appended to the beginning of the variable:

        ```sh
        echo $PATH
        ```

        returns:

        ```
        /opt:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
        ```

        Now you can try to execute it from the root directory (or any other):

        ```sh
        test_deseq2.R
        ```

    Instead of starting an interactive session with `/bin/bash` we can now more easily run the script non-interactively:

    ```sh
    docker run --rm deseq2 test_deseq2.R --rows 100
    ```

    Now it will directly print the output of `test_deseq2.R` to stdout. 

    When you pack your script inside a container, you are building a container specifically for your script, meaning you _almost_ want the container to behave as the program itself. In order to do that, you can use `ENTRYPOINT`. `ENTRYPOINT` is similar to `CMD`, but has two important differences:

    * `ENTRYPOINT` **can not be overwritten by positional arguments** (i.e. `docker run image [CMD]`), but has to be overwritten by `--entrypoint`.
    * The positional arguments (or `CMD`) are pasted to the `ENTRYPOINT` command. This means that you can use `ENTRYPOINT` as the executable and the positional arguments (or `CMD`) as the options. 

    Let's try it out:

    ```dockerfile
    FROM rocker/r2u:jammy

    RUN install2.r \
        DESeq2 \
        optparse \
        apeglm \
        IHW \
        limma \
        data.table \
        ggrepel \
        pheatmap \
        stringr

    COPY test_deseq2.R /opt 

    ENV PATH=/opt:$PATH

    # Note that if you want to be able to combine the two, both ENTRYPOINT and CMD need to written in the exec form
    ENTRYPOINT ["test_deseq2.R"]

    # Default option (if positional arguments are not specified)
    CMD ["--rows", "100"]
    ```

    **Exercise**: Re-build, and run the container non-interactively without any positional arguments. After that, try to pass a different number of rows to `--rows`. How do the commands look?

    ??? success "Answer"
        Just running the container non-interactively would be:

        ```sh
        docker run --rm deseq2
        ```

        Passing a different argument (i.e. overwriting `CMD`) would be:

        ```sh
        docker run --rm deseq2 --rows 200
        ```

        Here, the container behaves as the executable itself to which you can pass arguments. 

    !!! note 
        You can overwrite `ENTRYPOINT` with `--entrypoint` as an argument to `docker run`. 

    **Exercise**: Push the image to Docker Hub, so we can use it later with the Apptainer exercises.

    ??? success "Answer"

        Pushing it to Docker Hub:

        ```sh
        docker tag deseq2 [USER NAME]/deseq2:v1
        docker push [USER NAME]/deseq2:v1
        ```

    <h3> Extra: Get information on your image with `docker inspect` </h3>

    We have used `docker inspect` already in the previous chapter to find the default `Cmd` of the ubuntu image. However we can get more info on the image: e.g. the entrypoint, environmental variables, cmd, workingdir etc., you can use the `Config` record from the output of `docker inspect`. For our image this looks like:

    ```yaml
    [
        {
            "Id": "sha256:d19b445510187ef4d0ac77b1842ab906433fb07d8faba954372afa18ea8a51e2",
            "RepoTags": [
                "deseq2:latest"
            ],
            "RepoDigests": [],
            "Comment": "buildkit.dockerfile.v0",
            "Created": "2026-06-02T14:38:18.240243842+02:00",
            "Config": {
                "Env": [
                    "PATH=/opt:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin",
                    "LC_ALL=en_US.UTF-8",
                    "LANG=en_US.UTF-8",
                    "DEBIAN_FRONTEND=noninteractive",
                    "TZ=UTC"
                ],
                "Entrypoint": [
                    "test_deseq2.R"
                ],
                "Cmd": [
                    "--rows",
                    "100"
                ],
                "Labels": {
                    "maintainer": "Dirk Eddelbuettel <edd@debian.org>",
                    "org.label-schema.license": "GPL-2.0",
                    "org.label-schema.vcs-url": "https://github.com/rocker-org/",
                    "org.label-schema.vendor": "Rocker Project",
                    "org.opencontainers.image.version": "22.04"
                },
                "ArgsEscaped": true
            },
            "Architecture": "amd64",
            "Os": "linux",
            "Size": 1103316046,
            # ...
            "Metadata": {
                "LastTagTime": "2026-06-02T14:41:00.742210773+02:00"
            }
        }
    ]
    ```

    <h3> Extra: Adding metadata to your image </h3>

    You can annotate your `Dockerfile` and the image by using the instruction `LABEL`. You can give it any key and value with `<key>=<value>`. However, it is recommended to use the [Open Container Initiative (OCI) keys](https://github.com/opencontainers/image-spec/blob/v1.0.1/annotations.md).

    **Exercise**: Annotate our `Dockerfile` with the OCI keys on the creation date, author and description. After that, check whether this has been passed to the actual image with `docker inspect`. 

    !!! note
        You can type `LABEL` for each key-value pair, but you can also have it on one line by seperating the key-value pairs by a space, e.g.:

        ```dockerfile
        LABEL key_x="value_x" key_y="value_y"
        ```

    ??? success "Answer"

        The `Dockerfile` would look like:

        ```dockerfile
        FROM rocker/r2u:jammy

        LABEL org.opencontainers.image.created="2023-04-12" \
                org.opencontainers.image.authors="Geert van Geest" \
                org.opencontainers.image.description="Container with DESeq2 and friends"

        RUN install2.r \
            DESeq2 \
            optparse \
            apeglm \
            IHW \
            limma \
            data.table \
            ggrepel \
            pheatmap \
            stringr

        WORKDIR /opt

        COPY test_deseq2.R .

        ENV PATH=/opt:$PATH

        # Note that if you want to be able to combine the two, both ENTRYPOINT and CMD need to written in the exec form
        ENTRYPOINT ["test_deseq2.R"]

        # Default option (if positional arguments are not specified)
        CMD ["--rows", "100"]

        ```

        The `Config` record in the output of `docker inspect` was updated with:

        ```yaml
        "Labels": {
            "maintainer": "Dirk Eddelbuettel <edd@debian.org>",
            "org.label-schema.license": "GPL-2.0",
            "org.label-schema.vcs-url": "https://github.com/rocker-org/",
            "org.label-schema.vendor": "Rocker Project",
            "org.opencontainers.image.authors": "Geert van Geest",
            "org.opencontainers.image.created": "2024-10-12",
            "org.opencontainers.image.description": "Container with DESeq2 and friends",
            "org.opencontainers.image.version": "22.04"
        },
        ```

    <h3> Extra: Building an image with a browser interface </h3>

    In this exercise, we will use a different base image (`rocker/rstudio:4`), and we will install the same packages. [Rstudio server](https://posit.co/download/rstudio-server/) is a nice browser interface that you can use for a.o. (Aspect Oriented) programming in R. With the image we are creating, we will be able to run Rstudio server inside a container. Check out the `Dockerfile`:

    ```dockerfile
    FROM rocker/rstudio:4

    LABEL org.opencontainers.image.created="2024-10-12" \
        org.opencontainers.image.authors="Geert van Geest" \
        org.opencontainers.image.description="Container with DESeq2 and friends"

    RUN apt update && apt install -y \
        fonts-noto \
        libcurl4-openssl-dev \
        libfontconfig1-dev \
        libfreetype6-dev \
        libfribidi-dev \
        libharfbuzz-dev \
        libjpeg-dev \
        libpng-dev \
        libtiff5-dev \
        libxml2-dev \
        tzdata \
        && apt clean && apt autoremove && apt purge \
        && rm -rf /var/lib/apt/lists/*

    RUN install2.r \
        BiocManager

    RUN R -q -e "BiocManager::install(c( \
      'DESeq2', \
      'optparse', \
      'apeglm', \
      'IHW', \
      'limma', \
      'data.table', \
      'ggrepel', \
      'pheatmap', \
      'stringr' \
      ))"

    WORKDIR /opt

    COPY test_deseq2.R .

    ENV PATH=/opt:$PATH
    ```

    This will create an image from the existing `rstudio` image. It first installs required system libraries (`libcurl4-openssl-dev`, `libxml2-dev`...) via `apt`, which are needed to compile certain Bioconductor dependencies. It then installs `BiocManager` using `install2.r` (the helper script provided by rocker images) and next uses it to install several R packages (`DESeq2`, `optparse`, `apeglm`...) via an inline R command. Note that R packages themselves cannot be installed with `apt` in `rocker/rstudio` images, because R is compiled from source rather than installed from the system package manager. More information on how to install R packages in R containers in this [cheat sheet](https://raw.githubusercontent.com/sib-swiss/containers-introduction-training/main/r-docker-cheatsheet/r-docker-cheatsheet.pdf), or visit [rocker-project.org](https://rocker-project.org/).

    **Exercise:** Build an image based on this `Dockerfile` and give it a meaningful name.

    ??? success "Answer"
        === "x86_64 / AMD64"
            ```sh
            docker build -t rstudio-server .
            ```
        === "ARM64 (MacOS M1 chip)"
            ```sh
            export DOCKER_DEFAULT_PLATFORM=linux/amd64
            docker build -t rstudio-server .
            ```

    You can now run a container from the image. However, you will have to tell Docker where to publish port 8787 from the Docker container with `-p [HOSTPORT:CONTAINERPORT]`. We choose to publish it to the same port number:

    ```sh
    docker run --rm -it -p 8787:8787 rstudio-server
    ```

    !!! note "Networking"
        More info on Docker container networking [here](https://docs.docker.com/config/containers/container-networking/)

    By running the above command, a container will be started exposing Rstudio server at port 8787 on localhost. You can approach the instance of Rstudio by typing `localhost:8787` in your browser. You will be asked for a username (`rstudio`) and a password. You can find it in the terminal from which you have started the container.

    We can make this even more interesting by mounting a local directory to the container running the Rstudio image:

    ```sh
    docker run \
    -it \
    --rm \
    -p 8787:8787 \
    --mount type=bind,source=$(pwd),target=/working_dir/ \
    rstudio-server
    ```

    By doing this you have a completely isolated and shareable R environment running Rstudio server, but with your local files available to it. Pretty neat right? 

=== "Python"
    In the exercises, we will use a simple script called `test_pydeseq2.py`. You can download it [here](https://raw.githubusercontent.com/sib-swiss/containers-introduction-training/main/docker/exercise_python_script/test_pydeseq2.py), or copy-paste it:

    ```python title="test_pydeseq2.py"
    #!/usr/bin/env python3

    # Import libraries
    import sys

    # Load packages required for this script
    print('Loading packages required for this script', file=sys.stderr)

    # ruff: disable[E402]
    import argparse
    import numpy as np
    import pandas as pd
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    # ruff: enable[E402]

    # Create argument parser
    parser = argparse.ArgumentParser(description='Runs PyDESeq2 on dummy data')
    parser.add_argument(
        '--rows',
        type=int,
        default=100,
        help='Number of rows (genes) in dummy matrix [default: %(default)s]',
    )

    # Parse arguments
    print('Parsing arguments', file=sys.stderr)
    args = parser.parse_args()

    # Create a random dummy count matrix (samples x genes)
    print('Creating dummy count matrix', file=sys.stderr)
    counts_df = pd.DataFrame(
        np.random.negative_binomial(n=2, p=0.02, size=(10, args.rows)),
        index=[f'sample{i + 1}' for i in range(10)],
        columns=[f'gene{i + 1}' for i in range(args.rows)],
    )

    # Metadata / condition factor
    metadata = pd.DataFrame(
        {'condition': ['A'] * 5 + ['B'] * 5}, index=counts_df.index
    )

    # Run DESeq2 analysis
    print('Running PyDESeq2', file=sys.stderr)
    dds = DeseqDataSet(
        counts=counts_df, metadata=metadata, design_factors='condition'
    )
    dds.deseq2()

    # Get results
    print('Computing results...', file=sys.stderr)
    stat_res = DeseqStats(dds, contrast=["condition", "B", "A"])
    stat_res.summary()
    ```

    After you have downloaded it, make sure to set the permissions to executable:

    ```sh
    chmod +x test_pydeseq2.py
    ```

    It is a relatively simple script that runs PyDESeq2 on a dummy dataset. To execute it with default parameters, you can use:

    ```sh
    ./test_pydeseq2.py
    ```

    With this command, the dummy dataset will contain 100 rows (genes), the default value as written in the `# Create argument parser` code block of the script. If you want to run it on a different number of rows, you can use the `--rows` optional argument that specifies the number of rows generated in the input count matrix:

    ```sh
    ./test_pydeseq2.py --rows 75
    ```

    When running the script, it will return a bunch of messages and at the end an overview of differential gene expression analysis results:

    ```
    Log2 fold change & Wald test p-value: condition B vs A
               baseMean  log2FoldChange      lfcSE       stat     pvalue       padj
    gene1    116.049570       -0.329827   0.517247  -0.637659   0.523695   0.863125
    gene2     96.296548        1.081284   0.697623   1.549955   0.121152   0.583151
    gene3     62.700330       -0.336880   0.655994  -0.513541   0.607573   0.863779
    gene4     80.397130       -0.571927   0.792223  -0.721928   0.470339   0.863125
    gene5     79.193224        0.300298   0.565743   0.530802   0.595556   0.863125
    ...             ...             ...        ...        ...        ...        ...
    gene96    86.569540       -1.052608   0.604778  -1.740485   0.081774   0.583151
    gene97   100.885090       -1.025317   0.711971  -1.440109   0.149837   0.607309
    gene98   104.836128       -0.494590   0.766652  -0.645130   0.518843   0.863125
    gene99    72.800740       -1.096534   0.576206  -1.903025   0.057037   0.471511
    gene100   71.921340       -0.508122   0.630407  -0.806021   0.420231   0.863125
    ```

    From the script you can see it has `numpy`, `pandas` and `pydeseq2` as dependencies. If we want to run the script inside a container, we would have to install them. We do this in the `Dockerfile` below, with the following instructions:

    - Use the [official Python base image](https://hub.docker.com/_/python) version 3.13
    - Install the package `pydeseq2` with `pip` (this will also pull in `numpy` and `pandas` as dependencies)
    - Copy the script `test_pydeseq2.py` to `/opt` inside the container:

    ```dockerfile
    FROM python:3.13

    RUN pip install pydeseq2

    COPY test_pydeseq2.py /opt
    ```

    !!! note
        In order to use `COPY`, the file that needs to be copied should be in the same directory as the `Dockerfile` or one of its subdirectories.

    **Exercise:** Download the `test_pydeseq2.py` and build the image with `docker build`. Name the image `pydeseq2`. After that, start an interactive session and execute the script inside the container.

    !!! info
        Make an interactive session with the options `-i` and `-t` and use `/bin/bash` as the command. 

    ??? success "Answer"
        Build the container:

        === "x86_64 / AMD64"
            ```sh
            docker build -t pydeseq2 .
            ```
        === "ARM64 (MacOS M1 chip)"
            ```sh
            export DOCKER_DEFAULT_PLATFORM=linux/amd64
            docker build -t pydeseq2 .
            ```

        Run the container:

        ```sh
        docker run -it --rm pydeseq2 /bin/bash
        ```

        Inside the container we look up the script:

        ```sh
        cd /opt
        ls
        ```

        This should return `test_pydeseq2.py`.

        Now you can execute it from inside the container:

        ```sh
        ./test_pydeseq2.py --rows 100
        ```

    That's kind of nice. We can ship our Python script inside our container. However, we don't want to run it interactively every time. So let's make some changes to make it easy to run it as an executable. For example, we can add `/opt` to the global `$PATH` variable with `ENV`.

    !!! note "The `$PATH` variable"
        The path variable is a special variable that consists of a list of path seperated by colons (`:`). These paths are searched if you are trying to run an executable. More info this topic at e.g. [wikipedia](https://en.wikipedia.org/wiki/PATH_(variable)). 

    ```dockerfile
    FROM python:3.13

    RUN pip install pydeseq2

    COPY test_pydeseq2.py /opt

    ENV PATH=/opt:$PATH
    ```

    !!! note
        The `ENV` instruction can be used to set any variable. 

    **Exercise**: Rebuild the image and start an interactive bash session inside the new image. Is the path variable updated? (i.e. can we execute `test_pydeseq2.py` from anywhere?)

    ??? success "Answer"
        After re-building we start an interactive session:

        ```sh
        docker run -it --rm pydeseq2 /bin/bash
        ```

        The path is upated, `/opt` is appended to the beginning of the variable:

        ```sh
        echo $PATH
        ```

        returns:

        ```
        /opt:/usr/local/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
        ```

        Now you can try to execute it from the root directory (or any other):

        ```sh
        test_pydeseq2.py
        ```

    Instead of starting an interactive session with `/bin/bash` we can now more easily run the script non-interactively:

    ```sh
    docker run --rm pydeseq2 test_pydeseq2.py --rows 100
    ```

    Now it will directly print the output of `test_pydeseq2.py` to stdout.

    When you pack your script inside a container, you are building a container specifically for your script, meaning you _almost_ want the container to behave as the program itself. In order to do that, you can use `ENTRYPOINT`. `ENTRYPOINT` is similar to `CMD`, but has two important differences:

    * `ENTRYPOINT` **can not be overwritten by positional arguments** (i.e. `docker run image [CMD]`), but has to be overwritten by `--entrypoint`.
    * The positional arguments (or `CMD`) are pasted to the `ENTRYPOINT` command. This means that you can use `ENTRYPOINT` as the executable and the positional arguments (or `CMD`) as the options. 

    Let's try it out:

    ```dockerfile
    FROM python:3.13

    RUN pip install pydeseq2

    COPY test_pydeseq2.py /opt

    ENV PATH=/opt:$PATH

    # Note that if you want to be able to combine the two, both ENTRYPOINT and CMD need to written in the exec form
    ENTRYPOINT ["test_pydeseq2.py"]

    # Default option (if positional arguments are not specified)
    CMD ["--rows", "100"]
    ```

    **Exercise**: Re-build, and run the container non-interactively without any positional arguments. After that, try to pass a different number of rows to `--rows`. How do the commands look?

    ??? success "Answer"
        Just running the container non-interactively would be:

        ```sh
        docker run --rm pydeseq2
        ```

        Passing a different argument (i.e. overwriting `CMD`) would be:

        ```sh
        docker run --rm pydeseq2 --rows 200
        ```

        Here, the container behaves as the executable itself to which you can pass arguments. 

    !!! note 
        You can overwrite `ENTRYPOINT` with `--entrypoint` as an argument to `docker run`. 

    **Exercise**: Push the image to Docker Hub, so we can use it later with the Apptainer exercises.

    ??? success "Answer"

        Pushing it to Docker Hub:

        ```sh
        docker tag pydeseq2 [USER NAME]/pydeseq2:v1
        docker push [USER NAME]/pydeseq2:v1
        ```

    <h3> Extra: Get information on your image with `docker inspect` </h3>

    We have used `docker inspect` already in the previous chapter to find the default `Cmd` of the ubuntu image. However we can get more info on the image: e.g. the entrypoint, environmental variables, cmd, workingdir etc., you can use the `Config` record from the output of `docker inspect`. For our image this looks like:

    ```yaml
    [
        {
            "Id": "sha256:0051b8ea8110afef13bf4e5febea1dc890e7dab817842b39cfabb1c07221b47e",
            "RepoTags": [
                "pydeseq2:latest"
            ],
            "RepoDigests": [],
            "Comment": "buildkit.dockerfile.v0",
            "Created": "2026-06-02T13:52:10.621027918+02:00",
            "Config": {
                "Env": [
                    "PATH=/opt:/usr/local/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin",
                    "GPG_KEY=7169605F62C751356D054A26A821E680E5FA6305",
                    "PYTHON_VERSION=3.13.13",
                    "PYTHON_SHA256=2ab91ff401783ccca64f75d10c882e957bdfd60e2bf5a72f8421793729b78a71"
                ],
                "Entrypoint": [
                    "test_pydeseq2.py"
                ],
                "Cmd": [
                    "--rows",
                    "100"
                ],
                "ArgsEscaped": true
            },
            "Architecture": "amd64",
            "Os": "linux",
            "Size": 1721829929,
            # ...
            "Metadata": {
                "LastTagTime": "2026-06-02T13:53:14.094758235+02:00"
            }
        }
    ]
    ```

    <h3> Extra: Adding metadata to your image </h3>

    You can annotate your `Dockerfile` and the image by using the instruction `LABEL`. You can give it any key and value with `<key>=<value>`. However, it is recommended to use the [Open Container Initiative (OCI) keys](https://github.com/opencontainers/image-spec/blob/v1.0.1/annotations.md).

    **Exercise**: Annotate our `Dockerfile` with the OCI keys on the creation date, author and description. After that, check whether this has been passed to the actual image with `docker inspect`. 

    !!! note
        You can type `LABEL` for each key-value pair, but you can also have it on one line by seperating the key-value pairs by a space, e.g.:

        ```dockerfile
        LABEL key_x="value_x" key_y="value_y"
        ```

    ??? success "Answer"

        The `Dockerfile` would look like:

        ```dockerfile
        FROM python:3.13

        LABEL org.opencontainers.image.created="2026-06-02" \
            org.opencontainers.image.authors="Antonin Thiébaut" \
            org.opencontainers.image.description="Container with PyDESeq2 and friends"

        RUN pip install pydeseq2

        WORKDIR /opt

        COPY test_pydeseq2.py .

        ENV PATH=/opt:$PATH

        # Note that if you want to be able to combine the two, both ENTRYPOINT and CMD need to written in the exec form
        ENTRYPOINT ["test_pydeseq2.py"]

        # Default option (if positional arguments are not specified)
        CMD ["--rows", "100"]

        ```

        The `Config` record in the output of `docker inspect` was updated with:

        ```yaml
        "Labels": {
            "org.opencontainers.image.authors": "Antonin Thiébaut",
            "org.opencontainers.image.created": "2026-06-02",
            "org.opencontainers.image.description": "Container with PyDESeq2 and friends"
        },
        ```

    <h3> Extra: Building an image with a browser interface </h3>

    In this exercise, we will use a different base image ([Jupyter docker image stack](https://jupyter-docker-stacks.readthedocs.io/en/latest/using/common.html)). [JupyterLab](https://jupyter.org/) is a nice browser interface that you can use for a.o. (Aspect Oriented) programming in Python. With the image we are creating, we will be able to run JupyterLab inside a container. Check out the `Dockerfile`:

    ```dockerfile
    FROM jupyter/base-notebook:x86_64-ubuntu-22.04

    LABEL org.opencontainers.image.created="2026-06-02" \
        org.opencontainers.image.authors="Antonin Thiébaut" \
        org.opencontainers.image.description="Container with PyDESeq2 and JupyterLab"

    RUN pip install pydeseq2

    WORKDIR /opt

    COPY test_pydeseq2.py .

    ENV PATH=/opt:$PATH
    ```

    This will create an image from the existing `Python` image. It will also install `JupyterLab` with `pip`. As a default command it starts a Jupyter notebook at port 8888.

    !!! note "Ports"
        We have specified here that JupyterLab should use port 8888. However, this **inside** the container. We can not connect to it yet with our browser.

    **Exercise:** Build an image based on this `Dockerfile` and give it a meaningful name.

    ??? success "Answer"
        === "x86_64 / AMD64"
            ```sh
            docker build -t jupyter-lab .
            ```
        === "ARM64 (MacOS M1 chip)"
            ```sh
            export DOCKER_DEFAULT_PLATFORM=linux/amd64
            docker build -t jupyter-lab .
            ```

    You can now run a container from the image. However, you will have to tell Docker where to publish port 8888 from the Docker container with `-p [HOSTPORT:CONTAINERPORT]`. We choose to publish it to the same port number:

    ```sh
    docker run --rm -it -p 8888:8888 jupyter-lab
    ```

    !!! note "Networking"
        More info on Docker container networking [here](https://docs.docker.com/config/containers/container-networking/)

    By running the above command, a container will be started exposing JupyterHub at port 8888 on `localhost`. You can approach the instance of JupyterHub by typing `localhost:8888` in your browser. You will be asked for a token; you can find it in the terminal from which you have started the container. You can also find the direct URL to access the notebook without token in the same terminal; it will look like this:

    `http://127.0.0.1:8888/lab?token=e3434591f62012afb9de9d9b370f25540a8d859272dec55d`.

    We can make this even more interesting by mounting a local directory to the container running the JupyterLab image:

    ```sh
    docker run \
    -it \
    --rm \
    -p 8888:8888 \
    --mount type=bind,source=$(pwd),target=/working_dir/ \
    jupyter-lab
    ```

    By doing this you have a completely isolated and shareable Python environment running JupyterLab, but with your local files available to it. Pretty neat right?

    !!! note "Jupyter images"
        Jupyter has a wide range of pre-built images available [here](https://hub.docker.com/u/jupyter).
