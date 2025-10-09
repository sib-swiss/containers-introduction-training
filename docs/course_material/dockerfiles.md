## Learning outcomes

**After having completed this chapter you will be able to:**

* Build an image based on a dockerfile
* Use the basic dockerfile syntax
* Change the default command of an image and validate the change
* Map ports to a container to display interactive content through a browser

## Material

* [Official `Dockerfile` reference](https://docs.docker.com/engine/reference/builder/)
* [Ten simple rules for writing dockerfiles](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1008316)

## Exercises

To make your images shareable and adjustable, it's good practice to work with a `Dockerfile`. This is a script with a set of instructions to build your image from an existing image.

### Basic `Dockerfile`

You can generate an image from a `Dockerfile` using the command `docker build`. A `Dockerfile` has its own syntax for giving instructions. Luckily, they are rather simple. The script always contains a line starting with `FROM` that takes the image name from which the new image will be built. After that you usually want to run some commands to e.g. configure and/or install software. The instruction to run these commands during building starts with `RUN`.  In our `figlet` example that would be:

```dockerfile
FROM ubuntu:jammy-20250415.1
RUN apt-get update
RUN apt-get install figlet
```

!!! note "On writing reproducible `Dockerfiles`"
    At the `FROM` statement in the the above `Dockerfile` you see that we have added a specific tag to the image (i.e. `jammy-20250415.1`). We could also have written:

    ```dockerfile
    FROM ubuntu
    RUN apt-get update
    RUN apt-get install figlet
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
    The command `docker build` takes a directory as input (providing `.` means the current directory). This directory should contain the `Dockerfile`, but it can also contain more of the build context, e.g. (python, R, shell) scripts that are required to build the image.

What has happened? What is the name of the build image?

??? done "Answer"
    A new image was created based on the `Dockerfile`. You can check it with: `docker image ls`, which gives something like:

    ```
    REPOSITORY                        TAG       IMAGE ID       CREATED             SIZE
    <none>                            <none>    92c980b09aad   7 seconds ago       101MB
    ubuntu-figlet                     latest    e08b999c7978   About an hour ago   101MB
    ubuntu                            latest    f63181f19b2f   30 hours ago        72.9MB
    ```

    It has created an image without a name or tag. That's a bit inconvenient.

**Exercise:** Build a new image with a specific name. You can do that with adding the option `-t` to `docker build`. Before that, remove the nameless image.

!!! hint
    An image without a name is usually a "dangling image". You can remove those with `docker image prune`.

??? done "Answer"
    Remove the nameless image with `docker image prune`.

    After that, rebuild an image with a name:

    === "x86_64 / AMD64"
        ```sh
        docker build -t ubuntu-figlet:v2 .
        ```
    === "ARM (MacOS M1 chip)"
        ```sh
        DOCKER_DEFAULT_PLATFORM=linux/amd64
        docker build -t ubuntu-figlet:v2 .
        ```

### Using `CMD`

As you might remember the second positional argument of `docker run` is a command (i.e. `docker run IMAGE [CMD]`). If you leave it empty, it uses the default command. You can change the default command in the `Dockerfile` with an instruction starting with `CMD`. For example:

```dockerfile
FROM ubuntu:jammy-20250415.1
RUN apt-get update
RUN apt-get install figlet
CMD figlet My image works!
```

**Exercise:** Build a new image based on the above `Dockerfile`. Can you validate the change using `docker image inspect`? Can you overwrite this default with `docker run`?

??? done "Answer"
    Copy the new line to your `Dockerfile`, and build the new image like this:

    === "x86_64 / AMD64"
        ```sh
        docker build -t ubuntu-figlet:v3 .
        ```
    === "ARM64 (MacOS M1 chip)"
        ```sh
        DOCKER_DEFAULT_PLATFORM=linux/amd64
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
    RUN apt-get update
    RUN apt-get install figlet
    CMD figlet My image works!
    ```

    A `Dockerfile` with exec notation:

    ```dockerfile
    FROM ubuntu:jammy-20250415.1
    RUN apt-get update
    RUN apt-get install figlet
    CMD ["/bin/sh", "-c", "figlet My image works!"]
    ```

**Exercise:** Now push our created image (with a version tag) to docker hub. We will use it later for the [`apptainer` exercises](apptainer.md).

??? done "Answer"
    ```sh
    docker tag ubuntu-figlet:v3 [USER NAME]/ubuntu-figlet:v3
    docker push [USER NAME]/ubuntu-figlet:v3
    ```

### Build an image for your own script

Often containers are built for a specific purpose. For example, you can use a container to ship all dependencies together with your developed set of scripts/programs. For that you will need to add your scripts to the container. That is quite easily done with the instruction `COPY`. However, in order to make your container more user-friendly, there are several additional instructions that can come in useful. We will treat the most frequently used ones below. Depending on your preference, either choose **R** or **Python** below. 

=== "R"
    In the exercises will use a simple script called `test_deseq2.R`. You can download it [here](https://raw.githubusercontent.com/sib-swiss/containers-introduction-training/main/docker/exercise_r_script/test_deseq2.R), or copy-paste it:

    ```R title="test_deseq2.R"
    #!/usr/bin/env Rscript

    # load packages required for this script
    write("Loading packages required for this script", stderr())
    suppressPackageStartupMessages({
        library(DESeq2)
        library(optparse)
    })

    # workaround for issue 112: https://github.com/thelovelab/DESeq2/issues/112
    # this can probably be removed in the future
    setOldClass("ExpData")

    # load dependency packages for testing installations
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

    # parse options with optparse
    option_list <- list(
        make_option(c("--rows"),
            type = "integer",
            help = "Number of rows in dummy matrix [default = %default]",
            default = 100
        )
    )

    opt_parser <- OptionParser(
        option_list = option_list,
        description = "Runs DESeq2 on dummy data"
    )
    opt <- parse_args(opt_parser)

    # create a random dummy count matrix
    cnts <- matrix(rnbinom(n = opt$row * 10, mu = 100, size = 1 / 0.5), ncol = 10)
    cond <- factor(rep(1:2, each = 5))

    # object construction
    dds <- DESeqDataSetFromMatrix(cnts, DataFrame(cond), ~cond)

    # standard analysis
    dds <- DESeq(dds)
    res <- results(dds)

    # print results to stdout
    print(res)


    ```
    
    
    After you have downloaded it, make sure to set the permissions to executable:

    ```sh
    chmod +x test_deseq2.R
    ```
    It is a relatively simple script that runs DESeq2 on a dummy dataset. An example for execution would be:

    ```sh
    ./test_deseq2.R --rows 75
    ```

    Giving a list of results from DESeq2 on a dummy dataset with 75 rows.

    Here, `--rows` is an optional argument that specifies the number of rows generated in the input count matrix. When running the script, it will return a bunch of messages and at the end an overview of differential gene expression analysis results:

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

    From the script you can see it has `DESeq2` and `optparse` as dependencies. If we want to run the script inside a container, we would have to install them. We do this in the `Dockerfile` below. We give it the following instructions:

    - use the [r2u base image](https://hub.docker.com/r/rocker/r2u) version jammy
    - install the package `DESeq2`, `optparse` and some additional packages we will need later on. We perform the installations with `install2.r`, which is a helper command that is present inside most rocker images. More info [here](https://rocker-project.org/use/extending.html#install2.r). 
    - copy the script `test_deseq2.R` to `/opt` inside the container:

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
        In order to use `COPY`, the file that needs to be copied needs to be in the same directory as the `Dockerfile` or one of its subdirectories.

    !!! note "R image stack"
        The most used R image stack is from the [rocker project](https://rocker-project.org/). It contains many different base images (e.g. with shiny, Rstudio, tidyverse etc.). It depends on the type of image whether installations with `apt-get` or `install2.r` are possible. To understand more about how to install R packages in different containers, check it this [cheat sheet](https://raw.githubusercontent.com/sib-swiss/containers-introduction-training/main/r-docker-cheatsheet/r-docker-cheatsheet.pdf), or visit [rocker-project.org](https://rocker-project.org/).


    **Exercise:** Download the `test_deseq2.R` and build the image with `docker build`. Name the image `deseq2`. After that, start an interactive session and execute the script inside the container. 

    !!! hint
        Make an interactive session with the options `-i` and `-t` and use `/bin/bash` as the command. 


    ??? done "Answer"
        Build the container:

        === "x86_64 / AMD64"
            ```sh
            docker build -t deseq2 .
            ```
        === "ARM64 (MacOS M1 chip)"
            ```sh
            DOCKER_DEFAULT_PLATFORM=linux/amd64
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

    ??? done "Answer"
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
        /opt:/usr/local/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
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

    In the case you want to pack your script inside a container, you are building a container specifically for your script, meaning you almost want the container to behave as the program itself. In order to do that, you can use `ENTRYPOINT`. `ENTRYPOINT` is similar to `CMD`, but has two important differences:

    * `ENTRYPOINT` can not be overwritten by the positional arguments (i.e. `docker run image [CMD]`), but has to be overwritten by `--entrypoint`. 
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

    # note that if you want to be able to combine the two
    # both ENTRYPOINT and CMD need to written in the exec form
    ENTRYPOINT ["test_deseq2.R"]

    # default option (if positional arguments are not specified)
    CMD ["--rows", "100"]
    ```

    **Exercise**: Re-build, and run the container non-interactively without any positional arguments. After that, try to pass a different number of rows to `--rows`. How do the commands look?

    ??? done "Answer"
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

    **Exercise**: Push the image to dockerhub, so we can use it later with the apptainer exercises.



    ??? done "Answer"

        Pushing it to dockerhub: 

        ```sh
        docker tag deseq2 [USER NAME]/deseq2:v1
        docker push [USER NAME]/deseq2:v1
        ```

    <h3> Extra: Get information on your image with `docker inspect` </h3>

    We have used `docker inspect` already in the previous chapter to find the default `Cmd` of the ubuntu image. However we can get more info on the image: e.g. the entrypoint, environmental variables, cmd, workingdir etc., you can use the `Config` record from the output of `docker inspect`. For our image this looks like:

    ```yaml
    "Config": {
        "Hostname": "",
        "Domainname": "",
        "User": "",
        "AttachStdin": false,
        "AttachStdout": false,
        "AttachStderr": false,
        "Tty": false,
        "OpenStdin": false,
        "StdinOnce": false,
        "Env": [
            "PATH=/opt:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin",
            "LC_ALL=en_US.UTF-8",
            "LANG=en_US.UTF-8",
            "DEBIAN_FRONTEND=noninteractive",
            "TZ=UTC"
        ],
        "Cmd": [
            "--rows",
            "100"
        ],
        "ArgsEscaped": true,
        "Image": "",
        "Volumes": null,
        "WorkingDir": "/opt",
        "Entrypoint": [
            "test_deseq2.R"
        ],
        "OnBuild": null,
        "Labels": {
            "maintainer": "Dirk Eddelbuettel <edd@debian.org>",
            "org.label-schema.license": "GPL-2.0",
            "org.label-schema.vcs-url": "https://github.com/rocker-org/",
            "org.label-schema.vendor": "Rocker Project"
        }
    }
    ```

    You can annotate your `Dockerfile` and the image by using the instruction `LABEL`. You can give it any key and value with `<key>=<value>`. However, it is recommended to use the [Open Container Initiative (OCI) keys](https://github.com/opencontainers/image-spec/blob/v1.0.1/annotations.md).

    **Exercise**: Annotate our `Dockerfile` with the OCI keys on the creation date, author and description. After that, check whether this has been passed to the actual image with `docker inspect`. 

    !!! note
        You can type `LABEL` for each key-value pair, but you can also have it on one line by seperating the key-value pairs by a space, e.g.:

        ```dockerfile
        LABEL keyx="valuex" keyy="valuey"
        ```

    ??? done "Answer"

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

        # note that if you want to be able to combine the two
        # both ENTRYPOINT and CMD need to written in the exec form
        ENTRYPOINT ["test_deseq2.R"]

        # default option (if positional arguments are not specified)
        CMD ["--rows", "100"]

        ```

        The `Config` record in the output of `docker inspect` was updated with:

        ```yaml
            "Labels": {
                "org.opencontainers.image.authors": "Geert van Geest",
                "org.opencontainers.image.created": "2023-04-12",
                "org.opencontainers.image.description": "Container with DESeq2 and friends",
                "org.opencontainers.image.licenses": "GPL-2.0-or-later",
                "org.opencontainers.image.source": "https://github.com/rocker-org/rocker",
                "org.opencontainers.image.vendor": "Rocker Project"
            }
        ```

    <h3> Extra: Building an image with a browser interface </h3>

    In this exercise, we will use a different base image (`rocker/rstudio:4`), and we'll install the same packages. [Rstudio server](https://posit.co/download/rstudio-server/) is a nice browser interface that you can use for a.o. programming in R. With the image we are creating we will be able to run Rstudio server inside a container.  Check out the `Dockerfile`:

    ```dockerfile
    FROM rocker/rstudio:4

    RUN apt-get update && \
        apt-get install -y libz-dev

    RUN install2.r \
        optparse \
        BiocManager

    RUN R -q -e 'BiocManager::install("biomaRt")'
    ```

    This will create an image from the existing `rstudio` image. It will also install `libz-dev` with `apt-get`, `BiocManager` with `install2.r` and `biomaRt` with an R command. Despiste we're installing the same packages, the installation steps need to be different from the `r-base` image. This is because in the `rocker/rstudio` images R is installed from source, and therefore you can't install packages with `apt-get`. More information on how to install R packages in R containers in this [cheat sheet](https://raw.githubusercontent.com/sib-swiss/containers-introduction-training/main/r-docker-cheatsheet/r-docker-cheatsheet.pdf), or visit [rocker-project.org](https://rocker-project.org/).

    **Exercise:** Build an image based on this `Dockerfile` and give it a meaningful name.

    ??? done "Answer"
        === "x86_64 / AMD64"
            ```sh
            docker build -t rstudio-server .
            ```
        === "ARM64 (MacOS M1 chip)"
            ```sh
            export DOCKER_DEFAULT_PLATFORM=linux/amd64
            docker build -t rstudio-server .
            ```

    You can now run a container from the image. However, you will have to tell docker where to publish port 8787 from the docker container with `-p [HOSTPORT:CONTAINERPORT]`. We choose to publish it to the same port number:

    ```sh
    docker run --rm -it -p 8787:8787 rstudio-server
    ```

    !!! note "Networking"
        More info on docker container networking [here](https://docs.docker.com/config/containers/container-networking/)

    By running the above command, a container will be started exposing rstudio server at port 8787 at localhost. You can approach the instance of jupyterhub by typing `localhost:8787` in your browser. You will be asked for a password. You can find this password in the terminal from which you have started the container.

    We can make this even more interesting by mounting a local directory to the container running the jupyter-lab image:

    ```sh
    docker run \
    -it \
    --rm \
    -p 8787:8787 \
    --mount type=bind,source=/Users/myusername/working_dir,target=/home/rstudio/working_dir \
    rstudio-server
    ```

    By doing this you have a completely isolated and shareable R environment running Rstudio server, but with your local files available to it. Pretty neat right? 

=== "Python"
    In the exercises will use a simple script called `daterange.py`. You can download it [here](https://raw.githubusercontent.com/sib-swiss/containers-introduction-training/main/docker/exercise_python_script/daterange.py). Or copy-paste it from here:

    ```python title="daterange.py"
    #!/usr/bin/env python3

    import pandas as pd
    import argparse

    parser = argparse.ArgumentParser(description = "Get a daterange")

    parser.add_argument('-d', '--date', type=str, required=True, 
                        help='Date. Format: [YYYYMMDD]')

    args = parser.parse_args()

    dates = pd.date_range(args.date, periods=7)

    for d in dates:
        print(d)
    ```
    
    After you have downloaded it, make sure to set the permissions to executable:

    ```sh
    chmod +x daterange.py
    ```
    
    Have a look at `daterange.py`. It is a simple script that uses `pandas` and `argparse`. It takes a date (in the format `YYYYMMDD`) as provided by the option `--date`, and returns a list of all dates in the week starting from that date. An example for execution would be:

    ```sh
    ./daterange.py --date 20220226
    ```

    Giving a list of dates starting from 26-FEB-2022:

    ```
    2022-02-26 00:00:00
    2022-02-27 00:00:00
    2022-02-28 00:00:00
    2022-03-01 00:00:00
    2022-03-02 00:00:00
    2022-03-03 00:00:00
    2022-03-04 00:00:00
    ```

    From the script, you can see it has the dependecy `pandas`, which is not a [built-in module](https://docs.python.org/3/py-modindex.html). In the `Dockerfile` below we give the instruction to install `pandas` with `pip` and copy `daterange.py` to `/opt` inside the container:

    ```dockerfile
    FROM python:3.9.16

    RUN pip install pandas 

    COPY daterange.py /opt 
    ```

    !!! note
        In order to use `COPY`, the file that needs to be copied needs to be in the same directory as the `Dockerfile` or one of its subdirectories.

    **Exercise:** Download the `daterange.py` and build the image with `docker build`. After that, execute the script inside the container. 

    !!! hint
        Make an interactive session with the options `-i` and `-t` and use `/bin/bash` as the command. 

    ??? done "Answer"
        Build the container:

        === "x86_64 / AMD64"
            ```sh
            docker build -t daterange .
            ```
        === "ARM64 (MacOS M1 chip)"
            ```sh
            export DOCKER_DEFAULT_PLATFORM=linux/amd64
            docker build -t daterange .
            ```

        Run the container:

        ```sh
        docker run -it --rm daterange /bin/bash
        ```

        Inside the container we look up the script:

        ```sh
        cd /opt
        ls
        ```

        This should return `daterange.py`. 

        Now you can execute it from inside the container:

        ```sh
        ./daterange.py --date 20220226
        ```

    That's kind of nice. We can ship our python script inside our container. However, we don't want to run it interactively every time. So let's make some changes to make it easy to run it as an executable. For example, we can add `/opt` to the global `$PATH` variable with `ENV`. 

    !!! note "The `$PATH` variable"
        The path variable is a special variable that consists of a list of path seperated by colons (`:`). These paths are searched if you are trying to run an executable. More info this topic at e.g. [wikipedia](https://en.wikipedia.org/wiki/PATH_(variable)). 

    ```dockerfile
    FROM python:3.9.16

    RUN pip install pandas 

    COPY daterange.py /opt 

    ENV PATH=/opt:$PATH
    ```

    !!! note
        The `ENV` instruction can be used to set any variable. 

    **Exercise**: Start an interactive bash session inside the new container. Is the path variable updated? (i.e. can we execute `daterange.py` from anywhere?)

    ??? done "Answer"
        After re-building we start an interactive session:

        ```sh
        docker run -it --rm daterange /bin/bash
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
        daterange.py --date 20220226
        ```

    Instead of starting an interactive session with `/bin/bash` we can now more easily run the script non-interactively:

    ```sh
    docker run --rm daterange daterange.py --date 20220226
    ```

    Now it will directly print the output of `daterange.py` to stdout. 

    In the case you want to pack your script inside a container, you are building a container specifically for your script, meaning you almost want the container to behave as the program itself. In order to do that, you can use `ENTRYPOINT`. `ENTRYPOINT` is similar to `CMD`, but has two important differences:

    * `ENTRYPOINT` can not be overwritten by the positional arguments (i.e. `docker run image [CMD]`), but has to be overwritten by `--entrypoint`. 
    * The positional arguments (or `CMD`) are pasted to the `ENTRYPOINT` command. This means that you can use `ENTRYPOINT` as the executable and the positional arguments (or `CMD`) as the options. 

    Let's try it out:

    ```dockerfile
    FROM python:3.9.16

    RUN pip install pandas 

    COPY daterange.py /opt 

    ENV PATH=/opt:$PATH

    # note that if you want to be able to combine the two
    # both ENTRYPOINT and CMD need to written in the exec form
    ENTRYPOINT ["daterange.py"]

    # default option (if positional arguments are not specified)
    CMD ["--date", "20220226"]
    ```

    **Exercise**: Re-build, and run the container non-interactively without any positional arguments. After that, try to pass a different date to `--date`. How do the commands look?

    ??? done "Answer"
        Just running the container non-interactively would be:

        ```sh
        docker run --rm daterange
        ```

        Passing a different argument (i.e. overwriting `CMD`) would be:

        ```sh
        docker run --rm daterange --date 20210330
        ```

        Here, the container behaves as the executable itself to which you can pass arguments. 

    !!! note 
        You can overwrite `ENTRYPOINT` with `--entrypoint` as an argument to `docker run`. 

    **Exercise**: Push the image to dockerhub, so we can use it later with the apptainer exercises.



    ??? done "Answer"

        Pushing it to dockerhub: 

        ```sh
        docker tag daterange [USER NAME]/daterange:v1
        docker push [USER NAME]/daterange:v1
        ```

    <h3> Extra: Get information on your image with `docker inspect` </h3>

    We have used `docker inspect` already in the previous chapter to find the default `Cmd` of the ubuntu image. However we can get more info on the image: e.g. the entrypoint, environmental variables, cmd, workingdir etc., you can use the `Config` record from the output of `docker inspect`. For our image this looks like:

    ```yaml
    "Config": {
            "Hostname": "",
            "Domainname": "",
            "User": "",
            "AttachStdin": false,
            "AttachStdout": false,
            "AttachStderr": false,
            "Tty": false,
            "OpenStdin": false,
            "StdinOnce": false,
            "Env": [
                "PATH=/opt:/usr/local/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin",
                "LANG=C.UTF-8",
                "GPG_KEY=E3FF2839C048B25C084DEBE9B26995E310250568",
                "PYTHON_VERSION=3.9.4",
                "PYTHON_PIP_VERSION=21.1.1",
                "PYTHON_GET_PIP_URL=https://github.com/pypa/get-pip/raw/1954f15b3f102ace496a34a013ea76b061535bd2/public/get-pip.py",
                "PYTHON_GET_PIP_SHA256=f499d76e0149a673fb8246d88e116db589afbd291739bd84f2cd9a7bca7b6993"
            ],
            "Cmd": [
                "--date",
                "20220226"
            ],
            "ArgsEscaped": true,
            "Image": "",
            "Volumes": null,
            "WorkingDir": "/opt",
            "Entrypoint": [
                "daterange.py"
            ],
            "OnBuild": null,
            "Labels": null
        },
    ```

    <h3> Extra: Adding metadata to your image </h3>

    You can annotate your `Dockerfile` and the image by using the instruction `LABEL`. You can give it any key and value with `<key>=<value>`. However, it is recommended to use the [Open Container Initiative (OCI) keys](https://github.com/opencontainers/image-spec/blob/v1.0.1/annotations.md).

    **Exercise**: Annotate our `Dockerfile` with the OCI keys on the creation date, author and description. After that, check whether this has been passed to the actual image with `docker inspect`. 

    !!! note
        You can type `LABEL` for each key-value pair, but you can also have it on one line by seperating the key-value pairs by a space, e.g.:

        ```dockerfile
        LABEL keyx="valuex" keyy="valuey"
        ```

    ??? done "Answer"

        The `Dockerfile` would look like:

        ```dockerfile
        FROM python:3.9.16

        LABEL org.opencontainers.image.created="2022-04-12" \
            org.opencontainers.image.authors="Geert van Geest" \
            org.opencontainers.image.description="Great container for getting all dates in a week! \
            You will never use a calender again"

        RUN pip install pandas 

        WORKDIR /opt

        COPY daterange.py .

        ENV PATH=/opt:$PATH

        # note that if you want to be able to combine the two
        # both ENTRYPOINT and CMD need to written in the exec form
        ENTRYPOINT ["daterange.py"]

        # default option (if positional arguments are not specified)
        CMD ["--date", "20220226"]

        ```

        The `Config` record in the output of `docker inspect` was updated with:

        ```yaml
        "Labels": {
                    "org.opencontainers.image.authors": "Geert van Geest",
                    "org.opencontainers.image.created": "2022-04-12",
                    "org.opencontainers.image.description": "Great container for getting all dates in a week!     You will never use a calender again"
                }
        ```

    <h3> Extra: Building an image with a browser interface </h3>

    In this exercise, we will use a different base image from the [jupyter docker image stack](https://jupyter-docker-stacks.readthedocs.io/en/latest/using/common.html). [JupyterLab](https://jupyter.org/) is a nice browser interface that you can use for a.o. programming in python. With the image we are creating we will be able to run jupyter lab inside a container.  Check out the `Dockerfile`:

    ```dockerfile
    FROM jupyter/base-notebook:python-3.9

    RUN pip install pandas
    ```

    This will create an image from the existing `python` image. It will also install `jupyterlab` with `pip`. As a default command it starts a jupyter notebook at port 8888.

    !!! note "Ports"
        We have specified here that jupyter lab should use port 8888. However, this **inside** the container. We can not connect to it yet with our browser.

    **Exercise:** Build an image based on this `Dockerfile` and give it a meaningful name.

    ??? done "Answer"
        === "x86_64 / AMD64"
            ```sh
            docker build -t jupyter-lab .
            ```
        === "ARM64 (MacOS M1 chip)"
            ```sh
            export DOCKER_DEFAULT_PLATFORM=linux/amd64
            docker build -t jupyter-lab .
            ```

    You can now run a container from the image. However, you will have to tell docker where to publish port 8888 from the docker container with `-p [HOSTPORT:CONTAINERPORT]`. We choose to publish it to the same port number:

    ```sh
    docker run --rm -it -p 8888:8888 jupyter-lab
    ```

    !!! note "Networking"
        More info on docker container networking [here](https://docs.docker.com/config/containers/container-networking/)

    By running the above command, a container will be started exposing jupyterhub at port 8888 at localhost. You can approach the instance of jupyterhub by typing `localhost:8888` in your browser. You will be asked for a token. You can find this token in the terminal from which you have started the container.

    We can make this even more interesting by mounting a local directory to the container running the jupyter-lab image:

    ```sh
    docker run \
    -it \
    --rm \
    -p 8888:8888 \
    --mount type=bind,source=/Users/myusername/working_dir,target=/working_dir/ \
    jupyter-lab
    ```

    By doing this you have a completely isolated and shareable python environment running jupyter lab, but with your local files available to it. Pretty neat right? 

    !!! note
        Jupyter has a wide range of pre-built images available [here](https://jupyter-docker-stacks.readthedocs.io/en/latest/using/common.html).
