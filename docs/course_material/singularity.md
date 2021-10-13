## Learning outcomes

**After having completed this chapter you will be able to:**

* Login to a remote machine with `ssh`
* Use `singularity pull` to convert an image from dockerhub to the 'singularity image format' (`.sif`)
* Execute a singularity container
* Explain the difference in default mounting behaviour between `docker` and `singularity`
* Use `singularity shell` to generate an interactive shell inside a `.sif` image
* Search and use images with both `docker` and `singularity` from [biocontainers.pro](https://biocontainers.pro)


## Material

[:fontawesome-solid-file-pdf: Download the presentation](../assets/pdf/singularity.pdf){: .md-button }

* [Singularity documentation](https://sylabs.io/guides/3.7/user-guide/)
* [Singularity hub](https://singularity-hub.org/)
* [An article on Docker vs Singularity](https://pythonspeed.com/articles/containers-filesystem-data-processing/)
* [Using conda and containers with snakemake](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#containerization-of-conda-based-workflows)

## Exercises

### Login to remote

If you are enrolled in the course, you have received an e-mail with an IP, username, private key and password. To do the Singularity exercises we will login to a remote server. Below you can find instructions on how to login.

=== "mac OS/Linux"
    Open a terminal, and `cd` to the directory where you have stored your private key. After that, change the file permissions of the key:

    ```sh
    chmod 400 key_<username>.pem
    ```

    Then, login like this:

    ```sh
    ssh -i key_<username>.pem <username>@<IP>
    ```

=== "Windows"
    Below you can find video tutorials and information to log in with MobaXterm.

    MobaXterm is an SSH client for Windows. You can use it to connect to the remote host and edit remote scripts. With MobaXterm, you will automatically login to the remote server once you've started the SSH session. Set it up on your own computer using your own credentials and the video below.

    <iframe src="https://player.vimeo.com/video/473838657" width="640" height="360" frameborder="0" allow="autoplay; fullscreen" allowfullscreen></iframe>

### Pulling an image

Singularity can take several image formats (e.g. a `docker` image), and convert them into it's own `.sif` format. Unlike `docker` this image doesn't live in a local image cache, but it's stored as an actual file.

**Exercise:** On the remote server, pull the docker image that has the adjusted default `CMD` that we have pushed to dockerhub [in this exercise](../dockerfiles/#using-cmd) (`ubuntu-figlet-df:v2`) with `singularity pull`. The syntax is:

```sh
singularity pull docker://[USER NAME]/[IMAGE NAME]:[TAG]
```

??? done "Answer"
    ```sh
    singularity pull docker://[USER NAME]/ubuntu-figlet-df:v2
    ```
    This will result in a file called `ubuntu-figlet-df_v2.sif`

!!! note
    If you weren't able to push the image in the previous exercises to your docker hub, you can use `geertvangeest` as username to pull the image.

### Executing an image

These `.sif` files can be run as standalone executables:

```sh
./ubuntu-figlet-df_v2.sif
```

And you can overwrite the default command like this:

```sh
./[IMAGE NAME].sif [COMMAND]
```

**Exercise:** Run the `.sif` file without a command, and with a command that runs `figlet`. Do you get expected output?

??? done "Answer"
    Running it without a command (`./ubuntu-figlet-df_v2.sif`) should give:

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
    ./ubuntu-figlet-df_v2.sif figlet 'Something else'
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

### Mounting with Singularity

Singularity is also different from Docker in the way it handles mounting. By default, Singularity binds your home directory and a number of paths in the root directory to the container. This results in behaviour that is almost like if you are working on the directory structure of the host.  

Running the command `pwd` (full name of current working directory) will therefore result in a path on the host machine:

```sh
./ubuntu-figlet-df_v2.sif pwd
```

**Exercise:** Run the above command. What is the output? How would the output look like if you would run a similar command with Docker?

!!! hint
    A similar Docker command would look like (run this on your local computer):

    ```sh
    docker run --rm ubuntu-figlet-df:v2 pwd
    ```

??? done "Answer"
    The output of `./ubuntu-figlet-df_v2.sif pwd` is the current directory on the host: i.e. `/home/username` if you have it in your home directory. The output of `docker run --rm ubuntu-figlet-df:v2 pwd` (on the local host) would be `/`, which is the default workdir (root directory) of the container. As we did not mount any host directory, this directory exists only within the container (i.e. separated from the host).

### Interactive shell

If you want to debug or inspect an image, it can be helpful to have a shell inside the container. You can do that with `singularity shell`:

```sh
singularity shell ubuntu-figlet-df_v2.sif
```

!!! note
    To exit the shell type `exit`.

**Exercise:** Can you run `figlet` inside this shell?

??? done "Answer"
    Yes:
    ```
    Singularity> figlet test
     _            _
    | |_ ___  ___| |_
    | __/ _ \/ __| __|
    | ||  __/\__ \ |_
     \__\___||___/\__|

    ```

During the lecture you have learned that singularity takes over the user privileges of the user on the host. You can get user information with command like `whoami`, `id`, `groups` etc.

**Exercise:** Run the `figlet` container interactively. Do you have the same user privileges as if you were on the host? How is that with `docker`?

??? done "Answer"
    A command like `whoami` will result in your username printed at stdout:

    ```
    Singularity> whoami
    myusername
    Singularity> id
    uid=1030(myusername) gid=1031(myusername) groups=1031(myusername),1001(condausers)
    Singularity> groups
    myusername condausers
    ```

    With singularity, you have the same privileges inside the singularity container as on the host. If you do this in the docker container (based on the same image), you'll get output like this:

    ```
    root@a3d6e59dc19d:/# whoami
    root
    root@a3d6e59dc19d:/# groups
    root
    root@a3d6e59dc19d:/# id
    uid=0(root) gid=0(root) groups=0(root)
    ```

### A bioinformatics example (extra)

Check out [biocontainers.pro](https://biocontainers.pro). This registry contains a large collection of containerized bioinformatic tools. These are very powerful for pipeline development. You can pull one of their images from dockerhub using `singularity` like this:

```sh
singularity pull docker://biocontainers/fastqc:v0.11.9_cv7
```

!!! note "BioContainers and `singularity`"
     You can directly pull a `singularity` image like so (as shown on the `fastqc` [page](https://biocontainers.pro/tools/fastqc)):

    ```sh
    singularity pull fastqc_0.11.9.sif https://depot.galaxyproject.org/singularity/fastqc:0.11.9--0
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
./fastqc_v0.11.9_cv7.sif fastqc ./reads/ecoli_*.fastq.gz
```

This will result in `html` files in the directory `./reads`. These are quality reports for the sequence reads. If you'd like to view them, you can download them with `scp` or e.g. [FileZilla](https://filezilla-project.org/), and view them with your local browser.
