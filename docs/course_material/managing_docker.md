## Learning outcomes

**After having completed this chapter you will be able to:**

* Explain the concept of layers in the context of docker containers and images
* Explain the behaviour of `docker run` while creating a container from an image
* Use the command line to restart and re-attach to an exited container
* Create a new image with `docker commit`
* List locally available images with `docker image ls`
* Run a command inside a container non-interactively
* Find the default command of the image with `docker image inspect`
* Use the command line to prune dangling images and stopped containers
* Rename and tag a docker image
* Push a newly created image to dockerhub
* Use the option `--mount` to bind mount a host directory to a container

## Material

[:fontawesome-solid-file-pdf: Download the presentation](../assets/pdf/managing_docker.pdf){: .md-button }

* [Overview of how docker works](https://docs.docker.com/get-started/overview/)
* [More on bind mounts](https://docs.docker.com/storage/bind-mounts/)
* [Docker volumes in general](https://docs.docker.com/storage/volumes/)

## Exercises

### Restarting an exited container

If you would like to go back to your container with the `figlet` installation, you could try to run again:

```sh
docker run -it ubuntu
```

**Exercise:** Run the above command. Is your `figlet` installation still there? Why?

!!! hint
    Check the status of your containers:

    ```sh
    docker container ls -a
    ```

??? done "Answer"
    No, the installation is gone. Another container was created from the same ubuntu image, without the `figlet` installation.
    Running the command `docker container ls -a` results in:

    ```
    CONTAINER ID   IMAGE     COMMAND       CREATED              STATUS                     PORTS     NAMES
    8d7c4c611b70   ubuntu    "/bin/bash"   About a minute ago   Up About a minute                    kind_mendel
    27f7d11608de   ubuntu    "/bin/bash"   27 minutes ago       Exited (0) 2 minutes ago             great_moser
    ```

    In this case the container `great_moser` contains the `figlet` installation. But we have exited that container. We created a new container (`kind_mendel` in this case) with a fresh environment created from the original `ubuntu` image.

To restart your first created container, you'll have to look up its name. You can find it in the Docker dashboard, or with `docker container ls -a`.

!!! hint "Container names"
    The container name is the funny combination of two words separated by `_`, e.g.: `nifty_sinoussi`. Alternatively you can use the container ID (the first column of the output of `docker container ls`)

To restart a container you can use:

```sh
docker start [CONTAINER NAME]
```

And after that to re-attach to the shell:

```sh
docker attach [CONTAINER NAME]
```

And you're back in the container shell.

**Exercise:** Run the `docker start` and `docker attach` commands for the container that is supposed to contain the `figlet` installation. Is the installation of `figlet` still there?

??? done "Answer"
    yes:

    ```sh
    figlet 'try some more text!'
    ```

    Should give you output.

!!! note "`docker attach` and `docker exec`"
    In addition to `docker attach`, you can also "re-attach" a container with `docker exec`. However, these two are quite different. While `docker attach` gets you back to your stopped shell process, `docker exec` creates a new one (more information on [stackoverflow](https://stackoverflow.com/questions/30960686/difference-between-docker-attach-and-docker-exec)). The command `docker exec` enables you therefore to have multiple shells open in the same container. That can be convenient if you have one shell open with a program running in the foreground, and another one for e.g. monitoring. An example for using `docker exec` on a running container:

    ```sh
    docker exec -it [CONTAINER NAME] /bin/bash
    ```

    Note that  `docker exec` requires a CMD, it doesn't use the default.

### Creating a new image

You can store your changes and create a new image based on the `ubuntu` image like this:

```sh
docker commit [CONTAINER NAME] ubuntu-figlet
```

**Exercise:** Run the above command with the name of the container containing the `figlet` installation. Check out `docker image ls`. What have we just created?

??? done "Answer"
    A new image called `ubuntu-figlet` based on the status of the container.
    The output of `docker image ls` should look like:

    ```
    REPOSITORY                        TAG       IMAGE ID       CREATED         SIZE
    ubuntu-figlet                     latest    e08b999c7978   4 seconds ago   101MB
    ubuntu                            latest    f63181f19b2f   29 hours ago    72.9MB
    ```

Now you can generate a new container based on the new image:

```sh
docker run -it ubuntu-figlet
```

**Exercise:** Run the above command. Is the `figlet` installation in the created container?

??? done "Answer"
    yes

### Commands

The second positional argument of `docker run` can be a command followed by its arguments. So, we could run a container non-interactively (without `-it`), and just let it run a single command:

```sh
docker run ubuntu-figlet figlet 'non-interactive run'
```

Resulting in just the output of the `figlet` command.

In the previous exercises we have run containers without a command as positional argument. This doesn't mean that no command has been run, because the container would do nothing without a command. The default command stored in the image, and you can find it by `docker image inspect [IMAGE NAME]`.  

**Exercise:** What is the default command (`CMD`) of the ubuntu image?

??? done "Answer"
    Running `docker image inspect ubuntu` gives (amongst other information):
    ```
    "Cmd": [
               "/bin/sh",
               "-c",
               "#(nop) ",
               "CMD [\"/bin/bash\"]"
           ],
    ```
    The first part in the list following `"Cmd":` is the shell in which the command is executed (`/bin/sh -c`; i.e. *Bourne shell*), the second part, following `CMD`, is the default command. In the case of the ubuntu image this is `/bin/bash`, returning a shell in `bash` (i.e. *Bourne again shell* in stead of `sh`). Adding the options `-i` and `-t` (`-it`) to your `docker run` command will therefore result in an interactive `bash` shell. You can modify this default behaviour. More on that later, when we will work on [Dockerfiles](dockerfiles.md).

### Removing containers

In the meantime, with every call of `docker run` we have created a new container (check your containers with `docker container ls -a`). You probably don't want to remove those one-by-one. These two commands are very useful to clean up your Docker cache:

* `docker container prune`: removes stopped containers
* `docker image prune`: removes dangling images (i.e. images without a name)

So, remove your stopped containers with:

```sh
docker container prune
```

Unless you're developing further on a container, or you're using it for an analysis, you probably want to get rid of it once you have exited the container. You can do this with adding `--rm` to your `docker run` command, e.g.:

```sh
docker run --rm ubuntu-figlet figlet 'non-interactive run'
```

### Pushing to dockerhub

Now that we have created our first own docker image, we can store it and share it with the world on docker hub. Before we get there, we first have to (re)name and tag it.

Before pushing an image to dockerhub, `docker` has to know to which user and which repository the image should be added. That information should be in the name of the image, like this: `user/imagename`. We can rename an image with `docker tag` (which is a bit of misleading name for the command). So we could push to dockerhub like this:

```
docker tag ubuntu-figlet [USER NAME]/ubuntu-figlet
docker push [USER NAME]/ubuntu-figlet
```

!!! note "How docker makes money"
    All images pushed to dockerhub are open to the world. With a free account you can have one image on dockerhub that is private. Paid accounts can have more private images, and are therefore popular for commercial organisations. As an alternative to dockerhub, you can store images locally with [`docker save`](https://docs.docker.com/engine/reference/commandline/save/).

We didn't specify the tag for our new image. That's why `docker tag` gave it the default tag called `latest`. Pushing an image without a tag will overwrite the current image with the tag `latest` (more on (not) using `latest` [here](https://vsupalov.com/docker-latest-tag/)). If you want to maintain multiple versions of your image, you will have to add a tag, and push the image with that tag to dockerhub:

```
docker tag ubuntu-figlet [USER NAME]/ubuntu-figlet:v1
docker push [USER NAME]/ubuntu-figlet:v1
```

### Mounting a directory

For many analyses you do calculations with files or scripts that are on your host (local) computer. But how do you make them available to a docker container? You can do that in several ways, but here we will use bind-mount. You can bind-mount a directory with `-v` (`--volume`) or `--mount`. Most old-school `docker` users will use `-v`, but `--mount` syntax is easier to understand and now recommended, so we will use the latter here:

```sh
docker run \
--mount type=bind,source=/host/source/path,target=/path/in/container \
[IMAGE]
```

The target directory will be created if it does not yet exist. The source directory should exist.

!!! note "Using docker from Windows PowerShell"
    Most of the syntax for `docker` is the same for both PowerShell and UNIX-based systems. However, there are some differences, e.g. in Windows, directories in file paths are separated by `\` instead of `/`. Also, line breaks are not escaped by `\` but by `.

**Exercise:** Mount a host (local) directory to a target directory `/working_dir` in a container created from the `ubuntu-figlet` image and run it interactively. Check whether the target directory has been created.

??? done "Answer"
    e.g. on Mac OS this would be:

    ```sh
    docker run \
    -it \
    --mount type=bind,source=/Users/myusername/working_dir,target=/working_dir/ \
    ubuntu-figlet
    ```

    This creates a directory called `working_dir` in the root directory (`/`):

    ```
    root@8d80a8698865:/# ls
    bin   dev  home  lib32  libx32  mnt  proc  run   srv  tmp  var
    boot  etc  lib   lib64  media   opt  root  sbin  sys  usr  working_dir
    ```

This mounted directory is both available for the host (locally) and for the container. You can therefore e.g. copy files in there, and write output generated by the container.

**Exercise:** Write the output of `figlet "testing mounted dir"` to a file in `/working_dir`. Check whether it is available on the host (locally) in the source directory.

!!! hint
    You can write the output of `figlet` to a file like this:
    ```sh
    figlet 'some string' > file.txt
    ```

??? done "Answer"
    ```
    root@8d80a8698865:/# figlet "testing mounted dir" > /working_dir/figlet_output.txt
    ```

    This should create a file in both your host (local) source directory and the target directory in the container called `figlet_output.txt`.

!!! note "Using files on the host"
    This of course also works the other way around. If you would have a file on the host with e.g. a text, you can copy it into your mounted directory, and it will be available to the container.
