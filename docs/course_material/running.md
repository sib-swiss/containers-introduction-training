
## Restarting an exited container

If you would like to go back to your container with the `figlet` installation, you could try to run again:

```sh
docker run -it ubuntu
```

**Exercise:** Is your `figlet` installation still there? Why?

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

To restart your first created container, you'll have to look up it's name. You can find it in the Docker dashboard, or with `docker container ls -a`.

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
    figlet "try some more text!"
    ```

    Should give you output.

## Creating a new image

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

## Commands

The second positional argument of `docker run` can be a command followed by it's arguments. So, we could run a container non-interactively (without `-it`), and just let it run a single command:

```sh
docker run ubuntu-figlet figlet "non-interactive run"
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

## Removing containers

In the meantime, with every call of `docker run` we have created a new container (check your containers with `docker container ls -a`). You probably don't want to remove those one-by-one. The commands `docker container prune` and `docker image prune` removes stopped containers and dangling images respectively. So, remove your stopped containers with:

```sh
docker container prune
```

Unless you're developing further on a container, or you're using it for an analysis, you probably want to get rid of it once you have exited the container. You can do this with adding `--rm` to your `docker run` command, e.g.:

```sh
docker run --rm ubuntu-figlet figlet "non-interactive run"
```
