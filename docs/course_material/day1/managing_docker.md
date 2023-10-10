## Learning outcomes

**After having completed this chapter you will be able to:**

* Explain the concept of layers in the context of docker containers and images
* Use the command line to restart and re-attach to an exited container
* Create a new image with `docker commit`
* List locally available images with `docker image ls`
* Run a command inside a container non-interactively
* Use `docker image inspect` to get more information on an image
* Use the command line to prune dangling images and stopped containers
* Rename and tag a docker image
* Push a newly created image to dockerhub
* Use the option `--mount` to bind mount a host directory to a container

## Material

[:fontawesome-solid-file-pdf: Download the presentation](../../assets/pdf/managing_docker.pdf){: .md-button }

<iframe width="560" height="315" src="https://www.youtube.com/embed/y7eCm9NOHYg" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

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

In the previous exercises we have run containers without a command as positional argument. This doesn't mean that no command has been run, because the container would do nothing without a command. The default command is stored in the image, and you can find it by `docker image inspect [IMAGE NAME]`.  

**Exercise:** Have a look at the output of `docker image inspect`, particularly at `"Config"` (ignore `"ContainerConfig"` for now). What is the default command (`CMD`) of the ubuntu image?

??? done "Answer"
    Running `docker image inspect ubuntu` gives (amongst other information):
    
    ```yaml
    "Cmd": [
        "/bin/bash"
    ],
    ```

    In the case of the ubuntu the default command is `bash`, returning a shell in `bash` (i.e. *Bourne again shell*). Adding the options `-i` and `-t` (`-it`) to your `docker run` command will therefore result in an interactive `bash` shell. You can modify this default behaviour. More on that later, when we will work on [Dockerfiles](dockerfiles.md).

!!! note "The difference between `Config` and `ContainerConfig`"
    The configuration at `Config` represents the image, the configuration at `ContainerConfig` the last step during the build of the image, i.e. the last layer. More info e.g. at [this post at stackoverflow](https://stackoverflow.com/questions/36216220/what-is-different-of-config-and-containerconfig-of-docker-inspect).

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

!!! note "If on Linux"
    If you are on Linux and haven't connected to docker hub before, you will have login first. To do that, run:

    ```sh
    docker login
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

!!! note "MobaXterm users"
    You can specify your local path with the Windows syntax (e.g. `C:\Users\myusername`). However, you will have to use forward slashes (`/`) instead of backward slashes (`\`). Therefore, mounting a directory would look like:

    ```sh
    docker run \
    --mount type=bind,source=C:/Users/myusername,target=/path/in/container \
    [IMAGE]
    ```

    Do not use autocompletion or variable substitution (e.g. `$PWD`) in MobaXterm, since these point to 'emulated' paths, and are not passed properly to the docker command.

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
    root@8d80a8698865:/# figlet 'testing mounted dir' > /working_dir/figlet_output.txt
    ```

    This should create a file in both your host (local) source directory and the target directory in the container called `figlet_output.txt`.

!!! note "Using files on the host"
    This of course also works the other way around. If you would have a file on the host with e.g. a text, you can copy it into your mounted directory, and it will be available to the container.

### Managing permissions (extra)

Depending on your system, the user ID and group ID will be taken over from the user inside the container. If the user inside the container is root, this will be root. That's a bit inconvenient if you just want to run the container as a regular user (for example in certain circumstances your container could write in `/`). To do that, use the `-u` option, and specify the group ID and user ID like this:

```sh
docker run -u [uid]:[gid]
```

So, e.g.:

```sh
docker run \
-it \
-u 1000:1000 \
--mount type=bind,source=/Users/myusername/working_dir,target=/working_dir/ \
ubuntu-figlet
```

If you want docker to take over your current uid and gid, you can use:

```
docker run -u "$(id -u):$(id -g)"
```

!!! note "This behaviour is different on MacOS and MobaXterm"
    On MacOS and in the local shell of MobaXterm the uid and gid are taken over from the user running the container (even if you set `-u` as 0:0), i.e. your current ID. More info on [stackoverflow](https://stackoverflow.com/questions/43097341/docker-on-macosx-does-not-translate-file-ownership-correctly-in-volumes).

**Exercise:** Start an interactive container based on the `ubuntu-figlet` image, bind-mount a local directory and take over your current `uid` and `gid`. Write the output of a `figlet` command to a file in the mounted directory. Who and which group owns the file inside the container? And outside the container? Answer the same question but now run the container without setting `-u`.

??? done "Answer"
    === "Linux"

        **Running `ubuntu-figlet` interactively while taking over `uid` and `gid` and mounting my current directory:**

        ```sh
        docker run -it --mount type=bind,source=$PWD,target=/data -u "$(id -u):$(id -g)" ubuntu-figlet
        ```
        Inside container:

        ```
        I have no name!@e808d7c36e7c:/$ id
        uid=1000 gid=1000 groups=1000
        ```

        So, I have taken over uid 1000 and gid 1000.

        ```
        I have no name!@e808d7c36e7c:/$ cd /data
        I have no name!@e808d7c36e7c:/data$ figlet 'uid set' > uid_set.txt
        I have no name!@e808d7c36e7c:/data$ ls -lh
        -rw-r--r-- 1 1000 1000 0 Mar  400 13:37 uid_set.txt
        ```

        So the file belongs to user 1000, and group 1000.

        Outside container:

        ```
        ubuntu@ip-172-31-33-21:~$ ls -lh
        -rw-r--r-- 1 ubuntu ubuntu 400 Mar  5 13:37 uid_set.txt
        ```

        Which makes sense:

        ```
        ubuntu@ip-172-31-33-21:~$ id
        uid=1000(ubuntu) gid=1000(ubuntu) groups=1000(ubuntu)
        ```

        **Running `ubuntu-figlet` interactively without taking over `uid` and `gid`:**

        ```sh
        docker run -it --mount type=bind,source=$PWD,target=/data ubuntu-figlet
        ```
        Inside container:

        ```
        root@fface8afb220:/# id
        uid=0(root) gid=0(root) groups=0(root)
        ```

        So, uid and gid are `root`.

        ```
        root@fface8afb220:/# cd /data
        root@fface8afb220:/data# figlet 'uid unset' > uid_unset.txt
        root@fface8afb220:/data# ls -lh
        -rw-r--r-- 1 1000 1000 400 Mar  5 13:37 uid_set.txt
        -rw-r--r-- 1 root root 400 Mar  5 13:40 uid_unset.txt
        ```

        Outside container:

        ```
        ubuntu@ip-172-31-33-21:~$ ls -lh
        -rw-r--r-- 1 ubuntu ubuntu 0 Mar  5 13:37 uid_set.txt
        -rw-r--r-- 1 root   root   0 Mar  5 13:40 uid_unset.txt
        ```

        So, the uid and gid 0 (root:root) are taken over.

    === "MacOS"

        **Running `ubuntu-figlet` interactively while taking over `uid` and `gid` and mounting my current directory:**

        ```sh
        docker run -it --mount type=bind,source=$PWD,target=/data -u "$(id -u):$(id -g)" ubuntu-figlet
        ```
        Inside container:

        ```
        I have no name!@e808d7c36e7c:/$ id
        uid=503 gid=20(dialout) groups=20(dialout)
        ```
        So, the container has taken over uid 503 and group 20

        ```
        I have no name!@e808d7c36e7c:/$ cd /data
        I have no name!@e808d7c36e7c:/data$ figlet 'uid set' > uid_set.txt
        I have no name!@e808d7c36e7c:/data$ ls -lh
        -rw-r--r--  1 503 dialout    400 Mar  5 13:11 uid_set.txt
        ```

        So the file belongs to user 503, and the group `dialout`.

        Outside container:

        ```
        mac-34392:~ geertvangeest$ ls -lh
        -rw-r--r--   1 geertvangeest  staff     400B Mar  5 14:11 uid_set.txt
        ```

        Which are the same as inside the container:

        ```
        mac-34392:~ geertvangeest$ echo "$(id -u):$(id -g)"
        503:20
        ```

        The `uid` 503 was nameless in the docker container. However the group 20 already existed in the ubuntu container, and was named `dialout`.

        **Running `ubuntu-figlet` interactively without taking over `uid` and `gid`:**

        ```sh
        docker run -it --mount type=bind,source=$PWD,target=/data ubuntu-figlet
        ```
        Inside container:

        ```
        root@fface8afb220:/# id
        uid=0(root) gid=0(root) groups=0(root)
        ```

        So, inside the container I am `root`.
        Creating new files will lead to ownership of `root` inside the container:

        ```
        root@fface8afb220:/# cd /data
        root@fface8afb220:/data# figlet 'uid unset' > uid_unset.txt
        root@fface8afb220:/data# ls -lh
        -rw-r--r--  1 503 dialout    400 Mar  5 13:11 uid_set.txt
        -rw-r--r--  1 root root    400 Mar  5 13:25 uid_unset.txt
        ```

        Outside container:

        ```
        mac-34392:~ geertvangeest$ ls -lh
        -rw-r--r--   1 geertvangeest  staff     400B Mar  5 14:11 uid_set.txt
        -rw-r--r--   1 geertvangeest  staff     400B Mar  5 14:15 uid_unset.txt
        ```

        So, the uid and gid 0 (root:root) are not taken over. Instead, the uid and gid of the user running docker were used.

    === "MobaXterm"

        **Running `ubuntu-figlet` interactively while taking over `uid` and `gid` and mounting to a  specfied directory:**

        ```sh
        docker run -it --mount type=bind,source=C:/Users/geert/data,target=/data -u "$(id -u):$(id -g)" ubuntu-figlet
        ```
        Inside container:

        ```
        I have no name!@e808d7c36e7c:/$ id
        uid=1003 gid=513 groups=513
        ```
        So, the container has taken over uid 1003 and group 513

        ```
        I have no name!@e808d7c36e7c:/$ cd /data
        I have no name!@e808d7c36e7c:/data$ figlet 'uid set' > uid_set.txt
        I have no name!@e808d7c36e7c:/data$ ls -lh
        -rw-r--r--  1 1003 513    400 Mar  5 13:11 uid_set.txt
        ```

        So the file belongs to user 1003, and the group 513.

        Outside container:

        ```
        /home/mobaxterm/data$ ls -lh
        -rwx------   1 geert  UserGrp     400 Mar  5 14:11 uid_set.txt
        ```

        Which are the same as inside the container:

        ```
        /home/mobaxterm/data$ echo "$(id -u):$(id -g)"
        1003:513
        ```

        **Running `ubuntu-figlet` interactively without taking over `uid` and `gid`:**

        ```sh
        docker run -it --mount type=bind,source=C:/Users/geert/data,target=/data ubuntu-figlet
        ```
        Inside container:

        ```
        root@fface8afb220:/# id
        uid=0(root) gid=0(root) groups=0(root)
        ```

        So, inside the container I am `root`.
        Creating new files will lead to ownership of `root` inside the container:

        ```
        root@fface8afb220:/# cd /data
        root@fface8afb220:/data# figlet 'uid unset' > uid_unset.txt
        root@fface8afb220:/data# ls -lh
        -rw-r--r--  1 1003 503    400 Mar  5 13:11 uid_set.txt
        -rw-r--r--  1 root root    400 Mar  5 13:25 uid_unset.txt
        ```

        Outside container:

        ```
        /home/mobaxterm/data$ ls -lh
        -rwx------   1 geert  UserGrp     400 Mar  5 14:11 uid_set.txt
        -rwx------   1 geert  UserGrp     400 Mar  5 14:15 uid_unset.txt
        ```

        So, the uid and gid 0 (root:root) are not taken over. Instead, the uid and gid of the user running docker were used.
