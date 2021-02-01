## Material

## Exercises

Let's create our first container from an existing image. We do this with the image `ubuntu`, generating an environment with a minimal installation of ubuntu.  

```sh
docker run -it ubuntu
```

This will give you an interactive shell into the created container (this interactivity was invoked by the options `-i` and `-t`) .

**Exercise:** Check out the operating system of the container by typing `cat /etc/os-release` in the container's shell. Are we really in an ubuntu environment?

??? done "Answer"
    Yes:

    ```
    root@27f7d11608de:/# cat /etc/os-release
    NAME="Ubuntu"
    VERSION="20.04.1 LTS (Focal Fossa)"
    ID=ubuntu
    ID_LIKE=debian
    PRETTY_NAME="Ubuntu 20.04.1 LTS"
    VERSION_ID="20.04"
    HOME_URL="https://www.ubuntu.com/"
    SUPPORT_URL="https://help.ubuntu.com/"
    BUG_REPORT_URL="https://bugs.launchpad.net/ubuntu/"
    PRIVACY_POLICY_URL="https://www.ubuntu.com/legal/terms-and-policies/privacy-policy"

    VERSION_CODENAME=focal
    UBUNTU_CODENAME=focal
    ```

!!! note "Where does the image come from?"
    If the image `ubuntu` was not on your computer yet, `docker` will search and try to get them from [dockerhub](https://hub.docker.com/), and download it.

**Exercise:** Run the command `whoami` in the docker container. Who are you?

??? done "Answer"
    The command `whoami` returns the current user. In the container `whoami` will return `root`. This means you are the [`root` user](https://en.wikipedia.org/wiki/Superuser) i.e. within the container you are admin and can basically change anything.  

Check out the container panel at the Docker dashboard (the Docker gui) or open another host terminal and type:

```
docker container ls -a
```

**Exercise:** What is the container status?

??? done "Answer"
    In Docker dashboard you can see that the shell is running:

    <figure>
      <img src="../../assets/images/running_container_dashboard.png" width="300"/>
    </figure>

    The output of `docker container ls -a` is:

    ```
    CONTAINER ID   IMAGE     COMMAND       CREATED         STATUS         PORTS     NAMES
    27f7d11608de   ubuntu    "/bin/bash"   7 minutes ago   Up 6 minutes             great_moser
    ```

    Also showing you that the `STATUS` is `Up`.

Now let's install some software in our `ubuntu` environment. We'll install some simple software called [`figlet`](http://www.figlet.org/). Type into the container shell:

```sh
apt-get update
apt-get install figlet
```

!!! note "This will give some warnings"
    This installation will give some warnings. It's safe to ignore them.

Now let's try it out. Type into the container shell:

```sh
figlet 'SIB courses are great!'
```

Now you have installed and used software `figlet` in an `ubuntu` environment completely separated from your host computer. This already gives you an idea of the power of containerization.

Exit the shell by typing `exit`. Check out the container panel of Docker dashboard or type:

```sh
docker container ls -a
```

**Exercise:** What is the container status?

??? done "Answer"
    `docker container ls -a` gives:

    ```
    CONTAINER ID   IMAGE     COMMAND       CREATED          STATUS                     PORTS     NAMES
    27f7d11608de   ubuntu    "/bin/bash"   15 minutes ago   Exited (0) 8 seconds ago             great_moser
    ```

    Showing that the container has exited, meaning it's not running.
