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
FROM ubuntu:focal-20210401
RUN apt-get update
RUN apt-get install figlet
```

!!! note "On writing reproducible `Dockerfiles`"
    At the `FROM` statement in the the above `Dockerfile` you see that we have added a specific tag to the image (i.e. `focal-20210401`). We could also have written:

    ```dockerfile
    FROM ubuntu
    RUN apt-get update
    RUN apt-get install figlet
    ```

    This will automatically pull the image with the tag `latest`. However, if the maintainer of the `ubuntu` images decides to tag another `ubuntu` version as `latest`, rebuilding with the above `Dockerfile` will not give you the same result. Therefore it's always good practice to add the (stable) tag to the image in a `Dockerfile`. More rules on making your `Dockerfiles` more reproducible [here](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1008316).


**Exercise:** Create a file on your computer called `Dockerfile`, and paste the above instruction lines in that file. Make the directory containing the `Dockerfile` your current directory. Build a new image based on that `Dockerfile` with:

```sh
docker build .
```

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

    ```
    docker build -t ubuntu-figlet-df .
    ```

### Using `CMD`

As you might remember the second positional argument of `docker run` is a command (i.e. `docker run IMAGE [CMD]`). If you leave it empty, it uses the default command. You can change the default command in the `Dockerfile` with an instruction starting with `CMD`. For example:

```dockerfile
FROM ubuntu:focal-20210401
RUN apt-get update
RUN apt-get install figlet
CMD figlet 'My image works!'
```

**Exercise:** Build a new image based on the above `Dockerfile`. Can you validate the change using `docker image inspect`? Can you overwrite this default with `docker run`?

??? done "Answer"
    Copy the new line to your `Dockerfile`, and build the new image like this:

    ```sh
    docker build -t ubuntu-figlet-df:v2 .
    ```

    The command `docker inspect ubuntu-figlet-df:v2` will give:

    ```
    "Cmd": [
        "/bin/sh",
        "-c",
        "figlet 'My image works!'"
    ]
    ```

    So the default command (`/bin/bash`) has changed to `figlet 'My image works!'`

    Running the image (with clean-up (`--rm`)):

    ```sh
    docker run --rm ubuntu-figlet-df:v2
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
    docker run --rm ubuntu-figlet-df:v2 figlet 'another text'
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
    FROM ubuntu:focal-20210401
    RUN apt-get update
    RUN apt-get install figlet
    CMD figlet 'My image works!'
    ```

    A `Dockerfile` with exec notation:

    ```dockerfile
    FROM ubuntu:focal-20210401
    RUN apt-get update
    RUN apt-get install figlet
    CMD ["/bin/sh", "-c", "figlet 'My image works!'"]
    ```

**Exercise:** Now push our created image (with a version tag) to docker hub. We will use it later for the [`singularity` exercises](singularity.md).

??? done "Answer"
    ```sh
    docker tag ubuntu-figlet-df:v2 [USER NAME]/ubuntu-figlet-df:v2
    docker push [USER NAME]/ubuntu-figlet-df:v2
    ```

### A more real-world example (extra)

You might have gotten enough of `figlet`. Let's do something more fancy. Check out this `Dockerfile`:

```dockerfile
FROM python:3.9.4-buster

RUN pip install jupyterlab

CMD jupyter lab --ip=0.0.0.0 --port=8888 --allow-root
```

This will create an image from the existing `python` image. It will also install `jupyterlab` with `pip`. As a default command it starts a jupyter notebook at port 8888.

**Exercise:** Build an image based on this `Dockerfile` and give it a meaningful name.

??? done "Answer"
    ```sh
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

By doing this you have a completely isolated and shareable python environment running jupyter lab, but with your local files available to it.

!!! note
    Jupyter has a wide range of pre-built images available [here](https://jupyter-docker-stacks.readthedocs.io/en/latest/using/common.html). Example syntax with a pre-built jupyter image would look like:

    ```sh
    docker run \
    --rm \
    -e JUPYTER_ENABLE_LAB=yes \
    -p 8888:8888 \
    jupyter/base-notebook
    ```

    Using the above will also give you easier control over security, users and permissions.
