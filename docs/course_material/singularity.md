## Material

* [An article on Docker vs Singularity](https://pythonspeed.com/articles/containers-filesystem-data-processing/)

## Exercises

Singularity can take many several image formats (e.g. a `docker` image), and convert them into it's own `.sif` format. Unlike `docker` this image doesn't live in a local image cache, but it's stored as an actual file.

Pull the docker image that has the adjusted default `CMD` that we have pushed to dockerhub [in this exercise](../dockerfiles/#using-cmd) (`ubuntu-figlet-df:v2`) with `singularity pull`. The syntax is:

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
    ./ubuntu-figlet-df_v2.sif figlet "Something else"
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

```sh
singularity exec ubuntu-figlet-df_v2.sif figlet "now from singularity"
```

```sh
singularity exec ubuntu-figlet-df_v2.sif pwd
```

```sh
singularity shell ubuntu-figlet-df_v2.sif
```

```sh
which figlet
```

```sh
nano figlet_script.sh
```

```sh
#!/bin/bash

figlet "run from script"
```

```sh
singularity exec ubuntu-figlet-df_v2.sif ./figlet_script.sh
```
