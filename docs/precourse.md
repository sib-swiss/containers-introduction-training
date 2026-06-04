# Precourse instructions

## Background knowledge: UNIX

As stated in the course prerequisites at the [announcement web page](https://www.sib.swiss/training/course/20260609_DOCKR), we expect participants to have a basic understanding of working with the command line on UNIX-based systems. You can test your UNIX skills with a quiz [here](https://docs.google.com/forms/d/e/1FAIpQLSd2BEWeOKLbIRGBT_aDEGPce1FOaVYBbhBiaqcaHoBKNB27MQ/viewform?usp=sf_link). If you don't have experience with UNIX command line, or if you're unsure whether you meet the prerequisites, follow our [online UNIX tutorial](https://edu.sib.swiss/pluginfile.php/2878/mod_resource/content/4/couselab-html/content.html).

## Software

### OS and terminal

During the course exercises you will be mainly interacting with Docker through the command line. If you are using a UNIX or UNIX-like OS (_e.g_ MacOS), you already have a terminal readily usable for the course. If you are working with Windows, although Windows Powershell is suitable for that, we strongly recommend to install a UNIX or 'UNIX-like' terminal. You can get this by using [MobaXterm](https://mobaxterm.mobatek.net/ "get MobaXterm") or [WSL2](https://docs.microsoft.com/en-us/windows/wsl/install) (recommended solution).

### Code editor

You should have a modern code editor installed. During the course, we can give only limited help for installation and set-up issues, so we will only "officially" support [VS Code](https://code.visualstudio.com/download). If you are already very familiar with another combination of modern code editor/command line interface, feel free to use it, but know that support for this will be limited.

!!! info "Additional requirements"
    Besides VS Code, make sure you have completed the setup for the `Remote-SSH` extension:

    * [OpenSSH compatible client](https://code.visualstudio.com/docs/remote/troubleshooting#_installing-a-supported-ssh-client): A compatible SSH client is required and is usually pre-installed on most operating systems. To verify that it is available, open a terminal and run `ssh`.
    * [Remote-SSH extension](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-ssh): Install the Remote-SSH extension in VS Code. Open VS Code then open the **Extensions** view (the four-square icon in the Activity Bar), search for **Remote - SSH**, and click **Install**.

??? warning "Additional Windows requirement"
    If you are working on Windows, we recommend to install the [WSL extension for VS Code](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-wsl) to make file management easier.

### Docker

You will need to create a [Docker Hub](https://hub.docker.com/) account and install Docker on your local computer. This will allow you to build and push Docker images to Docker Hub, which is required for some of the exercises. You will also need Docker to run some of the exercises locally.

#### Docker Hub account creation

You can create an account [here](https://app.docker.com/signup). For this course, we recommend to create a **personal account**.

#### Docker CLI installation

The Docker CLI (`docker`) is the command-line tool used to interact with Docker daemons and manage containers, images, networks, and volumes.

Docker can be installed in two main ways:

* **Docker Desktop**: An all-in-one application for macOS and Windows that includes the Docker CLI, Docker Engine, a graphical interface, and additional developer tools
* **Docker Engine**: The container runtime itself. On Linux, Docker Engine and the Docker CLI are typically installed directly from the distribution packages or Docker repositories

For this course, you need to install:

* **Docker Engine** on Linux
* **Docker Desktop** on macOS and Windows

??? note "Docker and admin rights"
    Note that you need admin rights to install and use Docker, and if you are installing Docker on Windows, you need a recent Windows version.

##### Linux

Install Docker Engine by following the [official instructions](https://docs.docker.com/engine/install#installation-procedures-for-supported-platforms) for your Linux distribution. The installation will include both Docker Engine and the Docker CLI.

##### macOS

Install Docker Engine by following the [official instructions for MacOS](https://docs.docker.com/desktop/setup/install/mac-install/). Docker Desktop includes the Docker CLI and makes it available in the terminal.

##### Windows

!!! warning "WSL2 and Docker"
    If possible, we strongly recommend to use WSL2 on Windows and to install Docker Engine with the WSL2 backend. In this case, follow the Linux instructions above.

Install Docker Engine by following the [official instructions for Windows](https://docs.docker.com/desktop/setup/install/windows-install/). Docker Desktop includes the Docker CLI and makes it available for use in PowerShell, Command Prompt, and Windows Terminal.

##### Verify the installation

Open a terminal and run:

```bash
docker --version
```

You should see an output similar to:

```text
Docker version 29.5.2, build 79eb04c
```

You can also verify that the Docker daemon is reachable by running:

```bash
docker info
```

If Docker is correctly installed and running, information about the Docker server and client will be displayed:

```test
Client: Docker Engine - Community
 Version:    29.5.2
 Context:    default
 Debug Mode: false
 Plugins:
  buildx: Docker Buildx (Docker Inc.)
    Version:  v0.34.1
    Path:     /usr/libexec/docker/cli-plugins/docker-buildx
  compose: Docker Compose (Docker Inc.)
    Version:  v5.1.4
    Path:     /usr/libexec/docker/cli-plugins/docker-compose
...
```

### Apptainer

You don't need to install Apptainer on your local machine as you will be using a version we installed on the remote server. However, if you want to have a look before the course, you can find installation instructions for different OS [here](https://apptainer.org/docs/admin/main/installation.html). We recommend the `non-setuid` install.

### SSH connection to a server

In addition to your local computer, you will be working on an Amazon Web Services ([AWS](https://aws.amazon.com/)) Elastic Cloud (EC2) [server](https://aws.amazon.com/ec2/). This Ubuntu server behaves like a 'normal' remote server, and can be approached through [`ssh`](https://man7.org/linux/man-pages/man1/ssh.1.html). If you are enrolled in the course, you have access to a shared document containing instruction to retrieve your username and private ssh key, granting you access to a personal home directory on the server.

!!! danger "Server availability"
    Please note that, for cost reasons, the server will be **started at 8:00 on the morning of the course** and will be **stopped at 18:00 on the day after the course ends**, so you will not be able to connect before or after these dates. Likewise, if you need to retrieve data from the server, please do it before it is stopped.

!!! info "Help with SSH connection"
    The trainers will be available between **8:30 and 9:00 on the morning of the course** to help with technical set-up and troubleshooting.

??? tip "If you want to know more about `ssh`"
    If you are not familiar with `ssh`, you can check the [Heidelberg University tutorial](https://zah.uni-heidelberg.de/it-guide/ssh-tutorial-linux) for information on how to set it up and use it.

#### VS Code instructions

Here are instructions on how to use VS Code to connect with SSH to a remote server. First, place the `key_username.pem` file in the proper folder:

=== "Windows"
    Open a PowerShell terminal, `cd` to the directory where you have stored your private key (`key_username.pem` file) and move it to `~\.ssh`:
    ```powershell
    Move-Item -Path key_username.pem -Destination $HOME\.ssh
    ```

=== "macOS/Linux"
    Open a terminal, `cd` to the directory where you have stored your private key (`key_username.pem` file), change the permissions of the key file and move it to `~/.ssh`:
    ```sh
    chmod 400 key_username.pem
    mv key_username.pem ~/.ssh
    ```

Then:

* Open VS Code and click on the green or blue button in the bottom left corner
* Select `Connect to Host...` and then `Configure SSH Hosts...`
* Specify a location for the SSH config file (preferably the same directory as where your keys are stored): `~/.ssh/config`
* A skeleton config file will be provided. Edit it, so it looks like this (replace `username` with your username, and make sure the IP address in `HostName` match what the one given in the shared document):

    === "Windows"
        ```
        Host sib_course_remote
            User username
            HostName 52.57.40.7
            IdentityFile ~\.ssh\key_username.pem
        ```
        Note: if you are working with the Windows SSH executable (for example `C:\WINDOWS\System32\OpenSSH\ssh.exe`), you may have to use the full path of the key file instead of a relative one in `IdentityFile`:
        ```
            IdentityFile C:\Users\<windows_username>\.ssh\key_username.pem
        ```
    === "MacOS/Linux"
        ```
        Host sib_course_remote
            User username
            HostName 52.57.40.7
            IdentityFile ~/.ssh/key_username.pem
        ```

Finally:

* Save and close the config file
* Click again on the green or blue button in the bottom left corner
* Select `Connect to Host...`, and then `sib_course_remote`. You will be asked which operating system is used on the remote. Choose `Linux`

#### Video tutorial

You can find below a video tutorial showing how to configure SSH connections in VS Code:
<iframe width="560" height="315" src="https://www.youtube.com/embed/cOopQQIL8JU" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

??? info "If you are not working with VS Code"
    If you are not working with VS Code, you can login to the remote server with the following command in a terminal:
    ```sh
    ssh -i key_username.pem username@52.57.40.7
    ```
    If you want to edit files directly on the server, you can mount a directory with `sshfs`.
