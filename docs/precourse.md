
## UNIX

As is stated in the course prerequisites at the [announcement web page](https://www.sib.swiss/training/course/20251009_DOCK). We expect participants to have a basic understanding of working with the command line on UNIX-based systems. You can test your UNIX skills with a quiz [here](https://docs.google.com/forms/d/e/1FAIpQLSd2BEWeOKLbIRGBT_aDEGPce1FOaVYBbhBiaqcaHoBKNB27MQ/viewform?usp=sf_link). If you don't have experience with UNIX command line, or if you're unsure whether you meet the prerequisites, follow our [online UNIX tutorial](https://edu.sib.swiss/pluginfile.php/2878/mod_resource/content/4/couselab-html/content.html).

## Software

Install Docker on your local computer and create an account on [dockerhub](https://hub.docker.com/). You can find instructions [here](https://docs.docker.com/get-docker/). Note that you need admin rights to install and use Docker, and if you are installing Docker on Windows, you need a recent Windows version.

!!! note "If working with Windows"
    During the course exercises you will be mainly interacting with docker through the command line. Although windows powershell is suitable for that, it is easier to follow the exercises if you have UNIX or 'UNIX-like' terminal. You can get this by using [MobaXterm](https://mobaxterm.mobatek.net/ "get MobaXterm") or [WSL2](https://docs.microsoft.com/en-us/windows/wsl/install). Make sure you install the latest versions before installing docker. 

!!! note "If installing Docker is a problem"
    During the course, we can give only limited support for installation issues. If you do not manage to install Docker before the course, you can still do almost all exercises on [Play with Docker](https://labs.play-with-docker.com/). A Docker login is required.

### SSH connections

In addition to your local computer, you will be working on an Amazon Web Services ([AWS](https://aws.amazon.com/)) Elastic Cloud (EC2) [server](https://aws.amazon.com/ec2/). This Ubuntu server behaves like a 'normal' remote server, and can be approached through [`ssh`](https://man7.org/linux/man-pages/man1/ssh.1.html). If you are enrolled in the course, you have access to a shared document containing instruction to retrieve your username and private ssh key, granting you access to a personal home directory on the server.

??? tip "If you want to know more about `ssh`"
    If you are not familiar with `ssh`, you can check the [Heidelberg University tutorial](https://zah.uni-heidelberg.de/it-guide/ssh-tutorial-linux) for information on how to set it up and use it.

Here are instructions on how to use VScode to connect with SSH to a remote server. First, place the `key_username.pem` file in the proper folder:

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

* Open VScode and click on the green or blue button in the bottom left corner
* Select `Connect to Host...` and then `Configure SSH Hosts...`
* Specify a location for the SSH config file (preferably the same directory as where your keys are stored): `~/.ssh/config`
* A skeleton config file will be provided. Edit it, so it looks like this (replace `username` with your username, and make sure the IP address in `HostName` match what the one given in the shared document):

    === "Windows"
        ```
        Host sib_course_remote
            User username
            HostName 18.195.137.58
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
            HostName 18.195.137.58
            IdentityFile ~/.ssh/key_username.pem
        ```

Finally:

* Save and close the config file
* Click again on the green or blue button in the bottom left corner
* Select `Connect to Host...`, and then `sib_course_remote`. You will be asked which operating system is used on the remote. Choose `Linux`


!!! warning "Server availability"
    Please note that, for cost reasons, the server will be **started the morning of the course** (or the evening before, at the earliest), and will be **stopped the day after the course ends**, so you will not be able to connect before or after these dates. Likewise, if you need to retrieve data from the server, please do it before it is stopped.

#### Video tutorial

You can find a video tutorial to configure SSH connections in VScode below:
<iframe width="560" height="315" src="https://www.youtube.com/embed/cOopQQIL8JU" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

??? info "If you are not working with VScode"
    If you are not working with VScode, you can login to the remote server with the following command in a terminal:
    ```sh
    ssh -i key_username.pem username@18.195.137.58
    ```
    If you want to edit files directly on the server, you can mount a directory with `sshfs`.
