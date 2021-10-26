
## UNIX

As is stated in the course prerequisites at the [announcement web page](https://www.sib.swiss/training/course/20210413_DOCK). We expect participants to have a basic understanding of working with the command line on UNIX-based systems. You can test your UNIX skills with a quiz [here](https://docs.google.com/forms/d/e/1FAIpQLSd2BEWeOKLbIRGBT_aDEGPce1FOaVYBbhBiaqcaHoBKNB27MQ/viewform?usp=sf_link). If you don't have experience with UNIX command line, or if you're unsure whether you meet the prerequisites, follow our [online UNIX tutorial](https://edu.sib.swiss/pluginfile.php/2878/mod_resource/content/4/couselab-html/content.html).

## Software

Install Docker on your local computer and create an account on [dockerhub](https://hub.docker.com/). You can find instructions [here](https://docs.docker.com/get-docker/). Note that you need admin rights to install and use Docker, and if you are installing Docker on Windows, you need a recent Windows version.

!!! note "If working with Windows"
    During the course exercises you will be mainly interacting with docker through the command line. Although windows powershell is suitable for that, it is easier to follow the exercises if you have UNIX or 'UNIX-like' terminal. You can get this by using [MobaXterm](https://mobaxterm.mobatek.net/ "get MobaXterm") or [WSL2](https://docs.microsoft.com/en-us/windows/wsl/install). Make sure you install the latest versions before installing docker. 

!!! note "If installing Docker is a problem"
    During the course, we can give only limited support for installation issues. If you do not manage to install Docker before the course, you can still do almost all exercises on [Play with Docker](https://labs.play-with-docker.com/). A Docker login is required.

In addition to your local computer, we will be working on an Amazon Web Services ([AWS](https://aws.amazon.com/]))  Elastic Cloud (EC2) server. Our Ubuntu server behaves like a 'normal' remote server, and can be approached through `ssh` with a username, key and IP address. All participants will be granted access to a personal home directory.

Before the course, make sure you can work on a remote server. In MacOS and Linux you can use your terminal for this. On windows, you can use the terminal of the WSL (Windows Subsystem for Linux) or use [MobaXterm](https://mobaxterm.mobatek.net/ "get MobaXterm").
