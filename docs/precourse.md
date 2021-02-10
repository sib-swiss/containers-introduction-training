
## UNIX

As is stated in the course prerequisites at the [announcement web page](https://www.sib.swiss/training/course/20210413_DOCK). We expect participants to have a basic understanding of working with the command line on UNIX-based systems. You can test your UNIX skills with a quiz [here](https://docs.google.com/forms/d/e/1FAIpQLSd2BEWeOKLbIRGBT_aDEGPce1FOaVYBbhBiaqcaHoBKNB27MQ/viewform?usp=sf_link). If you don't have experience with UNIX command line, or if you're unsure whether you meet the prerequisites, follow our [online UNIX tutorial](https://edu.sib.swiss/pluginfile.php/2878/mod_resource/content/4/couselab-html/content.html).

## Software

Install `docker` on your local computer and create an account on [dockerhub](hub.docker.com/). You can find instructions [here](https://docs.docker.com/get-docker/). Note that you need admin rights to install and use `docker`, and if you are installing `docker` on Windows, you need a recent Windows version.

In addition to your local computer, we will be working on an Amazon Web Services ([AWS](https://aws.amazon.com/]))  Elastic Cloud (EC2) server. Our Ubuntu server behaves like a 'normal' remote server, and can be approached through `ssh` with a username, key and IP address. All participants will be granted access to a personal home directory.

Before the course, make sure you can comfortably work on a remote server. This means that you can approach it through the shell, modify scripts and transfer files. We can recommend `atom` for Linux and Mac, and MobaXterm for Windows. Therefore, install on your computer:

=== "mac OS/Linux"
    * SSH and scripting: [Atom](https://atom.io/) with packages like: [`terminus`](https://atom.io/packages/terminus) and [`ftp-remote-edit`](https://atom.io/packages/ftp-remote-edit)
    * Transferring files: [FileZilla](https://filezilla-project.org/)

=== "Windows"
    * SSH and scripting: [MobaXterm](https://mobaxterm.mobatek.net/ "get MobaXterm") and/or [Notepad++](https://notepad-plus-plus.org/downloads/) with the plugin `NppFTP`
    * Transferring files: [FileZilla](https://filezilla-project.org/)
