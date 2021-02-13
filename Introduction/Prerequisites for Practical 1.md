# Prerequisites for Practical 1

Assuming the issue with pinging `bactsrv` is resolved for next weeks tutorial, we will not have time to cover login troubleshooting.

***

### Logins

As Shane mentioned, your login information is available at the following link: https://docs.google.com/spreadsheets/d/1R51rkNyER1lTuXeccF7yQhZSTjgzrdVHBWT2WRLiFSQ/edit#gid=788379307.

To cover the login once more:

1. On your local machine, type `ssh <username>@bactsrv.nuigalway.ie`. This will prompt your password. If you have not yet changed it, the temporary password for `bactsrv` is `Galway2019!`

2. When you are on `bactsrv` prompt a password change by typing the command `passwd`. Follow the instructions on screen.
3. Log into `lugh` *from `bactsrv`* by typing `ssh <username>@lugh.nuigalway.ie` . Once again you will be prompted for a password, the temporary password is `Galway2020`.
4. Change your `lugh` password using `passwd`.

***

### `.ssh/config`

To cover the `~/.ssh/config` file, all we are doing is creating a shortcut to port to `<username>@bactsrv.nuigalway.ie`.

1. **On your local machine**, create a file `~/.ssh/config` using nano or vim. The contents should look like this:

   ```
   Host bactsrv
   	Hostname bactsrv.nuigalway.ie
   	User bdigby
   	ForwardX11 yes

   Host lugh
   	HostName lugh.nuigalway.ie
   	User bdigby
   	ForwardX11 yes

   Host *
   	ServerAliveInterval 180
   ```

   \* Substitute bdigby for your username.

Personally, I am able to port straight to `lugh` by typing `ssh lugh`. I am not sure if this will work on your machine, but the most important thing is that you are able to log in, even if it is via the scenic route.

Please get this working for next weeks tutorial. If you are having trouble with logging in from now until the next practical, email Shane as he is a system admin for `lugh`.

You are also assumed to be comfortable navigating the file system, so get back up to speed using `cd`, `pwd`, `ls`, `mv`, `cp` etc.

***

### Local installations

Week 1 tutorial assumes you have the following installed **on your local machine**:

1. Docker (https://docs.docker.com/engine/install/ubuntu/#install-using-the-repository)

2. Singularity (https://sylabs.io/guides/3.6/admin-guide/installation.html#id1)

    - Install Dependencies
    - Install GO
    - Install Stable Release (version 3.6.4)
    
    If you have installed singularity version 2.5.2 (check `singularity --version`) please delete the executable and install version 3.6.4 following the final step above. 
    ```bash
    ## delete old version 
    sudo rm -rf /usr/local/libexec/singularity
    ```

3. Anaconda (https://repo.anaconda.com/archive/Anaconda3-2020.11-Linux-x86_64.sh) - This will download an executable `.sh` file. Locate it on the command line (`~/Downloads`?) and type `bash Anaconda3-2020.11-Linux-x86_64.sh` and follow the prompts. You do not need to change the install directory, `/home/username/anaconda3` is fine.

4. Nextflow:

   ```
   curl -s https://get.nextflow.io | bash

   mv nextflow /usr/bin/
   ```

You must have Java 8 or higher installed for this to work:

```bash
sudo apt update
sudo apt install default-jre
java -version
```

You should see that Java 11 has been installed. Any Java version above 8 is fine. 

5. Please create a Dockerhub account. I would advise using a similar username to your one on lugh, and a similar/same password.

We will be creating docker containers locally and pushing them to dockerhub. This cannot be done on lugh because Docker requires sudo privileges which our accounts do not, and will never have. Once on lugh we can pull the images from dockerhub using singularity.

***

### Lugh
When you first log into lugh, you need to run the following commands to 1) Load EasyBuild 2) Set up Java via EasyBuild. 

```bash
module load EasyBuild
eb -r /data/MSc/2021/jdk-8u144-linux-x64.tar.gz   /data/MSc/2021/Java-1.8.0_144.eb
```

or 

```bash
module load EasyBuild
eb -r /data/MSc/2021/jdk-8u144-linux-x64.tar.gz
```

After running this, `module li` should return currently loaded modules:

```bash
Currently Loaded Modules:
 1) autotools  3) gnu7/7.2.0    5) ohpc       7) Java/1.8.0_144
 2) prun/1.2  4) openmpi3/3.0.0  6) EasyBuild/3.4.1
 ```

***

### `~/.bashrc`

We are going to configure your `~/.bashrc` file **on lugh** to automatically load modules (amongst other things) on startup so you dont have to do it manually every time. Open your `~/.bashrc` file using nano or vim and delete all of the contents, replacing them with the file below:

```
# .bashrc

# Source global definitions
if [ -f /etc/bashrc ]; then
	. /etc/bashrc
fi

# Load on LOGIN
module load EasyBuild/3.4.1
module load Java/1.8.0_144
module load singularity
source /home/mscstudent/anaconda3/bin/activate

# on startup, move to this directory
cd /data/MSc/2021

## Colour
case "$TERM" in
    xterm-color|*-256color) color_prompt=yes;;
esac

force_color_prompt=yes

if [ -n "$force_color_prompt" ]; then
    if [ -x /usr/bin/tput ] && tput setaf 1 >&/dev/null; then
        # We have color support; assume it's compliant with Ecma-48
        # (ISO/IEC-6429). (Lack of such support is extremely rare, and such
        # a case would tend to support setf rather than setaf.)
        color_prompt=yes
    else
        color_prompt=
    fi
fi

if [ "$color_prompt" = yes ]; then
    PS1='${debian_chroot:+($debian_chroot)}\[\033[01;32m\]\u@\h\[\033[00m\]:\[\033[01;34m\]\w\[\033[00m\]\$ '
else
    PS1='${debian_chroot:+($debian_chroot)}\u@\h:\w\$ '
fi
unset color_prompt force_color_prompt

if [ -x /usr/bin/dircolors ]; then
    test -r ~/.dircolors && eval "$(dircolors -b ~/.dircolors)" || eval "$(dircolors -b)"
    alias ls='ls --color=auto'
    #alias dir='dir --color=auto'
    #alias vdir='vdir --color=auto'

    alias grep='grep --color=auto'
    alias fgrep='fgrep --color=auto'
    alias egrep='egrep --color=auto'
fi

# bash alias
alias l="ls -lhg --color=auto"
alias la="ls -la --color=auto"
alias lt="ls -lcth --color=auto"
alias lth=" ls -lcth --color=auto | head"
alias ls='ls --color=auto'
alias tarzip="tar -cvzf"
alias tarunzip="tar -xvzf"
alias grep="grep --color=auto"
alias ht='ls -lht | head'
alias vbrc="vi ~/.bashrc"
alias sbrc="source ~/.bashrc"
alias sq="squeue"
alias squ="squeue -u bdigby"
```

Changes instances of `bdigby` to your username (just the last line).

To initialise the changes, type `source ~/.bashrc`.

`module li` should output:

```
Currently Loaded Modules:
 1) autotools  3) gnu7/7.2.0    5) ohpc       7) Java/1.8.0_144
 2) prun/1.2  4) openmpi3/3.0.0  6) EasyBuild/3.4.1  8) singularity/3.4.1`
```

***

### Nextflow on lugh

Similarly to the local install, run the command:

```
  curl -s https://get.nextflow.io | bash

  mv nextflow ~/bin
```

***

### Work Dirs

I have created directories for each of you under `/data/MSc/2021/<username>`

This is where you will run your analyses. Stay out of `/home/<username>` as the home directory is low on storage. Do not ever run an analysis there.

See you all next Thursday
