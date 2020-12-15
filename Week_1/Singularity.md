---
title:
layout: page
permalink: /Week_1/Singularity
---

<center>
<img src="https://raw.githubusercontent.com/BarryDigby/BarryDigby.github.io/master/_images/week1/singularity.png" width="100%" height="100%"/>
</center>


Jump to:
- [Singularity Pull](#pull)
- [Singularity Shell](#shell)

***

# Singularity Pull {#pull}
We have covered creating a docker container and pushing it to Dockerhub. Next, we can use `Singularity` to pull the image locally.

This is a slightly obtuse way of creating containers - utilising both Docker and Singularity - however in my opinion pushing containers to dockerhub is a great way to store containers vs. storing images locally, which can take up several GB of space depending on its contents.

To pull a docker container, use the following singularity command:

```bash
singularity pull --name week1.img docker://USERNAME/REPO:TAG
```

You will now have created an image called `week1.img`.

# Singularity shell {#shell}
To 'enter' the container and to utilise its tools, we use the `singularity shell` command:

```bash
singularity shell -B $(pwd) week1.img
```

The `-B` prefix specifies the *bind path*. Once inside the container, the container cannot "see" files that are in directories above the bind path specified. We used `$(pwd)` which is the same as the current directory path.

Enter your container and test that `fastqc`, `multiqc` and `bbduk.sh` all work (print help messages).
