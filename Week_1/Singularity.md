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

# Singularity Pull {#pull}
We have covered creating a docker image and pushing it to Dockerhub. Next, we can use `Singularity` to pull the image locally, creating an executable container.

This is a slightly obtuse way of creating containers - utilising both Docker and Singularity - however in my opinion pushing images to dockerhub comes with benefits:
1. Automatic linting on Dockerhub to test the stability of the image.
2. Link to github repository, can deploy automatic triggers when a commit is made to the repository.
3. Reduces the storage requirements of the container.

To pull a docker container, use the following singularity command specifying the name of the image to create and the Dockerhub repostiory to pull from:

```bash
singularity pull --name week1.img docker://USERNAME/REPO:TAG
```

You will now have created an image called `week1.img`.

# Singularity shell {#shell}
To 'enter' the container interactively and to utilise its tools, use the `singularity shell` command:

```bash
singularity shell -B $(pwd) week1.img
```

The `-B` prefix specifies the *bind path*. Once inside the container, the container cannot "see" files that are in directories above the bind path specified. We used `$(pwd)` which is the same as the current directory path.

### To Do:
1. Enter the container
2. Check that `fastqc`, `multiqc` and `bbduk.sh` all work.
3. Use `whereis` to print the path of each tool within the container. Can you see where conda installed the environment?
4. Exit the container. 
