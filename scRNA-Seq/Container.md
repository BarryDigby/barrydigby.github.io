---
title: Container
layout: page
permalink: /scRNA-Seq/Container
---

Developed for last years tutorial using a singularity definition file. As you can see using a conda yaml file with Docker is much cleaner.

```
Bootstrap: debootstrap
OSversion: bionic
MirrorURL: http://us.archive.ubuntu.com/ubuntu/

%labels
	scRNA-Seq tutorial
	Installs: Kallisto, Bustools and bamtofastq
	This image has been constructed to conduct
	an scRNA-Seq analysis using Kallisto pseudoaligner.  

%environment
	export LC_ALL=C.UTF-8
	export LANG=C.UTF-8
	PATH="/usr/bin/:$PATH"

%post
	apt-get update
        apt-get -y install \
                build-essential \
                wget \
                tar \
                unzip \
                sudo \
                vim \
		git \
                software-properties-common

	sudo add-apt-repository universe
        sudo add-apt-repository restricted
        sudo add-apt-repository multiverse

	mkdir -p repo/
	cd repo/

	wget --no-check-certificate https://github.com/pachterlab/kallisto/releases/download/v0.46.1/kallisto_linux-v0.46.1.tar.gz
	tar -xvzf kallisto_linux-v0.46.1.tar.gz
	sudo chmod 777 /repo/kallisto/kallisto
	sudo cp /repo/kallisto/kallisto /usr/bin/
	rm kallisto_linux-v0.46.1.tar.gz

	wget --no-check-certificate https://github.com/BUStools/bustools/releases/download/v0.39.3/bustools_linux-v0.39.3.tar.gz
	tar -xvzf bustools_linux-v0.39.3.tar.gz
	sudo chmod 777 /repo/bustools/bustools
	sudo cp /repo/bustools/bustools /usr/bin/
	rm bustools_linux-v0.39.3.tar.gz

	wget --no-check-certificate http://cf.10xgenomics.com/misc/bamtofastq-1.2.0
	mv bamtofastq-1.2.0 bamtofastq
	sudo chmod 777 bamtofastq
	sudo cp /repo/bamtofastq /usr/bin/

	wget --no-check-certificate https://www.python.org/ftp/python/3.6.5/    Python-3.6.5.tgz
        tar -xvzf Python-3.6.5.tgz
        cd Python-3.6.5/
        ./configure
        make
        make install

        git clone https://github.com/gpertea/gclib
        git clone https://github.com/gpertea/gffread
        cd gffread/
        make release
        cd ../
        sudo chmod 777 gffread/gffread
        sudo cp gffread/gffread /usr/bin/
```
