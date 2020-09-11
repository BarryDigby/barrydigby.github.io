This page offers a tutorial on how to construct a local ```conda environment``` which hosts a comprehensive toolkit to perform scRNA analysis in python. _This environment is not to be confused with the singularity container used to process the `fastq` reads on lugh._

Firstly, if you have not installed anaconda on your local machine, do so via the following link: https://www.anaconda.com/distribution/

![](https://github.com/BarryDigby/scRNA-Seq/blob/master/images/Screenshot%20from%202020-01-22%2016-48-45.png)

You may choose whichever version of python you want, we will be creating a separate environment with its own version of python regardless. 

Open a terminal and locate the ```Anaconda3-2019.10-Linux-x86_64.sh``` script. It will most likely be in your downloads.. 

``` bash
cd /home/barry/Downloads

bash Anaconda3-2019.10-Linux-x86_64.sh
```

Agree to all of the prompts. ```Anaconda3``` will be installed at your home directory. 


#### If Conda is already installed... 
Now we will create a ```conda environment``` to install all of the necessary tools without interfering with any of our local machines configurations. To do this, we will use the ```sc-tutorial.yml``` file provided in the **<> code** section of this repository. Download the file to your local machine and enter the following into your command line:
``` bash
conda env create -f sc-tutorial.yml
conda activate scRNA
```
You are now in the scRNA environment we created. To invoke a jupyter notebook, type `jupyter notebook`. To exit the environment, type ```conda deactivate```. Remember to activate this environment before attempting to conduct the scRNA analysis!

#### Jupyter notebook
This page is a handy reminder of the most used keyboard shortcuts for ```jupyter notebook``` http://maxmelnick.com/2016/04/19/python-beginner-tips-and-tricks.html
