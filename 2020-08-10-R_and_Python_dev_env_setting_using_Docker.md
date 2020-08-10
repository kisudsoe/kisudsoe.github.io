# R and Python Development Environment Setting Using Docker



# 1. Installs

## Download miniconda image

Ref: https://docs.anaconda.com/anaconda/user-guide/tasks/docker/

Install and Run Jupyter Notebook

```bash
docker run -it -p 8888:8888 continuumio/miniconda /bin/bash -c "/opt/conda/bin/conda install jupyter -y --quiet && mkdir /opt/notebooks && /opt/conda/bin/jupyter notebook --notebook-dir=/opt/notebooks --ip='*' --port=8888 --no-browser"
```



## Install R in the miniconda image

```bash
$ docker run -it continuumio/miniconda
container$ apt update
container$ apt -y upgrade
container$ apt -y install r-base
container$ R # Check R version
R console> q()
```





## Install Jupyter Lab

Ref: https://jupyterlab.readthedocs.io/en/stable/getting_started/installation.html

```bash
$ docker run -it continuumio/miniconda
container$ conda install -c conda-forge jupyterlab
container$ jupyter lab --version
```



Ref: https://jupyterlab.readthedocs.io/en/stable/getting_started/installation.html#docker

Ref: https://jupyter-docker-stacks.readthedocs.io/en/latest/

```bash
docker run -p 8888:8888 jupyter/scipy-notebook:17aba6048f44
```



## Install IRkernel

In R command of the miniconda image:

```R
install.packages(c('repr', 'IRdisplay', 'crayon', 'pbdZMQ', 'devtools','stringi','Rcpp'))
install.packages('IRkernel')
IRkernel::installspec(name='ir',displayname='R')
```

