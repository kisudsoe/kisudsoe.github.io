# R and Python Development Environment Setting Using Docker

By Seungsoo Kim, 2020-08-15-edit2



# 1. Installs

## Download my dev image

Ref: https://docs.anaconda.com/anaconda/user-guide/tasks/docker/

Install and Run Jupyter Notebook

```bash
docker pull kisudsoe/postgwas
```



## Install R

How to update R?

```bash
$ docker run -it kisudsoe/postgwas
container$ apt update && apt -y upgrade
container$ apt -y install r-base
container$ R # Check R version
R> q() # Exit R command
```



## Install Python 3

In container command:

```bash
apt update && apt -y upgrade
apt install software-properties-common -y
apt update && apt install python3.8 -y
python3 ––version
```



## Install pip

In container command:

```bash
apt install python3-pip -y
```



## Install Jupyter Notebook & Jupyter Lab

Ref: https://jupyterlab.readthedocs.io/en/stable/getting_started/installation.html

In container command:

```bash
pip3 install --upgrade pip
pip3 install jupyter
pip3 install jupyterlab
jupyter notebook --version # 6.1.1
jupyter lab --version # 2.2.4
```



## Install IRkernel in R

In R command of the miniconda image:

```R
#install.packages(c('repr', 'IRdisplay', 'crayon', 'pbdZMQ', 'devtools','stringi','Rcpp'))
install.packages('IRkernel')
IRkernel::installspec(name='ir',displayname='R')
```



## Save updated image

In CMD:

```CMD
docker ps #< Get container ID
docker stop 7685afdc7f4d
docker commit -a "jjy" 7685afdc7f4d kisudsoe/postgwas:latest
```





# 2. Run

Error, Ref: https://github.com/codenvy/codenvy/issues/2427

```bash
docker run -v "C:\working\directory\path:/data" -p 8888:8888 kisudsoe/postgwas jupyter lab --ip=0.0.0.0 --port=8888 --allow-root

docker run --rm -p 8888:8888 kisudsoe/postgwas jupyter lab --ip=0.0.0.0 --port=8888 --allow-root
```

Debugging:

```CMD
docker run -it kisudsoe/postgwas
```

