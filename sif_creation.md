# Using pseudoroot on Navin Lab GEO server

## Define the copykit image to be created

copykit.def
```bash
Bootstrap: docker
From: ubuntu:latest

%environment
	# set up all essential environment variables
	export LC_ALL=C
	export PATH=/opt/miniconda3/bin:$PATH
	export PYTHONPATH=/opt/miniconda3/lib/python3.9/:$PYTHONPATH
	export LC_ALL=C.UTF-8
%post
	# update and install essential dependencies
	apt-get -y update
	apt-get update && apt-get install -y automake \
	build-essential \
	bzip2 \
	wget \
	git \
	default-jre \
	unzip \
	zlib1g-dev \
	parallel \
	libglpk40 \
	gfortran

	# download, install, and update miniconda3
	wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
	bash Miniconda3-latest-Linux-x86_64.sh -b -f -p /opt/miniconda3/
	rm Miniconda3-latest-Linux-x86_64.sh

	# install dependencies via conda
	export PATH="/opt/miniconda3/bin:$PATH"
	conda install -y -c conda-forge mamba 
	mamba install -y -f bioconda::samtools #
	mamba install -y -f bioconda::bedtools #
	mamba install -y -f conda-forge::parallel #

	#install R packages
	mamba install -y -f conda-forge::r-base #=4.2
	mamba install -y -f conda-forge::r-devtools
	#mamba install -y -f bioconda::bioconductor-ensdb.hsapiens.v86

	R --slave -e 'devtools::install_github("navinlabcode/copykit")'
	conda install -y -f --no-deps conda-forge::r-igraph
	conda install -y -f --no-deps bioconda::bioconductor-bluster
	conda install -y -f --no-deps bioconda::bioconductor-copynumber
	conda install -y -f --no-deps bioconda::bioconductor-ggtree
	wget https://github.com/navinlabcode/copykit/releases/download/v.0.1.2/copykit_0.1.2.tar.gz
	R --slave -e 'install.packages("copykit_0.1.2.tar.gz", repos = NULL)' # the install_github is broken so pulling from archive

	R --slave -e 'install.packages("optparse", repos="http://cran.us.r-project.org")' #


%labels
	Author Ryan Mulqueen
	Version v0.1
	MyLabel Copykit 

```

```bash
singularity build --fakeroot copykit.sif copykit.def
```

## Bcl2Fastq Image

bcl2fastq.def
```bash
Bootstrap: docker
From: ubuntu:latest

%environment
	# set up all essential environment variables
	export LC_ALL=C
	export PATH=/opt/miniconda3/bin:$PATH
	export PYTHONPATH=/opt/miniconda3/lib/python3.9/:$PYTHONPATH
	export LC_ALL=C.UTF-8
%post
	# update and install essential dependencies
	apt-get -y update
	apt-get update && apt-get install -y automake \
	build-essential \
	bzip2 \
	wget \
	git \
	default-jre \
	unzip \
	zlib1g-dev \
	parallel \
	libglpk40 \
	gfortran

	# download, install, and update miniconda3
	wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
	bash Miniconda3-latest-Linux-x86_64.sh -b -f -p /opt/miniconda3/
	rm Miniconda3-latest-Linux-x86_64.sh

	# install dependencies via conda
	export PATH="/opt/miniconda3/bin:$PATH"
	conda install -y -c conda-forge mamba 
	mamba install -y -f bioconda::bwa #
	mamba install -y -f bioconda::samtools #
	mamba install -y -f bioconda::bedtools #
	mamba install -y -f conda-forge::parallel #

	#install R packages
	mamba install -y -f conda-forge::r-base=4.1
	mamba install -y -f conda-forge::r-devtools
	mamba install -y -f conda-forge::r-essentials
	mamba install -y -f bih-cubi::bcl2fastq2
	mamba install -y -f conda-forge::r-ggplot2 #
	mamba install -y -f bioconda::seqkit #

	R --slave -e 'install.packages("optparse", repos="http://cran.us.r-project.org")' #
	R --slave -e 'install.packages("patchwork", repos="http://cran.us.r-project.org")' #

	#mamba install -y -f conda-forge::biopython
	mamba install -y -f anaconda::pandas
	conda install -y -f conda-forge::biopython

%labels
	Author Ryan Mulqueen
	Version v0.1
	MyLabel Bcl2Fastq SIF

```

```bash
singularity build --fakeroot bcl2fastq2.sif bcl2fastq2.def

rsync -alPvz \
/home/ubuntu/copykit.sif \
mulqueen@qcprpgeo.mdanderson.edu

#sudo singularity build copykit.sif copykit.def 

#sudo singularity build copykit.sif copykit.def
#sudo singularity shell copykit.sif
```

## Generat scmet analysis sif

```bash
sudo singularity build --sandbox scmetR/ docker://ubuntu:latest
sudo singularity shell --writable scmetR/
#test r installs and stuff (just go line by line of def below)

```

scmetR.def
```bash
Bootstrap: docker
From: ubuntu:latest

%environment
# set up all essential environment variables
export LC_ALL=C
export PATH=/opt/miniconda3/bin:$PATH
export PYTHONPATH=/opt/miniconda3/lib/python3.12/:/opt/miniconda3/lib/python3.9/:$PYTHONPATH

%post
# update and install essential dependencies
apt-get -y update
apt-get update && apt-get install -y automake \
build-essential \
bzip2 \
wget \
git \
default-jre \
unzip \
zlib1g-dev \
parallel

# download, install, and update miniconda3
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -f -p /opt/miniconda3/
rm Miniconda3-latest-Linux-x86_64.sh

# install dependencies via conda
export PATH="/opt/miniconda3/bin:$PATH"
conda config --add channels bioconda
conda config --add channels conda-forge
conda install -y -c conda-forge mamba # general dependencies

#mamba installs
mamba install -y -f pip #
pip install pybedtools
pip install pandas
pip install scipy
mamba install -y -f numpy #
mamba install -y -f bioconda::bwa #
mamba install -y -f bioconda::samtools #
mamba install -y -f bioconda::bedtools #
mamba install -y -f conda-forge::parallel #

#install R packages
mamba install -y -f conda-forge::r-base #
mamba install -y -f conda-forge::r-devtools #
mamba install -y -f conda-forge::icu
mamba install -y -f bioconda::bioconductor-biocgenerics 
mamba install -y -f bioconda::bioconductor-sparsearray
mamba install -y -f bioconda::bioconductor-s4arrays
mamba install -y -f conda-forge::r-biocmanager #

#Funner stuff!
R --slave -e 'install.packages("Seurat", repos="http://cran.us.r-project.org")' #
R --slave -e 'devtools::install_github("stuart-lab/signac", "develop")' #
R --slave -e 'remotes::install_github("satijalab/seurat-wrappers")' #
R --slave -e 'install.packages(c("DescTools", "reshape2", "ggridges", "mice"), repos="http://cran.us.r-project.org")' #

mamba install -y -f conda-forge::r-rlang #
mamba install -y -f conda-forge::r-ggplot2 #
mamba install -y -f bioconda::bioconductor-dirichletmultinomial #
mamba install -y -f conda-forge::r-igraph #
mamba install -y -f conda-forge::r-rjags #
mamba install -y -f conda-forge::r-leiden #
mamba install -y -f conda-forge::r-hdf5r #
mamba install -y -f conda-forge::r-rmpfr #
mamba install -y -f conda-forge::r-ggraph #
mamba install -y -f conda-forge::r-nloptr #
mamba install -y -f conda-forge::r-jomo #

#R utility libraries
R --slave -e 'install.packages("remotes", repos="http://cran.us.r-project.org")' #
R --slave -e 'install.packages("circlize", repos="http://cran.us.r-project.org")' #
R --slave -e 'install.packages("optparse", repos="http://cran.us.r-project.org")' #
R --slave -e 'install.packages("patchwork", repos="http://cran.us.r-project.org")' #
R --slave -e 'install.packages("plyr", repos="http://cran.us.r-project.org")' #
R --slave -e 'install.packages("stringr", repos="http://cran.us.r-project.org")' #
R --slave -e 'install.packages("tidyverse", repos="http://cran.us.r-project.org")' #
R --slave -e 'install.packages("RColorBrewer", repos="http://cran.us.r-project.org")' #

#Bioconductor packages through conda
mamba install -y -f bioconda::bioconductor-biocparallel #
mamba install -y -f bioconda::bioconductor-bsgenome.hsapiens.ucsc.hg38 #
mamba install -y -f bioconda::bioconductor-ensdb.hsapiens.v86 #
mamba install -y -f bioconda::bioconductor-genomicranges
mamba install -y -f bioconda::bioconductor-jaspar2020 #
mamba install -y -f bioconda::bioconductor-org.hs.eg.db #
mamba install -y -f bioconda::bioconductor-tfbstools #
mamba install -y -f bioconda::bioconductor-txdb.hsapiens.ucsc.hg38.knowngene #
mamba install -y -f bioconda::bioconductor-universalmotif #
mamba install -y -f bioconda::bioconductor-chromvar #
mamba install -y -f bioconda::bioconductor-motifmatchr #
mamba install -y -f bioconda::bioconductor-scran #
mamba install -y -f bioconda::bioconductor-complexheatmap #
mamba install -y -f bioconda::bioconductor-biovizbase #

#correct cistopic install and matrix install
cd #
wget https://github.com/aertslab/cisTopic/archive/refs/tags/v2.1.0.tar.gz 
R --slave -e 'install.packages("v2.1.0.tar.gz", repos = NULL)' # the install_github is broken so pulling from archive

#Correct matrix version
mamba install -y -f conda-forge::r-matrix=1.6_1 # set this version
mamba install -y -f r::r-irlba # set this version
R --slave -e 'install.packages("Matrix", type = "source",repos="http://cran.us.r-project.org")' #reinstall from source
R --slave -e 'install.packages("irlba", type = "source",repos="http://cran.us.r-project.org")' #reinstall from source
R --slave -e 'install.packages("SeuratObject", type = "source",repos="http://cran.us.r-project.org")' #reinstall from source
R --slave -e 'oo <- options(repos = "https://cran.r-project.org/");
	tools::package_dependencies("Matrix", which = "LinkingTo", reverse = TRUE)[[1L]];
	install.packages("lme4", type = "source",repos = "http://cran.us.r-project.org");
	options(oo)'

%labels
    Author Ryan Mulqueen
    Version v0.0
    MyLabel Singlecell Methylation

```

```bash
sudo singularity build scmetR.sif scmetR.def
sudo singularity shell scmetR.sif
```

allcool.def
```bash
Bootstrap: docker
From: ubuntu:latest

%environment
# set up all essential environment variables
export LC_ALL=C
export PATH=/opt/miniconda3/bin:$PATH
#export PYTHONPATH=/opt/miniconda3/lib/python3.12/:$PYTHONPATH

%post
# update and install essential dependencies
apt-get -y update
apt-get update && apt-get install -y automake \
build-essential \
bzip2 \
wget \
git \
default-jre \
unzip \
zlib1g-dev \
parallel \
nano

# download, install, and update miniconda3
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -f -p /opt/miniconda3/
rm Miniconda3-latest-Linux-x86_64.sh

# install dependencies via conda
export PATH="/opt/miniconda3/bin:$PATH"
conda config --add channels bioconda
conda config --add channels conda-forge
conda install -y -c conda-forge mamba # general dependencies

#mamba installs
mamba install -y -f pip #
pip install pybedtools
pip install pandas
pip install scipy
mamba install -y -f numpy #
mamba install -y -f bioconda::bwa #
mamba install -y -f bioconda::samtools #
mamba install -y -f bioconda::bedtools #
mamba install -y -f conda-forge::parallel #

#from yap env yaml
mamba create -y -n bismark_env \
python=3.7 \
cutadapt=2.10 \
bismark \
picard \
subread=2.0 \
bowtie2=2.3 \
bowtie=1.3 \
matplotlib

#from allcool env yaml
mamba create -y -n allcool_env \
python=3.8 \
pip \
anndata \
biopython \
dask \
numba \
htslib>=1.9 \
jupyter \
jupyter_contrib_nbextensions \
leidenalg \
natsort \
netCDF4 \
networkx \
opentsne \
plotly \
pybedtools \
pyBigWig \
pynndescent \
pysam \
pytables \
scanpy \
scikit-learn=1.2.2 \
seaborn \
statsmodels \
xarray \
yaml \
zarr \
&& echo -e "#! /bin/bash\n\n# script to activate the conda environment" > ~/.bashrc \
&& conda init bash \
&& echo -e "\nsource activate allcool_env" >> ~/.bashrc \
&& conda clean -y -a \
&& mkdir -p /opt/etc/bashrc \
&& cp ~/.bashrc /opt/etc/bashrc

pip install zarr dask numba papermill imblearn matplotlib==3.6 scanpy opentsne leidenalg scikit-learn==1.2.2 allcools
pip install tables

%labels
    Author Ryan Mulqueen
    Version v0.0
    MyLabel AllCools SIF

```

```bash
singularity build --fakeroot allcool.sif allcool.def #trying on geo

singularity build --fakeroot --sandbox allcool/ docker://ubuntu:latest
singularity shell --fakeroot --writable allcool/

singularity shell allcool.sif
```

## Pull Images
Use sftp to get images off cluster. Move to local computer > geo > seadragon2

```bash

sftp -i ~/Downloads/newkey2.pem ubuntu@54.187.193.117
get copykit.sif

sftp mulqueen@qcprpgeo.mdanderson.edu
cd ./singularity

put copykit.sif
```

Now move to seadragon
```bash
ssh Rmulqueen@seadragon2
bsub -Is -W 6:00 -q transfer -n 1 -M 16 -R rusage[mem=16] /bin/bash #get interactive transfer node this has internet access for environment set up
lcd /rsrch4/home/genetics/rmulqueen/singularity
get scmetR.sif

```

```bash
#singularity shell --bind /rsrch4/home/genetics/rmulqueen/singularity/scmetR.sif
```


cellranger.def
Download cellranger and reference genome.
```bash
wget -O cellranger-8.0.1.tar.gz
wget "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz"

```
```bash
Bootstrap: docker
From: ubuntu:latest

%files
    /volumes/USR2/Ryan/tools/cellranger-8.0.1.tar.gz /opt

%environment
# set up all essential environment variables
export LC_ALL=C
export PATH=/opt/miniconda3/bin:$PATH
export PATH=/opt/cellranger-8.0.1:$PATH

#export PYTHONPATH=/opt/miniconda3/lib/python3.12/:$PYTHONPATH

%post
# update and install essential dependencies
apt-get -y update
apt-get update && apt-get install -y automake \
build-essential \
bzip2 \
wget \
git \
default-jre \
unzip \
zlib1g-dev \
parallel \
nano

# download, install, and update miniconda3
#wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
#bash Miniconda3-latest-Linux-x86_64.sh -b -f -p /opt/miniconda3/
#rm Miniconda3-latest-Linux-x86_64.sh

# install dependencies via conda
#export PATH="/opt/miniconda3/bin:$PATH"
#conda config --add channels bioconda
#conda config --add channels conda-forge
#conda install -y -c conda-forge mamba # general dependencies

#mamba installs
#mamba install -y -f pip #
#pip install pybedtools
#pip install pandas
#pip install scipy
#mamba install -y -f numpy #
#mamba install -y -f bioconda::bwa #
#mamba install -y -f bioconda::samtools #
#mamba install -y -f bioconda::bedtools #
#mamba install -y -f conda-forge::parallel #

tar -xzvf /opt/cellranger-8.0.1.tar.gz  -C /opt/

# #install R packages
# mamba install -y -f conda-forge::r-base #
# mamba install -y -f conda-forge::r-devtools #
# mamba install -y -f conda-forge::icu
# mamba install -y -f bioconda::bioconductor-biocgenerics 
# mamba install -y -f bioconda::bioconductor-sparsearray
# mamba install -y -f bioconda::bioconductor-s4arrays
# mamba install -y -f conda-forge::r-biocmanager #
# #Funner stuff!
# R --slave -e 'install.packages("Seurat", repos="http://cran.us.r-project.org")' #

%labels
    Author Ryan Mulqueen
    Version v0.0
    MyLabel Cellranger SIF

```

```bash
singularity build --fakeroot cellranger.sif cellranger.def #trying on geo

singularity build --fakeroot --sandbox allcool/ docker://ubuntu:latest
singularity shell --fakeroot --writable allcool/

singularity shell allcool.sif
```

