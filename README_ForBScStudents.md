# Setting up software and run jet analysis with JetToyHI framework

## Prerequisites

For the project you will need an envinronment that allows you to compile and run different types of software. In the student room there are linux machines that you can use. If you prefer to use your windows computer, you will need to run linux on a virtual machine or use docker (more modern).

### Docker
First you have to install docker on your laptop. For windows you can download it here: https://docs.docker.com/docker-for-windows/

Once you have installed docker you will have to install some libraries. Follow these instructions that have you can type into the terminal:
```sh
docker pull ubuntu:20.04 #this gets long-term support ubunut version 20.04

docker run -i -t ubuntu /bin/bash

apt update -y

apt install -y vim wget curl libcurl4-gnutls-dev build-essential gfortran cmake libmysqlclient-dev xorg-dev libglu1-mesa-dev libfftw3-dev libssl1.1 libxml2-dev git unzip python3-pip autoconf automake autopoint texinfo gettext libtool libtool-bin pkg-config bison flex libperl-dev libbz2-dev libboost-all-dev swig liblzma-dev libnanomsg-dev libyaml-cpp-dev rsync lsb-release unzip environment-modules

apt-get install xutils-dev libgsl23 libtbb-dev

apt-get install apt-utils libssl-dev
apt-get install libgsl0-dev

pip install matplotlib numpy certifi ipython==7.28.0 ipywidgets ipykernel notebook metakernel pyyaml
```

### Virtual Machine
If you prefer to use a virtual machine running ubuntu, use google to find one that you like. Afterwards you will have to install the same libraries as listed above for docker. Note that you will have to put `sudo` in front of all the `apt` commands.

### MacOS
On a mac there is no need for a virtual machine or docker. But you will need a C++ compiler that you can get by installing xcode to be found in the AppStore.

## JetToyHI installation

Start with creating a directory in which you want to install the software that you will use for this project. For example:
```sh
mkdir soft
```
Now go into the created directory:
```sh
cd soft
```

### Install ROOT
The easiest is to just grep a precompiled version from the root website https://root.cern/install/all_releases/ (take ROOT6). You can do this directly from the terminal:
```sh
wget https://root.cern.ch/download/root_v6.14.04.Linux-ubuntu18-x86_64-gcc7.3.tar.gz  #adjust this line with the appropriate version for you OS (see link above)
tar xvfz  root_v6.14.04.Linux-ubuntu18-x86_64-gcc7.3.tar.gz
rootsetup=$PWD/root/bin/thisroot.sh
. $rootsetup
echo $rootsetup >> ~/.bashrc
```

### Install PYTHIA8.3
```sh
wget http://home.thep.lu.se/~torbjorn/pythia8/pythia8303.tgz
tar xvfz pythia8303.tgz
cd pythia8303
./configure
make
PYTHIA=$PWD
cd ..
```

### Install fastjet

```sh
curl -O http://fastjet.fr/repo/fastjet-3.3.3.tar.gz 
tar zxvf fastjet-3.3.3.tar.gz
cd fastjet-3.3.3/

./configure --prefix=$PWD/../fastjet333-install
make
make check
make install
FASTJET=$PWD/../fastjet333-install
cd ..

export FJ_CONTRIB_VER=1.045 
curl -Lo source.tar.gz http://fastjet.hepforge.org/contrib/downloads/fjcontrib-"$FJ_CONTRIB_VER".tar.gz
tar xzf source.tar.gz
cd fjcontrib-"$FJ_CONTRIB_VER"
./configure --fastjet-config=$FASTJET/bin/fastjet-config --prefix=`$FASTJET/bin/fastjet-config --prefix`
make 
make install 
make fragile-shared #make shared library
make fragile-shared-install
cd ..
```

### Jet workshop software
```sh
git clone https://github.com/mverwe/JetToyHI.git
cd JetToyHI
git pull --rebase origin pyt83

echo `$FASTJET/bin/fastjet-config --prefix` > .fastjet
echo $PYTHIA > .pythia8
```

```sh
cd PU14
echo `$FASTJET/bin/fastjet-config --prefix` > .fastjet
./mkmk
make
cd ..

scripts/mkcxx.pl -f -s -1 -r -8 '-IPU14' -l '-LPU14 -lPU14 -lz'
make
```

Now you are done installing software. Let's generate 10 pythia events and run a simple jet analysis.
```sh
./runCreatePythiaEvents -nev 10 -pthat 120 -tune 14
./runSimpleJetAnalysis -hard PythiaEventsTune14PtHat120.pu14  -nev 10
```

You will have produced a root file with a tree. In this tree properties of jets are stored in std::vector format. To check what is inside do:
```
root JetToyHIResultSimpleJetAnalysis.root -l
TBrowser b
```
Click on `jetTree` and play around.

## Contribute
* If you want to contribute to this code you need to have a github account. Go here to do so: https://github.com/join.
* Fork the original repository. Go to: https://github.com/mverwe/JetToyHI and click 'Fork' in the upper right corner.
* Instead of cloning the original repository as shown above, clone your own.
* After committing your changes to your own branch, push them to your own fork. Don't know how to do this, ask your colleages or use google which might bring you here https://services.github.com/on-demand/downloads/github-git-cheat-sheet/
* Do a pull request once you have finished your developments.

## Running on quark cluster from home
```
ssh -Y [solisID]@gemini.uu.nl
ssh -Y quark.science.uu.nl
```
Now you are remotely logged in.

We will use the centrally installed software for ROOT, pythia and fastjet:
```
module load python/2.7
export PATH=/cm/local/apps/environment-modules/3.2.10/Modules/3.2.10/bin/:$PATH
export ALIBUILD_WORK_DIR=/data1/software/alisoft
alienv enter --shellrc VO_ALICE@pythia::v8243-3,VO_ALICE@ROOT::v6-18-04-alice1-2,VO_ALICE@fastjet::latest-v3.3.3-release,VO_ALICE@GSL::v1.16-4
```

Last step is to install JetToyHI. (note that a new branch was created to make it compatible with the fastjet version that is available on the quark cluster)
```
git clone https://github.com/mverwe/JetToyHI.git
cd JetToyHI
git pull --rebase origin forbsc2

pythia8-config --prefix > .pythia8
fastjet-config --prefix > .fastjet

cd PU14
cp ../.fastjet .
./mkmk
make
cd ..

scripts/mkcxx.pl -f -s -1 -r -8 '-IPU14' -l '-LPU14 -lPU14 -lz'
make
```

### Submitting job to quark cluster
A handy way is to make a submit script. Here I list an example and I name this file runJob.sh:
```
#!/bin/bash
mkdir gqoexjburg
cd gqoexjburg
cp /nethome/verwe121/soft2/JetToyHI/runCreatePythiaEvents .
./runCreatePythiaEvents -nev 10 -tune 14 -pthat 120
ls
```
To submit this job to the quark cluster you type:
```
qsub -V -cwd -N job1 runJob.sh
```
More info on the quark cluster can be found here: https://uugrasp.github.io/UUComputingDocs/quarkCluster.html

## Samples
Event samples can be found in the jet quenching CERNBOX:
* From lxplus (CERN account required): /eos/project/j/jetquenching/www
* Webbrowser CERNBOX (CERN account required): https://cernbox.cern.ch/index.php/s/kRy9M7NC9iilE9Z
* Webbrowser (publicly accessible): https://jetquenchingtools.web.cern.ch/JetQuenchingTools/ (You can use wget and curl on this)
* Mount eos on a laptop or local desktop (CERN account required): https://cern.service-now.com/service-portal/article.do?n=KB0003493 

You will find samples from various event generators. For underlying event we have: 'thermal' which is independent particle production using a Boltzmann distribution with a fixed multiplicity and mean p<sub>T</sub> (indicated in the file names). For the hard signal we have PYTHIA8 and JEWEL events with various p<sub>T,hat</sub> settings.

More details about the available samples can be found here: https://jetquenchingtools.github.io/ (public)
(old twiki at cern: https://twiki.cern.ch/twiki/bin/view/JetQuenchingTools/PU14Samples)


