# Setting up software and run jet analysis with JetToyHI framework with Docker

## Windows 
To run the JetToyHI software with Windows we need 3 pieces of software: WSL, Docker and Xming.

### WSL
Windows can not run the needed software natively. Therefore we will use some tricks to emulate a Linux environment. To do this we will use Windows Subsystem for Linux (WSL). You must be running Windows 10 version 2004 and higher (Build 19041 and higher) or Windows 11 to use WSL. 

To install WSL in windows open Powershell and type:
```sh
wsl --install
```

After WSL has finished installing you will need to reboot your pc. After you have rebooted, open WSL and follow the steps to setup your Linux environment inside of WSL. If the installation was succesful you will have a terminal that behaves as a Linux environment.

### Docker
Install docker using the instuctions on the website: https://docs.docker.com/desktop/windows/install/

After you have finished installing docker open WSL and type

```sh
docker run hello-world
```

After downloading the hello-world container docker should say 

```sh
Hello from Docker!                                                                                                                                           
This message shows that your installation appears to be working correctly. 
```

### Xming
In order to make docker work with your graphics adapter (so you can see nice plots) we need to use Xming. You can download Xming from: https://sourceforge.net/projects/xming/

To use Xming we need to know the relevant IP adress. After launching Xming there should be a Xming logo in your system tray. Rightclick the icon and press View Log.
You should now look for ```XdmcpRegisterConnection: newAddress``` in the 10th line. This is followed by an IP adress, we will need this.

Open WSL and type:

```sh
DISPLAY=ipadress:0.0
echo DISPLAY=$DISPLAY >> ~/.bashrc
```

where ipadress is the adress you found in the Xming log, so e.g. ```DISPLAY=1.2.3.4:0.0```

## Running the code
We will use a pre-configured docker container: https://hub.docker.com/repository/docker/bashofman/jettoyhi

To download this container do:

```sh
docker pull bashofman/jettoyhi:latest
```

After the container downloaded succesfully we will be cloning the JetToyHI code. Navigate to a directory into which you would like to download the code and do:

```sh
git clone https://github.com/mverwe/JetToyHI.git
cd JetToyHI
git pull --rebase origin pyt83
```

You should now be in the folder containing the JetToyHI code. From this folder we will be launching docker using:

```sh
docker run -it -v $PWD:/soft/JetToyHI -e DISPLAY=$DISPLAY bashofman/jettoyhi:latest
```

here ```-it``` runs the container interactively, so it stays open and you can work inside of it. ```-v $PWD:/soft/JetToyHI``` binds the folder your code is in (```$PWD```) to a folder inside of the container (```/soft/JetToyHI```). The ```-e``` flag sets the IP adress to send your graphical output to.

After running this command you are inside of the docker container. You can compile your code by using

```sh
compile
```

Now you are done installing software. Let's generate 10 pythia events and run a simple jet analysis.

```sh
./runCreatePythiaEvents -nev 10 -pthat 120 -tune 14
./runSimpleJetAnalysis -hard PythiaEventsTune14PtHat120.pu14  -nev 10
```

You will have produced a root file with a tree. In this tree properties of jets are stored in std::vector format. To check what is inside do:

```sh
root JetToyHIResultSimpleJetAnalysis.root -l
TBrowser b
```

A window should now pop up. Click ```jetTree``` and play around.

## Contribute
* If you want to contribute to this code you need to have a github account. Go here to do so: https://github.com/join.
* Fork the original repository. Go to: https://github.com/mverwe/JetToyHI and click 'Fork' in the upper right corner.
* Instead of cloning the original repository as shown above, clone your own.
* After committing your changes to your own branch, push them to your own fork. Don't know how to do this, ask your colleages or use google which might bring you here https://services.github.com/on-demand/downloads/github-git-cheat-sheet/
* Do a pull request once you have finished your developements.

## Samples
Event samples can be found in the jet quenching CERNBOX:
* From lxplus (CERN account required): /eos/project/j/jetquenching/www
* Webbrowser CERNBOX (CERN account required): https://cernbox.cern.ch/index.php/s/kRy9M7NC9iilE9Z
* Webbrowser (publicly accessible): http://jetquenchingtools.web.cern.ch/JetQuenchingTools/ (You can use wget and curl on this)
* Mount eos on a laptop or local desktop (CERN account required): https://cern.service-now.com/service-portal/article.do?n=KB0003493 

You will find samples from various event generators. For underlying event we have: 'thermal' which is independent particle production using a Boltzmann distribution with a fixed multiplicity and mean p<sub>T</sub> (indicated in the file names). For the hard signal we have PYTHIA8 and JEWEL events with various p<sub>T,hat</sub> settings.

More details about the available samples can be found here: https://jetquenchingtools.github.io/ (public)
(old twiki at cern: https://twiki.cern.ch/twiki/bin/view/JetQuenchingTools/PU14Samples)
