# JetToyHI

## Install on lxplus

First install fastjet and the contrib package (you only have to do this once)
```
source /afs/cern.ch/sw/lcg/contrib/gcc/4.8/x86_64-slc6-gcc48-opt/setup.sh

source /afs/cern.ch/sw/lcg/app/releases/ROOT/6.06.08/x86_64-slc6-gcc48-opt/root/bin/thisroot.sh

#Install fastjet (you just have to do this once)

curl -O http://fastjet.fr/repo/fastjet-3.3.0.tar.gz 
tar zxvf fastjet-3.3.0.tar.gz
cd fastjet-3.3.0/

./configure --prefix=$PWD/../fastjet-install
make
make check
make install
FASTJET=$PWD
cd ..

export FJ_CONTRIB_VER=1.026 
curl -Lo source.tar.gz http://fastjet.hepforge.org/contrib/downloads/fjcontrib-"$FJ_CONTRIB_VER".tar.gz
tar xzf source.tar.gz
cd fjcontrib-"$FJ_CONTRIB_VER"
./configure --fastjet-config=$FASTJET/fastjet-config --prefix=`fastjet-config --prefix`
make 
make install 
make fragile-shared #make shared library
make fragile-shared-install
cd ..
```

Next steps
```
git clone git@github.com:mverwe/JetToyHI.git
cd JetToyHI
. setup.sh
echo $FASTJET > .fastjet
```

```
cd PU14
echo $FASTJET > .fastjet
./mkmk
make
cd ..

scripts/mkcxx.pl -f -s -1 -r -8 '-IPU14' -l '-LPU14 -lPU14 -lz'
make
./runFromFile -hard samples/PythiaEventsTune14PtHat120.pu14 -pileup samples/ThermalEventsMult12000PtAv0.70.pu14 -upu 1 -nev 2
```
