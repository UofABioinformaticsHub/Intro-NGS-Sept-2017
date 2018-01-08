# Ubuntu installation

To set up your own computer for today's session, follow these instructions.
This is a good example of how software installation is often performed on an ubuntu system.
Copying and pasting the given code may be the easiest way to make sure everything works.

- Open a Terminal and enter these lines one at a time. Enter `y` where required.
- The prefix `sudo` means *super user do* and gives the system authorisation to install software. You'll need to enter your password the first time, but this will only be required once, or until the session times out. If you see a message saying `sudo: unable to resolve host...` just ignore it. (No-one knows why all the VMs do this.)

- This first line just updates the locations for all of the software repositories

```
sudo apt-get update
```

- Now we'll use the command `install` from the `apt-get` set of tools to perform installation fromthe default set of software repositories.
This is similar to `brew install` which may be familiar to OSX users.

```
sudo apt-get install cmake
sudo apt-get install bwa
sudo apt-get install samtools
sudo apt-get install igv
sudo apt-get install picard-tools
```

- To install `FastQC`, we'll place this in the folder `/opt`, which is a common location for *optional* software.
Installing things here can be helpful as you won't accidentally overwrite any important system tools.

```
cd /opt
sudo wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip
sudo unzip fastqc_v0.11.5.zip
echo 'export PATH="/opt/FastQC:$PATH"' >> ~/.bashrc
source ~/.bashrc
```

- To install `bamtools` we need to clone the software using `git`, then use `cmake` & `make`.
(The following can be copied and pasted as a single command.)

```
cd /opt
sudo git clone https://github.com/pezmaster31/bamtools
cd bamtools
sudo mkdir build
cd build
sudo cmake ..
sudo make
sudo make install
echo 'export PATH="/opt/bamtools/lib:$PATH"' >> ~/.bashrc
source ~/.bashrc
```

- To install `bcftools` we also need to clone two git repos before the install.

```
mkdir ~/Downloads
cd ~/Downloads
git clone git://github.com/samtools/htslib.git
git clone git://github.com/samtools/bcftools.git
cd bcftools
make
sudo make install
```

- Next we'll need to install `freebayes` which we'll use later for variant calling

```
cd /opt
sudo git clone --recursive git://github.com/ekg/freebayes.git
cd freebayes
sudo make
echo 'export PATH="/opt/freebayes/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc
```

- We'll also use `sabre` for demultiplexing

```
cd /opt
sudo git clone https://github.com/najoshi/sabre.git
cd sabre
sudo make
echo 'export PATH="/opt/sabre:$PATH"' >> ~/.bashrc
source ~/.bashrc
```

Finally, to install cutadapt, we need to use the `pip` installer instead of `apt-get` which is based on python

```
cd
sudo apt-get install python-pip python-dev build-essential
sudo pip install --upgrade pip
sudo pip install --upgrade virtualenv
sudo pip install --upgrade cutadapt
```

[Home](../)
