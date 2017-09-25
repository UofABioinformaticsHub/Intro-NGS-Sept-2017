# OSX installation

To set up your own computer for today's session, follow these instructions.
Copying and pasting the given code may be the easiest way to make sure everything works.

1. Open a Terminal `Applications -> Utilities -> Terminal`.
2. If you don't have `homebrew` installed, type the following to your terminal. If you *do* have `homebrew` installed, please skip to the next step

```
/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
```

- If you haven't installed the commands `sed` and `wget` in the `Bash` workshop last week, then install them now:

```
brew install gnu-sed --default-names
brew install wget
echo "export PATH=\${PATH}:/usr/local/Cellar/gnu-sed/4.4/bin:/usr/local/Cellar/wget/1.19.1_1/bin" >> ~/.bash_profile
```

- Install NGS tools: `fastqc`, `cutadapt`, `bwa`, `samtools`, `bcftools`, `freebayes`, `sabre`, `IGV`, `picard`

```
brew update
brew upgrade
brew tap homebrew/science
brew install fastqc
brew install cutadapt
brew install bwa
brew install samtools
brew install bcftools
brew install freebayes
brew install igv
brew install picard-tools
```

- Install `sabre`

```
brew install git
git clone https://github.com/najoshi/sabre
cd sabre
make
SABRE_HOME=`pwd`
echo "export PATH=\${PATH}:${SABRE_HOME}" >> ~/.bash_profile
```

After completing the above steps, close the terminal then open it again to start.

[Home](../)
