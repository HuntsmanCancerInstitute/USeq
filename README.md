# USeq
Applications for analyzing next generation sequencing data from Illumina, SOLiD, and 454 platforms. 
ChIP-seq, RNA-seq, Bis-seq, re-sequencing, SNP INDELs, capture array design tools, IGB/ DAS2/IGV/UCSC 
graph manipulation tools.... GUI and cmd line interface.

## Installing the latest version
1. Download the USeq_xxx.zip file from https://github.com/HuntsmanCancerInstitute/USeq/releases
2. Uncompress
3. Execute an app with no parameters to pull the help menue, e.g. java -jar USeq_9.1.3/Apps/BamConcordance

## Installation in Eclipse for developers
1. Clone the github repo.
2. Download Apache ANT libraries (apache-ant-1.9.6-bin.tar.gz) from http://ant.apache.org/bindownload.cgi and unzip.
3. Set ANT_HOME in your bash profile to point to the unzipped ANT directory. Don't forget to source the file to implement the changes.
4. Download the ant-contrib-0.6 (ant-contrib-0.6-bin.zip) from http://sourceforge.net/projects/ant-contrib/files/ant-contrib/. 
Unzip and move the lib/ant-contrib-0.6.jar to $ANT_HOME/lib directory.
5. In Eclipse, set the Ant HOME: `Window > Preferences > Ant Home … > $ANT_HOME`
6. Refresh the project in Eclipse.

## Compiling USeq
To install USeq into the Releases/ subfolder:

`Right click on the "build.xml" > “Run as … “ > 1 Ant Build.`

This will create an unzipped and zipped archive.
