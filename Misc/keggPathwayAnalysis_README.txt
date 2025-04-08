# Test example for running a combine differential gene expression and differential gene variant KEGG pathway analysis
# 24 March 2025, david.nix@hci.utah.edu, https://github.com/HuntsmanCancerInstitute/USeq

# Download the latest KeggResourcesXXX.zip file and uncompress, see https://github.com/HuntsmanCancerInstitute/USeq/blob/master/Misc/KeggResourcesFeb2025.zip
# If needed, install Java version 11 or higher
# If needed, install base R
# Download and uncompress the latest USeq release, see https://github.com/HuntsmanCancerInstitute/USeq/releases

# Set the paths to match your installation
keggResourcesDir=/Users/u0028003/Code/USeq/Misc/KeggResourcesFeb2025
testData=$keggResourcesDir/TestDatasets
keggGenes=$keggResourcesDir/GeneSymbolLookup/geneSymbol2KeggGeneInfo1Feb2025.txt
keggNetworks=$keggResourcesDir/Networks
useqApps=/Users/u0028003/Code/USeq/Releases/USeq_9.3.6/Apps
rApp=/usr/local/bin/R
results=$keggResourcesDir/TestResults

# Execute the following
set -e
ls $keggResourcesDir $testData $keggGenes $keggNetworks $useqApps $rApp &> /dev/null || echo Missing Required Resources
java -version &> /dev/null || echo Missing Java
$rApp --version &> /dev/null || echo Missing R
mkdir -p $results
cd $keggResourcesDir
set +e

java -jar -Xmx1G $useqApps/KeggGeneAndVariantPathwayAnalyzer \
  -a $testData/varSeq_EO.txt \
  -b $testData/varSeq_TLO.txt \
  -g $testData/rnaSeqDiffExp_EOvsTLO.txt \
  -i $testData/rnaSeqAllInterrogatedGenes.txt \
  -r $results \
  -k $keggGenes -n $keggNetworks \
  -t 'Pathogen,Ev factor' -e $rApp -o

# If you have just gene expression or just variant analysis use the USeq KeggGenePathwayAnalyzer and KeggVariantPathwayAnalyzer apps respectively.

# To build the KeggResources, run the USeq KeggGeneSymbolIdExtractor and KeggResourceExtractor apps.  KEGG is continuously updating their databases so refresh these resources every 6 months.