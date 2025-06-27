# Quickstart tutorial for running a combine differential gene expression and differential gene variant KEGG pathway analysis
# 27 June 2025, david.nix@hci.utah.edu, https://github.com/HuntsmanCancerInstitute/USeq

# Download the latest KeggResourcesXXX.zip file and uncompress, see https://github.com/HuntsmanCancerInstitute/USeq/blob/master/Misc/KeggResourcesFeb2025.zip
# If needed, install Java version 11 or higher
# If needed, install base R
# Download and uncompress the latest USeq release, see https://github.com/HuntsmanCancerInstitute/USeq/releases

# Set the paths to match your installation
keggResourcesDir=/Users/u0028003/Downloads/KEGGTest/KeggResourcesFeb2025
testData=$keggResourcesDir/TestDatasets
keggGenes=$keggResourcesDir/GeneSymbolLookup/geneSymbol2KeggGeneInfo1Feb2025.txt
keggNetworks=$keggResourcesDir/Networks
useqApps=/Users/u0028003/Code/USeq/Releases/USeq_9.3.7.beta/Apps
rApp=/usr/local/bin/R
results=$keggResourcesDir/TestResults

# Execute the following to check your paths and resources
ls $keggResourcesDir $testData $keggGenes $keggNetworks $useqApps $rApp &> /dev/null || echo STOP! You are missing required resources
java -version &> /dev/null || echo STOP! Missing Java
$rApp --version &> /dev/null || echo STOP! Missing R
mkdir -p $results
cd $keggResourcesDir

# List the help menu and options
java -jar -Xmx1G $useqApps/KeggGeneAndVariantPathwayAnalyzer

## CmdLine output:
**************************************************************************************
**                 Kegg Gene and Variant Pathway Analyzer : March 2025              **
**************************************************************************************
Runs both the KeggGenePathwayAnalyzer and KeggVariantPathwayAnalyzer applications.
Combines the results at the KEGG Pathway level, coloring each gene for interacive
exploration in the KEGG Pathway Viewer (https://www.kegg.jp). Calculates
a combine Pathway p-value (Fisher's method) and FDR (Benjamini-Hochberg) for Pathways
with multiple significant Networks. Best to use this tool if you are comparing two
large cohorts with both differential gene expression and somatic mutation datasets.

Related Applications:

KeggGeneSymbolIdExtractor - Uses the Kegg API to look up and parse Kegg Gene Ids that
   match each of the provided HUGO Gene Symbols. Use this tool to generate the
   required '-k KEGG gene Id HUGO gene symbol lookup file'.
KeggResourceExtractor - Pulls the latest list of KEGG Networks, their associated
   genes, and pathways. Use this to generate the required '-n KEGG network dir'.
KeggGenePathwayAnalyzer - Differential gene expression KEGG Network and Pathway
   analysis.
KeggVariantPathwayAnalyzer - Differential gene mutation cohort KEGG Network and
   Pathway analysis.
AnnotatedVcfParser - App to select high impact, gain/loss of function gene mutations.
DESeq2, edgeR - R packages for selecting differentially expressed gene sets.

KEGG Gene Color Key:

   Negative Diff Exp LgRto - Light Blue - #87CEEB
   Positive Diff Exp LgRto - Light Red, Pink - #FFB6C1
   Negative Mutation Freq LgRto - Light Violet - #D6B4FC
   Positive Mutation Freq LgRto - Light Orange - #FFD580
   Abs(Mut Freq) < 1.25x - Light Yellow - #FFFFAC
   Genes in red text are disease associated

App Parameters:

-i File containing all of the interrogated genes in your gene expression study, one
     gene per line. Typically > 15K genes
-g File containing the selected genes of interest with their log2Rto, one per line,
     tab delimited, from your study, e.g. differentially expressed, typically <2K
-e File path to the R executable

-a File containing cohort A gene sets, each line represents a subject's genes of
     interest (e.g. those with HIGH impact mutations), tab delimited, the first cell
     is the subject ID, subsequent cells are the gene names.
-b File containing cohort B gene sets, ditto.
-o (Optional) Add one to zero count A or B fractions when calculating the
     log2Rto(fracA/fracB)

-r Directory to save the spreadsheet results
-k KEGG gene Id HUGO gene symbol lookup file, see above
-n KEGG network directory, see above
-t (Optional) Network TYPEs to exclude from testing, comma delimited no spaces.
-m (Optional) Minimum number of interrogated genes in a network for analysis,
     defaults to 4
-x (Optional) Maximum FDR for including networks into the combine Kegg Pathway
     spreadsheet, defaults to 0.15

Example:

java -Xmx1G -jar pathTo/USeq/Apps/KeggGeneAndVariantPathwayAnalyzer
   -i allTestedGenes.txt -g diffExpGenes.txt -a earlyCRC.txt -b lateCRC.txt 
   -k keggIdLookup.txt -n KeggNetworks/ -o -r CombinePathwayAnalyzerResults -t
   'Pathogen,Ev factor'

**************************************************************************************"

# Run the joint analyzer
## If you have just gene expression or just variant analysis use the USeq KeggGenePathwayAnalyzer and KeggVariantPathwayAnalyzer apps respectively.

java -jar -Xmx1G $useqApps/KeggGeneAndVariantPathwayAnalyzer \
  -a $testData/varSeq_EO.txt \
  -b $testData/varSeq_TLO.txt \
  -g $testData/rnaSeqDiffExp_EOvsTLO.txt \
  -i $testData/rnaSeqAllInterrogatedGenes.txt \
  -r $results \
  -k $keggGenes -n $keggNetworks \
  -t 'Pathogen,Ev factor' -e $rApp -o

## CmdLine output:
  
[27 June 2025 14:29] USeq_9.3.7.beta Arguments: -a /Users/u0028003/Downloads/KEGGTest/KeggResourcesFeb2025/TestDatasets/varSeq_EO.txt -b /Users/u0028003/Downloads/KEGGTest/KeggResourcesFeb2025/TestDatasets/varSeq_TLO.txt -g /Users/u0028003/Downloads/KEGGTest/KeggResourcesFeb2025/TestDatasets/rnaSeqDiffExp_EOvsTLO.txt -i /Users/u0028003/Downloads/KEGGTest/KeggResourcesFeb2025/TestDatasets/rnaSeqAllInterrogatedGenes.txt -r /Users/u0028003/Downloads/KEGGTest/KeggResourcesFeb2025/TestResults -k /Users/u0028003/Downloads/KEGGTest/KeggResourcesFeb2025/GeneSymbolLookup/geneSymbol2KeggGeneInfo1Feb2025.txt -n /Users/u0028003/Downloads/KEGGTest/KeggResourcesFeb2025/Networks -t Pathogen,Ev factor -e /usr/local/bin/R -o
Checking required files and directories...
Launching gene and variant pathway analysis...
Done - Gene Pathway Analysis! 4 seconds
Done - Variant Pathway Analysis! 13 seconds
Saving combine pathway results...
Done! 13 seconds

# Examine the results
cd TestResults

## Check the logs for issues
less TestResults/geneAnalysisLog.txt
less variantAnalysisLog.txt

## Open the following spreadsheets in Excel and use the provided KeggLink URLs to launch the KEGG website and explore the results. 
## Keep in mind that KEGG Pathways are an organized collection of KEGG Networks (or sub pathways) and often have a clickable interaction map.
## See the color key above to interpret the highlighted gene boxes in the KEGG Pathway maps

combineGeneVariantPathwaysMinGen4MaxFdr0.15.xls	- combine gene expression and variant results at the KEGG Pathway level
geneNetworksMinGen4.xls - detailed differential gene expression Network results
genePathwaysMinGen4MaxFdr0.15.xls - detailed differential gene expression Pathway results
variantNetworksMinGen4.xls - ditto but for differential variant mutation frequencies
variantPathwaysMinGen4MaxFdr0.15.xls - ditto but for differential variant mutation frequencies


# To build the KeggResource files, run the USeq KeggGeneSymbolIdExtractor and KeggResourceExtractor apps.  KEGG is continuously updating their databases so refresh these resources every 6 months.