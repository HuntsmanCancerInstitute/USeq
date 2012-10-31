#set root and change into BioTools directory and launch this shell script to test xxx.useq formats
root=/Users/davidnix/Code/BioTools/USeqTestBinary

testFile=$root/useqTestFile.bed
useqFile=$root/UseqTestFile.useq
useqDir=$root/Useq
realDir=$root/CorrectUSeqOutput
mkdir $useqDir

#Position
java edu/utah/seq/useq/apps/Text2USeq -g H_sapiens_Mar_2006 -c 0 -b 1 -f $testFile
mv $useqFile $useqDir"/Position.useq"

#Position Score
java edu/utah/seq/useq/apps/Text2USeq -g H_sapiens_Mar_2006 -c 0 -b 1 -v 4 -f $testFile
mv $useqFile $useqDir"/PositionScore.useq"

#Position Score Text
java edu/utah/seq/useq/apps/Text2USeq -g H_sapiens_Mar_2006 -c 0 -b 1 -v 4 -t 3 -f $testFile
mv $useqFile $useqDir"/PositionScoreText.useq"

#Position Text
java edu/utah/seq/useq/apps/Text2USeq -g H_sapiens_Mar_2006 -c 0 -b 1 -t 3 -f $testFile
mv $useqFile $useqDir"/PositionText.useq"

#Position Strand
java edu/utah/seq/useq/apps/Text2USeq -g H_sapiens_Mar_2006 -c 0 -b 1 -s 5 -f $testFile
mv $useqFile $useqDir"/PositionStrand.useq"

#Position Score Strand 
java edu/utah/seq/useq/apps/Text2USeq -g H_sapiens_Mar_2006 -c 0 -b 1 -v 4 -s 5 -f $testFile
mv $useqFile $useqDir"/PositionScoreStrand.useq"

#Position Score Text Strand 
java edu/utah/seq/useq/apps/Text2USeq -g H_sapiens_Mar_2006 -c 0 -b 1 -v 4 -t 3 -s 5 -f $testFile
mv $useqFile $useqDir"/PositionScoreTextStrand.useq"

#Position Text Strand 
java edu/utah/seq/useq/apps/Text2USeq -g H_sapiens_Mar_2006 -c 0 -b 1 -t 3 -s 5 -f $testFile
mv $useqFile $useqDir"/PositionTextStrand.useq"


#Region
java edu/utah/seq/useq/apps/Text2USeq -g H_sapiens_Mar_2006 -c 0 -b 1 -e 2 -f $testFile
mv $useqFile $useqDir"/Region.useq"

#Region Score
java edu/utah/seq/useq/apps/Text2USeq -g H_sapiens_Mar_2006 -c 0 -b 1 -v 4 -e 2 -f $testFile
mv $useqFile $useqDir"/RegionScore.useq"

#Region Score Text
java edu/utah/seq/useq/apps/Text2USeq -g H_sapiens_Mar_2006 -c 0 -b 1 -v 4 -t 3 -e 2 -f $testFile
mv $useqFile $useqDir"/RegionScoreText.useq"

#Region Text
java edu/utah/seq/useq/apps/Text2USeq -g H_sapiens_Mar_2006 -c 0 -b 1 -t 3 -e 2 -f $testFile
mv $useqFile $useqDir"/RegionText.useq"

#Region Strand
java edu/utah/seq/useq/apps/Text2USeq -g H_sapiens_Mar_2006 -c 0 -b 1 -s 5 -e 2 -f $testFile
mv $useqFile $useqDir"/RegionStrand.useq"

#Region Score Strand 
java edu/utah/seq/useq/apps/Text2USeq -g H_sapiens_Mar_2006 -c 0 -b 1 -v 4 -s 5 -e 2 -f $testFile
mv $useqFile $useqDir"/RegionScoreStrand.useq"

#Region Score Text Strand 
java edu/utah/seq/useq/apps/Text2USeq -g H_sapiens_Mar_2006 -c 0 -b 1 -v 4 -t 3 -s 5 -e 2 -f $testFile
mv $useqFile $useqDir"/RegionScoreTextStrand.useq"

#Region Text Strand 
java edu/utah/seq/useq/apps/Text2USeq -g H_sapiens_Mar_2006 -c 0 -b 1 -t 3 -s 5 -e 2 -f $testFile
mv $useqFile $useqDir"/RegionTextStrand.useq"


#Convert all back to text
bedDir=$useqDir"/BackToText"
mkdir $bedDir
java edu/utah/seq/useq/apps/USeq2Text -f $useqDir
mv $useqDir"/"*bed $bedDir

#Strip comment lines and sort
cd $bedDir
for x in $(ls *bed)
do
cat $x | grep chr > $x"_Filt"
sort $x"_Filt" > $x"_Filt_Sort"
done

#Compare with real
echo " "
echo "####### Only the file name should print, nothing more ##########"
for x in $(ls *bed_Filt_Sort)
do
echo $x
diff $x $realDir/$x
done

echo " "


