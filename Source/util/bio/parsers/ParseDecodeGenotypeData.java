package util.bio.parsers;

import util.gen.*;
import java.io.*;

/**Parses a directory of joined decode data files, ie 'chr10	358624	3	0	0.992874	cnv1071p9' chr,pos,#A,#B,1-pval,SNPName.
 * Tosses snps with a 1-pval < 0.95.*/
public class ParseDecodeGenotypeData {

	public static void main(String[] args) {

		try{
			File directory = new File(args[0]);	
			File[] files = IO.extractFiles(directory,".txt");
		
			System.out.println("Name\t#Pass\t#Fail\tA\tB\tAB\t>2\t>3\t>4\t>5\tPVal");
			
			//for each file/ sample
			for (int i=0; i< files.length; i++){
				//make save Directory
				File saveDir = new File(directory, Misc.removeExtension(files[i].getName())+"_FilteredSGRs");
				saveDir.mkdir();
				BufferedReader in = new BufferedReader(new FileReader (files[i]));
				//files
				File sumFile = new File(saveDir,"Sum.sgr");
				File aFile = new File(saveDir,"A.sgr");
				File bFile = new File(saveDir,"B.sgr");
				File pvalFile = new File(saveDir,"PVal.sgr");
				//PWs
				PrintWriter sumPW = new PrintWriter(new FileWriter(sumFile));
				PrintWriter aPW = new PrintWriter(new FileWriter(aFile));
				PrintWriter bPW = new PrintWriter(new FileWriter(bFile));
				PrintWriter pvalPW = new PrintWriter(new FileWriter(pvalFile));
				
				String line;
				// 0       1        2   3       4          5
				//chr5	85597382	0	2	0.999942	cnv9762p1
				double numberLinesPass = 0;
				double numberLinesFail = 0;
				double totalA = 0;
				double totalB = 0;
				double totalPVal = 0;
				double number2CNVs = 0;
				double number3CNVs = 0;
				double number4CNVs = 0;
				double number5CNVs = 0;
				while ((line=in.readLine())!=null){
					//parse
					String[] bits = line.split("\\s+");
					String chrom = bits[0];
					String pos = bits[1];
					int a = Integer.parseInt(bits[2]);
					int b = Integer.parseInt(bits[3]);
					int ab = a+b;
					float pval = Float.parseFloat(bits[4]);

					//double p = Double.parseDouble(bits[4]);
					//p = -Num.log10(1-p)*10;
					//print lines if pval > 0.95
					if (pval >= 0.95){
						sumPW.println(chrom+"\t"+pos+"\t"+ab);
						aPW.println(chrom+"\t"+pos+"\t"+a);
						bPW.println(chrom+"\t"+pos+"\t"+b);
						pvalPW.println(chrom+"\t"+pos+"\t"+bits[4]);
						numberLinesPass++;
						//increment
						totalA += a;
						totalB += b;
						totalPVal += pval;
						double cnvs = a+b;
						if ( cnvs > 2) {
							number2CNVs++;
							if (cnvs >3){
								number3CNVs++;
								if (cnvs >4){
									number4CNVs++;
									if (cnvs >5) number5CNVs++;
								}
							}
						}
					}
					else numberLinesFail++;
				}
				in.close();
				sumPW.close();
				aPW.close();
				bPW.close();
				pvalPW.close();
				//zip
				IO.zipAndDelete(sumFile);
				IO.zipAndDelete(aFile);
				IO.zipAndDelete(bFile);
				IO.zipAndDelete(pvalFile);
				//make averages
				double ab = (totalA+totalB)/numberLinesPass;
				totalA = totalA/numberLinesPass;
				totalB = totalB/numberLinesPass;
				totalPVal = totalPVal/numberLinesPass;
				//summary line
				System.out.println(files[i].getName()+"\t"+numberLinesPass+"\t"+numberLinesFail+"\t"+totalA+"\t"+totalB+"\t"+ab+"\t"+number2CNVs+"\t"+number3CNVs+"\t"+number4CNVs+"\t"+number5CNVs+"\t"+totalPVal);
			}
			
			
		}catch (Exception e){
			e.printStackTrace();
		}
	}



}
