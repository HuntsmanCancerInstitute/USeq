
package edu.utah.seq.parsers;
import java.io.*;
import java.util.regex.*;

import util.gen.*;

import java.util.*;

import edu.utah.seq.data.sam.*;

/**
 * Parses a dbSNP vcf file for particular rs IDs
 * @author david.nix@hci.utah.edu 
 **/
public class DbSNPExtractor{
	
	//constructors
	public DbSNPExtractor(String[] args){
		try {
			if (args.length==0) System.out.println("dbSNP.vcf.gz   rsIDs.txt   output.vcf.gz");
			else {
				int numFound = 0;
				HashSet<String> rsIDs = IO.loadFileIntoHashSet(new File(args[1]));
				Gzipper outVcf = new Gzipper(new File(args[2]));
				BufferedReader in = IO.fetchBufferedReader(new File(args[0]));
				String line;
				String[] tokens;
				while ((line = in.readLine()) != null){
					if (line.startsWith("#"))outVcf.println(line);
					else {
						tokens = Misc.TAB.split(line);
						if (rsIDs.contains(tokens[2])) {
							numFound++;
							outVcf.println(line);
						}
					}
				}
				System.out.println("Total IDs\t"+rsIDs.size());
				System.out.println("Num found\t"+numFound);
				in.close();
				outVcf.close();
			}

		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public DbSNPExtractor(){
		File t1 = new File("/Users/u0028003/Desktop/ToNormalize/block_norm_clinvar_common_and_clinical_20150603.vcf");
		File t2 = new File("/Users/u0028003/Desktop/ToNormalize/block_norm_clinvar_rsIDsFromSpreadSheet.vcf");
		File t3 = new File("/Users/u0028003/Desktop/ToNormalize/block_norm_clinvar_NoIDsFromSpreadSheet.vcf");
		File t4 = new File("/Users/u0028003/Desktop/ToNormalize/block_norm_dbSNP_Pathogenic_ID.vcf");
		LinkedHashMap<String, String> uniqueVCF = new LinkedHashMap<String, String>();
		addRecords(uniqueVCF, t1, "ClinVarVCF");
		System.out.println (uniqueVCF.size());
		addRecords(uniqueVCF, t2, "ClinVarSpreadSheetRs");
		System.out.println (uniqueVCF.size());
		addRecords(uniqueVCF, t3, "ClinVarSpreadSheetNoID");
		System.out.println (uniqueVCF.size());
		addRecords(uniqueVCF, t4, "DbSNPPathogenic");
		System.out.println (uniqueVCF.size());
		File out = new File ("/Users/u0028003/Desktop/ToNormalize/merge.xls");
		IO.writeHashMap(uniqueVCF, out);
	}
	
	

	private void addRecords(LinkedHashMap<String, String> uniqueVCF, File file, String tier) {
		try {
			BufferedReader in = IO.fetchBufferedReader(file);
			String line;
			String[] tokens;
			while ((line = in.readLine())!= null){
				if (line.startsWith("#")) continue;
				tokens = Misc.TAB.split(line);
				//#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
				String key = tokens[0]+"_"+tokens[1];
				if (uniqueVCF.containsKey(key) == false) uniqueVCF.put(key, "SRC="+tier+"\t"+line);
			}
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	public static void main(String[] args) {
		//new DbSNPExtractor(args);
		new DbSNPExtractor();
	}		


	
}
