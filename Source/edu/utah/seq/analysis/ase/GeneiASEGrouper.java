package edu.utah.seq.analysis.ase;

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.TreeMap;
import java.util.TreeSet;

import util.gen.IO;
import util.gen.Misc;

public class GeneiASEGrouper {

	public static void main(String[] args) {
		
		new GeneiASEGrouper();
		
	}
	
	public GeneiASEGrouper(){
		File f1 = new File("/Users/u0028003/HCI/Labs/DeAngelis/ASE_2.0/CombineAnalysis/MergedResData/AREDS3_RPE");
		File f2 = new File("/Users/u0028003/HCI/Labs/DeAngelis/ASE_2.0/CombineAnalysis/MergedResData/AREDS3_Retina");
		File f3 = new File("/Users/u0028003/HCI/Labs/DeAngelis/ASE_2.0/CombineAnalysis/MergedResData/NeoAMD_RPE");
		File f4 = new File("/Users/u0028003/HCI/Labs/DeAngelis/ASE_2.0/CombineAnalysis/MergedResData/NeoAMD_Retina");
		File f5 = new File("/Users/u0028003/HCI/Labs/DeAngelis/ASE_2.0/CombineAnalysis/MergedResData/Normal_RPE");
		File f6 = new File("/Users/u0028003/HCI/Labs/DeAngelis/ASE_2.0/CombineAnalysis/MergedResData/Normal_Retina");
		File[] dirs = new File[]{f1, f2, f3, f4, f5, f6};
		
		//create groups and id all genes that pass
		ASEGroup[] groups = new ASEGroup[dirs.length];
		TreeSet<String> allGenes = new TreeSet<String>();
		for (int i=0; i< dirs.length; i++) {
			groups[i] = new ASEGroup(dirs[i]);
			allGenes.addAll(groups[i].allGenes);
			System.out.println("\n"+dirs[i]);
			loadSingle(dirs[i]);
		}
		
		//for each gene that passed, count how many times it was seen
		System.out.print("\nGene/Locus");
		for (ASEGroup g : groups){
			System.out.print("\t");
			System.out.print("#Pass "+g.dir.getName());
			System.out.print("\t");
			System.out.print("#Fail "+g.dir.getName());
		}
		System.out.println();
		
		
		for (String pGene: allGenes){
			System.out.print(pGene);
			//for each group
			for (ASEGroup g : groups){
				System.out.print("\t");
				System.out.print(g.getPassFail(pGene));
			}
			System.out.println();
		}
		
	}
	
	private class ASEGroup {
		
		File dir;
		
		//Map of all passing genes observed
		TreeSet<String> allGenes = new TreeSet<String>();
		
		//Map of source and genes that pass and fail
		TreeMap<String, HashSet<String>[]> fileNamePassingGenes = new TreeMap<String, HashSet<String>[]>();
		
		private ASEGroup(File f){
			dir = f;
			
			File[] files = IO.extractFiles(f, ".xls");
			
			for (int i=0; i< files.length; i++){
				String[] lines = IO.loadFile(files[i]);
				HashSet<String> passingGenes = new HashSet<String>();
				HashSet<String> failingGenes = new HashSet<String>();
				
				for (int j=1; j< lines.length; j++){
					String[] t = Misc.TAB.split(lines[j]);
					double fdr = Double.parseDouble(t[4]);
					if (fdr <= 0.1) passingGenes.add(t[0]);
					else failingGenes.add(t[0]);
					
				}
				allGenes.addAll(passingGenes);
				fileNamePassingGenes.put(Misc.removeExtension(files[i].getName()), new HashSet[]{passingGenes, failingGenes});
			}
		}

		public String getPassFail(String pGene) {
			//scan the files
			int pass = 0;
			int fail = 0;
			
			for (String fName: fileNamePassingGenes.keySet()){
				HashSet<String>[] passFail = fileNamePassingGenes.get(fName);
				if (passFail[0].contains(pGene)) pass++;
				if (passFail[1].contains(pGene)) fail++;
			}
			return pass+"\t"+fail;
		}
	}
	
	private void loadSingle(File f){
		//load the files of interest
				File[] files = IO.extractFiles(f, ".xls");
				
				//Map of all passing genes observed
				TreeSet<String> allGenes = new TreeSet<String>();
				
				//Map of source and genes that pass and fail
				TreeMap<String, HashSet<String>[]> fileNamePassingGenes = new TreeMap<String, HashSet<String>[]>();
				
				for (int i=0; i< files.length; i++){
					String[] lines = IO.loadFile(files[i]);
					HashSet<String> passingGenes = new HashSet<String>();
					HashSet<String> failingGenes = new HashSet<String>();
					
					for (int j=1; j< lines.length; j++){
						String[] t = Misc.TAB.split(lines[j]);
						double fdr = Double.parseDouble(t[4]);
						if (fdr <= 0.1) passingGenes.add(t[0]);
						else failingGenes.add(t[0]);
						
					}
					allGenes.addAll(passingGenes);
					fileNamePassingGenes.put(Misc.removeExtension(files[i].getName()), new HashSet[]{passingGenes, failingGenes});
				}
				
				//for each gene that passed, count how many times it was seen
				System.out.println("Gene/Locus\t#Pass\t#Fail\tPassing\tFailing");
				for (String pGene: allGenes){
					System.out.print(pGene +"\t");
					//scan the files
					ArrayList<String> pass = new ArrayList<String>();
					ArrayList<String> fail = new ArrayList<String>();
					
					for (String fName: fileNamePassingGenes.keySet()){
						HashSet<String>[] passFail = fileNamePassingGenes.get(fName);
						if (passFail[0].contains(pGene)) pass.add(fName);
						if (passFail[1].contains(pGene)) fail.add(fName);
					}
					System.out.println (pass.size()+"\t"+fail.size()+"\t"+Misc.stringArrayListToString(pass, ",")+"\t"+Misc.stringArrayListToString(fail, ","));
				}

	}

}
