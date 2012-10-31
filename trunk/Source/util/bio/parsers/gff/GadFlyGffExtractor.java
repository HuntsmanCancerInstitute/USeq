package util.bio.parsers.gff;

import java.util.regex.*;
import java.util.*;
import java.io.*;

import util.bio.annotation.*;
import util.gen.*;
/**
 Class for building annotation objects from GffFeatures, basically a converter, specific to the old dmel GadFly annotation.
 GffFeature object are built from the following type of gff file.
 Assumes a specific order of incoming features:
 	(exon(s), translation(optional), transcript) many then a gene/transposon/rRNA etc
 	example gff gadfly format
 	
 	2R	gadfly	exon	67118	67499	.	+	.	genegrp=CG8416; transgrp=CG8416-RD; text=CG8416:8
	2R	gadfly	exon	68758	68907	.	+	.	genegrp=CG8416; transgrp=CG8416-RD; text=CG8416:3
	2R	gadfly	exon	69442	69522	.	+	.	genegrp=CG8416; transgrp=CG8416-RD; text=CG8416:4
	2R	gadfly	exon	69585	70954	.	+	.	genegrp=CG8416; transgrp=CG8416-RD; text=CG8416:6
	2R	gadfly	translation	67344	69773	.	+	.	genegrp=CG8416; transgrp=CG8416-RD
	2R	gadfly	transcript	67118	70954	.	+	.	genegrp=CG8416; transgrp=CG8416-RD; text=CG8416-RD
	2R	gadfly	gene	66411	70954	.	+	.	genegrp=CG8416; text=CG8416; dbxref=GO:0016318; dbxref=GO:0007391; dbxref=GO:0007391; dbxref=GO:0007391; dbxref=GO:0007010; dbxref=GO:0003931; dbxref=GO:0016318; dbxref=GO:0007369; dbxref=GO:0007254; dbxref=GO:0007164; dbxref=GO:0003931; dbxref=GO:0003924; dbxref=GO:0003931; dbxref=GO:0007405; dbxref=GO:0030239; dbxref=GO:0016203; symbol=Rho1; dbxref=FlyBase:FBgn0014020; cytorange=52E3-52E4; cdna_clone=GH20776; cdna_clone=LD03419
	2L	gadfly	exon	47514	52519	.	+	.	genegrp=TE19092; transgrp=TE19092-RA; text=TE19092:1
	2L	gadfly	transcript	47514	52519	.	+	.	genegrp=TE19092; transgrp=TE19092-RA; text=TE19092-RA
	2L	gadfly	transposable_element	47514	52519	.	+	.	genegrp=TE19092; text=TE19092; symbol=jockey{}277; dbxref=FlyBase:FBti0019092; cytorange=21A3-21A3
	2L	gadfly	exon	1791017	1792026	.	-	.	genegrp=CR31930; transgrp=CR31930-RA; text=CR31930:1
	2L	gadfly	exon	1790806	1790958	.	-	.	genegrp=CR31930; transgrp=CR31930-RA; text=CR31930:2
	2L	gadfly	transcript	1790806	1792026	.	-	.	genegrp=CR31930; transgrp=CR31930-RA; text=CR31930-RA
	2L	gadfly	pseudogene	1790806	1792026	.	-	.	genegrp=CR31930; text=CR31930; symbol=Gr22d; dbxref=FlyBase:FBgn0045498; cytorange=22B2-22B2
	
Assumes that there is at least one exon and transcript per gene/transposon/rRNA etc.		
Assumes that start is always less than stop/stop.  Use orientation to get 1 or -1, (+ or -)
 */
public class GadFlyGffExtractor {
	//fields
	private ArrayList geneGrps;
	private ArrayList exons;
	private Translation translation;
	private Transcript transcript;
	private ArrayList transGrps;
	private LinkedHashSet genericFeaturesHash = new LinkedHashSet(); //used to keep track of the generic features
	private ArrayList genericFeaturesAL = new ArrayList();
	private boolean genericFeaturesFound = false;
	/*
	public static void main(String[] args) {
		GadFlyGFFExtractor x = new GadFlyGFFExtractor();
		System.out.println(x);
	}
	*/
	public String toString(){
		StringBuffer x = new StringBuffer();
		for (int i=0; i<geneGrps.size(); i++){
			x.append(((GeneGroup)geneGrps.get(i)).toString()+"\n");
		}
		return x.toString();
	}
	public GadFlyGffExtractor(File gffFile, int startNt, int stopNt) {
		//initialize object arrays
			geneGrps = new ArrayList();
			exons = new ArrayList();
			transGrps = new ArrayList();
			
		//parse gff file
		GffParser anno = new GffParser(gffFile, startNt, stopNt);
		GffFeature[] features = anno.getGffFeatures();
		
		//build annotation objects, note everything is compared to lower case names
		geneGrps = new ArrayList();
		Pattern geneTypes = Pattern.compile("gene|rna|transpos|misc"); //add new gene type items here!
		//loop thru all features
		for (int i = 0; i < features.length; i++) {
			String n = features[i].getFeature();
			//this is actually an exon, gene, transcript, translation, rna etc
			//check if empty
			if (Misc.isEmpty(n)) {
				System.err.println(
					"Fatal error with GFF extraction, empty feature found in GadFlyGFFExtractor");
				System.exit(1);
			}
			n = n.toLowerCase();			
			if (n.equals("exon"))
				buildExon(features[i]);
			else if (n.equals("translation"))
				buildTranslation(features[i]);	
			else if (n.equals("transcript")){
				buildTranscript(features[i]);
				buildTransGrp();
			}
			else if (geneTypes.matcher(n).find())
				buildGeneGrp(features[i]);
			else {
				buildGenericFeature(features[i]);
				genericFeaturesHash.add(n);  //add text to hash to get list 
			} 
		}
		//set whether generic features were found
		if (genericFeaturesHash.size()>0) genericFeaturesFound = true;
	}
	/**Contains all the feature types that were not recognized by the GadFlyGFFExtractor.
	 * Should be user added items like CRMs, enhancers, etc.*/
	public LinkedHashSet getGenericFeatureHash(){
		return genericFeaturesHash;
	}
	public ArrayList getGenericFeatures(){
		return genericFeaturesAL;
	}
	public void buildGenericFeature(GffFeature f){
		genericFeaturesAL.add(new GenericFeature(f));
	}
	
	public void buildExon(GffFeature f) {
		HashMap hash = f.getAttsHash();
		ExonIntron ex = new ExonIntron(f.getStart(), f.getEnd(),(String)hash.get("text"));
		exons.add(ex);
	}
	public void buildTranslation(GffFeature f) {
		HashMap hash = f.getAttsHash();
		translation = new Translation(f.getStart(), f.getEnd(),GeneGroup.convertPlusToNumOrientation(f.getStrand()));
	}
	public void buildTranscript(GffFeature f) {
		HashMap hash = f.getAttsHash();
		transcript = new Transcript(f.getStart(), f.getEnd(),(String)hash.get("text"),GeneGroup.convertPlusToNumOrientation(f.getStrand()));
	}
	public void buildTransGrp(){
		ExonIntron[] x = new ExonIntron[exons.size()];  //note, translation may be null
		transGrps.add(new TransGroup( transcript,  translation, (ExonIntron[])exons.toArray(x)));
		//wipeout exon arraylist, transcript and translation
		exons.clear();
		transcript=null;
		translation=null;
	}
	public void buildGeneGrp(GffFeature f){
		HashMap hash = f.getAttsHash();
		String chrom = f.getSeqName();
		String type = f.getFeature();
		String name = (String)hash.get("text");
		TransGroup[] x = new TransGroup[transGrps.size()];
		String attributes = f.getAttsString();
		geneGrps.add(new GeneGroup(name, chrom, f.getStart(), f.getEnd(), type, (TransGroup[])transGrps.toArray(x),GeneGroup.convertPlusToNumOrientation(f.getStrand()),attributes));
		transGrps.clear();
	}
	public ArrayList getGeneGrps(){
		return geneGrps;
	}
	public boolean isGenericFeaturesFound() {
		return genericFeaturesFound;
	}

}
