package edu.utah.kegg.api;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import util.gen.IO;
import util.gen.Misc;

public class KeggApiNetwork {

	//fields
	private String networkId;
	private String networkName;
	private String networkType;
	private KeggApiClass[] parentClasses;
	private KeggApiPathway[] pathways = null;
	private KeggApiGene[] genes;
	private HashMap<String, KeggApiGene> geneNameKeggApiGene = null;
	
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append("Id\t"); sb.append(networkId); sb.append("\n");
		sb.append("Name\t"); sb.append(networkName); sb.append("\n");
		sb.append("Type\t"); sb.append(networkType); sb.append("\n");
		if (parentClasses !=null) {
			sb.append("Classes\t"); 
			sb.append(KeggApiClass.toString(parentClasses, ", "));
			sb.append("\n");
		}
		if (pathways !=null) {
			sb.append("Pathways\t"); 
			sb.append(KeggApiPathway.toString(pathways, ", ")); 
			sb.append("\n");
		}
		if (genes !=null) {
			sb.append("Genes\t"); 
			sb.append(KeggApiGene.toString(genes, ", "));
		}
		return sb.toString();
	}

	public KeggApiNetwork (File txt, String networkId) throws IOException {
		this.networkId = networkId;

		String[] lines = IO.loadFileNoTrim(txt);
		for (int i=0; i< lines.length; i++) {
			String l = lines[i];
			//ENTRY, check 
			if (l.startsWith("ENTRY")) {
				if (l.contains(networkId) == false) throw new IOException("ERROR: the network id "+networkId+" doesn't match "+lines[i]);
			}
			//NAME
			else if (l.startsWith("NAME")) {
				String[] f = Misc.WHITESPACE.split(l);
				StringBuilder sb = new StringBuilder(f[1]);
				for (int j=2; j< f.length; j++) {
					sb.append(" ");
					sb.append(f[j]);
				}
				networkName = sb.toString();
			}
			//NAME
			else if (l.startsWith("TYPE")) {
				String[] f = Misc.WHITESPACE.split(l);
				networkType = f[1];
			}
			//PATHWAY
			else if(l.startsWith("PATHWAY")) {
				ArrayList<KeggApiPathway> pIds = new ArrayList<KeggApiPathway>();
				//PATHWAY     hsa00380  Tryptophan metabolism
				String[] f = Misc.WHITESPACE.split(l);
				pIds.add(new KeggApiPathway(f[1], mergeStringArray(f, 2)));
				i++;
				for (; i<lines.length; i++) {
					if (lines[i].startsWith(" ")) {
						//hsa00760  Nicotinate and nicotinamide metabolism
						f = Misc.WHITESPACE.split(lines[i].trim());
						pIds.add(new KeggApiPathway(f[0], mergeStringArray(f, 1)));
					}
					else {
						pathways = new KeggApiPathway[pIds.size()];
						pIds.toArray(pathways);
						i--;
						break;
					}
				}
			}
			//CLASS es
			else if(l.startsWith("CLASS")) {
				ArrayList<KeggApiClass> pIds = new ArrayList<KeggApiClass>();
				//CLASS       nt06513 Complement cascade
				String[] f = Misc.WHITESPACE.split(l);
				pIds.add(new KeggApiClass(f[1], mergeStringArray(f, 2)));
				i++;
				for (; i<lines.length; i++) {
					if (lines[i].startsWith(" ")) {
						//nt06164 Kaposi sarcoma-associated herpesvirus (KSHV)
						f = Misc.WHITESPACE.split(lines[i].trim());
						pIds.add(new KeggApiClass(f[0], mergeStringArray(f, 1)));
					}
					else {
						parentClasses = new KeggApiClass[pIds.size()];
						pIds.toArray(parentClasses);
						i--;
						break;
					}
				}
			}

			//GENE   
			else if (l.startsWith("GENE")) {
				geneNameKeggApiGene = new HashMap<String, KeggApiGene>();
				ArrayList<KeggApiGene> genesAL = new ArrayList<KeggApiGene>();
				//GENE        6999  TDO2; tryptophan 2,3-dioxygenase
				String[] f = Misc.WHITESPACE.split(l);
				KeggApiGene kg = new KeggApiGene(f[2].substring(0, f[2].length()-1), f[1]);
				genesAL.add(kg);
				geneNameKeggApiGene.put(kg.getName(), kg);
				i++;
				for (; i<lines.length; i++) {			
					if (lines[i].startsWith(" ")) {
						//169355  IDO2; indoleamine 2,3-dioxygenase 2
						f = Misc.WHITESPACE.split(lines[i].trim());
						KeggApiGene kg2 = new KeggApiGene(f[1].substring(0, f[1].length()-1), f[0]);
						genesAL.add(kg2);
						geneNameKeggApiGene.put(kg2.getName(), kg2);
					}
					else {
						genes = new KeggApiGene[genesAL.size()];
						genesAL.toArray(genes);
						i--;
						break;
					}
				}
			}
		}
	}

	public static String mergeStringArray(String[] s, int startingWithIndex) {
		if (s.length == 1) return "";
		StringBuilder sb = new StringBuilder(s[startingWithIndex]);
		for (int i=(startingWithIndex+1); i< s.length; i++) {
			sb.append(" ");
			sb.append(s[i]);
		}
		return sb.toString();
	}
	public String getNetworkIdName() {
		return networkId+" : "+networkName;
	}
	public String getNetworkId() {
		return networkId;
	}

	public String getNetworkName() {
		return networkName;
	}

	public KeggApiPathway[] getPathways() {
		return pathways;
	}

	public KeggApiGene[] getGenes() {
		return genes;
	}
	
	public KeggApiClass[] getClasses() {
		return parentClasses;
	}

	public String getPathwayLinks() {
		if (pathways == null) return networkId+"\tNADA\tNADA\n";
		
		StringBuilder sb = new StringBuilder();
		for (KeggApiPathway kp: pathways) {
			//N00003	hsa05221+N00003	Acute myeloid leukemia
			sb.append(networkId); sb.append("\t");
			sb.append(kp.getId()); sb.append("+");
			sb.append(networkId); sb.append("\t");
			sb.append(networkName); sb.append("\n");
		}
		return sb.toString();
	}
	
	public String getPathwayGenes() {
		//REFERENCE_WNT_SIGNALING_PATHWAY	N00056	APC	AXIN1	AXIN2	CCND1	CTNNB1	DVL1...
		//  type           name                id     genes....
		StringBuilder sb = new StringBuilder(networkType);
		sb.append("_");
		sb.append(Misc.WHITESPACE.matcher(networkName).replaceAll("_"));
		String fullName = sb.toString().toUpperCase();
		
		sb = new StringBuilder(fullName);
		sb.append("\t");
		sb.append(networkId);
		for (KeggApiGene kg: genes) {
			sb.append("\t");
			sb.append(kg.getName());
		}
		
		return sb.toString();
		
	}

	public HashMap<String, KeggApiGene> getGeneNameKeggApiGene() {
		return geneNameKeggApiGene;
	}

	public String getNetworkType() {
		return networkType;
	}
	
/* https://www.genome.jp/kegg/document/help_network.html for symbols key
ENTRY       N01723                      Network
NAME        NAD biosynthesis
DEFINITION  Trp -- (TDO2,IDO1/2) >> AFMID -> Kyn -- KMO >> KYNU >> HAAO >> QPRT >> NMNAT >> NADSYN1 -> NAD+
  EXPANDED  C00078 -- (6999,169355,3620) >> 125061 -> C00328 -- 8564 >> 8942 >> 23498 >> 23475 >> (23057,349565,64802) >> 55191 -> C00003
CLASS       nt06513 Complement cascade
            nt06164 Kaposi sarcoma-associated herpesvirus (KSHV)
TYPE        Reference
PATHWAY     hsa00380  Tryptophan metabolism
            hsa00760  Nicotinate and nicotinamide metabolism
GENE        6999  TDO2; tryptophan 2,3-dioxygenase
            169355  IDO2; indoleamine 2,3-dioxygenase 2
            3620  IDO1; indoleamine 2,3-dioxygenase 1
            125061  AFMID; arylformamidase
            8564  KMO; kynurenine 3-monooxygenase
            8942  KYNU; kynureninase
            23498  HAAO; 3-hydroxyanthranilate 3,4-dioxygenase
            23475  QPRT; quinolinate phosphoribosyltransferase
            23057  NMNAT2; nicotinamide nucleotide adenylyltransferase 2
            349565  NMNAT3; nicotinamide nucleotide adenylyltransferase 3
            64802  NMNAT1; nicotinamide nucleotide adenylyltransferase 1
            55191  NADSYN1; NAD synthetase 1
METABOLITE  C00078  L-Tryptophan
            C00328  L-Kynurenine
            C00003  NAD+
REFERENCE   PMID:31954874
  AUTHORS   Castro-Portuguez R, Sutphin GL
  TITLE     Kynurenine pathway, NAD(+) synthesis, and mitochondrial function: Targeting tryptophan metabolism to promote longevity and healthspan.
  JOURNAL   Exp Gerontol 132:110841 (2020)
            DOI:10.1016/j.exger.2020.110841
///
*/

}
