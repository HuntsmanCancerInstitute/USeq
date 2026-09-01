package edu.utah.kohli;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.regex.Pattern;

import org.freehep.graphicsio.swf.SWFAction.Call;

import edu.utah.kohli.KohliGene.CopyTestResult;
import util.bio.annotation.Bed;
import util.gen.IO;
import util.gen.Misc;
import util.gen.Num;

public class CompareGATKDatasets {
	
	private TreeMap<String, KohliPatient> patients = new TreeMap<String,KohliPatient>();
	private TreeMap<String, KohliSample> samples = new TreeMap<String,KohliSample>();
	private LinkedHashMap<String, KohliGene> genes = null;
	private File resultsDirectory = null;
	//Just the panel genes
	//private String[] testedGenes = new String[]{"AR-Enh","ARID1A","PTEN","CCND1","ATM","ZBTB16","CDKN1B","KMT2D","CDK4","MDM2","BRCA2","RB1","FOXA1","AKT1","ZFHX3","TP53","NCOR1","CDK12","BRCA1","SPOP","RNF43","MSH2","ERG","TMPRSS2","CHEK2","CTNNB1","FOXP1","PIK3CB","PIK3CA","PIK3R1","CHD1","APC","CDK6","BRAF","KMT2C","NKX3-1","NCOA2","MYC","COL22A1","CDKN2A","NOTCH1","AR","OPHN1"};
	//Just the Zaki set
	//private String[] testedGenes = new String[]{"AR","BRCA2","CHEK2","MYC","NKX3-1","OPHN1","PIK3CA","PIK3CB","TP53","ZBTB16"};
	//1235 OncoKB genes + Panel genes
	private String[] testedGenes = new String[] {"A1CF", "AAMP", "ABCB1", "ABCC3", "ABI1", "ABL1", "ABL2", "ABRAXAS1", "ACACA", "ACKR3", "ACP3", "ACSL3", "ACSL4", "ACSL6", "ACTB", "ACTG1", "ACVR1", "ACVR1B", "ACVR2A", "ADAR", "ADARB2", "ADGRA2", "ADGRG4", "ADHFE1", "AFDN", "AFF1", "AFF3", "AFF4", "AGGF1", "AGK", "AGO1", "AGO2", "AIM2", "AIP", "AJUBA", "AKT1", "AKT1S1", "AKT2", "AKT3", "ALB", "ALDH1A3", "ALDH1L2", "ALDH2", "ALK", "ALOX12B", "ALOX15B", "ALOX5", "AMER1", "ANK1", "ANKRD11", "ANKRD26", "APC", "APCDD1", "APEX1", "APH1A", "APLNR", "APOBEC3B", "AR", "AR-Enh", "ARAF", "ARFRP1", "ARHGAP26", "ARHGAP35", "ARHGEF12", "ARHGEF28", "ARID1A", "ARID1B", "ARID2", "ARID3A", "ARID3B", "ARID3C", "ARID4A", "ARID4B", "ARID5A", "ARID5B", "ARNT", "ASCL1", "ASMTL", "ASPSCR1", "ASXL1", "ASXL2", "ASXL3", "ATF1", "ATG5", "ATIC", "ATM", "ATMIN", "ATP1A1", "ATP2B3", "ATP6AP1", "ATP6V1B2", "ATR", "ATRIP", "ATRX", "ATXN2", "ATXN7", "AURKA", "AURKB", "AURKC", "AXIN1", "AXIN2", "AXL", "B2M", "BAALC", "BABAM1", "BABAM2", "BACH1", "BACH2", "BAP1", "BARD1", "BAX", "BBC3", "BCL10", "BCL11A", "BCL11B", "BCL2", "BCL2L1", "BCL2L11", "BCL2L2", "BCL3", "BCL6", "BCL7A", "BCL9", "BCL9L", "BCOR", "BCORL1", "BCR", "BIRC3", "BIRC5", "BLM", "BMPR1A", "BRAF", "BRCA1", "BRCA2", "BRD2", "BRD3", "BRD4", "BRD9", "BRIP1", "BRSK1", "BTG1", "BTG2", "BTK", "BTLA", "BUB1B", "CACNA1D", "CAD", "CALR", "CAMTA1", "CANT1", "CARD11", "CARM1", "CARS1", "CASP8", "CASR", "CAV1", "CBFA2T3", "CBFB", "CBL", "CBLB", "CBLC", "CCAR1", "CCDC6", "CCN6", "CCNB1IP1", "CCNB3", "CCND1", "CCND2", "CCND3", "CCNE1", "CCNQ", "CCT6B", "CD19", "CD22", "CD274", "CD276", "CD28", "CD36", "CD58", "CD70", "CD74", "CD79A", "CD79B", "CDC42", "CDC73", "CDH1", "CDH11", "CDH2", "CDH4", "CDK12", "CDK4", "CDK6", "CDK8", "CDKN1A", "CDKN1B", "CDKN1C", "CDKN2A", "CDKN2B", "CDKN2C", "CDX2", "CEBPA", "CENPA", "CEP43", "CHCHD7", "CHD1", "CHD2", "CHD4", "CHEK1", "CHEK2", "CHIC2", "CHN1", "CHTF8", "CHUK", "CIC", "CIITA", "CILK1", "CKS1B", "CLIP1", "CLP1", "CLTC", "CLTCL1", "CMTR2", "CNBP", "CNOT3", "CNTRL", "COL1A1", "COL22A1", "COL2A1", "COL5A1", "COP1", "CPS1", "CRBN", "CREB1", "CREB3L1", "CREB3L2", "CREBBP", "CREM", "CRKL", "CRLF2", "CRTC1", "CRTC3", "CSDE1", "CSF1", "CSF1R", "CSF3R", "CTC1", "CTCF", "CTDNEP1", "CTLA4", "CTNNA1", "CTNNB1", "CTR9", "CUL3", "CUL4A", "CUX1", "CXCR4", "CYLD", "CYP17A1", "CYP19A1", "CYSLTR2", "D2HGDH", "DAXX", "DAZAP1", "DCTN1", "DCUN1D1", "DDB2", "DDIT3", "DDR1", "DDR2", "DDX10", "DDX3X", "DDX4", "DDX41", "DDX5", "DDX6", "DEK", "DHX15", "DICER1", "DIS3", "DIS3L2", "DKK1", "DKK2", "DKK3", "DKK4", "DLL3", "DNAJB1", "DNM2", "DNMT1", "DNMT3A", "DNMT3B", "DOT1L", "DPYD", "DROSHA", "DTX1", "DUSP2", "DUSP22", "DUSP4", "DUSP9", "E2F3", "EBF1", "ECSIT", "ECT2L", "EED", "EGFL7", "EGFR", "EGR1", "EGR2", "EIF1AX", "EIF2B1", "EIF3E", "EIF4A2", "EIF4E", "ELF3", "ELF4", "ELK4", "ELL", "ELL2", "ELN", "ELOC", "ELP2", "EML4", "EMSY", "EP300", "EP400", "EPAS1", "EPCAM", "EPHA3", "EPHA5", "EPHA7", "EPHB1", "EPHB4", "EPOR", "EPS15", "ERBB2", "ERBB3", "ERBB4", "ERC1", "ERCC1", "ERCC2", "ERCC3", "ERCC4", "ERCC5", "ERCC6", "ERF", "ERG", "ERRFI1", "ESCO1", "ESCO2", "ESR1", "ETAA1", "ETNK1", "ETS1", "ETV1", "ETV4", "ETV5", "ETV6", "EWSR1", "EXOSC6", "EXT1", "EXT2", "EZH1", "EZH2", "EZHIP", "EZR", "FAF1", "FANCA", "FANCB", "FANCC", "FANCD2", "FANCE", "FANCF", "FANCG", "FANCI", "FANCL", "FANCM", "FAS", "FAT1", "FAT4", "FBXO11", "FBXO31", "FBXW2", "FBXW7", "FCGR2B", "FCRL4", "FES", "FEV", "FGF1", "FGF10", "FGF12", "FGF14", "FGF19", "FGF2", "FGF23", "FGF3", "FGF4", "FGF5", "FGF6", "FGF7", "FGF8", "FGF9", "FGFR1", "FGFR2", "FGFR3", "FGFR4", "FH", "FHIT", "FIP1L1", "FLCN", "FLI1", "FLT1", "FLT3", "FLT4", "FLYWCH1", "FNBP1", "FOLH1", "FOLR1", "FOXA1", "FOXF1", "FOXL2", "FOXN4", "FOXO1", "FOXO3", "FOXO4", "FOXP1", "FRS2", "FSTL1", "FSTL3", "FUBP1", "FURIN", "FUS", "FYN", "FZR1", "GAB1", "GAB2", "GABRA6", "GADD45B", "GAS7", "GATA1", "GATA2", "GATA3", "GATA4", "GATA6", "GEN1", "GID4", "GLI1", "GMPS", "GNA11", "GNA12", "GNA13", "GNAQ", "GNAS", "GNB1", "GOLGA5", "GOPC", "GPC3", "GPHN", "GPS2", "GRB7", "GREM1", "GRIN2A", "GRM3", "GSK3B", "GSTO1", "GSTP1", "GTF2I", "GTSE1", "H1-2", "H1-3", "H1-4", "H1-5", "H2AC11", "H2AC16", "H2AC17", "H2AC6", "H2BC11", "H2BC12", "H2BC17", "H2BC4", "H2BC5", "H2BC8", "H3-3A", "H3-3B", "H3-4", "H3-5", "H3C1", "H3C10", "H3C11", "H3C12", "H3C13", "H3C14", "H3C15", "H3C2", "H3C3", "H3C4", "H3C6", "H3C7", "H3C8", "H3P6", "H4C6", "H4C9", "HDAC1", "HDAC2", "HDAC4", "HDAC7", "HERPUD1", "HEY1", "HFE", "HGF", "HIF1A", "HIP1", "HIRA", "HLA-A", "HLA-B", "HLA-C", "HLF", "HMGA1", "HMGA2", "HNF1A", "HNF1B", "HNRNPA2B1", "HOOK3", "HOXA11", "HOXA13", "HOXA3", "HOXA9", "HOXB13", "HOXC11", "HOXC13", "HOXD11", "HOXD13", "HRAS", "HSD17B2", "HSD3B1", "HSP90AA1", "HSP90AB1", "HTATIP2", "ICOSLG", "ID1", "ID3", "IDH1", "IDH2", "IFNAR1", "IFNGR1", "IGF1", "IGF1R", "IGF2", "IGH", "IGK", "IGL", "IKBKB", "IKBKE", "IKZF1", "IKZF2", "IKZF3", "IL10", "IL2", "IL21R", "IL3", "IL6ST", "IL7R", "ING1", "INHA", "INHBA", "INPP4A", "INPP4B", "INPP5D", "INPPL1", "INSR", "INTS6", "IQGAP1", "IRF1", "IRF2", "IRF4", "IRF8", "IRS1", "IRS2", "IRS4", "ITK", "ITPKB", "JAK1", "JAK2", "JAK3", "JARID2", "JAZF1", "JUN", "KAT6A", "KAT6B", "KAT7", "KBTBD4", "KCNJ5", "KDM2B", "KDM4C", "KDM5A", "KDM5C", "KDM5D", "KDM6A", "KDR", "KDSR", "KEAP1", "KEL", "KIAA1549", "KIF5B", "KIT", "KLF2", "KLF3", "KLF4", "KLF5", "KLF6", "KLHL6", "KLK2", "KMT2A", "KMT2B", "KMT2C", "KMT2D", "KMT5A", "KNL1", "KNSTRN", "KRAS", "KSR2", "KTN1", "L1TD1", "LARP4B", "LASP1", "LATS1", "LATS2", "LCK", "LCP1", "LEF1", "LGR5", "LIFR", "LMNA", "LMO1", "LMO2", "LPP", "LRIG3", "LRP1B", "LRP5", "LRP6", "LRRK2", "LTB", "LTK", "LYL1", "LYN", "LZTR1", "LZTS1", "MAD2L2", "MAF", "MAFB", "MAGED1", "MAGI2", "MAL2", "MALT1", "MAML2", "MAP2K1", "MAP2K2", "MAP2K4", "MAP3K1", "MAP3K13", "MAP3K14", "MAP3K21", "MAP3K6", "MAP3K7", "MAP4K4", "MAPK1", "MAPK3", "MAPKAP1", "MASTL", "MAX", "MBD4", "MBD6", "MCL1", "MDC1", "MDH2", "MDM2", "MDM4", "MDS2", "MECOM", "MED12", "MEF2B", "MEF2C", "MEF2D", "MEN1", "MERTK", "MET", "MGA", "MGAM", "MIB1", "MIDEAS", "MITF", "MKI67", "MKNK1", "MLF1", "MLH1", "MLH3", "MLLT1", "MLLT10", "MLLT11", "MLLT3", "MLLT6", "MN1", "MNX1", "MOB3B", "MPEG1", "MPL", "MRE11", "MRTFA", "MS4A1", "MSH2", "MSH3", "MSH6", "MSI1", "MSI2", "MSN", "MST1", "MST1R", "MTAP", "MTCP1", "MTHFD2", "MTHFR", "MTOR", "MUC1", "MUTYH", "MYB", "MYBL1", "MYC", "MYCL", "MYCN", "MYD88", "MYH11", "MYH9", "MYO18A", "MYO5A", "MYOD1", "NAB2", "NACA", "NADK", "NBEAP1", "NBN", "NCOA1", "NCOA2", "NCOA3", "NCOA4", "NCOR1", "NCOR2", "NCSTN", "NDRG1", "NEGR1", "NF1", "NF2", "NFATC2", "NFE2", "NFE2L2", "NFIB", "NFKB2", "NFKBIA", "NFKBIE", "NHERF1", "NIN", "NKX2-1", "NKX3-1", "NOD1", "NONO", "NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4", "NPM1", "NQO1", "NR4A3", "NRAS", "NRG1", "NSD1", "NSD2", "NSD3", "NT5C2", "NTHL1", "NTRK1", "NTRK2", "NTRK3", "NUF2", "NUMA1", "NUP214", "NUP93", "NUP98", "NUTM1", "NUTM2A", "NUTM2B", "NUTM2D", "OGA", "OGT", "OLIG2", "OMD", "ONECUT2", "OPHN1", "P2RY8", "PAFAH1B2", "PAG1", "PAK1", "PAK3", "PAK5", "PALB2", "PARP1", "PARP2", "PARP3", "PASK", "PATZ1", "PAX3", "PAX5", "PAX7", "PAX8", "PBRM1", "PBX1", "PC", "PCBP1", "PCLO", "PCM1", "PCSK7", "PDCD1", "PDCD11", "PDCD1LG2", "PDE4DIP", "PDGFB", "PDGFRA", "PDGFRB", "PDK1", "PDPK1", "PDS5B", "PER1", "PGBD5", "PGR", "PHF1", "PHF19", "PHF6", "PHLPP1", "PHLPP2", "PHOX2B", "PICALM", "PIGA", "PIK3C2B", "PIK3C2G", "PIK3C3", "PIK3CA", "PIK3CB", "PIK3CD", "PIK3CG", "PIK3R1", "PIK3R2", "PIK3R3", "PIM1", "PLAG1", "PLCG1", "PLCG2", "PLK2", "PMAIP1", "PML", "PMS1", "PMS2", "PNRC1", "POLD1", "POLE", "POLG", "POLH", "POLQ", "POT1", "POU2AF1", "POU2F2", "POU3F2", "POU3F4", "POU5F1", "PPARG", "PPFIBP1", "PPM1D", "PPP1CB", "PPP2R1A", "PPP2R2A", "PPP4R2", "PPP6C", "PRCC", "PRDM1", "PRDM14", "PRDM16", "PREX2", "PRF1", "PRKACA", "PRKAR1A", "PRKCB", "PRKCI", "PRKD1", "PRKDC", "PRKN", "PRPF8", "PRRX1", "PRSS1", "PRSS8", "PSIP1", "PSMB2", "PTCH1", "PTEN", "PTK6", "PTK7", "PTP4A1", "PTPN1", "PTPN11", "PTPN13", "PTPN14", "PTPN2", "PTPN6", "PTPRB", "PTPRC", "PTPRD", "PTPRK", "PTPRO", "PTPRS", "PTPRT", "PUM1", "QKI", "RAB35", "RABEP1", "RAC1", "RAC2", "RAD17", "RAD21", "RAD50", "RAD51", "RAD51B", "RAD51C", "RAD51D", "RAD52", "RAD54L", "RAF1", "RALGDS", "RANBP2", "RAP1GDS1", "RARA", "RASA1", "RASGEF1A", "RB1", "RBM10", "RBM15", "RECQL", "RECQL4", "REL", "RELA", "RELN", "REST", "RET", "REV3L", "RHEB", "RHOA", "RHOH", "RICTOR", "RIOK2", "RIT1", "RMI2", "RNASEH2A", "RNASEH2B", "RNF213", "RNF217-AS1", "RNF43", "ROBO1", "ROS1", "RPL10", "RPL22", "RPL5", "RPN1", "RPS15", "RPS6KA4", "RPS6KB1", "RPS6KB2", "RPTOR", "RRAGC", "RRAS", "RRAS2", "RSPO2", "RSPO3", "RTEL1", "RUNX1", "RUNX1T1", "RUNX2", "RXRA", "RYBP", "S1PR2", "SALL4", "SAMD9", "SAMD9L", "SAMHD1", "SBDS", "SCG5", "SDC4", "SDHA", "SDHAF2", "SDHB", "SDHC", "SDHD", "SEC31A", "SEPTIN5", "SEPTIN6", "SEPTIN9", "SERP2", "SERPINB3", "SERPINB4", "SESN1", "SESN2", "SESN3", "SET", "SETBP1", "SETD1A", "SETD1B", "SETD2", "SETD3", "SETD4", "SETD5", "SETD6", "SETD7", "SETDB1", "SETDB2", "SF3B1", "SF3B2", "SFPQ", "SFRP1", "SFRP2", "SFRP4", "SGK1", "SH2B3", "SH2D1A", "SH3GL1", "SHOC2", "SHQ1", "SIX1", "SLC1A2", "SLC34A2", "SLC45A3", "SLFN11", "SLIT2", "SLIT3", "SLX4", "SMAD2", "SMAD3", "SMAD4", "SMARCA1", "SMARCA2", "SMARCA4", "SMARCB1", "SMARCD1", "SMARCE1", "SMC1A", "SMC3", "SMG1", "SMO", "SMYD3", "SNCAIP", "SND1", "SNX29", "SOCS1", "SOCS2", "SOCS3", "SOS1", "SOX10", "SOX17", "SOX2", "SOX9", "SP140", "SPEN", "SPOP", "SPRED1", "SPRTN", "SQSTM1", "SRC", "SRP72", "SRSF2", "SRSF3", "SS18", "SS18L1", "SSX1", "SSX2", "SSX4", "STAG1", "STAG2", "STAT1", "STAT2", "STAT3", "STAT4", "STAT5A", "STAT5B", "STAT6", "STIL", "STK11", "STK19", "STK40", "STRN", "SUFU", "SUZ12", "SYK", "SZT2", "TACSTD2", "TAF1", "TAF15", "TAL1", "TAL2", "TAP1", "TAP2", "TBL1XR1", "TBX3", "TCEA1", "TCF12", "TCF3", "TCF7L2", "TCL1A", "TCL1B", "TEC", "TEK", "TENT5C", "TERC", "TERT", "TET1", "TET2", "TET3", "TFE3", "TFEB", "TFG", "TFPT", "TFRC", "TGFBR1", "TGFBR2", "TIGAR", "TIPARP", "TLE1", "TLE2", "TLE3", "TLE4", "TLL2", "TLX1", "TLX3", "TMEM127", "TMEM30A", "TMPRSS2", "TMSB4XP8", "TNFAIP3", "TNFRSF11A", "TNFRSF14", "TNFRSF17", "TNFRSF9", "TNFSF13", "TONSL", "TOP1", "TOP2A", "TOX", "TP53", "TP53BP1", "TP63", "TPM3", "TPM4", "TPMT", "TPR", "TRA", "TRAF2", "TRAF3", "TRAF5", "TRAF7", "TRB", "TRD", "TRG", "TRIB3", "TRIM24", "TRIM27", "TRIM33", "TRIP11", "TRIP13", "TRRAP", "TSC1", "TSC2", "TSHR", "TTL", "TUSC3", "TYK2", "TYMS", "TYRO3", "U2AF1", "U2AF2", "UBA1", "UBE2A", "UBR5", "UBTF", "UCHL1", "UPF1", "USP1", "USP6", "USP8", "VAV1", "VAV2", "VEGFA", "VEGFB", "VHL", "VTCN1", "WAS", "WDCP", "WDR90", "WEE1", "WIF1", "WRN", "WT1", "WWP1", "WWTR1", "XBP1", "XIAP", "XPA", "XPC", "XPO1", "XRCC1", "XRCC2", "XRCC3", "YAP1", "YES1", "YPEL5", "YWHAE", "YY1", "YY1AP1", "ZBTB16", "ZBTB20", "ZBTB7A", "ZFHX3", "ZFP36L1", "ZFP36L2", "ZMYM2", "ZMYM3", "ZNF217", "ZNF24", "ZNF292", "ZNF331", "ZNF384", "ZNF521", "ZNF703", "ZNF750", "ZNRF3", "ZRSR2"};	private HashMap<String, String> geneNameInfo = null;
	private Bed arEnhancer = new Bed("chrX", 66899317, 66907698);
	private HashMap<String, String> vcfKeyRecord = null;
	private File fullPathToR = new File("/usr/local/bin/R");
	private File saveDir = null;

	public static void main(String[] args) throws IOException {
		File patientSampleInfo = new File ("/Users/u0028003/HCI/Labs/Kohli_Manish/PriorTo2026/ProstateCfDNAAnalysis/FinalAggregateAnalysis/sampleMatchingPersonGNomIDs_13June2024.txt");
		File geneInfo = new File ("/Users/u0028003/HCI/Labs/Kohli_Manish/PriorTo2026/ProstateCfDNADesign/CBioReviewProstateCancer23July2024/targetGeneSummary.txt");
		//File gatkResultsDir = new File ("/Users/u0028003/HCI/Labs/Kohli_Manish/PriorTo2026/ProstateCfDNAAnalysis/WgsAnalysis/GATKForWGS/SubSamplingComparision/SubSampled");
		File gatkResultsDir = new File ("/Users/u0028003/HCI/Labs/Kohli_Manish/2026/WgsCopyAnalysis/Subsampling/PassingBeds/MergedWindowSets/Splits");
		//File gatkKeyResultsDir = new File("/Users/u0028003/HCI/Labs/Kohli_Manish/PriorTo2026/ProstateCfDNAAnalysis/WgsAnalysis/GATKForWGS/SubSamplingComparision/S0_100_PassBed");
		File gatkKeyResultsDir = new File("/Users/u0028003/HCI/Labs/Kohli_Manish/2026/WgsCopyAnalysis/Subsampling/PassingBeds/MergedWindowSets/S0");
		new CompareGATKDatasets(patientSampleInfo, geneInfo, gatkResultsDir, gatkKeyResultsDir);
	}
	
	public CompareGATKDatasets(File patientSampleInfo, File geneInfo, File gatkResultsDir, File gatkKeyResultsDir) throws IOException {

		//load patients and samples
		loadPatients(patientSampleInfo);
		for (KohliPatient kp: patients.values()) IO.pl(kp);

		//load gene info
		loadGeneInfo(geneInfo);

		//load the key
		loadGatkCopyRatioCalls(gatkKeyResultsDir, "Key");

		//load dirs of cnv results
		File[] splitDirs = IO.extractOnlyDirectories(gatkResultsDir);
		for (File sd: splitDirs) loadGatkCopyRatioCalls(sd, sd.getName());

		compareResultsWithMergedGatkWgsConfusionMatrix();

		IO.pl("\nCOMPLETE!");

	}
	
	




	private void loadGeneInfo(File geneInfo) {
		String[] lines = IO.loadFileIntoStringArray(geneInfo);
		geneNameInfo = new HashMap<String,String>();
		for (String l: lines) {
			if (l.startsWith("#") || l.length() ==0) continue;
			String[] f = Misc.TAB.split(l);
			IO.pl(l);
			geneNameInfo.put(f[0], f[1]+" -> "+f[2]);
		}
		
	}

	private void compareResultsWithMergedGatkWgsConfusionMatrix() {
		IO.pl("\nComparing copy number results...");
		IO.pl("Header:\t"+ConfusionMatrix.toStringHeader());
		//for each patient
		for (KohliPatient kp: patients.values()) {

			ArrayList<KohliSample> nonGermline = kp.getNonGermlineSamples();
			
			//for each cfDNA sample
			for (KohliSample ks: nonGermline) {
				IO.p(kp.getHciPatientId()+"\t"+ks.getSampleId());
				//find the key
				CnvCallSet key = null;
				
				for (CnvCallSet cc: ks.getCnvCallSets()) {
					if (cc.getMethodName().equals("Key")) {
						key = cc;
						break;
					}
				}
				if (key == null) {
					IO.pl("\tKey not found....");
					continue;
				}

				String keyGenes = key.getCopyAlteredGeneString(testedGenes);
				TreeSet<String> keyGeneCalls = mergeGeneStrings(keyGenes);
				IO.p("\tKey: "+keyGenes);
				
				for (CnvCallSet cc: ks.getCnvCallSets()) {
					if (cc.getMethodName().equals("Key")) continue;
					IO.p("\t"+cc.getMethodName()+": "+cc.getCopyAlteredGeneString(testedGenes));
					ArrayList<String> genesDirection = cc.getCopyAlteredGenes(testedGenes);
					
					//create confusion matrix
					HashSet<String> testGeneCalls = new HashSet<String>();
					for (String s: genesDirection) testGeneCalls.add(s);
					ConfusionMatrix cm = new ConfusionMatrix(testedGenes, keyGeneCalls, testGeneCalls);
					
					IO.p("\t"+ cm.toString());
					
				}
				IO.pl();
			}	
		}
	}


	public String formatFraction(int a, int b) {
		double f = (double)a / (double)b;
		return a+"/"+b+"("+Num.formatNumber(f, 4)+")";
	}



	public TreeSet<String> mergeGeneStrings(String one) {
		TreeSet<String> merge = new TreeSet<String>();
		String[] oneSplit = Misc.COMMA.split(one);
		for (String x: oneSplit) merge.add(x);
		if (merge.size()>1 && merge.contains("None")) merge.remove("None");
		return merge;
	}

	private void loadGatkCopyRatioCalls(File gatkResultsDir, String name) {
		IO.pl("\nLoading GATK "+name+" results...");
		IO.pl("\tDataset\tNumBedLines\tNumInterrogatedGenesCalled");
		File[] beds = IO.extractFiles(gatkResultsDir, ".bed.gz");
		if (beds.length==0) beds = IO.extractFiles(gatkResultsDir, ".bed");
		Pattern genePattern = Pattern.compile(";genes=");
		//for each results bed file results set
		for (File b: beds) {
			//1110039_22712G29_23802G10_10_Hg38.called.seg.pass.bed
			//   0       1         2     3   4
			String[] patient_cfDNA_germline = Misc.UNDERSCORE.split(b.getName().replaceAll("G", "X"));
			
			//pull KohliSample
			KohliSample ks = samples.get(patient_cfDNA_germline[1]);
			if (ks == null) Misc.printErrAndExit("\nFailed to find the sample associated with "+b);
			
			//load all of the bed results
			Bed[] bedResults = Bed.parseFile(b, 0, 0);
			int numCalls = bedResults.length;
			TreeMap<String, GeneCallResult> affectedGeneStats = new TreeMap<String, GeneCallResult>();
			
			//for each bed line
			for (Bed bed: bedResults) {
				//name    numOb=28;lg2Tum=-0.191;lg2Norm=0;genes=RALGAPA1P1,SNORD121B,LOC101928775
				String[] f = genePattern.split(bed.getName());
				if (f.length!=2 || f[0].startsWith("numOb=")==false) Misc.printErrAndExit("Failed to split on genes the bed result: "+bed.toString());
				String[] genesAffected = Misc.COMMA.split(f[1]);
				boolean containsMinus = f[0].contains("lg2Tum=-");
				HashSet<String> gas = Misc.loadHashSet(genesAffected);
				//for each interrogated gene
				for (String intGene: testedGenes) {
					if (gas.contains(intGene)) {
						affectedGeneStats.put(intGene, new GeneCallResult(true, containsMinus==false, f[0]));
					}
				}
				//does it intersect the AR Enhancer?
				if (bed.intersects(arEnhancer)) {
					affectedGeneStats.put("AR-Enh", new GeneCallResult(true, containsMinus==false, f[0]));
				}
			}
			IO.pl("\t"+b.getName()+"\t"+numCalls+"\t"+affectedGeneStats.size()+"\t"+affectedGeneStats.keySet());
			ks.getCnvCallSets().add(new CnvCallSet(name, ks.getSampleId(), numCalls, affectedGeneStats));
		}
	}
	
	private void loadPatients(File anno) {
		IO.pl("Loading patient and sample meta data...");
		String[] lines = IO.loadFileIntoStringArray(anno);
		if (lines[0].equals("#HCIPersonID\tSampleID\tType\tDateDrawn") == false) Misc.printErrAndExit("First line in the anno file isn't '#HCIPersonID SampleID Type	DateDrawn'? -> "+lines[0]);
		for (int i=1; i< lines.length; i++) {
			String[] tokens = Misc.TAB.split(lines[i]);
			KohliPatient kp = patients.get(tokens[0]);
			if (kp == null) {
				kp = new KohliPatient(tokens);
				patients.put(tokens[0], kp);
			}
			else {
				boolean isGermline = false;
				if (tokens[2].toLowerCase().contains("germline")) isGermline = true;
				kp.getSamples().add(new KohliSample(kp, tokens[1], isGermline, tokens[3]));
			}
		}
		//load the sample treemap
		for (KohliPatient kp: patients.values()) {
			for (KohliSample ks: kp.getSamples()) {
				if (samples.containsKey(ks.getSampleId())) Misc.printErrAndExit("Duplicate sample found! "+ks.getSampleId());
				samples.put(ks.getSampleId(), ks);
			}
		}
	}

}
