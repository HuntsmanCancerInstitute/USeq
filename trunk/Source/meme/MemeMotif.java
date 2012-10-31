package meme;


/**This contains info generated from parsing a meme report, each found motif is represented by this class.*/
public class MemeMotif {
    
    private int numberName; //text of motif, 1st found by meme, 2nd found by meme...
    private String motifSumLn;  //one line summary
    private String motifDesc;   //multiline visual descript of the motif
    private String motifSites;  //multiline visual descript of the sites found containing this motif
    private String[] hits;      //the actual sites parsed from motifSites
    private double[][] memePSPM;   //meme's position-specific probability matrix
    private double[][] add1PSPM;   //an add one PSPM calculated from the hits
    private int[][] freqMatrix; 	//a frequency Matrix generated from hits;
    private MemeParser mpRef;   //reference to the MemeParser that generated this MemeMotif
    
    //getter methods
    public int getNumberName(){return numberName;}
    public String getMotifSumLn(){return motifSumLn;}
    public String getMotifDesc(){return motifDesc;}
    public String getMotifSites(){return motifSites;}
    public String[] getHits(){return hits;}
    public double[][] getMemePSPM(){return memePSPM;}
    public double[][] getLLMemePSPM(double pA, double pC, double pG, double pT){
        //pN is background probability of a particular nucleotide (ie 0.25)
        return makeLLPSPM(memePSPM, pA, pC, pG, pT);
    }
    public double[][] getAdd1PSPM(){return add1PSPM;}
    public double[][] getLLAdd1PSPM(double pA, double pC, double pG, double pT){
        //pN is background probability of a particular nucleotide (ie 0.25)
        return makeLLPSPM(add1PSPM, pA, pC, pG, pT);
    }
    public MemeParser getMemeParser(){return mpRef;}    
    
    //setter methods
    public void setNumberName(int x){numberName = x;}
    public void setMotifSumLn(String x){motifSumLn =x;}
    public void setMotifDesc(String x){motifDesc =x;}
    public void setMotifSites(String x){motifSites =x;}
    public void setHits(String[] x){
        hits =x;
        add1PSPM = makeAddOnePSPM(hits);
        freqMatrix = makeFrequencyMatrix(hits);
    }
    public void setMemePSPM(double[][] x){memePSPM =x;}
    public void setMemeParserRef(MemeParser x){mpRef =x;}
    
    //primary methods
	/**takes a double[][] where array[0] refers to the first position in a motif
	 * array[0][0 to 3] = probability of A, C, G, or T nucleotide at that position.
	 * MakeLLPSPM converts each nucleotide probablility to a log likelihood value by
	 * taking the natural log of the probability ration ie array[0][0]/ pA  .
	 * pN is background probability of a particular nucleotide (ie 0.25)*/    
    public static double[][] makeLLPSPM(double[][] PSPM, double pA, double pC, double pG, double pT){
        int len = PSPM.length;
        double[][] LL = new double[len][4];
        for (int i=0; i<len; i++){
            LL[i][0] = Math.log(PSPM[i][0]/pA);
            LL[i][1] = Math.log(PSPM[i][1]/pC);
            LL[i][2] = Math.log(PSPM[i][2]/pG);
            LL[i][3] = Math.log(PSPM[i][3]/pT);
        }
        return LL;
    }
    
    public String toString(){
        StringBuffer dataDump = new StringBuffer(
        "Motif Summary: "+motifSumLn+"\n\n"+
        motifDesc+
        motifSites+"\n\n"+
        "MEME Position Specific Probability Matrix - ACGT\n");
        
        for (int x=0; x<memePSPM.length; x++){
            dataDump.append(memePSPM[x][0]+" \t"+memePSPM[x][1]+" \t"+memePSPM[x][2]+" \t"+memePSPM[x][3]+"\n");
        }
        
        dataDump.append("\n\nLog-Likelihood MEME Position Specific Probability Matrix- ACGT (pA,pC,pG,pT = 0.25)\n");
        double[][] LL = getLLMemePSPM(0.25,0.25,0.25,0.25);
        for (int x=0; x<LL.length; x++){
            dataDump.append(LL[x][0]+" \t"+LL[x][1]+" \t"+LL[x][2]+" \t"+LL[x][3]+"\n");
        }
        
        dataDump.append("\n\nAdd One Position Specific Probability Matrix- ACGT\n");
        
        for (int x=0; x<add1PSPM.length; x++){
            dataDump.append(add1PSPM[x][0]+" \t"+add1PSPM[x][1]+" \t"+add1PSPM[x][2]+" \t"+add1PSPM[x][3]+"\n");
        }
        
        dataDump.append("\n\nLog-Likelihood Add One Position Specific Probability Matrix- ACGT (pA,pC,pG,pT = 0.25)\n");
        LL = getLLAdd1PSPM(0.25,0.25,0.25,0.25);
        for (int x=0; x<LL.length; x++){
            dataDump.append(LL[x][0]+" \t"+LL[x][1]+" \t"+LL[x][2]+" \t"+LL[x][3]+"\n");
        }
        
		dataDump.append("\n\nBase Count Frequency Matrix (# ACGT's observed at each position in motif oligos)\n");
		
		for (int x=0; x<freqMatrix.length; x++){
			dataDump.append(freqMatrix[x][0]+" \t"+freqMatrix[x][1]+" \t"+freqMatrix[x][2]+" \t"+freqMatrix[x][3]+"\n");
		}        
        
        dataDump.append("\n************************************************************************************************\n\n");
        return dataDump.toString();
    }

	/**This method returns a position specific probability weight matrix given an array
	 * of Strings, each of the same size, containing GATC, no gaps, no N's.
	 * The weight matrix is an Add One matrix where a single psuedo count is added to
	 * all base counts when a zero is found in a particular base bin.
	 * The order of the array is ACGT for each position in a single hit.*/    
    public static double[][] makeAddOnePSPM(String[] hits){
        int lenOfAHit = hits[0].length();
        int len = hits.length;
        double[][] PSPM = new double[lenOfAHit][4];
        for (int i=0; i<lenOfAHit; i++){
            //set base counters
            double As =0;
            double Cs =0;
            double Gs =0;
            double Ts =0;
            //run through each hit looking at a particular base
            for (int j=0; j<len; j++){
                char test = hits[j].charAt(i);
                switch (test){
                    case 'A': As++; break;
                    case 'C': Cs++; break;
                    case 'G': Gs++; break;
                    case 'T': Ts++; break;
                    default: System.out.println("\n\nFatal error: odd base in makeAddOnePSPM hits -> "+test+"\n\n"); System.exit(0);
                }
            }
            if (As==0 || Cs==0 || Gs==0 || Ts==0){
                As++; Cs++; Gs++; Ts++;
            }
            double total = As+Cs+Gs+Ts;
            double[] Ps = {As/total, Cs/total, Gs/total, Ts/total};
            PSPM[i] = Ps;
        }
        return PSPM;
    }
	/**This method returns a position specific probability weight matrix given an array
	 * of Strings, each of the same size, containing GATC, no gaps, no N's.
	 * The order of the array is ACGT for each position in a single hit.*/    
	public static double[][] makePSPM(String[] hits){
		int lenOfAHit = hits[0].length();
		int len = hits.length;
		double[][] PSPM = new double[lenOfAHit][4];
		for (int i=0; i<lenOfAHit; i++){
			//set base counters
			double As =0;
			double Cs =0;
			double Gs =0;
			double Ts =0;
			//run through each hit looking at a particular base
			for (int j=0; j<len; j++){
				char test = hits[j].charAt(i);
				switch (test){
					case 'A': As++; break;
					case 'C': Cs++; break;
					case 'G': Gs++; break;
					case 'T': Ts++; break;
					default: System.out.println("\n\nFatal error: odd base in makeAddOnePSPM hits -> "+test+"\n\n"); System.exit(0);
				}
			}
			double total = As+Cs+Gs+Ts;
			double[] Ps = {As/total, Cs/total, Gs/total, Ts/total};
			PSPM[i] = Ps;
		}
		return PSPM;
	}  
	/**Generates a matrix of the number of As Cc Gg Ts (top) by 1,2,3,4...positions in the
	 * motif (side) observed
	 * in all the Strings of the String[].  The String[] should contain equal length Strings,
	 * comprised entirely of GATC.*/
	public static int[][] makeFrequencyMatrix(String[] hits){
		int lenOfAHit = hits[0].length();
		int len = hits.length;
		int[][] matrix = new int[lenOfAHit][4];
		//System.out.println("      \tAs\tCs\tGs\tTs");
		for (int i=0; i<lenOfAHit; i++){
			//set base counters
			int As =0;
			int Cs =0;
			int Gs =0;
			int Ts =0;
			//run through each hit looking at a particular base
			for (int j=0; j<len; j++){
				char test = hits[j].charAt(i);
				switch (test){
					case 'a':
					case 'A': As++; break;
					case 'c':
					case 'C': Cs++; break;
					case 'g':
					case 'G': Gs++; break;
					case 't':
					case 'T': Ts++; break;
					default: System.out.println("\n\nFatal error: odd base in makeFrequencyMatrix hits -> "+test+"\n\n"); System.exit(0);
				}
			}
			//System.out.println("Pos: "+(i+1)+"\t"+As+"\t"+Cs+"\t"+Gs+"\t"+Ts);
			int[] counts = {As, Cs, Gs, Ts};
			matrix[i] = counts;
		}
		return matrix;
	}
	public int[][] getFreqMatrix() {
		return freqMatrix;
	}

}
