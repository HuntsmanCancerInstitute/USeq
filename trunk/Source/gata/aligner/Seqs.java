package gata.aligner;

import java.util.*;

import util.bio.parsers.*;

/**
 * @author Nix
 * Object representing the two user input sequences with methods to fire the multiFastaParsers.
 *  * */
public class Seqs {
    private String refSeq;
    private String compSeq;
    AlignParams AP;
    
    public Seqs(AlignParams ap, MultiFastaParser ref, MultiFastaParser comp, int index){
        AP = ap;
        refSeq = ref.getSeqs()[0];
        compSeq = comp.getSeqs()[index];
        String refSeqName = ref.getNames()[0];
        String compSeqName = comp.getNames()[index];
        
        //set parameters in AP
        AP.setLENGTH_REFSEQ(refSeq.length());
        AP.setLENGTH_COMPSEQ(compSeq.length());
        AP.setNAME_REFSEQ(refSeqName);
        AP.setNAME_COMPSEQ(compSeqName);
    }
    public LocalAlignment[] makeLocalAlignments(ArrayList data){
        // takes cycles of the following and makes local alignments: 
        //      0=score, 1=start, 2=stop, 3=seq, 4=start, 5=stop, 6=seq, 7=orientation (0 for +/+, 1 for +/-);
        int len = data.size();
        LocalAlignment[] locs = new LocalAlignment[len/8];
        int counter =0;
        int realStartRefSeq = AP.getSTART_INDEX_REFSEQ()-1;
        int realStartCompSeq = AP.getSTART_INDEX_COMPSEQ()-1;
        for (int i=0;i<len;i+=8){
            //constructor for loc align (String seq_1, int startSeq_1, int stopSeq_1, String seq_2, int startSeq_2,int stopSeq_2, int ori_, AlignParams ap_, int score_)
           locs[counter] = new LocalAlignment(
            (String)data.get(i+3), ((Integer)data.get(i+1)).intValue()+realStartRefSeq, ((Integer)data.get(i+2)).intValue()+realStartRefSeq,
            (String)data.get(i+6), ((Integer)data.get(i+4)).intValue()+realStartCompSeq, ((Integer)data.get(i+5)).intValue()+realStartCompSeq,
            ((Integer)data.get(i+7)).intValue(), AP, ((Integer)data.get(i)).intValue());
           
            counter++;
        }
        return locs;
    }

}