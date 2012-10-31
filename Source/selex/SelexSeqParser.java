package selex;
import java.util.*;

/**
 * Parses out oligos from multiple rounds of SELEX.
 * */
public class SelexSeqParser {
    
    public static void main(String[] args) {
        //make SelexParams object
        SelexParams sp = new SelexParams();
            sp.processArgs(args);
        
        //print header info and save to SelexParams running StringBuffer field
        sp.printSave("************************* New Selex Sequence Parser Run *****************************\n"+
        "Date: " + new Date()+"\n\n");
        
        //Create a new sequence file
        int len = (sp.getSeqFiles()).length;
        for (int i=0; i<len; i++){
            System.out.println("");
            SeqFiles sf = new SeqFiles(sp.getSeqFile(i), sp.getQualFile(i), sp);
        }
        //tell SelexParam object how many files are being processed
        sp.incNumFiles(len);
        
        //print final report and histogram
        sp.printSave(sp.makeFinalReport());
        sp.printSave(sp.plotHistogram());
        
        //write files
            //for each file write files using sp.names[] stuff
        sp.writeSubSeqFasta();
        sp.writeReport();
        
        
        System.out.println("\n\nDone!");
    }
    
}

