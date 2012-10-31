package selex;
import java.util.*;
import java.text.*;

/**Utilites for the Selex package.*/
public class UtilSelex {
	
	/**Averages an ArrayList of Integer objects.*/
	public static int averageIntegerArrayList(ArrayList Integers){
		int len = Integers.size();
		if (len==0) return 0;
		int total = 0;
		for (int i=0; i<len; i++){
			int num = ((Integer)Integers.get(i)).intValue();
			total += num;
			}
		double ave = (double)total/(double)len;
		return (int)Math.round(ave);
		}
	/**Averages an ArrayList of Double objects.*/
		public static String averageDoubleArrayList(ArrayList Doubles){
			int len = Doubles.size();
			if (len==0) return "0";
			double total = 0;
			for (int i=0; i<len; i++){
				double num = ((Double)Doubles.get(i)).doubleValue();
				total += num;
				}
			double ave = total/(double)len;
			NumberFormat formatter = NumberFormat.getNumberInstance();
			formatter.setMinimumFractionDigits(2);
			return formatter.format(ave);
			}	
    
    public static void printDoc(){
        System.out.print( "\n***************************************************************************************\n"+
                            "*                       Selex Sequence Parser (0.2):                                  *\n"+
                            "***************************************************************************************\n\n"+
                         
                         "This program extracts SELEX 'oligos' from concatenated SELEX experiment sequences\n"+
                         "  given an identifying, intervening motif (i.e. a restriction site)\n\n"+
                         
                         "Two files are returned: a FASTA file containing all the SELEX 'oligos' that pass both\n"+
                         "  quality and length filters, as well as a results file with process information, \n"+
                         "  summary details, and an 'oligo' length distribution histogram.\n\n"+
                         
                         "To run this program you need Java 1.4+ and two files: a sequence FASTA file and a \n"+
                         "  FASTA-like quality file. Use the command line parameters to, at minimum, indicate\n"+
                         "  which pair of files should be parsed. For example, on unix based machines, use the\n"+
                         "  terminal to get into the directory containing the SelexSeqParser and type:\n\n"+
                         
                         "  java SelexSeqParser -s /my/seqfile/ishere/test1.seq -q /my/qualfile/ishere/test1.qual\n\n"+
                         
                         "Command line parameters:\n"+
                         "***************************************************************************************\n"+
                         "-s  required  Sequence file (ex: -s /my/seq/file/is/here/test1.seq )\n"+
                         "-q  required  Quality file (ex: -q /my/qual/file/is/here/test1.qual )\n"+
                         "-c  optional  Quality score cut off, defaults to 20. (ex: -q 18 )\n"+
                         "-i  optional  Minimal length cut off, defaults to 15 bases. (ex: -i 10 )\n"+
                         "-a  optional  Maximal length cut off, defaults to 17 bases. (ex: -a 18 )\n"+
						 "-e  optional  Expected parsed seq length (w/o the motif), defaults to 16 nt (ex: -e 18)\n"+
                         "-m  optional  Parsing motif/ restriction site, defaults to GGATCC. (ex: -m GGGCCC )\n"+
                         "-p  optional  Path to where you would like the two result files placed, defaults to the\n"+
                         "                sequence file directory. (ex: -p /data/goes/here/ )\n"+
                         "-n  optional  Prefix text for result files, defaults to the text of the first sequence\n"+
                         "                file. (ex: -n hunchback26Sept03 )\n"+
                         "-h  optional  Print documentation. (ex: -h )\n"+
                         "***************************************************************************************\n\n"+
                         
                         "Multiple files can be processed and results pooled.  Separate each with a comma and use \n"+
                         "  the same ordering. (ex: -s /data/t1.seq,/data/t2.seq -q /data/t1.qual,/data/t2.qual )\n\n"+
                         
                         "Q's, P's, and S's? nix@uclink.berkeley.edu or mbeisen@lbl.gov Version: 0.2 (1Dec03)\n\n");
    }
    
}

