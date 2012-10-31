package meme;

import java.io.*;
import java.util.*;
import java.util.regex.*;

    /** A catch all class to hold static methods and parameters set by different classes
     *associated by MemeR.*/
public class UtilMeme {
 

    
    public static String makeFullPathName(String resultsDirectory, String fileName){
        File dump = new File(resultsDirectory, fileName);
        String fullPathName = "";
        try{
            fullPathName = dump.getCanonicalPath();
        }
        catch (IOException e){
            System.err.println("\n\nProblem building text to write files!\n\n");
            e.printStackTrace();
            System.exit(1);
        }
        return fullPathName;
    } 
    public static String extractFileName(String data){
        Pattern pat = Pattern.compile ("meme\\s(.+?)\\s");
        Matcher mat = pat.matcher(data);
        if (mat.lookingAt()) {
            File file = new File(mat.group(1));
            return file.getName();
        }
        else return "problem with parsing file text from commandline";
    }
    public static ArrayList createStringAL (String[] data){
        int len = data.length;
        ArrayList x = new ArrayList(len);
        for (int i=0; i<len; i++)x.add(data[i]);
        return x;
    }
    public static boolean checkFile(String par){
        File f = new File(par);
        if (f.exists()) return true;
        return false;
    }
    public static void printDocs(){
        System.out.println(UtilMeme.getDocs());
    }
    public static String getDocs(){
        return "\n\n"+
        "********************************************************************\n"+
        "MemeR: a Java Wrapper for MEME, version 0.2, Nix@uclink.berkeley.edu\n"+
        "********************************************************************\n"+
        "To use MemeR you will need to install MEME and Java 1.4+\n"+
        "   ftp://ftp.sdsc.edu/pub/sdsc/biology/meme/\n"+
        "   type 'java -version' to obtain your local java stats\n\n"+
        "You will also need one or more multi-FASTA files containing DNA seqs.\n"+
        "   MemeR runs MEME on each FASTA file to identify sequence motifs.\n"+
        "   These motifs are then used to search for matches in each FASTA\n"+
        "   file.\n\n"+
        "MemeR returns lots of information:\n"+
        "   Parsed MEME reports for each motif...\n"+
        "   Log Likelihood Position Specific Probability Matrices for both\n"+
        "       MEME's PSPM and an Add One Pseudo Count PSPM...\n"+
        "   Motif matches to sequences in each FASTA file...\n"+
        "   Scores for the number of sequences found to contain the\n"+
        "       MEME and Add One defined LLPSPM motifs.  (Use this as an \n"+
        "       estimation of SELEX oligo enrichment. Add One motifs are a\n"+
        "       little more relaxed.)\n"+
        "   Two summary tables detailing the number of motifs found and the\n"+
        "       number of sequences with the motif in each file.\n\n"+
        "Three result files are written to disk:\n"+
        "   yourResults.screenDump = the terminal screen contents\n"+
        "   yourResults.motifs = motif reports for each motif\n"+
        "   yourResults.sumResults = summary tables for the motifs, motif\n"+
        "       hit statistics, and a graph plotting the percentage of seqs\n"+
        "       in each file that contain each motif\n\n"+
        "TO RUN MEMER type 'java MemeR run' on the command line and fill out\n"+
        "   the five dialog boxes.\n\n";  
    }
}
