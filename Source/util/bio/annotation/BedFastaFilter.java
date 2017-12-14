package util.bio.annotation;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import util.gen.Gzipper;
import util.gen.IO;
import util.gen.Misc;

public class BedFastaFilter {


	public static void main(String[] args) throws IOException {
		if (args.length == 0) Misc.printErrAndExit("\nPlease provide the file path to a bed file or directory of such containing fasta sequence in the name field.\n");

		File[] files = IO.extractFiles(new File(args[0]));

		//Pattern pat = Pattern.compile("[AG]AT{3,6}G", Pattern.CASE_INSENSITIVE);

		Pattern[] pats = new Pattern[]{
				Pattern.compile("[AG]AT{3,6}G", Pattern.CASE_INSENSITIVE),
				Pattern.compile("AT{3,6}G", Pattern.CASE_INSENSITIVE),
				Pattern.compile("[AG]AT{3,6}", Pattern.CASE_INSENSITIVE),
				Pattern.compile("[AG].T{3,6}G", Pattern.CASE_INSENSITIVE),

				Pattern.compile("[AG]A.TTG", Pattern.CASE_INSENSITIVE),
				Pattern.compile("[AG]AT.TG", Pattern.CASE_INSENSITIVE),
				Pattern.compile("[AG]ATT.G", Pattern.CASE_INSENSITIVE),

				Pattern.compile("[AG]A.TTTG", Pattern.CASE_INSENSITIVE),
				Pattern.compile("[AG]AT.TTG", Pattern.CASE_INSENSITIVE),
				Pattern.compile("[AG]ATT.TG", Pattern.CASE_INSENSITIVE),
				Pattern.compile("[AG]ATTT.G", Pattern.CASE_INSENSITIVE),

				Pattern.compile("[AG]A.TTTTG", Pattern.CASE_INSENSITIVE),
				Pattern.compile("[AG]AT.TTTG", Pattern.CASE_INSENSITIVE),
				Pattern.compile("[AG]ATT.TTG", Pattern.CASE_INSENSITIVE),
				Pattern.compile("[AG]ATTT.TG", Pattern.CASE_INSENSITIVE),
				Pattern.compile("[AG]ATTTT.G", Pattern.CASE_INSENSITIVE),

				Pattern.compile("[AG]A.TTTTTG", Pattern.CASE_INSENSITIVE),
				Pattern.compile("[AG]AT.TTTTG", Pattern.CASE_INSENSITIVE),
				Pattern.compile("[AG]ATT.TTTG", Pattern.CASE_INSENSITIVE),
				Pattern.compile("[AG]ATTT.TTG", Pattern.CASE_INSENSITIVE),
				Pattern.compile("[AG]ATTTT.TG", Pattern.CASE_INSENSITIVE),
				Pattern.compile("[AG]ATTTTT.G", Pattern.CASE_INSENSITIVE)
		};

		for (File f: files){
			BufferedReader in = IO.fetchBufferedReader(f);
			System.out.println(f.getName());
			Gzipper out = new Gzipper(new File(f.getParentFile(), Misc.removeExtension(f.getName())+"filtered.bed.gz"));
			String line = null;
			String[] tokens = null;
			ArrayList<String> matches = new ArrayList<String>();
			while ((line= in.readLine()) != null){
				tokens = Misc.TAB.split(line);
				//enough reads?
				double numReads = Double.parseDouble(tokens[4]);
				if (numReads < 10.0) continue;
				matches.clear();
				//loop through each pattern
				for (Pattern pat: pats){
					Matcher mat = pat.matcher(tokens[3]);
					if (mat.find()) {
						matches.add(mat.group());
					}
				}
				if (matches.size() != 0){
					 out.println(line+"\t"+ Misc.stringArrayListToString(matches, ","));
				}
			}
			out.close();
			in.close();
		}

	}

}
