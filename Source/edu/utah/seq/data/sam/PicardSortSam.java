package edu.utah.seq.data.sam;

/*
 * The MIT License
 *
 * Copyright (c) 2009 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 * 
 * @original author alecw@broadinstitute.org 
 * 
 * Modification of their SortSam.java class to accommodate the USeq RNASeq wrapper app.
 * 
 * @author davidnix
 * 
 */
import picard.sam.*;
import picard.cmdline.CommandLineProgram;
import util.gen.Misc;
import util.gen.NullPrintStream;
import htsjdk.samtools.*;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Iterator;


public class PicardSortSam extends CommandLineProgram {

	public File INPUT;
	public File OUTPUT;
	public SAMFileHeader.SortOrder SORT_ORDER;

	protected int doWork() {

		SamReader reader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(INPUT);
		reader.getFileHeader().setSortOrder(SORT_ORDER);
		SAMFileWriter writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(reader.getFileHeader(), false, OUTPUT);

		Iterator<SAMRecord> iterator = reader.iterator();
		while (iterator.hasNext()) writer.addAlignment(iterator.next());

		try { reader.close(); } catch (IOException e) { e.printStackTrace();}
		writer.close();

		return 0;
	}
	
	/**Makes a bam based on the name of the inputSamFile.*/
	public PicardSortSam(File inputSamFile){
		//make bam File
		String name = Misc.removeExtension(inputSamFile.getName());
		File bam = new File (inputSamFile.getParentFile(), name+".bam");
		new PicardSortSam (inputSamFile, bam, true);
	}

	/** Launches SortSam suppressing all output except warnings and errors
	 * Temp files are written to the parent of the outputBamFile
	 */
	 public PicardSortSam(File inputSamFile, File outputBamFile, boolean byCoordinate){
		 File realOutputFile;
		 try {
			 realOutputFile = outputBamFile.getCanonicalFile();
			 String[] argv = null;
			 if (byCoordinate) {
				 argv = new String[] { "I="+inputSamFile,
					 "O="+outputBamFile,
					 "SO=coordinate",
					 "CREATE_INDEX=true",
					 "TMP_DIR="+realOutputFile.getParent(),
					 "QUIET=true",
					 "USE_JDK_DEFLATER=true",
					 "USE_JDK_INFLATER=true",
					 "VERBOSITY=ERROR",
					 "VALIDATION_STRINGENCY=SILENT"};
			 }
			 else {
				 argv = new String[] { "I="+inputSamFile,
						 "O="+outputBamFile,
						 "SO=queryname",
						 "TMP_DIR="+realOutputFile.getParent(),
						 "QUIET=true",
						 "USE_JDK_DEFLATER=true",
						 "USE_JDK_INFLATER=true",
						 "VERBOSITY=ERROR",
						 "VALIDATION_STRINGENCY=SILENT"};
 
			 }
			 
			 new SortSam().instanceMain(argv);

		 } catch (IOException e) {
			 e.printStackTrace();
		 }
	 }
	 
		/** Launches SortSam suppressing all output except warnings and errors
		 * Temp files are written to the parent of the outputBamFile
		 
		 public PicardSortSam(File inputSamFile, File outputBamFile, boolean byCoordinate){
			 File realOutputFile;
			 try {
				 realOutputFile = outputBamFile.getCanonicalFile();
				 String[] argv = null;
				 if (byCoordinate) {
					 argv = new String[] { "-I", inputSamFile.getCanonicalPath(),
						 "-O", outputBamFile.getCanonicalPath(),
						 "-SO","coordinate",
						 "-CREATE_INDEX", "true",
						 "-TMP_DIR", realOutputFile.getParentFile().getCanonicalPath(),
						 "-QUIET", "true",
						 "-USE_JDK_DEFLATER","true",
						 "-USE_JDK_INFLATER","true",
						 "-VERBOSITY","ERROR",
						 "-VALIDATION_STRINGENCY","SILENT"};
				 }
				 else {
					 argv = new String[] { "-I", inputSamFile.getCanonicalPath(),
							 "-O", outputBamFile.getCanonicalPath(),
							 "-SO","queryName",
							 "-TMP_DIR", realOutputFile.getParentFile().getCanonicalPath(),
							 "-QUIET", "true",
							 "-USE_JDK_DEFLATER","true",
							 "-USE_JDK_INFLATER","true",
							 "-VERBOSITY","ERROR",
							 "-VALIDATION_STRINGENCY","SILENT"};
				 }

				 new SortSam().instanceMain(argv);

			 } catch (IOException e) {
				 e.printStackTrace();
			 }
		 }*/
		 
	 /*
	 public static void main(String[] args) {
		 if (args.length == 3) {
		 File sam = new File("/Users/u0028003/Downloads/BamPileup/snv0.05.bam");
		 File samSorted = new File("/Users/u0028003/Downloads/delme.snv0.05.bam");
		 new PicardSortSam(sam, samSorted);
		 }
	 }*/

}
