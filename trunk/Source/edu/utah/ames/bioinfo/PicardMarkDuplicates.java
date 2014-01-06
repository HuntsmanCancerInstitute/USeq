package edu.utah.ames.bioinfo;

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
 * Modification of their MarkDuplicates class to accommodate the USeq NovoalignBisulfiteParser app.
 * 
 * @author darren.ames@hci.utah.edu
 * 
 */

import net.sf.picard.sam.*;
import net.sf.picard.cmdline.*;
import net.sf.picard.io.IoUtil;
import net.sf.samtools.*;
import java.io.File;
import java.io.IOException;
import java.util.Iterator;


public class PicardMarkDuplicates extends CommandLineProgram {

	public File INPUT;
	public File OUTPUT;
	public File METRICS;

	protected int doWork() {
		
		SAMFileReader reader = new SAMFileReader(IoUtil.openFileForReading(INPUT));
		SAMFileWriter writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(reader.getFileHeader(), false, OUTPUT);
		
		//iterate over the records
		Iterator<SAMRecord> iterator = reader.iterator();
		while (iterator.hasNext()) writer.addAlignment(iterator.next());

		reader.close();
		writer.close();

		return 0;
	}

	/** Launches MarkDuplicates suppressing all output except warnings and errors
	 * Temp files are written to the parent of the outputBamFile
	 */
	 public PicardMarkDuplicates(File INPUT, File OUTPUT, File METRICS){
		 File realOutputFile;
		 try {
			 realOutputFile = OUTPUT.getCanonicalFile();
			 String[] argv = {
					 "I="+INPUT,
					 "O="+OUTPUT,
					 "M="+METRICS,
					 "REMOVE_DUPLICATES=true",
					 "TMP_DIR="+realOutputFile.getParent(),
					 "QUIET=true",
					 "VERBOSITY=WARNING",
					 "VALIDATION_STRINGENCY=SILENT"
			 };

			 new MarkDuplicates().instanceMain(argv);

		 } catch (IOException e) {
			 e.printStackTrace();
		 }
	 }

}
