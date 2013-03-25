package edu.utah.ames.bioinfo;

import java.io.File;

/**
 * Recursively finds regex-defined novoindices in user directories
 * 
 * @author darren.ames@hci.utah.edu
 *
 */
public class FileFindRegex {

	private File fileObject;
	private Sample sample;

	public static void main(String[] args) {
		String path = "/Users/darren/Desktop/novoalignerTestDir/";
		Sample sample = new Sample(args);
		FileFindRegex ffr = new FileFindRegex(new File(path), sample);
		ffr.find();
	}

	//constructor
	public FileFindRegex(File fileObject, Sample sample) {
		this.fileObject = fileObject;
		this.sample = sample;
	}

	public void find() {
		recursiveFind(fileObject);
	}

	//recursively search dir
	public void recursiveFind(File fileObject) {
		//System.out.println(fileObject + " RT");
		if (fileObject.isDirectory()) {
			//System.out.println(indent + fileObject.getName());
			File allFiles[] = fileObject.listFiles();
			for (File aFile : allFiles) {
				recursiveFind(aFile);
			}
		} 
		//TODO set pattern matching stuff here 
		if (fileObject.getName().equalsIgnoreCase("zv9EnsTransRad46Num100kMin10SplicesChrPhiXAdaptr.nov.illumina.nix")) {
			//sample.getNovoindex().add(fileObject);
			System.out.println(fileObject);
		}
	}
}
