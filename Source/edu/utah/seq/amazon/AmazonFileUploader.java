package edu.utah.seq.amazon;

import java.io.File;
import java.io.IOException;
import java.util.UUID;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import com.amazonaws.AmazonClientException;
import com.amazonaws.AmazonServiceException;
import com.amazonaws.auth.profile.ProfileCredentialsProvider;
import com.amazonaws.services.s3.transfer.TransferManager;
import com.amazonaws.services.s3.transfer.Upload;

import edu.utah.seq.vcf.VCFMerger;
import util.gen.IO;
import util.gen.Misc;

public class AmazonFileUploader {
	
	
	//hmm this really isn't needed, use jetS3
	
	
	//user defined fields
	private String amazonPathToResultsDir = "results";
	private Integer year = 0; //2016, 2017, ....
	private String serverType = null; //dev, cert, prod
	private String uuid = null; //optional, will be generated if null
	private String amazonBucket = "aruplab";
	private File directory = null; //15329403054-670197
	private File amazonProfileCredentials = null;
	private String amazonUser = null; //david.nix@aruplab.com
	
	//internal fields
	private File[] filesToUpload = null;
	private String amazonPath = null;
	
	
	/**Generates the path sans bucket name,
	 * amazonPathToResultsDir/serverType/year/directorySampleName/uuid/
	 * e.g. results/dev/2016/15329403054-670197/067e6162-3b6f-4ae2-a171-2470b63dff00/  */
	private void buildAmazonPath() {
		StringBuilder sb = new StringBuilder();
		sb.append(amazonPathToResultsDir); sb.append("/");
		sb.append(serverType); sb.append("/");
		sb.append(year.toString()); sb.append("/");
		sb.append(directory.getName()); sb.append("/");
		sb.append(uuid); sb.append("/");
		amazonPath = sb.toString();
	}
	
	private void printArgs(){
		//System.out.println(x);
	}
	
	public AmazonFileUploader (String[] args){
		try {
			
		processArgs(args);
		
		buildAmazonPath();

			
			//String keyName            = "NixTests/block_norm_calls.vcf.gz";
			String keyName            = "results/dev/2016/12-123-123456/block_norm_calls.vcf.gz";
			String filePath           = "/Users/DavidNix/Annotation/Hg19/Nist2.19/block_norm_calls.vcf.gz";  
 
			/* Create a ProfileCredentialsProvider pointing to a txt file with the following
			[userName]
			aws_access_key_id={YOUR_ACCESS_KEY_ID}
			aws_secret_access_key={YOUR_SECRET_ACCESS_KEY}
			* be sure to chmod 600 this file! */
			ProfileCredentialsProvider pcp = new ProfileCredentialsProvider(amazonProfileCredentials.getCanonicalPath(), amazonUser);
			
			//Use a TransferManager to chunk large files and upload each over many streams
			//From the ClientConfiguration for the TransferManager : "The default configuration is to use HTTPS for all requests for increased security."
			System.out.println("Initializing...");
			TransferManager tm = new TransferManager(pcp);
			
			for (int i=0; i< filesToUpload.length; i++){
				//Create an Uploader for each file
				Upload upload = tm.upload(amazonBucket, keyName, new File(filePath));

				//Start transfer
				System.out.println("Starting transfer...");
				upload.waitForCompletion();
				System.out.println("Upload complete! ");
			}
			
			tm.shutdownNow();

		} catch (AmazonServiceException ase) {
			System.err.println("\nUnable to upload file! ASE");
			System.out.println("Error Message:    " + ase.getMessage());
			System.out.println("HTTP Status Code: " + ase.getStatusCode());
			System.out.println("AWS Error Code:   " + ase.getErrorCode());
			System.out.println("Error Type:       " + ase.getErrorType());
			System.out.println("Request ID:       " + ase.getRequestId());
			ase.printStackTrace();
		} catch (Exception e) {
			e.printStackTrace();
		} 
	}


	




	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new AmazonFileUploader(args);
	}		

	/**	private String resultsDir = "results";
	private int year = 0; //2016, 2017, ....
	private String serverType = null; //dev, cert, prod
	private String resultsName = null; //15329403054-670197
	private String uuid = null; //optional, will be generated if null
	private String existingBucketName = "aruplab";
	private File[] filesToUpload = null;/

	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'p': amazonPathToResultsDir = args[++i]; break;
					case 'c': amazonProfileCredentials = new File(args[++i]); break;
					case 'b': amazonBucket = args[++i]; break;
					case 'd': directory = new File (args[++i]); break;
					case 'y': year = Integer.parseInt(args[++i]); break;
					case 's': serverType = args[++i]; break;
					case 'u': uuid = args[++i]; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		
		//pull files
		if (directory == null || directory.canRead() == false) Misc.printErrAndExit("\nError: please provide a file or directory of files to upload with -d.\n");
		if (directory.isFile()) Misc.printErrAndExit("\nError: please provide a directory containing files to upload, see -d.\n");
		filesToUpload = IO.fetchFilesRecursively(directory);
		if (filesToUpload == null || filesToUpload.length == 0) Misc.printErrAndExit("\nError: no files were found to upload? See -d and -> "+ directory+ "\n");
		
		//check path, shouldn't be set, keep hidden
		if (amazonPathToResultsDir == null) Misc.printErrAndExit("\nError: please provide a path within your bucket to save the results. See -p. e.g. results\n");
		amazonPathToResultsDir = Misc.removeLeadingTrailingForwardSlashes(amazonPathToResultsDir);
		
		//check bucket
		if (amazonBucket == null) Misc.printErrAndExit("\nError: please provide the name of an existing bucket to save the files. See -b\n");
		amazonBucket = Misc.removeLeadingTrailingForwardSlashes(amazonBucket);
		
		//check year
		if (year < 2010) Misc.printErrAndExit("\nError: please provide the year of this result (e.g. 2016, 2017, etc). See -y\n");
		
		//check serverType
		if (serverType == null) Misc.printErrAndExit("\nError: please provide a server type, either dev, cert, or prod. See -s\n");
		serverType = serverType.toLowerCase();
		if (serverType.equals("dev") == false && serverType.equals("cert") == false && serverType.equals("prod") == false) Misc.printErrAndExit("\nError: please provide a server type, either dev, cert, or prod. See -s\n");
		
		//create uuid?
		if (uuid == null || uuid.length()==0) uuid = UUID.randomUUID().toString();
		
	}	
	

	
	

	
	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                              VCF Merger : Feb 2016                               **\n" +
				"**************************************************************************************\n" +
				"Merges VCF files with the same samples. Collapses the headers with a simple hash. Will\n"+
				"not work well with downstream apps that cannot process mixed INFO and FORMAT records.\n" +

				"\nRequired:\n"+
				"-v Full path to a vcf file (xxx.vcf(.gz/.zip OK)) or directory containing such. Note,\n"+
				"       Java often fails to parse tabix compressed vcf files.  Best to uncompress.\n\n"+
								
				"Optional:\n" +
				"-o Full path to an output vcf file, defaults to merged.vcf.gz in parent -v dir.\n" +

				"\n"+

				"Example: java -Xmx4G -jar pathTo/USeq/Apps/VCFMerger -v /CancerSamples/\n\n"+

		"**************************************************************************************\n");

	}
}
