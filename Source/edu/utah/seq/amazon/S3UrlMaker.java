package edu.utah.seq.amazon;

import java.io.File;
import java.net.URL;
import java.util.ArrayList;
import java.util.Date;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import com.amazonaws.HttpMethod;
import com.amazonaws.auth.profile.ProfileCredentialsProvider;
import com.amazonaws.services.s3.*;
import com.amazonaws.services.s3.model.GeneratePresignedUrlRequest;
import com.amazonaws.services.s3.model.ListObjectsRequest;
import com.amazonaws.services.s3.model.ObjectListing;
import com.amazonaws.services.s3.model.S3ObjectSummary;
import util.gen.IO;
import util.gen.Misc;
import util.gen.Num;

/**
 * Generates signed timed urls for amazon s3 objects.
 * e.g. https://aruplab.s3.amazonaws.com/results/dev/2016/12-123-123456/block_norm_calls.vcf.gz?AWSAccessKeyId=AKIAIU6AVRF7G3Z2UXRA&Expires=1460134006&Signature=jb6Mk2y2xBYXl6K1MO8SZBxLoG0%3D
 * @author DavidNix*/
public class S3UrlMaker {

	//user defined fields
	private String bucketName;
	private String amazonPath;
	private File amazonProfileCredentials = null;
	private String amazonUser = null; //david.nix@aruplab.com
	private long hoursUntilExpiration = 72l;
	private boolean verbose = true;

	//internal
	private AmazonS3 s3client;
	private ArrayList<String> keysToFetch = new ArrayList<String>();
	private String[] signedUrls;
	private Date expiration;

	public S3UrlMaker(String[] args) {

		//parse and check args
		processArgs(args);

		//print fields
		if (verbose) printArgs();

		//set up log4j to dump warn and error messages to the console
		BasicConfigurator.configure();
		Logger.getRootLogger().setLevel(Level.WARN);

		//create a credentialed s3 client worker
		ProfileCredentialsProvider pcp = new ProfileCredentialsProvider(amazonProfileCredentials.toString(), amazonUser);
		s3client = new AmazonS3Client(pcp); 

		//lookup objects to make urls 
		lookupObjects();

		if (keysToFetch.size() ==0) System.err.println("\nError: no S3 objects found?! Aborting.\n");
		else {
			//when expire
			setExpirationDate();

			//make URLS
			makeUrls();
			
			for (String u: signedUrls) System.out.println(u);
		}
	}

	private void makeUrls() {
		if (verbose) System.out.println("\nSigned Urls:");
		signedUrls = new String[keysToFetch.size()];
		for (int i=0; i< signedUrls.length; i++){
			GeneratePresignedUrlRequest gpur = new GeneratePresignedUrlRequest(bucketName, keysToFetch.get(i));
			gpur.setMethod(HttpMethod.GET);
			gpur.setExpiration(expiration);
			URL s = s3client.generatePresignedUrl(gpur); 
			signedUrls[i] = s.toString();
		}
	}

	private void setExpirationDate() {
		expiration = new Date();
		long msec = expiration.getTime();
		msec += 1000 * 60 * 60 * hoursUntilExpiration;
		expiration.setTime(msec);
	}

	private void lookupObjects() {
		if (verbose) System.out.println("\nLooking up S3 Objects:");
		ListObjectsRequest listObjectsRequest = new ListObjectsRequest().withBucketName(bucketName).withPrefix(amazonPath);
		ObjectListing objectListing;

		do {
			objectListing = s3client.listObjects(listObjectsRequest);
			for (S3ObjectSummary objectSummary : objectListing.getObjectSummaries()) {
				if (objectSummary.getSize() !=0){
					if (verbose) System.out.println(objectSummary.getKey() + "\t" + Num.formatBytesToHumanReadableSize (objectSummary.getSize()));
					keysToFetch.add(objectSummary.getKey());
				}
			}
			listObjectsRequest.setMarker(objectListing.getNextMarker());
		} while (objectListing.isTruncated());
	}



	private void printArgs() {
		StringBuilder sb = new StringBuilder();
		sb.append("Params:");
		sb.append("\nBucket\t"); sb.append(bucketName);
		sb.append("\nPathToFetchFiles\t"); sb.append(amazonPath);
		sb.append("\nUser\t"); sb.append(amazonUser);
		sb.append("\nCredentials\t"); sb.append(amazonProfileCredentials.toString());
		sb.append("\nHrs2Expiration\t"); sb.append(hoursUntilExpiration);
		System.out.println(sb);
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new S3UrlMaker(args);
	}		

	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'p': amazonPath = args[++i]; break;
					case 'c': amazonProfileCredentials = new File(args[++i]); break;
					case 'u': amazonUser = args[++i]; break;
					case 'b': bucketName = args[++i]; break;
					case 't': hoursUntilExpiration = Long.parseLong(args[++i]); break;
					case 's': verbose = false; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}

		if (verbose) System.out.println(IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");

		//check path, shouldn't be set, keep hidden
		if (amazonPath == null) Misc.printErrAndExit("\nError: please provide a path within your bucket to the files you which to make signed URLs. "
				+ "See -p. e.g. results/prod/2016/15329403054-670197\n");
		amazonPath = Misc.removeLeadingTrailingForwardSlashes(amazonPath);

		//check bucket
		if (bucketName == null) Misc.printErrAndExit("\nError: please provide the name of an existing bucket to containing your path. See -b\n");
		bucketName = Misc.removeLeadingTrailingForwardSlashes(bucketName);

		//check credentials
		if (amazonProfileCredentials == null || amazonProfileCredentials.canRead() == false) Misc.printErrAndExit("\nError: can't find or read your amazon credentials. See -c "+ amazonProfileCredentials+ "\n");

		//check user name
		if (amazonUser == null || amazonUser.length()==0) Misc.printErrAndExit("\nError: please provide the name of the user specified in the amazon profile credentials. See -u "+ amazonUser+ "\n");

	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                               S3UrlMaker : April 2016                            **\n" +
				"**************************************************************************************\n" +
				"Generates Amazon S3 signed and timed URLs.\n" +

				"\nRequired:\n"+
				"-b Amazon bucket name.\n"+ 
				"-p Relative 'path' to a 'directory' or 'file' within the bucket to fetch S3 objects,\n"+
				"     recursive.\n"+
				"-c Path to an Amazon profile credentials text file containing:\n"+
				"     [USER_NAME]\n"+
				"     aws_access_key_id=YOUR_ACCESS_KEY_ID\n"+
				"     aws_secret_access_key=YOUR_SECRET_ACCESS_KEY\n"+
				"-u Amazon USER_NAME in credentials file.\n\n"+

				"Optional:\n" +
				"-t Hours until URLs expire, defaults to 72.\n"+
				"-s Silence non error messages.\n"+

				"\nExample: java -Xmx4G -jar pathTo/USeq/Apps/S3UrlMaker -b aruplab -t 24 \n"+
				"     -u gip@gmail.com -p results/dev/2016/12-123-123456 -c ~/Amazon/cred.txt \n\n"+

				"**************************************************************************************\n");

	}

	public String getBucketName() {
		return bucketName;
	}

	public String getAmazonPath() {
		return amazonPath;
	}

	public long getHoursUntilExpiration() {
		return hoursUntilExpiration;
	}

	public ArrayList<String> getKeysToFetch() {
		return keysToFetch;
	}

	public String[] getSignedUrls() {
		return signedUrls;
	}

}
