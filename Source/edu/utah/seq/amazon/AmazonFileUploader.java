package edu.utah.seq.amazon;

import java.io.File;
import java.io.IOException;

import com.amazonaws.AmazonClientException;
import com.amazonaws.AmazonServiceException;
import com.amazonaws.auth.profile.ProfileCredentialsProvider;
import com.amazonaws.services.s3.transfer.TransferManager;
import com.amazonaws.services.s3.transfer.Upload;

public class AmazonFileUploader {

	public static void main(String[] args) throws IOException {
		try {
			String existingBucketName = "pipeline2-test";
			String keyName            = "NixTests/block_norm_calls.vcf.gz";
			String filePath           = "/Users/DavidNix/Annotation/Hg19/Nist2.19/block_norm_calls.vcf.gz";  
 
			//Create a ProfileCredentialsProvider pointing to a txt file with the following
			//[userName]
			//aws_access_key_id={YOUR_ACCESS_KEY_ID}
			//aws_secret_access_key={YOUR_SECRET_ACCESS_KEY}
			//be sure to chmod 600 this file!
			ProfileCredentialsProvider pcp = new ProfileCredentialsProvider("/Users/DavidNix/Arup/Amazon/KeysCredentials/credentials.txt", "david.nix@aruplab.com");
			
			//Use a TransferManager to chunk large files and upload each over many streams
			//From the ClientConfiguration for the TransferManager : "The default configuration is to use HTTPS for all requests for increased security."
			System.out.println("Initializing...");
			TransferManager tm = new TransferManager(pcp);
			
			
			//Create an Uploader for each file
			Upload upload = tm.upload(existingBucketName, keyName, new File(filePath));

			//Start transfer
			System.out.println("Starting transfer...");
			upload.waitForCompletion();
			System.out.println("Upload complete! ");

			tm.shutdownNow();

			//for debugging
		} catch (AmazonServiceException ase) {
			System.err.println("\nUnable to upload file! ASE");
			System.out.println("Error Message:    " + ase.getMessage());
			System.out.println("HTTP Status Code: " + ase.getStatusCode());
			System.out.println("AWS Error Code:   " + ase.getErrorCode());
			System.out.println("Error Type:       " + ase.getErrorType());
			System.out.println("Request ID:       " + ase.getRequestId());
			ase.printStackTrace();

		} catch (AmazonClientException ace) {
			ace.printStackTrace();

		} catch (InterruptedException ie) {
			ie.printStackTrace();
		}


	}
}
