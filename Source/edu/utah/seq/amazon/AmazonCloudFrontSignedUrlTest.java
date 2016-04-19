package edu.utah.seq.amazon;
import com.amazonaws.services.cloudfront.CloudFrontUrlSigner.Protocol;
import com.amazonaws.util.DateUtils;
import java.io.File;
import java.util.Date;
import com.amazonaws.services.cloudfront.*;

public class AmazonCloudFrontSignedUrlTest {

	public static void main(String[] args) {
		try {
			//Constants, modify the location of your private key file and keyPairId
			String distributionDomain = "d3adh8xh2m5mz3.cloudfront.net";
			File privateKeyFile = new File("/Users/DavidNix/Arup/Amazon/KeysCredentials/pk-APKAIKVQM23ZJHPC4PEQ.pem");
			String keyPairId = "APKAIKVQM23ZJHPC4PEQ";

			//Amazon S3 file to make a signed URL
			//Note, don't include the bucket name (aruplab), the Distribution prepends this
			//String s3PathAndFileName = "aruplab/results/dev/2016/NixTest/FullRunCNV36B_final.bai";
			String s3PathAndFileName = "aruplab/results/dev/2016/NixTest/CompV1Run4_NA12878_1_final.bam";

			//Set a date for when the URL will expire
			Date urlLinkExpirationDate = DateUtils.parseISO8601Date("2017-03-14T00:00:00.000Z");

			String url1 = CloudFrontUrlSigner.getSignedURLWithCannedPolicy(
					Protocol.https, distributionDomain, privateKeyFile,
					s3PathAndFileName, keyPairId, urlLinkExpirationDate);

			System.out.println("wget -O test "+url1);
			
			
			//Note this URL will throw an access denied for 5+min until it becomes available, so patience!

		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
