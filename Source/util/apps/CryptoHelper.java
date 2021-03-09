package util.apps;

import java.io.File;
import java.security.Key;
import java.security.SecureRandom;

import javax.crypto.Cipher;
import javax.crypto.KeyGenerator;
import javax.crypto.SecretKey;
import javax.crypto.spec.IvParameterSpec;

import org.apache.commons.codec.binary.Base64;

import util.gen.IO;

public class CryptoHelper {
	
	//fields
	private Key key;
	
	//constructor
	public CryptoHelper( Key key ) {
		this.key = key;
	}
	
	public static void main( String [] args ) throws Exception {
		Key key = generateSymmetricKey();

		String plaintext = "My secret.";
		System.out.println( plaintext );

		byte[] iv = generateIV();
		String ciphertext = encrypt(iv, plaintext, key );
		System.out.println( ciphertext );

		String decrypted = decrypt( ciphertext, key );
		System.out.println( decrypted );


	}


	public static String encrypt( byte [] iv, String plaintext, Key key ) throws Exception {
		byte [] decrypted = plaintext.getBytes();
		byte [] encrypted = encrypt( iv, decrypted, key );
		StringBuilder ciphertext = new StringBuilder();
		ciphertext.append( Base64.encodeBase64String( iv ) );
		ciphertext.append( ":" );
		ciphertext.append( Base64.encodeBase64String( encrypted ) );
		return ciphertext.toString();
	}

	public static String decrypt( String ciphertext, Key key ) throws Exception {
		String [] parts = ciphertext.split( ":" );
		byte [] iv = Base64.decodeBase64( parts[0] );
		byte [] encrypted = Base64.decodeBase64( parts[1] );
		byte [] decrypted = decrypt( iv, encrypted, key);
		return new String( decrypted );
	}
	public static byte [] generateIV() {
		SecureRandom random = new SecureRandom();
		byte [] iv = new byte [16];
		random.nextBytes( iv );
		return iv;
	}
	public static Key generateSymmetricKey() throws Exception {
		KeyGenerator generator = KeyGenerator.getInstance( "AES" );
		SecretKey key = generator.generateKey();
		return key;
	}
	public static byte [] encrypt( byte [] iv, byte [] plaintext, Key key ) throws Exception {
		Cipher cipher = Cipher.getInstance( key.getAlgorithm() + "/CBC/PKCS5Padding" );
		cipher.init( Cipher.ENCRYPT_MODE, key, new IvParameterSpec( iv ) );
		return cipher.doFinal( plaintext );
	}
	public static byte [] decrypt( byte [] iv, byte [] ciphertext, Key key ) throws Exception {
		Cipher cipher = Cipher.getInstance( key.getAlgorithm() + "/CBC/PKCS5Padding" );
		cipher.init( Cipher.DECRYPT_MODE, key, new IvParameterSpec( iv ) );
		return cipher.doFinal( ciphertext );
	}

}