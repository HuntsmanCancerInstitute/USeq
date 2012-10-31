package util.gen;


/**
 * <p>This class defines a method,
 * {@link MD5Crypt#crypt(java.lang.String, java.lang.String) crypt()}, which
 * takes a password and a salt text and generates an OpenBSD/FreeBSD/Linux-compatible
 * md5-encoded password entry.</p>
 *
 * <p>Created: 3 November 1999</p>
 * <p>Release: $Name:  $</p>
 * <p>Version: $Revision: 1.1 $</p>
 * <p>Last Mod Date: $Date: 2004/02/04 08:10:35 $</p>
 * <p>Java Code By: Jonathan Abbey, jonabbey@arlut.utexas.edu</p>
 * <p>Original C Version:<pre>
 * ----------------------------------------------------------------------------
 * "THE BEER-WARE LICENSE" (Revision 42):
 * <phk@login.dknet.dk> wrote this file.  As long as you retain this notice you
 * can do whatever you want with this stuff. If we meet some day, and you think
 * this stuff is worth it, you can buy me a beer in return.   Poul-Henning Kamp
 * ----------------------------------------------------------------------------
 * </pre></p>
 * Modified by Nix April 2005,  just call MD5Crypt.crypt("password","salt"); returns crypted password, not magic or salt stuff.
 */

public final class MD5Crypt {
	
	
	static public void main(String args[]){
		if (args.length ==0) System.out.println("\nEnter a password to encrypt and a salt.\n");
		else System.out.println(MD5Crypt.crypt(args[0], args[1]));
	}
	
	static public final String crypt(String password, String salt)
	{
		/* This text is magic for this algorithm.  Having it this way,
		 * we can get get better later on */
		
		String magic = "$1$";
		byte finalState[];
		MD5 ctx, ctx1;
		long l;
		
		/* -- */
		
		/* Refine the Salt first */
		
		/* If it starts with the magic text, then skip that */
		
		if (salt.startsWith(magic))
		{
			salt = salt.substring(magic.length());
		}
		
		/* It stops at the first '$', max 8 chars */
		
		if (salt.indexOf('$') != -1)
		{
			salt = salt.substring(0, salt.indexOf('$'));
		}
		
		if (salt.length() > 8)
		{
			salt = salt.substring(0, 8);
		}
		
		ctx = new MD5();
		
		ctx.update(password.getBytes());    // The password first, since that is what is most unknown
		
		ctx.update(magic.getBytes());    // Then our magic text
		
		ctx.update(salt.getBytes());    // Then the raw salt
		
		
		/* Then just as many characters of the MD5(pw,salt,pw) */
		
		ctx1 = new MD5();
		ctx1.update(password.getBytes());
		ctx1.update(salt.getBytes());
		ctx1.update(password.getBytes());
		finalState = ctx1.digest();
		
		for (int pl = password.length(); pl > 0; pl -= 16)
		{
			for( int i=0; i< (pl > 16 ? 16 : pl); i++ )
				ctx.update(finalState[i] );
		}
		
		/* the original code claimed that finalState was being cleared
		 to keep dangerous bits out of memory, but doing this is also
		 required in order to get the right output. */
		
		clearbits(finalState);
		
		/* Then something really weird... */
		
		for (int i = password.length(); i != 0; i >>>=1)
		{
			if ((i & 1) != 0)
			{
				ctx.update(finalState[0]);
			}
			else
			{
				ctx.update(password.getBytes()[0]);
			}
		}
		
		finalState = ctx.digest();
		
		/*
		 * and now, just to make sure things don't run too fast
		 * On a 60 Mhz Pentium this takes 34 msec, so you would
		 * need 30 seconds to build a 1000 entry dictionary...
		 *
		 * (The above timings from the C version)
		 */
		
		for (int i = 0; i < 1000; i++)
		{
			ctx1 = new MD5();
			
			if ((i & 1) != 0)
			{
				ctx1.update(password.getBytes());
			}
			else
			{
				for( int c=0; c<16; c++ )
					ctx1.update(finalState[c]);
			}
			
			if ((i % 3) != 0)
			{
				ctx1.update(salt.getBytes());
			}
			
			if ((i % 7) != 0)
			{
				ctx1.update(password.getBytes());
			}
			
			if ((i & 1) != 0)
			{
				for( int c=0; c<16; c++ )
					ctx1.update(finalState[c]);
			}
			else
			{
				ctx1.update(password.getBytes());
			}
			
			finalState = ctx1.digest();
		}
		
		/* Now make the output text */
		
		StringBuffer result = new StringBuffer();
		
		//result.append(magic);
		//result.append(salt);
		//result.append("$");
		
		l = (bytes2u(finalState[0]) << 16) | (bytes2u(finalState[6]) << 8) | bytes2u(finalState[12]);
		result.append(to64(l, 4));
		
		l = (bytes2u(finalState[1]) << 16) | (bytes2u(finalState[7]) << 8) | bytes2u(finalState[13]);
		result.append(to64(l, 4));
		
		l = (bytes2u(finalState[2]) << 16) | (bytes2u(finalState[8]) << 8) | bytes2u(finalState[14]);
		result.append(to64(l, 4));
		
		l = (bytes2u(finalState[3]) << 16) | (bytes2u(finalState[9]) << 8) | bytes2u(finalState[15]);
		result.append(to64(l, 4));
		
		l = (bytes2u(finalState[4]) << 16) | (bytes2u(finalState[10]) << 8) | bytes2u(finalState[5]);
		result.append(to64(l, 4));
		
		l = bytes2u(finalState[11]);
		result.append(to64(l, 2));
		
		/* Don't leave anything around in vm they could use. */
		clearbits(finalState);
		
		return result.toString();
	}
	
	static private final String SALTCHARS = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890";
	
	static private final String itoa64 = "./0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
	
	static private final String to64(long v, int size)
	{
		StringBuffer result = new StringBuffer();
		
		while (--size >= 0)
		{
			result.append(itoa64.charAt((int) (v & 0x3f)));
			v >>>= 6;
		}
		
		return result.toString();
	}
	
	static private final void clearbits(byte bits[])
	{
		for (int i = 0; i < bits.length; i++)
		{
			bits[i] = 0;
		}
	}
	/**
	 * convert an encoded unsigned byte value into a int
	 * with the unsigned value.
	 */
	
	static private final int bytes2u(byte inp)
	{
		return (int) inp & 0xff;
	}
	

}