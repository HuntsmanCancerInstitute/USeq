package edu.utah.ames.bioinfo;

/**
 * 
 * @author darren
 *
 */

public class SequencingApplicationCodeEnumerator {

	//current sequencing application codes used in our pipeline
	private static enum Apps {
		APP2, APP3, MRNASEQ, APP1, EXCAPSSC, TDNASEQ, SMRNASEQ, APP4, DNASEQ, 
		APP5, APP6, CHIPSEQ, APP33, APP29, APP30, APP26, APP20, APP21, APP22, 
		APP23, APP24, APP25, APP8, APP11, APP12, APP13, APP27, APP9;
	}

	static String value = "EXCAPNIM"; // assume input
	static Apps apps = Apps.valueOf(value); // surround with try/catch

	public static void main(String[] args) {
		caseSwitchEnum();
	}

	/**
	 * Method to perform case/switch on enumerated sequencing application codes 
	 * to determine application-specific alignment parameters
	 */
	public static void caseSwitchEnum() {
		try {
			switch(apps) {
			case APP2: doStuff1(); break;
			case CHIPSEQ: doStuff1(); break;
			case APP1: doStuff2(); break;
			case MRNASEQ: doStuff4(); break;
			case APP5: doStuff3(); break;
			case APP3: doStuff1(); break;
			case EXCAPSSC: doStuff2(); break;
			case TDNASEQ: doStuff3(); break;
			case SMRNASEQ: doStuff4(); break;
			case APP4: doStuff2(); break;
			case DNASEQ: doStuff3(); break;
			case APP6: doStuff4(); break;
			case APP33: doStuff1(); break;
			case APP29: doStuff2(); break;
			case APP30: doStuff4(); break;
			case APP26: doStuff3(); break;
			case APP20: doStuff1(); break;
			case APP21: doStuff2(); break;
			case APP22: doStuff3(); break;
			case APP23: doStuff4(); break;
			case APP24: doStuff2(); break;
			case APP25: doStuff3(); break;
			case APP8: doStuff4(); break;
			case APP11: doStuff1(); break;
			case APP12: doStuff2(); break;
			case APP13: doStuff4(); break;
			case APP27: doStuff3(); break;
			case APP9: doStuff1(); break;
			default: break;
			}
		}
		catch (Exception eof) {
			System.out.println("Error. Not a recognized sequencing application code!");
		}
	}

	public static void doStuff1() {
		System.out.println("1");
	}

	public static void doStuff2() {
		System.out.println("2");
	}

	public static void doStuff3() {
		System.out.println("3");
	}

	public static void doStuff4() {
		System.out.println("4");
	}
}


