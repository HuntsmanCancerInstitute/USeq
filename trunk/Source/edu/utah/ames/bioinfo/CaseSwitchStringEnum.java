package edu.utah.ames.bioinfo;

/**
 * 
 * @author darren.ames@hci.utah.edu
 *
 */

public class CaseSwitchStringEnum {
	
	private static enum Apps {
		APP2, APP3, MRNASEQ, APP1, EXCAPNIM, TDNASEQ, SMRNASEQ, APP4, DNASEQ, APP5, APP6;
	}
	static String value = "APP2"; // assume input
	static Apps apps = Apps.valueOf(value); // surround with try/catch
	
	public static void main(String[] args) {
		figureThisOut();
	}
	
	public static void figureThisOut() {
		 try {
			 switch(apps) {
			    case APP2: doStuff1(); 
			    doStuff2();
			    doStuff3();
			    break;
			    case APP1: doStuff2(); break;
			    case MRNASEQ: doStuff4(); break;
			    case APP5: doStuff3(); break;
			    case APP3: doStuff1(); break;
			    case EXCAPNIM: doStuff2(); break;
			    case TDNASEQ: doStuff3(); break;
			    case SMRNASEQ: doStuff4(); break;
			    case APP4: doStuff2(); break;
			    case DNASEQ: doStuff3(); break;
			    case APP6: doStuff4(); break;
			    default: break;
				}
			}
		 catch (Exception eof) {
				System.out.println("error");
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
