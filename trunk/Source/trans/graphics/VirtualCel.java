package trans.graphics;
import javax.swing.JFrame;
import util.gen.*;
import java.io.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;


/**Builds virtual chips from float[][] files and saves each as a PNG image file.*/
public class VirtualCel {
	//fields
	private float[][] intensityValues;
	private File[] celFiles;
	private int maxIntensityValueToPlot = 20000;
	
	//constructors
	/**For generating png files.*/
	public VirtualCel(String[] args){
		//kill the need to display
		System.setProperty("java.awt.headless","true");
		processArgs(args);
		try {
			//for each file make a panel and save
			VirtualCelPanel panel;
			for (int i=0; i<celFiles.length; i++){
				System.out.println("\nLoading -> "+celFiles[i].getName());
				intensityValues = (float[][])IO.fetchObject(celFiles[i]);
				String fileName = celFiles[i].getName().replaceFirst(".cela",".png");
				panel = new VirtualCelPanel(intensityValues, maxIntensityValueToPlot, new File( celFiles[i].getParent()+File.separator+ fileName));	
			}
		}catch (Exception e){
			System.out.println("\nAn error occured while trying to read one of your cela files?!\n");
			e.printStackTrace();
		}
		
	}
	
	/**For integration with CelMasker.*/
	public VirtualCel(float[][] intensityValues, CelMasker celMasker){
		VirtualCelFrame frame = new VirtualCelFrame(intensityValues, celMasker.getMaxIntensityValueToPlot(), celMasker.getMaskValue());
		frame.setCelMasker(celMasker);
		frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		frame.show();
	}
	
	/**To write a matrix to png*/
	public VirtualCel(float[][] intensityValues, File file, int maxIntensityValueToPlot){
		//kill the need to display
		System.setProperty("java.awt.headless","true");
		VirtualCelFrame frame = new VirtualCelFrame();
		frame.makePNGPanel(intensityValues, maxIntensityValueToPlot, file);
	}
	
	public static void main(String[] args) {
		if (args.length == 0) {
			printDocs();
			System.exit(0);
		}
		new VirtualCel(args);
	}
	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                           Virtual Cel: August 2005                               **\n" +
				"**************************************************************************************\n" +
				"VC builds virtual chips from converted cel files (see CovertCelFiles.java) and saves\n" +
				"each as a PNG image file. Examine these images for inconsistencies that need to be\n" +
				"masked using the Cel Masker application. Use the following options when running VC:\n\n" +
				
				"-c Full path file text for the directory containing converted 'xxx.cela' float[][] cel\n" +
				"     files.\n" +
				"-m Maximum intensity value to color scale, defaults to 20000.\n"+
				"\n" +
				
				"Example: java -Xmx256M -jar pathTo/T2/Apps//VirtualCel -c /affy/chips/\n\n"+

		"**************************************************************************************\n");
	}
	/**This method will process each argument and assign any new varibles*/
	public void processArgs(String[] args){
		File dirFile = null;
		Pattern pat = Pattern.compile("-[a-z]");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'c': dirFile = new File(args[i+1]); i++; break;
					case 'm': maxIntensityValueToPlot =Integer.parseInt(args[i+1]); i++; break;
					case 'h': printDocs(); System.exit(0);
					default: {
						System.out.println("\nProblem, unknown option! " + mat.group()+"\n");
						System.exit(0);
					}
					}
				}
				catch (Exception e){
					System.out.print("\nSorry, something doesn't look right with this parameter request: -"+test);
					System.out.println();
					System.exit(0);
				}
			}
		}
		//extract cel files
		if (dirFile != null){
			celFiles = IO.extractFiles(dirFile, ".cela");
			if (celFiles.length == 0){
				System.out.println("\nError: No, 'xxx.cela' files were found?!\n");
				System.exit(0);
			}
		}
		else {
			System.out.println("\nError: Please enter a directory containing 'xxx.cela' files?!\n");
			System.exit(0);
		}
	}
}