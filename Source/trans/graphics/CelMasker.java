package trans.graphics;
import java.io.*;
import java.util.regex.*;
import util.bio.parsers.*;
import trans.main.*;
import trans.misc.Util;
import util.gen.*;
import java.util.*;

/**
 * Application to draw and mask a raw text cel file.
 */
public class CelMasker {
	//fields
	private int rows = 2560;
	private File celFile;
	private float[][] intensities;
	private int maxIntensityValueToPlot = 20000;
	private float maskValue = -1;
	
	public CelMasker (String[] args) {
		processArgs(args);
		System.out.println("\tLoading data from "+celFile+" ...");
		
		//load virtual slide with intensities
		if (celFile.getName().endsWith("cela")) intensities = (float[][])IO.fetchObject(celFile);
		else intensities = Util.createVirtualCel (celFile);
		
		//check if mask value is set
		if (maskValue == -1){
			//not set therefore set to median
			float[] collapsed = Num.collapseFloatArray(intensities);
			Arrays.sort(collapsed);
			double median = Num.median(collapsed);
			maskValue = new Double(median).floatValue();
			System.out.println("\tSetting maskValue to the chip's median intensity -> "+maskValue);
		}
		
		//make plot
		new VirtualCel(intensities, this);	
		
	}
	
	
	public void writeCelFile(){
		System.out.println("\tWriting masked cel file...");
		//from a converted cela file?
		if (celFile.getName().endsWith("cela")){
			File masked = new File (celFile.getParent()+File.separator+"masked"+celFile.getName());
			IO.saveObject(masked, intensities);
		}
		//no from a text cel file
		else {
			//load virtual slide with .cel file intensities
			try{
				BufferedReader in = new BufferedReader(new FileReader(celFile));
				PrintWriter out = new PrintWriter(new FileWriter(celFile.getParent()+File.separator+"masked"+celFile.getName()));
				//skip header
				boolean inHeader = true;
				String line;
				int counter = 0;
				while (inHeader) {
					line = in.readLine();
					out.println(line);
					if (line.indexOf("CellHeader") != -1) inHeader = false;
				}
				//parse xy coords and intensity and add to float array
				String[] tokens;
				int x;
				int y;
				while ((line = in.readLine()) !=null){
					counter++;
					tokens = line.trim().split("\\s+");
					//check that a data line was found and not stop of file or junk
					if (tokens.length!=5) {
						out.println(line);
						continue;
					}
					//parse coords and reassign value
					x = Integer.parseInt(tokens[0]);
					y = Integer.parseInt(tokens[1]);
					tokens[2] = Float.toString(intensities[x][y]);
					line = Misc.stringArrayToString(tokens,"\t");
					out.println(line);
				}
				out.close();
				in.close();
			} catch (IOException e){
				e.printStackTrace();
			}
		}
	}
	
	
	/**This method will process each argument and assign any new varibles*/
	public void processArgs(String[] args){
		
		Pattern pat = Pattern.compile("-[a-z]");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'c': celFile = new File(args[i+1]); i++; break;
					case 'r': rows=Integer.parseInt(args[i+1]); i++; break;
					case 'v': maskValue=Float.parseFloat(args[i+1]); i++; break;
					case 'm': maxIntensityValueToPlot=(int)Double.parseDouble(args[i+1]); i++; break;
					case 'h': printDocs(); System.exit(0);
					default: System.out.println("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					System.out.print("\nSorry, something doesn't look right with this parameter request: -"+test);
					System.out.println();
					System.exit(0);
				}
			}
		}
		if (celFile == null || celFile.canRead() == false){
			System.out.println("\nError: Cannot find or read you cel file?!\n");
			System.exit(0);
		}
		
	}
	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                            Cel Masker: Feb 2006                                  **\n" +
				"**************************************************************************************\n" +
				"CM builds and displays a virtual chip from a text version cel file.  Problem areas on\n" +
				"the chip can be circled and assigned a different intensity. Use the median chip value\n" +
				"to neutrally mask the area.\n" +
				
				"-c Full path file text for a text version cel file or converted xxx.cela file\n" +
				"-r The number of rows on the chip, defaults to 2560, only needed for text cel files.\n" +
				"-m Maximum intensity value to color scale, defaults to 20000.\n"+
				"-v Value to assign to the masked region. Defaults to the median chip intensity.\n"+
				"\n" +
				"   Mouse click once and drag to create a selection ellipse.\n"+
				"   Use the arrow keys for fine adjustment.\n"+
				"   Use shift+ arrow keys to modify the ellipse.\n"+
				"   Double click the mouse to set the ellipse and create another.\n" +
				"   Use the i or o keys to zoom in and out.\n" +
				"   Press the m key to apply the masks and  redraw the image.\n" +
				"   Press the q key to apply the masks and write the modified cel file to disk.\n\n"+
				
				"Example: java -Xmx256M -jar pathTo/T2/Apps/CelMasker -c /affy/badChip.cel -r 1280\n" +
				"      -m 250\n\n"+

		"**************************************************************************************\n");
	}	
	public static void main(String[] args) {	
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new CelMasker(args);
	}
	
	
	public float getMaskValue() {
		return maskValue;
	}
	
	
	public int getMaxIntensityValueToPlot() {
		return maxIntensityValueToPlot;
	}
	
}

