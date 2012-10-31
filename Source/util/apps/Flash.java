package util.apps;
import java.io.*;
import util.gen.*;
import java.util.*;
import javax.swing.JOptionPane;

public class Flash {
	private File file;
	private String trimmedFileName;
	private FlashCard[] flashCards;
	private boolean keysFirst = true;
	private boolean queryTilCorrect = true;
	private String answer = "";
	private TextFrame textFrame;
	
	public Flash (String[] args){
		//read in file
		
		file = new File (args[0]);
		flashCards = parseFile(file);
		trimmedFileName = Misc.removeExtension(file.getName());
		
		//make text window
		//make text frame for notes
		textFrame = new TextFrame(0, 0, 500, 500, "Notes");
		//textFrame.getTextArea().setFont(new Font("Courier"));
		textFrame.setVisible(true);
		
		//get order and shuffle
		queryOrder();
		queryShuffle();
		
		//test
		interrogate();
	}
	
	public void interrogate(){
		boolean go = true;
		while (go){
			//query
			boolean complete = true;
			for (int i=0; i< flashCards.length; i++){
					//what is the question?
					String question;
					if (keysFirst) question = flashCards[i].key;
					else question = flashCards[i].value;
					//ask it
					Object ob =JOptionPane.showInputDialog(null,question,trimmedFileName,JOptionPane.PLAIN_MESSAGE,null,null,answer);
					if (ob == null) {
						textFrame.append("\nCiao!\n");
						Misc.printExit("\nCiao!\n");
					}
					String response = ((String)ob).trim();
					if (Misc.isEmpty(response)) i--;
					else {
						if (keysFirst && response.toLowerCase().equalsIgnoreCase(flashCards[i].value)) {
							textFrame.append("Giusto! "+flashCards[i]+"\n");
							answer = "";
						}
						else if (keysFirst == false && response.toLowerCase().equalsIgnoreCase(flashCards[i].key)) {
							textFrame.append("Giusto! "+flashCards[i]+"\n");
							answer = "";
						}
						else {
							answer = response;
							StringBuffer sb = new StringBuffer("Sbagliato :( ");
							if (flashCards[i].numFailed == 0) addFlashCard(flashCards[i]);
							if (flashCards[i].numFailed > 0) {
								sb.append("\tHint: "+flashCards[i].misc);
								if (keysFirst) sb.append(" - "+flashCards[i].value);
								else sb.append(" - "+flashCards[i].key);
							}
							sb.append("\n");
							textFrame.append(sb.toString());
							flashCards[i].numFailed++;
							if (queryTilCorrect) i--;
							else complete = false;
						}
					}
			}
			if (complete){
				textFrame.append("\nComplete! Resetting.\n");
				TreeSet incorrectAnswers = new TreeSet();
				//reset
				for (int i=0; i< flashCards.length; i++) {
					if (flashCards[i].numFailed !=0){
						flashCards[i].numFailed = 0;
						incorrectAnswers.add(flashCards[i].toString());
					}
				}
				//print incorrects
				if (incorrectAnswers.size() !=0){
					textFrame.append("\nWatch out for:\n");
					Iterator it = incorrectAnswers.iterator();
					while (it.hasNext()){
						textFrame.append((String)it.next() +"\n");
					}
				}
				textFrame.append("\n");
				queryOrder();
				queryShuffle();
			}
		}
	}
	
	public void queryOrder(){
		try {
			System.out.print("Would you like to be shown the keys? (yes or no) ");
			BufferedReader stdin = new BufferedReader (new InputStreamReader(System.in));
			String response = stdin.readLine();
			if (response !=null && response.matches("[nN].*")){
				keysFirst = false;
			}
			else keysFirst = true;
		} catch (IOException e){
			e.printStackTrace();
		}
	}
	
	public void queryShuffle(){
		try {
			System.out.print("Shuffle the cards? (yes or no) ");
			BufferedReader stdin = new BufferedReader (new InputStreamReader(System.in));
			String response = stdin.readLine();
			if (response !=null && response.matches("[yY].*"))shuffleCards();
		} catch (IOException e){
			e.printStackTrace();
		}
	}
	
	/**Adds flash card to the stop of the array.*/
	public void addFlashCard(FlashCard card){	
		FlashCard[] fc = new FlashCard[flashCards.length+1];
		fc[flashCards.length] = card;
		System.arraycopy(flashCards,0,fc,0,flashCards.length);
		flashCards = fc;			
	}
	
	/**Randomly permutates the FlashCard[].*/
	public void shuffleCards (){
		int len = flashCards.length;
		FlashCard current;
		FlashCard random;
		int index;
		Random generator = new Random();
		for (int i=0; i<len; i++){
			index = generator.nextInt(len);
			current = flashCards[i];
			random = flashCards[index];
			flashCards[i] = random;
			flashCards[index] = current;
		}
	}
	
	public static void main(String[] args) {
		if (args.length==0) Misc.printExit("\nEnter the full path text for a tab delimited text file " +
		"containing columns of key, value, and optionally notes or hints.\n");
		new Flash(args);
	}
	
	private class FlashCard{
		String key;
		String value;
		String misc = "";
		int numFailed = 0;
		
		public FlashCard (String key, String value, String misc){
			this.key = key;
			this.value = value;
			this.misc = misc;
		}
		
		public String toString(){
			StringBuffer sb = new StringBuffer(key);
			sb.append("  \t");
			sb.append(value);
			sb.append("  \t");
			sb.append(misc);
			return sb.toString();
		}
	}
	
	/**Strips out " that stupid excel adds in. Skips lines beginning with #*/
	public FlashCard[] parseFile(File file){
		ArrayList al = new ArrayList();
		try {
			BufferedReader in = new BufferedReader( new FileReader(file) );
			String line;
			while ((line = in.readLine()) != null){
				line = line.trim();
				if (line.length() == 0 || line.startsWith("#")) continue;
				line = line.replaceAll("\"", "");
				String[] tokens = line.split("\\t");
				if (tokens.length == 3) al.add(new FlashCard(tokens[0].trim(), tokens[1].trim(), tokens[2].trim()));
				else if (tokens.length == 2) al.add(new FlashCard(tokens[0].trim(), tokens[1].trim(), ""));
			}
			
		} catch (IOException e){
			e.printStackTrace();
		}
		FlashCard[] fcs = new FlashCard[al.size()];
		al.toArray(fcs);
		return fcs;
	}
	
}
