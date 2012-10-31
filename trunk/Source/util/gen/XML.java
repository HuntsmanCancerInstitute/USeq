package util.gen;
import javax.xml.parsers.*;
import org.xml.sax.*;
/**
 Utilities for XML parsing.
 */
public class XML {

	/** create a new XML reader, returns null if a problem is encountered. */
	public static XMLReader makeXMLReader() { 
		try{
			SAXParserFactory saxParserFactory   = SAXParserFactory.newInstance();
			SAXParser        saxParser = saxParserFactory.newSAXParser();
			XMLReader        parser    = saxParser.getXMLReader();
			return parser; 
		} catch (Exception e){
			e.printStackTrace();
		}
		return null;
	}
}

