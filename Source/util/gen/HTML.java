package util.gen;
import java.util.regex.*;
import org.apache.commons.beanutils.BeanUtils;
import org.apache.log4j.Logger;
import javax.servlet.http.*;
import java.util.*;
import org.apache.commons.fileupload.*;


/**
 * Static methods for HTML related tasks.
 */
public class HTML {
	
	/**Parses a multipart request returning a HashMap of fieldName: String, 
	 * if the same fieldName is found more than once, a String[] is created. 
	 * Files are referenced by their fieldName : FileItem obj. Call FileItem.write(File) to save 
	 * it to disk.  Assumes file field names are always unique. */
	public static HashMap parseMultiPartRequest(HttpServletRequest request){
		//Process multipart params, look for multiple markers
		HashMap nameValue = new HashMap();
		try{
			DiskFileUpload upload = new DiskFileUpload();
			List items = upload.parseRequest(request);
			Iterator iter = items.iterator();
			String fieldName;
			ArrayList arrayList = new ArrayList();
			FileItem item;
			boolean arrayListsMade = false;
			
			while (iter.hasNext()) {
				item = (FileItem) iter.next();	
				fieldName = item.getFieldName();
				//System.out.println("fieldName: "+fieldName+" value: "+item.getString());
				//check if it exists, will only exist if it is a FormField, assumes file field names will always be different.
				if (nameValue.containsKey(fieldName)){
					//text exits
					Object obj = nameValue.get(fieldName);
					//is it an ArrayList?
					if (obj.getClass().isInstance(arrayList)){
						//it is an ArrayList so add new item String
						((ArrayList)obj).add(item.getString());
					}
					else{
						//is not an ArrayList, make new and add both
						ArrayList al = new ArrayList();
						al.add(obj);	//old
						al.add(item.getString()); //new
						nameValue.put(fieldName, al);	//replace entry
						//set flag
						arrayListsMade = true;
					}
				}
				else{
					//doesn't exits
					if (item.isFormField()) nameValue.put(fieldName, item.getString());
					else nameValue.put(fieldName,item);
				}
			}
			//run thru converting ArrayLists to String[]
			if (arrayListsMade){
				iter = nameValue.keySet().iterator();
				Object obj;
				Object name;
				while (iter.hasNext()){
					name = iter.next();
					obj = nameValue.get(name);
					if (obj.getClass().isInstance(arrayList)){
						nameValue.put(name, Misc.stringArrayListToStringArray((ArrayList)obj));
					}
					
				}
			}
		}catch (FileUploadException e){
			System.err.println("Problem with parseMultiPartRequest()");
			e.printStackTrace();
		}
		return nameValue;
	}

	
	/** For extracting either a String[] or making a new String[] by casting the object as a String.*/
	public static String[] parseStringArray(Object obj){
		if (obj == null) return null;
			if (obj.getClass().getName().equals("java.lang.String")){
				return new String[]{(String)obj};
			}
			return (String[])obj;
	}
	

	/**Filters a String for bad characters replacing them with html equivalents.
	 * Truncates the text to a lengthCutOff*/
	public static String filterField(String field, int lengthCutOff){
		if (Misc.isEmpty(field)) return "";
		field = field.replaceAll("<", "&lt;");
		field = field.replaceAll(">", "&gt;");
		field = field.replaceAll("\"", "&#34;");
		field = field.replaceAll("'", "&#39;");  
		field = field.replaceAll("`", "&#39;");
		field = field.replaceAll("Õ", "&#39;");
		if (field.length()>lengthCutOff) return field.substring(0,lengthCutOff);
		return field;
	}
	
	/**Filters a String[] for bad characters replacing them with html equivalents.
	 * Truncates the Strings to a lengthCutOff. Modifies original String[]*/
	public static String[] filterFields(String[] fields, int lengthCutOff){
		int num = fields.length;
		for (int i=0; i< num; i++){
			fields[i] = filterField(fields[i], lengthCutOff);
		}
		return fields;
	}
	
	
	/**Returns an upDown arrow ala Mac look and feel that will set a formParm to a given value and
	 * submit the form provided the form is named "form" (ie <form text= "form" action="/ etc >).
	 * If more than one arrow is going to be put on a page then provide different buttonIncrements.
	 * ie buttonIncrement =1 then 2 then 3...etc.
	 * Hint, use a hidden value to hold the new value.*/
	public static String fetchArrowButton (String formParamToSet, String value, String pathToArrowGifs, int buttonIncrement){
		StringBuffer sb = new StringBuffer("<A HREF=\"javascript:document.form.submit()\" onmouseover=\"document.form.sub_but");
		sb.append(buttonIncrement);			
		sb.append(".src='"); 
		sb.append(pathToArrowGifs);
		sb.append("PrArrow.gif'\" onmouseout=\"document.form.sub_but");
		sb.append(buttonIncrement);			
		sb.append(".src='"); 
		sb.append(pathToArrowGifs);
		sb.append("unPrArrow.gif'\" onclick=\"document.form.");
		sb.append(formParamToSet);				
		sb.append(".value='");
		sb.append(value);
		sb.append("'; return val_form_this_page()\"><IMG SRC=\"");
		sb.append(pathToArrowGifs);
		sb.append("unPrArrow.gif\" width='16' height='22' BORDER='0' align='middle' NAME='sub_but"); 
		sb.append(buttonIncrement);
		sb.append("'></A>");				
		return sb.toString();
	}
	/**Example table row with arrow buttons.*/
	public static String makeArrowedTableRow(String[] headings, String[] values){
		StringBuffer sb = new StringBuffer("<tr>");
		int num = headings.length;
		for (int i=0; i<num; i++){
			sb.append("<td font color=blue>");
			sb.append(headings[i]);
			sb.append("&nbsp;");
			sb.append(HTML.fetchArrowButton("switchOrder",values[i], "/Images", i));
		}
		sb.append("</tr>");
		return sb.toString();
	}	
	
	/**Returns a selected option group with the toSelect selected if exactly matches one of the String[] options*/
	public static String fetchOptionGroup(String toSelect, String[] options){
		if (toSelect == null || options== null) return "";
		StringBuffer sb = new StringBuffer();
		int len = options.length;
		for (int i=0; i<len; i++){
			if (toSelect.equals(options[i])) sb.append("<option selected>"+options[i]+"</option>");
			else sb.append("<option>"+options[i]+"</option>");
		}
		return sb.toString();
	}
	
	/**Returns a selected option group with the toSelect selected if at all present in the String[] options.
	 * Case-insensitive.*/
	public static String fetchOptionGroupLoose(String toSelect, String[] options){
		//System.out.println("toSelect: '"+toSelect+"'");
		Pattern pat = Pattern.compile(".*"+toSelect.trim()+".*", Pattern.CASE_INSENSITIVE);
		Matcher mat;
		StringBuffer sb = new StringBuffer();
		int len = options.length;
		for (int i=0; i<len; i++){
			mat = pat.matcher(options[i]);
			if (mat.matches()) sb.append("<option selected>"+options[i]+"</option>");
			else sb.append("<option>"+options[i]+"</option>");
		}
		return sb.toString();
	}	
	
	
	/**Returns a selected option group with the toSelect selected if present in the String[].
	 * This method assumes that an id is at the beginning of some options (ie options[i] = 
	 * "7: Human: modern").  The method will strip off the "7:" and display "Human: modern"
	 * but set the return value as the original options[i].  This is just a good way to hide the
	 * return id number from the user.
	 * */
	public static String fetchValueOptionGroup(String toSelect, String[] options){	
		StringBuffer sb = new StringBuffer();
		int len = options.length;
		String selected = "";
		int firstIndex = 0;
		for (int i=0; i<len; i++){
			if (toSelect.equals(options[i])) selected = " SELECTED"; 
			else selected = "";
			sb.append ("<OPTION VALUE=\"");
			sb.append (options[i]);
			sb.append ("\"");
			sb.append (selected);
			sb.append (">");
			//create clipped option (ie remove first "anything...:") to display in browser, otherwise show options[i]
			firstIndex = options[i].indexOf(":");
			if (firstIndex!=-1) sb.append(options[i].substring(firstIndex+1).trim());
			else sb.append (options[i]);
			sb.append ("</OPTION>");
		}
		return sb.toString();
	}	
	
	
	/**Returns a selected option group with the toSelect selected if present in the String[]
	 * of options.  Wiggle room allowed in toSelect, checks for blanks too.  
	 * Set toSelect ==null if you just want to list the options.*/
	public static String fetchMultipleOptionGroup(String[] toSelect, String[] options) {
		StringBuffer sb = new StringBuffer();
		int len = options.length;
		if (toSelect==null) {
			for (int i = 0; i < len; i++) {
				if (Misc.isEmpty(options[i])) continue; //check for blanks
				sb.append("<option>" + options[i] + "</option>");
			}
		}
		else {
			String toSelectString = Misc.stringArrayToString(toSelect,":");
			for (int i = 0; i < len; i++) {
				if (Misc.isEmpty(options[i])) continue; //check for blanks
				if (toSelectString.indexOf(options[i])!=-1)
					sb.append("<option selected>" + options[i] + "</option>");
				else
					sb.append("<option>" + options[i] + "</option>");
			}
		}
		return sb.toString();
	}	
	
	
	/**Returns a selected check box group with the toSelect selected if present in the String[]*/
	public static String fetchCheckGroup(String name, String toCheck, String[] boxes){	
		StringBuffer sb = new StringBuffer();
		int len = boxes.length;
		for (int i=0; i<len; i++){
			if (Pattern.matches(".*"+boxes[i]+".*",toCheck)) sb.append(
					"<INPUT TYPE='checkbox' NAME='"+name+"' VALUE='"+ boxes[i]+"' CHECKED>"+boxes[i]+"&nbsp;&nbsp;");
			else sb.append("<INPUT TYPE='checkbox' NAME='"+name+"' VALUE='"+ boxes[i]+"'>"+boxes[i]+"&nbsp;&nbsp;");
		}
		return sb.toString();
	}
	
	
	
	/**Returns a selected radio group with the toSelect selected if present in the String[]*/
	public static String fetchRadioGroup(String name, String toCheck, String[] boxes){	
		StringBuffer sb = new StringBuffer();
		int len = boxes.length;
		for (int i=0; i<len; i++){
			if (Pattern.matches(".*"+boxes[i]+".*",toCheck)) sb.append(
					"<INPUT TYPE='radio' NAME='"+name+"' VALUE='"+ boxes[i]+"' CHECKED>"+boxes[i]+"&nbsp;&nbsp;");
			else sb.append("<INPUT TYPE='radio' NAME='"+name+"' VALUE='"+ boxes[i]+"'>"+boxes[i]+"&nbsp;&nbsp;");
		}
		return sb.toString();
	}	
	
	
	/**Returns a radio input line with CHECKED if the value is equal to the potential*/
	public static String fetchRadio(String name, String value, String potential){	
		if (potential.equals(value))  return "<INPUT TYPE='radio' NAME='"+name+"' VALUE='"+ value+"' CHECKED>";
		return "<INPUT TYPE='radio' NAME='"+name+"' VALUE='"+ value+"'>";
	}
	
	
	
	/**Returns a selected radio group with the toSelect selected if present in the String[]
	 * Also throws in a little javascript to submit the form if the button is clicked.
	 * Add a text = 'form' attribute to the <form .... tag.
	 * Use the nbsp String to add multiple &nbsp; concats to vary space between buttons.*/
	public static String fetchRadioGroupSubmit(String name, String toCheck, String[] boxes, String nbsp){	
		StringBuffer sb = new StringBuffer();
		int len = boxes.length;
		for (int i=0; i<len; i++){
			if (Pattern.matches(".*"+boxes[i]+".*",toCheck)) sb.append(
					"<INPUT TYPE='radio' NAME='"+name+"' VALUE='"+ boxes[i]+"' onClick=\"document.form.submit()\" CHECKED >"+boxes[i]+nbsp);
			else sb.append("<INPUT TYPE='radio' NAME='"+name+"' VALUE='"+ boxes[i]+"' onClick=\"document.form.submit()\">"+boxes[i]+nbsp);
		}
		return sb.toString();
	}	
	
	/**Fills in all the params of a bean.*/
	public static void populateBean(Object bean, HttpServletRequest request) { 
		try {
			BeanUtils.populate(bean, request.getParameterMap());
		}
		catch (Exception e) {
			System.err.println("Problem with populateBean()");
			e.printStackTrace();
		}
	}	

	
	
}
