package util.gen;

import org.json.JSONException;
import org.json.JSONObject;

public class Json {

	/**Returns null if not found.*/
	public static String getStringAttribute(JSONObject obj, String key) throws JSONException {
		if (obj.has(key) && obj.isNull(key) == false) {
			String trimmed = obj.getString(key).trim();
			if (trimmed.length() ==0) return null;
			return trimmed;
		}
		return null;
	}
	
	/**Returns null if not found.*/
	public static Integer getIntegerAttribute(JSONObject obj, String key) throws JSONException {
		if (obj.has(key) && obj.isNull(key) == false) return obj.getInt(key);
		return null;
	}
	
	/**Returns null if not found.*/
	public static Double getDoubleAttribute(JSONObject obj, String key) throws JSONException {
		if (obj.has(key) && obj.isNull(key) == false) return obj.getDouble(key);
		return null;
	}
	
	/**Returns null if not found, is a JSON null, or length of String rep is zero*/
	public static String forceGetString(JSONObject obj, String key) throws JSONException {
		if (obj.has(key) == false || obj.isNull(key) == true) return null;
		Object o = obj.get(key);
		String s = o.toString().trim();
		if (s.length() == 0) return null;
		return s;
	}
}
