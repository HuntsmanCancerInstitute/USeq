package util.gen;
import java.util.*;
import java.sql.*;


/**
 * Utility class to hold database and relevant fields and methods, some from Hall and Brown,
 * most from BioRoot. The following is a collection of stuff for interacting with a MySQL database.
 * @author nix
 */
public class SQL {
		
	//fields
	private String dbURL = "jdbc:mysql://localhost:3306/";
	private String dbDriver = "com.mysql.jdbc.Driver";
	private Connection connection;
	private Statement statement;
	private ResultSet results;
	
	//constructors
	public SQL(String database, String user, String password, String dbURL, String dbDriver){
		this.dbURL = dbURL;
		this.dbDriver = dbDriver;
		connect2Db(database, user, password);
	}
	public SQL(){}
	
	//methods
	
	/**Returns a concatinate of the results where each String in the String[] represents a results
	 * table row where each column entry is separated by the divider String.  If a particular column
	 * has no value, nothing is added, no double dividers.  Use the lengthCutOff to limit the size of
	 * the Strings; a ... is added if it is trunkated.  Returns new String[]{""} if no
	 * entries are present.*/
	public String[] fetchResultConcat(String sql, int numColumns, String divider, int lengthCutOff){
		ResultSet results = executeSQL(sql);
		ArrayList al = new ArrayList();
		String[] headers = {""};
		int realNum = numColumns+1;
		try {
			while(results.next()){
				StringBuffer sb = new StringBuffer ();
				for (int i=1; i<realNum; i++){
					String res = results.getString(i);
					if (Misc.isEmpty(res)==false){
						sb.append(res);
						sb.append(divider);
					}
				}
				//delete the last divider
				int len = sb.length();
				int lenDivider = divider.length();
				if (len>lenDivider)sb.delete(len-lenDivider, len);
				len = sb.length();
				if (len>lengthCutOff){
					sb.delete(lengthCutOff-3, len);
					sb.append("...");
				}
				al.add(sb.toString());
			}
		}
		catch (SQLException e) {
			e.printStackTrace();
		}
		int len = al.size();
		if (len !=0) {
			headers = new String[len];
			al.toArray(headers);		
		}
		return headers;
	}

	/**Returns a concatinate of the results where each String in the String[] represents a results
	 * table row where each column entry is separated by the divider String.  If a particular column
	 * has no value, nothing is added, no double dividers.  Use the lengthCutOff to limit the size of
	 * the Strings; a ... is added if it is trunkated.  Returns new String[]{""} if no
	 * entries are present.  The first query is assumed to be the id value and is always added to the stop, 
	 * even if the text was truncated.  Thus the sql should look like "SELECT id, bla, bla FROM Table.."
	 * The result will look like bla,bla,id or bla, bla...,id    
	 * This is a way of moving the ID to the stop of the result yet trunkating in needed.*/
	public String[] fetchConcatMakeIdLast(String sql, int numColumns, String divider, int lengthCutOff){
		ResultSet results = executeSQL(sql);
		ArrayList al = new ArrayList();
		String[] headers = {""};
		int realNum = numColumns+1;
		try {
			while(results.next()){
				StringBuffer sb = new StringBuffer ();
				String id = results.getString(1);
				for (int i=2; i<realNum; i++){
					String res = results.getString(i);
					if (Misc.isNotEmpty(res)){
						sb.append(res);
						sb.append(divider);
					}
				}
				//check the length
				int len = sb.length();
				if (len>lengthCutOff){
					sb.delete(lengthCutOff-(5+id.length()), len);
					sb.append("...: ");
				}
				//add on id
				sb.append(id);
				
				al.add(sb.toString());
			}
		}
		catch (SQLException e) {
			e.printStackTrace();
		}
		int len = al.size();
		if (len !=0) {
			headers = new String[len];
			al.toArray(headers);		
		}
		return headers;
	}


	/**Returns a String[] of the items in a particular column.
	 * Don't call on a table without the designated heading!  
	 * Returns a new String[]{""} if no entries are present.*/
	public String[] fetchSingleColumn(String tableName, String columnName){
		String statement = "SELECT "+columnName+" FROM "+ tableName;
		executeSQL(statement);
		ArrayList al = new ArrayList();
		String[] headers = {""};
		try {
			while(results.next()){
				al.add(results.getString(1));
			}
		}
		catch (SQLException e) {
			System.err.println("Problem fetchSingleColumn( "+tableName+", "+columnName+" )");
			e.printStackTrace();
		}
		int len = al.size();
		if (len !=0) {
			headers = new String[len];
			al.toArray(headers);		
		}
		return headers;
	}
	
	/**Returns a String[] of the results.  
	 * Returns a new String[]{""} if no entries are present.*/
	public String[] fetchMultipleCells(String sql){
		executeSQL(sql);
		ArrayList al = new ArrayList();
		String[] headers = {""};
		try {
			while(results.next()){
				al.add(results.getString(1));
			}
		}
		catch (SQLException e) {
			System.err.println("Problem fetchMultipleCells( "+sql+" )");
			e.printStackTrace();
		}
		int len = al.size();
		if (len !=0) {
			headers = new String[len];
			al.toArray(headers);		
		}
		return headers;
	}

	
	/**Closes the Connection and subsequently the Statement and ResultSet objects*/
	public void closeConnection() {
		if (connection != null) {
			try {
				connection.close();
			}
			catch (SQLException e) {
				System.err.println("Problem with closing Connection or related objects!");
				e.printStackTrace();
			}
		}
	}


	/**Method to connect to the database and create a statement object*/
	public boolean connect2Db(String database, String userName, String password) {
		//check to see if a connection has already been made
		if (connection != null)
			return true;
		String url = dbURL + database; 
		//lower case only
		try {
			//register driver		
			Class.forName(dbDriver);
			//make connection
			connection = DriverManager.getConnection(url, userName, password);
			//make statement
			statement = connection.createStatement();
		}
		catch (ClassNotFoundException e) {
			System.err.println("Can't register driver!");
			e.printStackTrace();
			return false;
		}
		catch (SQLException e) {
			System.err.println("Problem with making a connection to the database!");
			e.printStackTrace();
			return false;
		}
		return true;
	}
	
	
	/**Use to avoid a null return when attempting to fetch a ResultSet.getString().*/
	public String getResultSetString(String columnName){
		try{
			String resultString = results.getString(columnName);
			if (resultString == null || resultString.equals("null")) return "";
			else return resultString;
			
		}catch (SQLException e) {
			System.err.println("Problem with getResultSetString( "+columnName+")");
			e.printStackTrace();
		}
		return "";		
	}

	/**Use to avoid a null return when attempting to fetch a ResultSet.getInt(). Sorta stupid*/
	public int getResultSetInt(String columnName){
		try{
			int result = results.getInt(columnName);
			if (result != 0) return result;
			
		}catch (SQLException e) {
			System.err.println("Problem with getResultSetInt( "+columnName+" )");
			e.printStackTrace();
		}
		return 0;		
	}
	
	/**Returns the last inserted row auto increment number.*/
	public int getLastInsertId(String table){
		String SQL = "SELECT LAST_INSERT_ID() FROM "+table;
		executeSQLAdvance(SQL);	
		return getResultSetInt("last_insert_id()");
	}

	
	public ResultSet executeSQL(String sqlQuery) {
		try {
			results = statement.executeQuery(sqlQuery);
		}
		catch (SQLException e) {
			System.err.println("Problem with executing an SQL statement!\nSQL: "+sqlQuery);
			e.printStackTrace();
		}
		return results;
	}
	/**Execute and advance the ResultSet.  Good for retrieving a single table row.*/
	public ResultSet executeSQLAdvance(String sqlQuery) {
		try {
			results = statement.executeQuery(sqlQuery);
			results.next();
		}
		catch (SQLException e) {
			System.err.println("Problem with executing an SQL statement!\nSQL: "+sqlQuery);
			e.printStackTrace();
		}
		return results;
	}
	
	/**Fires SQL Update statement returning true if no problems were encountered.*/
	public boolean executeSQLUpdate(String sqlQuery) {
		try {
			statement.executeUpdate(sqlQuery);
			return true;
		}
		catch (SQLException e) {
			System.err.println("Problem with executing SQL statement!\n" + sqlQuery);
			e.printStackTrace();
			return false;
		}
	}
	
	public Connection getConnection() {
		return connection;
	}	
	

	/**Looks for a particular text in a table with a Name column.  Throws
	 * error if no text column.*/
	public boolean isNameUnique(String tableName, String name){
		if (Misc.isEmpty(name)) return false;
		String sql = "Select id FROM "+ tableName+" WHERE text = '"+name+"'";
		String result = getCell(sql);
		return result.equals(""); 
	}
	
	/**Returns the first cell from an sql statement.*/
	public String getCell(String sql) {
		executeSQL(sql);
		try {
			if (results.next()) return results.getString(1);
			else return "";
		}
		catch (SQLException e) {
			System.err.println("Problem with getCell( " + sql+" )");
			e.printStackTrace();
			return "";
		}
	}
	public ResultSet getResults() {
		return results;
	}
}	
	

