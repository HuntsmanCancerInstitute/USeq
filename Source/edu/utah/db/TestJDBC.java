package edu.utah.db;

import java.sql.*;

public class TestJDBC {

	//were running Microsoft SQL Server 2014 (SP2) (KB3171021) - 12.0.5000.0 (X64)
	public static void main (String[] args){
		new TestJDBC().test();
	}
	
	public void test(){

		System.out.println("Attempting Connection....");
		
		try {
			Driver d = (Driver) Class.forName("com.microsoft.sqlserver.jdbc.SQLServerDriver").newInstance();
			System.out.println("Driver "+d);
		} catch (Exception e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		String connectionUrl = "jdbc:sqlserver://localhost:1433;databaseName=gnomex;user=pipeline;password=yuk0nJ@ck";


		Connection con = null;
		Statement stmt = null;
		ResultSet rs = null;

		try {

			//establish connection
			Class.forName("com.microsoft.sqlserver.jdbc.SQLServerDriver");
			con = DriverManager.getConnection(connectionUrl);
			
			System.out.println("Attempting query....");
			
			String SQL = "select idLab,firstName,lastName from Lab where lastName='Moon'";
			stmt = con.createStatement();
			rs = stmt.executeQuery(SQL);

			while (rs.next()) {
				System.out.println(rs.getString(1) + " " +rs.getString(2));
			}

		}catch(Exception e){

			e.printStackTrace();

		}finally{

			if(rs != null) try { rs.close(); } catch(Exception e) {}
			if(stmt != null) try { stmt.close(); } catch(Exception e) {}
			if(con != null) try { con.close(); } catch(Exception e) {}

		}

	}

}



