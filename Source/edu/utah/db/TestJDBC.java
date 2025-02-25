package edu.utah.db;

import java.sql.*;
import java.util.HashSet;
import util.gen.IO;

public class TestJDBC {

	//The standard GNomEx instance is running Microsoft SQL Server, it might be 2014 (SP2) (KB3171021) - 12.0.5000.0 (X64)
	public static void main (String[] args){
		new TestJDBC().test();
	}
	
	public void test(){

		
		
		try {
			IO.pl("Instantiating a driver...");
			Driver d = (Driver) Class.forName("com.microsoft.sqlserver.jdbc.SQLServerDriver").newInstance();
			IO.pl("Driver "+d);
		} catch (Exception e1) {
			// TODO Auto-generated catch block
			IO.el("Failed driver instantiation");
			e1.printStackTrace();
		}
		String connectionUrl = "jdbc:sqlserver://hci-db.hci.utah.edu:1433;databaseName=gnomex;user=pipeline;password=xxxxxxxxx;encrypt=true;trustServerCertificate=true";
		//replace xxxxx pwd from https://ri-confluence.hci.utah.edu/pages/viewpage.action?pageId=38076459


		Connection con = null;
		Statement stmt = null;
		ResultSet rs = null;

		try {

			//establish connection
			IO.pl("Attempting to make a connection...");
			Class.forName("com.microsoft.sqlserver.jdbc.SQLServerDriver");
			con = DriverManager.getConnection(connectionUrl);
			
			IO.pl("Attempting query....");
			
			//String SQL = "select a.number, a.name from Lab JOIN analysis a on Lab.idLab=a.idLab where Lab.lastName='Snyder'";
			/*String SQL = "SELECT "+
					"request.number,  "+
					"request.name,  "+
					"request.createDate,  "+
					"project.name,  "+
					"appuser.email, "+
					"appuser.firstname,  "+
					"appuser.lastname,  "+
					"lab.firstname,  "+
					"lab.lastname,  "+
					"application.application, "+
					"request.analysisInstructions "+
					"FROM request  "+
					"join project on project.idproject = request.idproject  "+
					"join lab on lab.idlab = request.idlab  "+
					"join appuser on appuser.idappuser = request.idappuser  "+
					"join application on application.codeapplication = request.codeapplication "+
					"WHERE request.bioInformaticsAssist = 'Y' AND request.codeRequestStatus = 'COMPLETE'"+
					"ORDER BY request.createDate; ";
			int numReturnValues = 11;
			*/
			
			/*
			String SQL = "SELECT "+
					"lab.firstname,  "+
					"lab.lastname, "+
					"appuser.firstname,  "+
					"appuser.lastname,  "+
					"appuser.email, "+
					"request.number, "+
					"project.name, "+
					"organism.organism, "+
					"genomebuild.genomebuildname, "+
					"request.analysisInstructions, "+
					"seqlibprotocolapplication.codeapplication, "+
					"application.application, "+
					"request.createdate, "+
					"numbersequencingcycles.numbersequencingcycles, "+
					"seqruntype.seqruntype "+
					"FROM lab "+
					"join request on request.idlab = lab.idlab "+
					"join appuser on appuser.idappuser = request.idappuser "+
					"join project on request.idproject = project.idproject "+
					"join sample on sample.idrequest = request.idrequest "+
					"join organism on sample.idorganism = organism.idorganism "+
					"join sequencelane on sequencelane.idsample = sample.idsample "+
					"left outer join genomebuild on sequencelane.idgenomebuildalignto = genomebuild.idgenomebuild "+
					"join seqruntype on sequencelane.idseqruntype = seqruntype.idseqruntype "+
					"join numbersequencingcycles on sequencelane.idnumbersequencingcycles = numbersequencingcycles.idnumbersequencingcycles "+
					"join seqlibprotocol on sample.idseqlibprotocol = seqlibprotocol.idseqlibprotocol "+
					"join seqlibprotocolapplication on seqlibprotocol.idseqlibprotocol = seqlibprotocolapplication.idseqlibprotocol "+
					"join application on seqlibprotocolapplication.codeapplication = application.codeapplication "+
					"WHERE request.bioInformaticsAssist = 'Y' AND request.codeRequestStatus = 'COMPLETE' ";
			int numReturnValues = 15; 
			*/
			
			/*
			 * From Tim M:
			 * SELECT DISTINCT
			       request.number, 
			       request.createDate, 
			       appuser.email,
			       organism.organism,
			       application.application
			       from request
			       join project on project.idproject = request.idproject 
			       join sample on sample.idrequest = request.idrequest
			       join organism on sample.idorganism = organism.idorganism
			       join appuser on appuser.idappuser = request.idappuser 
			       join application on application.codeapplication = request.codeapplication
			       join sequencelane on sequencelane.idrequest = request.idrequest
			       WHERE request.codeRequestStatus = 'COMPLETE' and sequencelane.idGenomeBuildAlignTo is not null
			       ORDER BY request.createDate;
			 */
			
			String SQL = "SELECT DISTINCT "+
					"request.number,  "+
					"request.createDate,  "+
					"appuser.email, "+
					"lab.lastname,  "+
					"lab.firstname,  "+
					"organism.organism, "+
					"genomebuild.genomebuildname, "+
					"application.application, "+
					"request.analysisInstructions, "+
					"request.alignToGenomeBuild "+
					"FROM request  "+
					"join project on project.idproject = request.idproject  "+
					"join lab on lab.idlab = request.idlab  "+
					"join sample on sample.idrequest = request.idrequest "+
					"join organism on sample.idorganism = organism.idorganism "+
					"join sequencelane on sequencelane.idsample = sample.idsample "+
					"left outer join genomebuild on sequencelane.idgenomebuildalignto = genomebuild.idgenomebuild "+
					"join appuser on appuser.idappuser = request.idappuser  "+
					"join application on application.codeapplication = request.codeapplication "+
					"WHERE request.codeRequestStatus = 'COMPLETE' AND (request.bioInformaticsAssist = 'Y' OR sequencelane.idGenomeBuildAlignTo IS NOT NULL) "+
					"ORDER BY request.createDate; ";
			int numReturnValues = 10;
			
			//SQL = "exec sp_columns request;";
			SQL = "SELECT request.alignToGenomeBuild, request.bioInformaticsAssist, request.analysisInstructions FROM request ORDER BY request.createDate;";
			numReturnValues = 3;
			
			//alignToGenomeBuild is NA, N, or Y
			//bioInformaticsAssist is N or Y
			//analysisInstructions is NA or freeform txt
			
			SQL = "SELECT DISTINCT "+
					"request.number,  "+				//0
					"request.createDate,  "+			//1
					"appuser.email, "+					//2
					"lab.lastname,  "+					//3
					"lab.firstname,  "+					//4
					"organism.organism, "+				//5
					"application.application, "+		//7
					"request.analysisInstructions, "+	//8   NA or freeform txt
					"request.alignToGenomeBuild, "+     //9  NA, N, or Y
					"request.bioInformaticsAssist "+    //10  N or Y
					"FROM request  "+
					"join project on project.idproject = request.idproject  "+
					"join lab on lab.idlab = request.idlab  "+
					"join sample on sample.idrequest = request.idrequest "+
					"join organism on sample.idorganism = organism.idorganism "+
					"join appuser on appuser.idappuser = request.idappuser  "+
					"join application on application.codeapplication = request.codeapplication "+
					"WHERE request.createDate > (select dateadd(month, -5, getdate())) "+
					"AND (request.bioInformaticsAssist = 'Y' OR request.alignToGenomeBuild = 'Y') "+
					"ORDER BY request.createDate; ";
			numReturnValues = 10;
			
			stmt = con.createStatement();
			rs = stmt.executeQuery(SQL);
			IO.pl("Results....\n");
			int numResults = 0;
			HashSet<String> uniqueResults = new HashSet<String>();
			while (rs.next()) {
				StringBuilder sb = new StringBuilder();
				for (int i=1; i<numReturnValues+1; i++) {
					String val = rs.getString(i);
					if (val != null) sb.append(val.replaceAll("\\s+", " "));
					else sb.append("NA");
					sb.append("\t");
				}
				IO.pl(sb);
				uniqueResults.add(sb.toString());
				numResults++;
			}
			IO.pl("\nNumRes\t"+numResults);
			IO.pl("NumUni\t"+uniqueResults.size());

		}catch(Exception e){
			IO.el("Failed either connecting or attempting a query or reading the results");
			e.printStackTrace();

		}finally{

			if(rs != null) try { rs.close(); } catch(Exception e) {}
			if(stmt != null) try { stmt.close(); } catch(Exception e) {}
			if(con != null) try { con.close(); } catch(Exception e) {}

		}

	}

}



