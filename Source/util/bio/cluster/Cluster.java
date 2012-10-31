package util.bio.cluster;

import java.io.File;
import javax.swing.*;
import java.awt.geom.*;
import util.gen.IO;
import util.gen.Num;
import java.io.*;

public class Cluster implements Serializable{ 
	
	private double correlationCoefficient = -2; //-2 = not set
	private Cluster parentCluster = null;
	private Cluster topCluster = null;
	private Cluster bottomCluster = null;
	private double totalClusters =1;
	private String name;
	private File floatArrayFile;
	private JLabel nameLabel;
	private Line2D line;
	private Point2D parentPoint;
	private double nameLabelCenterYCoordinate;
	
	/**Initializing Constructor*/
	public Cluster(File floatArrayFile, String name){
		this.floatArrayFile = floatArrayFile;
		this.name = name;
	}
	
	/**Secondary Constructor*/
	public Cluster(Cluster one, Cluster two, double correlationCoefficient, StringBuffer results){
		this.correlationCoefficient = correlationCoefficient;
		//set parent in children
		one.setParentCluster(this);
		two.setParentCluster(this);
		//set top and bottom
		if (one.totalClusters < two.getTotalClusters()) {
			topCluster = two;
			bottomCluster = one;
		}
		else {
			topCluster = one;
			bottomCluster = two;
		}
		//set total number of clusters
		totalClusters = one.totalClusters+two.getTotalClusters();
		//average their float[]s, save, and set the file
		merge();
		//save results summary
		String cc = Num.formatNumber((correlationCoefficient * correlationCoefficient * 100), 1);
		results.append(cc+"\t"+topCluster.getFloatArrayFile().getName()+"\t"+
				bottomCluster.getFloatArrayFile().getName()+ "\t"+ floatArrayFile.getName()+"\n");
	}
	
	public float[] getValues(){
		return (float[])IO.fetchObject(floatArrayFile);
	}

	/**Average float arrays from top and bottom clusters, same the float[], set the file*/
	private void merge(){
		float[] one = topCluster.getValues();
		float[] two = bottomCluster.getValues();
		float[] merged = Num.mean(one, two);
		File file = new File(topCluster.getFloatArrayFile().getParentFile(),topCluster.getFloatArrayFile().getName()+"_"+bottomCluster.getFloatArrayFile().getName());
		IO.saveObject(file, merged);
		floatArrayFile = file;
		one = null;
		two = null;
		merged = null;
	}

	public Cluster getBottomCluster() {
		return bottomCluster;
	}

	public void setBottomCluster(Cluster bottomCluster) {
		this.bottomCluster = bottomCluster;
	}

	public double getCorrelationCoefficient() {
		return correlationCoefficient;
	}

	public void setCorrelationCoefficient(double correlationCoefficient) {
		this.correlationCoefficient = correlationCoefficient;
	}

	public File getFloatArrayFile() {
		return floatArrayFile;
	}

	public void setFloatArrayFile(File floatArrayFile) {
		this.floatArrayFile = floatArrayFile;
	}

	public String getName() {
		return name;
	}

	public void setName(String name) {
		this.name = name;
	}

	public Cluster getParentCluster() {
		return parentCluster;
	}

	public void setParentCluster(Cluster parentCluster) {
		this.parentCluster = parentCluster;
	}

	public Cluster getTopCluster() {
		return topCluster;
	}

	public void setTopCluster(Cluster topCluster) {
		this.topCluster = topCluster;
	}

	public double getTotalClusters() {
		return totalClusters;
	}

	public JLabel getNameLabel() {
		return nameLabel;
	}

	public void setNameLabel(JLabel nameLabel) {
		this.nameLabel = nameLabel;
	}

	public void setTotalClusters(double totalClusters) {
		this.totalClusters = totalClusters;
	}

	public Line2D getLine() {
		return line;
	}

	public void setLine(Line2D line) {
		this.line = line;
	}

	public double getNameLabelCenterYCoordinate() {
		return nameLabelCenterYCoordinate;
	}

	public void setNameLabelCenterYCoordinate(double centerYCoordinate) {
		this.nameLabelCenterYCoordinate = centerYCoordinate;
	}

	public Point2D getParentPoint() {
		return parentPoint;
	}

	public void setParentPoint(Point2D parentPoint) {
		this.parentPoint = parentPoint;
	}

}