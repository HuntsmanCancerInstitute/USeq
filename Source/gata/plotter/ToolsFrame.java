package gata.plotter;

import gata.aligner.*;
import gata.main.*;

import java.awt.*;
import javax.swing.*;
import java.awt.event.*;
import java.text.*;

/**
 * 
 * @author nix
 * Frame to hold the tools for manipulating the alignment image, sliders, zoom buttons, position reporters, etc.
 */
public class ToolsFrame extends JFrame {
	
	private AlignParams AP;
	private GATAParams params;
	private AlignPanel alignPanel;
	private JPanel panel = new JPanel();
	private JLabel ntRefSeq;
	private JLabel ntCompSeq;
	private JTextField minBit;
	private JTextField minE;
	private JTextField maxBit;
	private JTextField maxE;
	private double lambdaNats;
	private double K;
	private double effectiveN;
	private double effectiveM;
	private double magicNum;
	private DecimalFormat decimalFormat = new DecimalFormat("0E0");
	private NumberFormat numberFormat = NumberFormat.getNumberInstance();
	
	
	public ToolsFrame (int locX, int locY, GATAParams params){
		this.params = params;
		AP = params.getAlignParams();
		alignPanel= params.getAlignPanel();
		//set blast constants
		lambdaNats = AP.getLAMBDA();
		K = AP.getK();
		effectiveN = AP.getEFF_N();
		effectiveM = AP.getEFF_M();
		magicNum = K * effectiveN * effectiveM;
		numberFormat.setMaximumFractionDigits(1);
		//modify JFrame
		setResizable(true);
		setTitle("Tools");
		setLocation(locX, locY);
		
		//pass ToolsFrame ref to plotter params
		params.setToolsFrameRef(this);
		
		//make widgets
		Sliders two = new Sliders(AP, params);
		JSlider[] slider = two.getSliders();
		JLabel min = new JLabel("min");
		JLabel max = new JLabel("max");
		JButton zoomIn = makeButton("Zoom In");
		zoomIn = addAction(zoomIn, 25);
		JButton zoomOut = makeButton("Zoom Out");
		zoomOut = addAction (zoomOut, -25);
		JButton noZoom = makeButton("No Zoom");
		noZoom = addAction(noZoom, 0);
		
		ntRefSeq = makeLabel("Ref:  bp        ");
		ntCompSeq = makeLabel("Cmp:  bp        ");
		
		//set blank panel, needed to color the JFrame white
		panel.setBackground(Color.WHITE);
		panel.setLayout(null);
		
		//add widgets
		panel.add(slider[0]);
		panel.add(slider[1]);
		panel.add(min);
		panel.add(max);
		panel.add(ntRefSeq);
		
		//my very own layout manager!
		//set slider 1
		Dimension dimSliderMin = slider[0].getPreferredSize();
		slider[0].setBounds(0,0,(int)dimSliderMin.getWidth(), (int)dimSliderMin.getHeight());
		
		//set slider 2
		Dimension dimSliderMax = slider[0].getPreferredSize();
		slider[1].setBounds((int)dimSliderMin.getWidth(),0,(int)dimSliderMax.getWidth(),(int)dimSliderMax.getHeight());
		
		//set ZoomIn and ZoomOut
		Dimension dimZoomOut = zoomOut.getPreferredSize();
		Dimension dimZoomIn = zoomIn.getPreferredSize();
		int widthOfSliders = (int)(dimSliderMin.getWidth()*2);
		int vertDist =8;
		zoomIn.setBounds((int)(10+(dimZoomOut.getWidth()-dimZoomIn.getWidth())/2)+widthOfSliders, vertDist, (int)dimZoomIn.getWidth(),(int)dimZoomIn.getHeight());
		vertDist += (int)dimZoomIn.getHeight();
		zoomOut.setBounds(10+widthOfSliders,vertDist,(int)dimZoomOut.getWidth(),(int)dimZoomOut.getHeight());
		
		//set NoZoom
		Dimension dimNoZoom = noZoom.getPreferredSize();
		vertDist += (int)dimZoomOut.getHeight();
		noZoom.setBounds((int)(10+(dimZoomOut.getWidth()-dimNoZoom.getWidth())/2)+widthOfSliders, vertDist, (int)dimNoZoom.getWidth(),(int)dimNoZoom.getHeight());
		
		//set text fields for locations
		Dimension dimNtRefSeq = ntRefSeq.getPreferredSize();
		vertDist+= (int)dimNoZoom.getHeight()+10;
		ntRefSeq.setBounds(widthOfSliders, vertDist, (int)dimNtRefSeq.getWidth()+40, (int)dimNtRefSeq.getHeight());
		
		Dimension dimNtCompSeq = ntCompSeq.getPreferredSize();
		vertDist += (int)dimNtRefSeq.getHeight();
		ntCompSeq.setBounds(widthOfSliders, vertDist, (int)dimNtCompSeq.getWidth()+40, (int)dimNtCompSeq.getHeight());
		
		//create score boxes for min and max
		Font smallFont = new Font("Dialog",1,10);
		int x = widthOfSliders;
		int y = vertDist + (int)dimNtCompSeq.getHeight()+10;
		int xAdj = x;
		JLabel maxScore = makeLabel ("Max:");
		Dimension maxDim = maxScore.getPreferredSize();
		JLabel minScore = makeLabel ("Min:");
		minBit = new JTextField();
		minBit.setFont(smallFont);
		minE = new JTextField();
		minE.setFont(smallFont);
		int[] minScoreSize = GATAUtil.setLabel(minScore,x,y);
		xAdj+=(int)maxDim.getWidth();
		int[] minBitSize = GATAUtil.setField(minBit,xAdj, y, 40);
		xAdj+=minBitSize[0]+5;
		int[] minESize = GATAUtil.setField(minE,xAdj,y,40);
		xAdj+=minESize[0]+10;
		panel.add(minBit);
		panel.add(minE);
		
		y+=minBitSize[1]+5;
		xAdj = x;
		maxBit = new JTextField();
		maxBit.setFont(smallFont);
		maxE = new JTextField();
		maxE.setFont(smallFont);
		int[] maxScoreSize = GATAUtil.setLabel(maxScore,x,y);
		xAdj+=maxScoreSize[0];
		int[] maxBitSize = GATAUtil.setField(maxBit,xAdj, y, 40);
		xAdj+=maxBitSize[0]+5;
		int[] maxESize = GATAUtil.setField(maxE,xAdj,y,40);
		xAdj+=maxESize[0]+10;
		panel.add(maxBit);
		panel.add(maxE);
		y+=maxBitSize[1];
		
		JLabel bits = makeLabel("(bits)");
		JLabel expect = makeLabel("(Expect)");
		Dimension bitsDim = bits.getPreferredSize();
		Dimension expectDim = expect.getPreferredSize();
		x = widthOfSliders + (int)maxDim.getWidth(); 
		int sizeTextFields = maxBitSize[0];
		x+=(sizeTextFields-bitsDim.getWidth())/2;
		GATAUtil.setLabel(bits,x,y);
		x = widthOfSliders + (int)maxDim.getWidth() +sizeTextFields + 5 + (int)(sizeTextFields-expectDim.getWidth())/2;
		
		GATAUtil.setLabel(expect,x,y);
		
		//set min and max labels under sliders
		Dimension dimMin = min.getPreferredSize();
		min.setBounds((int)((dimSliderMin.getWidth()-dimMin.getWidth())/2), (int)dimSliderMin.getHeight(), (int)dimMin.getWidth(),(int)dimMin.getHeight());
		//set max
		Dimension dimMax = max.getPreferredSize();
		max.setBounds((int)((dimSliderMax.getWidth()-dimMax.getWidth())/2)+(int)dimSliderMin.getWidth(), (int)dimSliderMax.getHeight(),(int)dimMax.getWidth(), (int)dimMax.getHeight());
		
		//calculate smallest dimensions for Frame
		int smallestWidth = xAdj;//10+ 10 + widthOfSliders + (int)dimZoomOut.getWidth();
		int smallestHeight = 30+(int)dimSliderMin.getHeight() + (int)dimMin.getHeight();
		setSize(smallestWidth, smallestHeight);
		
		//set initial scores
		setScores(AP.getMIN_SCORE(),true);
		setScores(AP.getMATCH() * AP.getWIN_SIZE(),false);
		
		//add panel to JFrame and make visible
		Container contentPane = getContentPane();
		contentPane.add(panel);
		setVisible(true);
	}
	public void setScores(double rawS, boolean minOrMax){
		if (minOrMax){
			minBit.setText(convertRawScoreToBit(rawS));
			minE.setText(calculateEValue(rawS));
			return;
		}
		maxBit.setText(convertRawScoreToBit(rawS));
		maxE.setText(calculateEValue(rawS));
	}
	/**Converts a raw S score to a bit score given lambda (in nats) and K*/
	public String convertRawScoreToBit(
			double rawS) {
		double bs = (lambdaNats * rawS - Math.log(K)) / Math.log(2);
		return numberFormat.format(bs);
	}
	public String calculateEValue(double rawS) {
		double num = magicNum
		* Math.exp(-1 * lambdaNats * (double) rawS);
		String numS = decimalFormat.format(num);
		if (numS.equals("0E0"))
			return "<1E-99";
		return numS;
		
	}		
	public JLabel makeLabel (String text){
		JLabel lab = new JLabel(text);
		lab.setFont(new Font("Dialog",1,10));
		panel.add(lab);
		return lab;
	}
	
	public void setNtRefSeq(int num){ntRefSeq.setText("Ref: "+Integer.toString(num)+" bp");}
	public void setNtCompSeq(int num){ntCompSeq.setText("Cmp: "+Integer.toString(num)+" bp");}
	
	public JButton makeButton (String text){
		JButton but = new JButton(text);
		but.setFont(new Font("Dialog", 1, 10));
		but.setBackground(Color.WHITE);
		panel.add(but);
		return but;
	}
	
	public JButton addAction (JButton but, int scalar){
		ButtonZoom zoom = new ButtonZoom(scalar);
		but.addActionListener(zoom);
		return but;
	}
	
	private class ButtonZoom implements ActionListener {
		private int zoomMultiplier=1;
		public ButtonZoom (int num){
			zoomMultiplier = num;
		}
		public void actionPerformed(ActionEvent event){
			alignPanel.zoom(zoomMultiplier);
		} 
	}
	public JTextField getMaxBit() {
		return maxBit;
	}
	public JTextField getMaxE() {
		return maxE;
	}
	public JTextField getMinBit() {
		return minBit;
	}
	public JTextField getMinE() {
		return minE;
	}
	
}

