package trans.main;
import java.util.*;

/**
 * Finds peaks in {@link Interval}s, both flat topped or peaks with sloping sides.
 * */

public class PeakPicker {
	//slope picking
	private int minRun = 5;			//minimum number of oligos to include in look ahead window
	private double minPercent = .51;	//minimum fraction in look ahead window that must be positive or negative to keep growing, will attempt to continue growing one oligo at a time if look a head fails.
	private int minBpSides = 60;		//minimum bp width to call a peak side, total width must be 2* to be called a peak
	private double scoreCutOff = 1.3;	//only consider peaks >= to this cut off, this is set internally when looking for flattop peaks
	private int maxNumPeaks =4 ;		//maximum number of peaks to save in each Interval, they will be ranked by score.
	//flat top picking
	private int maxBpGap = 75;		//max bp gap between oligos
	private double maxScoreFractionDrop = 0.80;	//expanded score cannot be less than this * original sub window score
	private boolean makeFlattops = false;
	//internal do not touch
	private Oligo[] oligos;
	private Interval interval;
	private int numOligos;
	private int peakIndex;
	private int leftIndex;
	private int rightIndex;
	private int peakBp;
	
	//methods
	
	/**Attempts to expand the leftIndex and rightIndex*/
	public boolean expandFlatTopPeak (){
		boolean expanded = false;
		//look left
		while (true){
			
			int testIndex = leftIndex -1;
			System.out.println("LookingLeft "+testIndex+" right "+rightIndex);
			//is there a left score
			if (testIndex < 0) break;
			//check gap
			if (oligos[leftIndex].getStart() - oligos[testIndex].getStart() > maxBpGap) break;
			//check score
			double testScore = averageWindowScores(testIndex, rightIndex);
			System.out.println("\tTestScore "+testScore);
			if (testScore < scoreCutOff) break;
			//all ok so expand
			leftIndex--;
			expanded = true;
		}
		System.out.println("LookingRight");
		//look right
		while (true){
			int testIndex = rightIndex +1;
			//is there a left score
			if (testIndex >= oligos.length) break;
			//check gap
			if (oligos[testIndex].getStart() - oligos[rightIndex].getStart() > maxBpGap) break;
			//check score
			double testScore = averageWindowScores(leftIndex, testIndex);
			if (testScore < scoreCutOff) break;
			//all ok so expand
			rightIndex++;
			expanded = true;
		}
		return expanded;
	}
	
	/**Averages the smoothed oligo scores, stop is inclusive.*/
	public double averageWindowScores(int start, int stop){
		double length = stop - start + 1;
		double total = 0;
		for (int i=start; i<=stop; i++){
			total += oligos[i].getSmoothedScore();
		}
		return total/length;
	}
	
	/**Averages the smoothed scores in the best sub window. and sets this * maxScoreFractionDrop as the scoreCutOff.*/
	public boolean setSubWindowScoreAsMinScore(){
		SubWindow sw = interval.getBestSubWindow();
		if (sw == null) return false;
		Oligo[] swOligos = sw.getOligos();
		double length = swOligos.length;
		double total = 0;
		for (int x = 0; x < length; x++) total += swOligos[x].getSmoothedScore();
		double score = total/length;
		scoreCutOff = score * maxScoreFractionDrop;
System.out.println("Score "+score +" minScore "+scoreCutOff);		
		return true;
	}

	public void pickPeaks(Interval[] intervals){	
		int numIntervals = intervals.length;
		for (int x=0; x<numIntervals; x++){
			interval = intervals[x];
			oligos = interval.getOligos();
			numOligos = oligos.length;
			//Slope oligos
			slopeOligos();
			//find peaks
			ArrayList bindingPeaksAL = new ArrayList();
			int numBPs = 0;
			boolean go = true;
			while (go){
				int flag = -1;
				//find a flat top?
				if (makeFlattops){						
					//set score cutOff
					if (setSubWindowScoreAsMinScore() == false) break;
					flag = findFlatTopPeak();
				}
				//find a peak!
				else flag = findPeak();

				//check flag and if good make a peak
				if ( flag == -1) go = false;
				else if (flag == 1){
					//make and add binding region
					BindingPeak bp = new BindingPeak(peakIndex, leftIndex, rightIndex, peakBp, oligos[peakIndex].getSmoothedScore());
					bindingPeaksAL.add(bp);
					numBPs ++;
					if (numBPs == maxNumPeaks) go = false;
				}
			}
			//found any peaks, if so set in interval
			numBPs = bindingPeaksAL.size();			
			if (numBPs !=0 ){
				BindingPeak[] bps = new BindingPeak[bindingPeaksAL.size()];
				bindingPeaksAL.toArray(bps);
				interval.setBindingPeaks(bps);
				/*System.out.println("Score\tLeft\tPeak\tRight");
				for (int i=0; i<bps.length; i++){
					System.out.println(bps[i].toString(oligos, interval.getSizeOfOligo()));
				}*/
			}
			else interval.setBindingPeaks(null);
		}
	}
	
	public void mask(int start, int stop){
		int realStop = stop+1;
		for (int i=start; i<realStop; i++){
			oligos[i].setInPeak(true);
		}
	}
	
	public int findFlatTopPeak(){
		//find index of the highest availible oligo smoothedScore
		peakIndex = findHighestScore();  //returns the first highest score might be more of the same to the right
		//will return -1 if score is less than minimum
		if (peakIndex == -1) return -1;
		peakBp = findPeakBp(peakIndex);
		//expand
		leftIndex = peakIndex;
		rightIndex = peakIndex;
		if (expandFlatTopPeak() == false){
			oligos[peakIndex].setSkip(true);
			return 0;
		}
		//mask
		mask(leftIndex,rightIndex);
		return 1;
		
	}
	
	public int findPeak(){
		//find index of the highest availible oligo smoothedScore
		peakIndex = findHighestScore();  //returns the first highest score might be more of the same to the right
		if (peakIndex == -1) return -1;
		peakBp = findPeakBp(peakIndex);
		//growleft
		leftIndex = growLeft(peakIndex);
		//System.out.println("\tleftIndex "+leftIndex);
		if (leftIndex == -1) {
			oligos[peakIndex].setSkip(true);
			return 0;
		}
		//growright
		rightIndex = growRight(peakIndex);
		//System.out.println("\trightIndex "+rightIndex);
		if (rightIndex == -1) {
			oligos[peakIndex].setSkip(true);
			return 0;
		}
		//check width and mask if peak is OK
		if (checkWidth(leftIndex, peakIndex, rightIndex)) {
			mask(leftIndex,rightIndex);
			return 1;
		}
		//otherwise set high oligo to skip
		oligos[peakIndex].setSkip(true);
		return 0;
	}
	
	public boolean checkWidth(int leftIndex, int peakIndex, int rightIndex){
		int peakBp = oligos[peakIndex].getStart();
		int leftBp= oligos[leftIndex].getStart();
		if ((peakBp-leftBp) < minBpSides) return false;
		int rightBp = oligos[rightIndex].getStart();
		if ((rightBp-peakBp) < minBpSides) return false;
		return true;
	}
	
	public int findPeakBp(int peakIndex){
		//look right
		int right= peakIndex;
		int test = peakIndex+1;
		//look right
		while (test<numOligos){
			if (oligos[test].getSmoothedScore() == oligos[peakIndex].getSmoothedScore()) {
				right = test;
				test++;
			}
			else break;
		}
		double startBp = oligos[peakIndex].getStart();
		double stopBp = oligos[right].getStart() + interval.getSizeOfOligoMinusOne();
		return (int)Math.round((stopBp - startBp)/2) + (int)startBp;
	}
	
	/**Uses the idea of a local window, counts the number of neg slopes within the window, if > minPercent, keeps advancing.*/
	public int growRight(int index){
		int rightIndex = -1;
		//run thru oligos or until a break is called or a inPeak Oligo is encountered
		for (int i=index; i<numOligos; i++){
			//are there minRun oligos to the right?
			int rightStop = i+minRun -1;
			//System.out.print("\t"+i+" :i then Right Stop: "+rightStop);
			if ( rightStop < numOligos){
				//are their combine slope % > minPercent
				//calc % watching for in Peak oligos
				double numNeg = 0;
				rightStop++;
				for (int j=rightStop-minRun; j<rightStop; j++){
					//System.out.print(" "+j);
					if (oligos[j].isInPeak()) {
						break;
					}
					if (oligos[j].getRightSlope()<=0) numNeg++; //include zeros
				}
				//System.out.println("  slope: "+(numNeg/minRun) +" "+numNeg+" "+minRun);
				//check %
				if ( (numNeg/minRun) < minPercent ) {
					break;
				}
				else rightIndex = i;
			}
			//can't advance, no more windows of oligos
			else {
				break;
			}
		}
		//System.out.println("\nnoMoreOligos "+noMoreOligos+" inPeakFound "+inPeakFound+" badSlope "+badSlope);
		//System.out.println("rI prior to grow and trim: "+rightIndex);
		if (rightIndex != -1){
			//attempt single oligo advance til pos slope encountered
			for (int i=rightIndex+1; i<numOligos; i++){
				if (oligos[i].isInPeak()==false && oligos[i].getLeftSlope()<=0) rightIndex++;
				else break;
			}
			//System.out.println("after single oligo grow: "+rightIndex);
			//trim back til neg slope encountered, remove any zeros
			for (int i=rightIndex; i>index; i--){
				if (oligos[i].getLeftSlope()>=0) rightIndex--;
				else break;
			}	
		}
		return rightIndex;
	}
	
	/**Uses the idea of a local window, counts the number of pos slopes within the window, if > minPercent, keeps advancing.*/
	public int growLeft(int index){
		int leftIndex = -1;
		//run thru oligos or until a break is called or a inPeak Oligo is encountered
		for (int i=index; i>=0; i--){
			//are there minRun oligos to the right?
			int leftStop = i-minRun+1;
			//System.out.print("\t"+i+" :i then Left Stop: "+leftStop);
			if ( leftStop >= 0){
				//are their combine slope % > minPercent
				//calc % watching for in Peak oligos
				double numPos = 0;
				leftStop--;
				for (int j=i; j>leftStop; j--){
					//System.out.print(" "+j);
					if (oligos[j].isInPeak()) {
						break;
					}
					if (oligos[j].getLeftSlope()>=0) numPos++; //include zeros
				}
				//System.out.println("  slope: "+(numPos/minRun) +" "+numPos+" "+minRun);
				
				//check %
				if ( (numPos/minRun) < minPercent ) {
					break;
				}
				else leftIndex = i;
			}
			//can't advance, no more windows of oligos
			else {
				break;
			}
		}
		//System.out.println("\nnoMoreOligos "+noMoreOligos+" inPeakFound "+inPeakFound+" badSlope "+badSlope);
		//System.out.println("rI prior to grow and trim: "+leftIndex);
		if (leftIndex != -1){
			//attempt single oligo advance til neg slope encountered
			for (int i=leftIndex-1; i>=0; i--){
				if (oligos[i].isInPeak()==false && oligos[i].getRightSlope()>=0) leftIndex--;
				else break;
			}		
			//System.out.println("after single oligo grow: "+leftIndex);
			
			//trim back til 1st good pos slope encountered to remove any zeros
			for (int i=leftIndex; i<index; i++){
				if (oligos[i].getRightSlope()<=0) leftIndex++;
				else break;
			}	
		}
		return leftIndex;
	}
	
	public int findHighestScore(){
		int maxIndex = -1;
		double maxScore = 0;		
		for (int i=0; i<numOligos; i++){			
			if (oligos[i].isInPeak()==false && oligos[i].skip()==false && oligos[i].getSmoothedScore()>maxScore){
				maxScore = oligos[i].getSmoothedScore();
				maxIndex = i;
			}
		}
		if (maxScore >= scoreCutOff) return maxIndex;
		return -1;
	}
	
	public void slopeOligos(){
		//calc median Interval to use to adjust cut off
		//reset skip and inPeak in case these intervals are being rerun
		double[] smoothed = new double[numOligos];
		for (int i=0; i<numOligos; i++){
			smoothed[i] = oligos[i].getSmoothedScore();
			oligos[i].setInPeak(false);
			oligos[i].setSkip(false);
		}
		Arrays.sort(smoothed);
		double median;
		if (interval.getBestSubWindow()!=null) median = interval.getBestSubWindow().getMedianRatio();
		else median = 3;
		double cutOff;
		if (median>5) cutOff = 0.25;
		else if (median>4) cutOff = 0.2;
		else if (median>3.5) cutOff = 0.175;
		else if (median>3) cutOff = 0.15;
		else if (median>2.75) cutOff = 0.125;
		else if (median>2.5) cutOff = 0.1;
		else if (median>2.25) cutOff = 0.075;
		else if (median>2.0) cutOff = 0.05;
		else if (median>1.75) cutOff = 0.025;
		else if (median>1.5) cutOff = 0.02;
		else cutOff = 0.015;
		
		for (int i=0; i<numOligos; i++){
			//is there another to the right?
			if ((i+1)<numOligos){
				double slope;
				double deltaY = oligos[i+1].getSmoothedScore()-oligos[i].getSmoothedScore();
				if (Math.abs(deltaY)< cutOff) slope = 0;
				else {
					double deltaX = oligos[i+1].getStart()-oligos[i].getStart();
					slope = deltaY/deltaX;
				}
				oligos[i].setRightSlope(slope);
				oligos[i+1].setLeftSlope(slope);
			}
		}
	}
	public Interval getInterval() {
		return interval;
	}
	public void setInterval(Interval interval) {
		this.interval = interval;
	}
	public int getLeftIndex() {
		return leftIndex;
	}
	public void setLeftIndex(int leftIndex) {
		this.leftIndex = leftIndex;
	}
	public int getMinBpSides() {
		return minBpSides;
	}
	public void setMinBpSides(int minBpSides) {
		this.minBpSides = minBpSides;
	}
	public double getMinPercent() {
		return minPercent;
	}
	public void setMinPercent(double minPercent) {
		this.minPercent = minPercent;
	}
	public int getMinRun() {
		return minRun;
	}
	public void setMinRun(int minRun) {
		this.minRun = minRun;
	}
	public int getNumOligos() {
		return numOligos;
	}
	public void setNumOligos(int numOligos) {
		this.numOligos = numOligos;
	}
	public Oligo[] getOligos() {
		return oligos;
	}
	public void setOligos(Oligo[] oligos) {
		this.oligos = oligos;
	}
	public int getPeakBp() {
		return peakBp;
	}
	public void setPeakBp(int peakBp) {
		this.peakBp = peakBp;
	}
	public int getPeakIndex() {
		return peakIndex;
	}
	public void setPeakIndex(int peakIndex) {
		this.peakIndex = peakIndex;
	}
	public int getRightIndex() {
		return rightIndex;
	}
	public void setRightIndex(int rightIndex) {
		this.rightIndex = rightIndex;
	}
	public double getScoreCutOff() {
		return scoreCutOff;
	}
	public void setScoreCutOff(double scoreCutOff) {
		this.scoreCutOff = scoreCutOff;
	}
	public int getMaxNumPeaks() {
		return maxNumPeaks;
	}
	public void setMaxNumPeaks(int maxNumPeaks) {
		this.maxNumPeaks = maxNumPeaks;
	}

	public void setMakeFlattops(boolean makeFlattops) {
		this.makeFlattops = makeFlattops;
	}

	public void setMaxBpGap(int maxBpGap) {
		this.maxBpGap = maxBpGap;
	}

	public void setMaxScoreFractionDrop(double maxScoreFractionDrop) {
		this.maxScoreFractionDrop = maxScoreFractionDrop;
	}
}
