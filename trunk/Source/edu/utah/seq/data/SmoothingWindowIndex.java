package edu.utah.seq.data;

public class SmoothingWindowIndex {
	
	private int index;
	private SmoothingWindow smoothingWindow;
	
	public SmoothingWindowIndex (int index, SmoothingWindow smoothingWindow){
		this.index = index;
		this.smoothingWindow = smoothingWindow;
	}

	public static SmoothingWindowIndex[] makeSmoothingWindowIndex(SmoothingWindow[] smoothingWindow){
		SmoothingWindowIndex[] swi = new SmoothingWindowIndex[smoothingWindow.length];
		for (int i=0; i< swi.length; i++) swi[i] = new SmoothingWindowIndex(i, smoothingWindow[i]);
		return swi;
	}
	
	public int getIndex() {
		return index;
	}

	public void setIndex(int index) {
		this.index = index;
	}

	public SmoothingWindow getSmoothingWindow() {
		return smoothingWindow;
	}

	public void setSmoothingWindow(SmoothingWindow smoothingWindow) {
		this.smoothingWindow = smoothingWindow;
	}
}
