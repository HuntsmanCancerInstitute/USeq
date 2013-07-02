package edu.utah.seq.data;

import java.util.ArrayList;
import util.gen.Num;

/**Container for paired beta values 0-1 off Illumina methylation data.  Zero values are excluded.
 * @author davidnix*/
public class MethylationArraySamplePair{

	private float[] treatment;
	private float[] control;
	private boolean performPairedAnalysis;

	public MethylationArraySamplePair(float[] treatment, float[] control, boolean performPairedAnalysis) {
		this.treatment = treatment;
		this.control = control;
		this.performPairedAnalysis = performPairedAnalysis;
	}

	public void randomize() {
		Num.randomizePairedValues(treatment, control, System.currentTimeMillis());
	}

	/**StopIndex is included!*/
	public void fetchScoresByIndex(int startIndex, int stopIndex, ArrayList<Float> treatmentAL, ArrayList<Float> controlAL){
		if (performPairedAnalysis){
			for (int i=startIndex; i<= stopIndex; i++){
				if (treatment[i] !=0.0f && control[i] != 0.0f){
					treatmentAL.add(treatment[i]);
					controlAL.add(control[i]);
				}
			}
		}
		else{
			if (treatment !=null){
				for (int i=startIndex; i<= stopIndex; i++){
					if (treatment[i] !=0.0f ) treatmentAL.add(treatment[i]);
				}
			}
			if (control !=null){
				for (int i=startIndex; i<= stopIndex; i++){
					if (control[i] != 0.0f) controlAL.add(control[i]);
				}
			}
		}
	}
}