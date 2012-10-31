package trans.cel;
import util.gen.*;
import java.util.*;

/**Generates multiplicative scalars based on a set of points using a combination of a cubic spline for values within the provided points 
 * and linear regression for points outside the value range. Lagrange Interpolation might be something possibly more appropriate.
 * Assumes that each sucessive point is >= than last.*/
public class SplineScalar {
	
	//fields
	private double[] modelPoints;
	private int numModPtsMinusOne;
	private CubicSpline spline;
	private LinearRegression leftLR;
	private LinearRegression rightLR;
	private double backOff = 0.01;	//fraction to back into spline to get point for linear regression
	
	//constructor
	/** Assumes modelPoints are sorted min to max.
	 * @param modelPoints - the points that represent the average/ model you want to generate scalars for.*/
	public SplineScalar (double[] modelPoints) {
		this.modelPoints = modelPoints;
		//make CubicSpline obj for covered values
		spline = new CubicSpline(modelPoints);
		//make LinearRegression objects for left and right sides
		double[] x = {0,1};
		//left
		double[] ls = {modelPoints[0], spline.interpolate(0,backOff)};
		leftLR = new LinearRegression (x, ls);
		//right
		double fracDistToEnd = 1-backOff;
		numModPtsMinusOne = modelPoints.length -1;
		int nextTwoLastIndex = numModPtsMinusOne -1;
		double[] rs = {spline.interpolate(nextTwoLastIndex, fracDistToEnd), modelPoints[0]};
		rightLR = new LinearRegression (x, rs);
	}
	
	/**Use to calculate a scaler that when multiplied by the testValue brings the testValue up to the model.
	 * @param testValue a non zero double
	 * @return scaler x testValue = model.*/
	public double scalar(double testValue, SplineScalar testScalar){
		double[] testModelPoints = testScalar.getModelPoints();
		double scalar = 0;
		//is it inside?
		if (testValue >= testModelPoints[0] && testValue <= testModelPoints[numModPtsMinusOne]){
			//System.out.println("Inside");
			//start index
			int startIndex = Arrays.binarySearch(testModelPoints,testValue);
			//interpolation?
			if (startIndex < 0) {
				startIndex = -1 * (startIndex+2);
				//System.out.println("Interpolating\nstartIndex "+startIndex);
				//distance
				double range = testModelPoints[startIndex+1] - testModelPoints[startIndex];
				double diff = testValue - testModelPoints[startIndex];
				double fractionDistance = diff/ range;
				//System.out.println("fractionDistance "+fractionDistance);
				//find modelValue from this CubicSpline using startIndex and distance
				double modelValue = spline.interpolate(startIndex, fractionDistance);
				//System.out.println("modelValue "+modelValue);
				//calculate scalar
				scalar = modelValue/testValue;
			}
			//read directly
			else {
				//System.out.println("Reading Directly");
				scalar = modelPoints[startIndex]/ testModelPoints[startIndex];
			}
		}
		//is it left side
		else if (testValue < testModelPoints[0]){
			//System.out.println("Left Side");
			//estimate distance from testCelScalar's left side LinearRegression
			double xTest = testScalar.getLeftLR().calculateX(testValue);
			double modelValue = leftLR.calculateY(xTest);
			scalar = modelValue/testValue;
			
		}
		//right side by default
		else {
			//System.out.println("Right Side");			
			//estimate distance from testCelScalar's right side LinearRegression
			double xTest = testScalar.getRightLR().calculateX(testValue);
			double modelValue = rightLR.calculateY(xTest);
			scalar = modelValue/testValue;
		}
		return scalar;
	}
	
	
	
	public static void main (String[] test){
		//10th percentile, neg controls, gapdh, actinB, 90th percentile
		double[] modelValues = new double[]{2,3,4.5,5.5,5.75};
		SplineScalar s = new SplineScalar(modelValues);

		double[] testValues = new double[]{2.1,3.1,5,6,7};
		SplineScalar testSC = new SplineScalar(testValues);
		
		double[] toTest = {1.9, 2.2, 3.05, 5.1, 5.9,6.9, 7.2, 7.4, 7.6};
		for (int i=0; i<toTest.length; i++){
			double multiplier = s.scalar(toTest[i], testSC);
			System.out.print("\t"+multiplier*toTest[i]);
		}
	}

	public LinearRegression getLeftLR() {
		return leftLR;
	}

	public double[] getModelPoints() {
		return modelPoints;
	}

	public LinearRegression getRightLR() {
		return rightLR;
	}

	public CubicSpline getSpline() {
		return spline;
	}
	
	
	
}
