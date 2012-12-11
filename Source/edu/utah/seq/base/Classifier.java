package edu.utah.seq.base;

import org.apache.commons.math3.util.FastMath;
import util.gen.Num;

public class Classifier {
	
	//float[4][nMerLength * 4 channels]
	float[][] meansACGT;
	float[][] stdsACGT;
    private static final double SQRT2PI = FastMath.sqrt(2 * FastMath.PI);
    
	//constructor
	public Classifier (float[][] meansACGT, float[][] stdsACGT){
		this.meansACGT = meansACGT;
		this.stdsACGT = stdsACGT;
	}
	
	//methods
	public double[] calculateLikelihoods(double[] obsInts){
		
		//calculate summed log densities for each case
		double[] acgtSums = new double[meansACGT.length];
		double totalSums = 0;
		for (int i=0; i< meansACGT.length; i++){
			//for each of the nmer bases
			for (int j=0; j< obsInts.length; j++){
				double d = density(obsInts[j], meansACGT[i][j], stdsACGT[i][j]);
				//if (d <=0) System.err.println("Negative/ zero density observed, skipping "+d+" "+obsInts[j]+" "+meansACGT[i][j]+" "+stdsACGT[i][j]);
				if (d>0) acgtSums[i] += Math.log(d);
			}
			//delog it
			acgtSums[i] = Num.antiLog(acgtSums[i], Math.E);
			totalSums+= acgtSums[i];
		}

		//calculate likelihoods
		double[] actgLikelihoods = new double[meansACGT.length];
		for (int i=0; i< meansACGT.length; i++){
			actgLikelihoods[i] = acgtSums[i] / totalSums;
		}

		return actgLikelihoods;
	}
	
	public double[] calculateLikelihoodsTest(double[] obsInts){
		
		//calculate summed log densities for each case
		double[] acgtSums = new double[meansACGT.length];
		double totalSums = 0;
		for (int i=0; i< meansACGT.length; i++){
			System.out.println("\nTest case "+i);
			//for each of the nmer bases
			for (int j=0; j< obsInts.length; j++){
				double d = density(obsInts[j], meansACGT[i][j], stdsACGT[i][j]);
				System.out.println("\tDensity for "+j+"\t"+d+"  Logged density "+Math.log(d));
				if (d <=0) System.err.println("Negative/ zero density observed, skipping "+obsInts[j]+" "+meansACGT[i][j]+" "+stdsACGT[i][j]);
				else acgtSums[i] += Math.log(d);
			}
			//delog it
			acgtSums[i] = Num.antiLog(acgtSums[i], Math.E);
			System.out.println("Final delogged sum for "+i+"\t"+acgtSums[i]);
			totalSums+= acgtSums[i];
		}
		
		//calculate average of sums

		//calculate likelihoods
		System.out.println("Likelihoods for each test case");
		double[] actgLikelihoods = new double[meansACGT.length];
		for (int i=0; i< meansACGT.length; i++){
			actgLikelihoods[i] = acgtSums[i] / totalSums;
			System.out.println(i+" case  "+actgLikelihoods[i] +" = "+acgtSums[i]+" / "+totalSums);
		}

		return actgLikelihoods;
	}
	
	/**From apache.commons.math3*/
    public static double density(double x, double mean, double standardDeviation) {
        final double x0 = x - mean;
        final double x1 = x0 / standardDeviation;
        return FastMath.exp(-0.5 * x1 * x1) / (standardDeviation * SQRT2PI);
    }
    
    /*For testing.*/
    public static void main (String[] args){
    	//all are ACGT vs xxNxx
    	float[][] means =
    	{
    	{127.2825203f, 843.7784553f, 2514.060976f, 2335.296748f, 238.1382114f, 330.5914634f, 3052.806911f, 3224.786585f, 2994.587398f, 386.9207317f, 191.0792683f, 62.24390244f, 53.05894309f, 199.5447154f, 2728.028455f, 5157.414634f, 396.2276423f, 85.93699187f, 121.6707317f, 1600.107724f },
    	{102.4744608f, 734.4267877f, 839.415437f, 2257.863791f, 192.0737798f, 321.6049943f, 3007.700341f, 3069.659478f, 2972.011351f, 355.0681044f, 183.1418842f, 50.52894438f, 57.05675369f, 197.8274688f, 2681.215664f, 5058.802497f, 331.2088536f, 291.8513053f, 119.3200908f, 1530.944381f },
    	{113.4363636f, 720.9090909f, 224.5181818f, 2004.763636f, 180.5272727f, 355.7818182f, 2894.927273f, 426.2f, 2547.154545f, 332.3181818f, 196.9181818f, 203.0363636f, 2678.763636f, 299.9181818f, 2653.036364f, 4918.763636f, 342.7f, 1531.881818f, 177.6f, 1482.936364f },
    	{121.1476323f, 663.2896936f, 206.1448468f, 2258.167131f, 217.4011142f, 341.1754875f, 2870.412256f, 402.3593315f, 2888.587744f, 369.6350975f, 183.4373259f, 56.51253482f, 160.9777159f, 197.1197772f, 2717.910864f, 5256.924791f, 661.3927577f, 5331.740947f, 344.8718663f, 1675.571031f },
    	};
    	
    	float[][] stds =
    	{
    	{151.1290282f, 179.5676912f, 399.5074671f, 387.8890142f, 166.3095328f, 246.9477747f, 447.7579282f, 485.0016579f, 474.3115007f, 272.6513403f, 144.7154375f, 154.4982598f, 158.5477967f, 158.2012435f, 483.8462448f, 1007.22716f, 356.7293054f, 330.1752329f, 317.5954101f, 390.453737f },
    	{148.0493379f, 164.4975253f, 162.8982252f, 352.0819049f, 161.7689858f, 236.0648085f, 428.2122441f, 431.3686725f, 449.4663299f, 254.1513684f, 145.6924815f, 150.2229839f, 158.934017f, 151.4582592f, 448.020326f, 957.6247709f, 336.4352085f, 337.1800361f, 298.2333355f, 382.7656628f },
    	{154.6703228f, 171.0274695f, 133.5057111f, 359.3732803f, 166.8299529f, 225.1590827f, 431.2097493f, 229.8018213f, 428.754676f, 225.2686442f, 133.841213f, 124.0715928f, 462.0044864f, 138.0790315f, 445.7498358f, 983.3051542f, 307.4750844f, 369.9467785f, 314.3810053f, 349.1201596f},
    	{172.6364829f, 166.5932546f, 166.4235398f, 377.2022188f, 180.4364219f, 293.8123995f, 429.7374052f, 282.9276841f, 460.5800781f, 298.2088485f, 155.53593f, 150.8248429f, 151.8438243f, 148.5565544f, 489.9388491f, 994.6335966f, 398.605676f, 876.7700165f, 374.9295357f, 405.9969932f },
    	};
    	
    	Classifier c = new Classifier(means, stds);
    	
    	double[] obsInts = {57, 544, 859, 2344, 130, 97, 2808, 2995, 2662, 248, -160, -21, -436, 199, 2552, 4016, -97, -45, -117, 807};
     	
    	c.calculateLikelihoodsTest(obsInts);
    	
    }
}
