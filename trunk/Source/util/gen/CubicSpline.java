package util.gen;

/**For calculating cubic splines of incremental data, interpolation, from http://www.frank-buss.de/CubicSpline.java */
public class CubicSpline {
	private double y[];
	private double y2[];
	
	/**
	 * The constructor calculates the second derivatives of the interpolating function
	 * at the tabulated points xi, with xi = (i, y[i]).
	 * Based on numerical recipes in C, http://www.library.cornell.edu/nr/bookcpdf/c3-3.pdf .
	 * @param y Array of y coordinates for cubic-spline interpolation.
	 */
	public CubicSpline(double y[]) {
		this.y = y;
		int n = y.length;
		y2 = new double[n];
		double u[] = new double[n];
		for (int i = 1; i < n - 1; i++) {
			y2[i] = -1.0 / (4.0 + y2[i - 1]);
			u[i] = (6.0 * (y[i + 1] - 2.0 * y[i] + y[i - 1]) - u[i - 1]) / (4.0 + y2[i - 1]);
		}
		for (int i = n - 2; i >= 0; i--) {
			y2[i] = y2[i] * y2[i + 1] + u[i];
		}
	}

	/**
	 * Returns a cubic-spline interpolated value y for the point between
	 * point (n, y[n]) and (n+1, y[n+1), with t ranging from 0 for (n, y[n])
	 * to 1 for (n+1, y[n+1]).  
	 * @param startIndex The start index/ point.
	 * @param dist2NextPoint The distance to the next point (0..1).
	 * @return A cubic-spline interpolated y value.
	 */
	public double interpolate(int startIndex, double dist2NextPoint) {
		return dist2NextPoint * y[startIndex + 1] - ((dist2NextPoint - 1.0) * dist2NextPoint * ((dist2NextPoint - 2.0) * y2[startIndex] - (dist2NextPoint + 1.0) * y2[startIndex + 1])) / 6.0 + y[startIndex] - dist2NextPoint * y[startIndex];
	}
	
	public static void main (String[] args){
		double[] d = {30,60,90};//{18,38,77,172,310,498,676,945,1061,1183,1249,1064,897,690,460,301,170,110,33,28,7};
		
		CubicSpline s = new CubicSpline(d);
		System.out.println(s.interpolate(1,0.5));

	}
}

