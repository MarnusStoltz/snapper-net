package snapper;

import java.util.Arrays;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.transform.DctNormalization;
import org.apache.commons.math3.transform.FastCosineTransformer;
import org.apache.commons.math3.transform.TransformType;
import org.apache.commons.math3.transform.FastFourierTransformer;
import org.apache.commons.math3.transform.DftNormalization;


public class ChebyshevPolynomial {	
	static FastCosineTransformer transformer = null;
	static FastFourierTransformer FFTtransformer = null;
	
	public double [] f; // function values
	public double [] a; // Chebyshev polynomial coefficients

	public ChebyshevPolynomial(int N) {
		f = new double[N];
		a = new double[N];
		if (transformer == null) {
			transformer = new org.apache.commons.math3.transform.FastCosineTransformer(DctNormalization.STANDARD_DCT_I);
		}
		if (FFTtransformer == null) {
			FFTtransformer = new org.apache.commons.math3.transform.FastFourierTransformer(DftNormalization.STANDARD);
		}
	}
	
	/**
	 * initialises f and a to approximate a binomial(r,n) distribution 
	 * @param r number of red lineages
	 * @param n total number of lineages
	 */
	public void init(int r, int n) {
		double c = logBinom(r, n);
		int N = f.length;
//		double [] f = new double[n]; 
		int nf = f.length;


		
		if (n == 0) { 
			// deal with missing data
			Arrays.fill(a, 0);
			a[0] = 1;
			aToF();
		} else {
			// there is some data
	    	// first, set f[0] and f[N-1];
	    	if (r == 0) {
	    		f[0] = 1;
	    		f[nf-1] = 0;            		
	    	} else if (r == n) {
	    		f[0] = 0;
	    		f[nf-1] = 1;
	    	} else {
	    		f[0] = 0;
	    		f[nf-1] = 0;
	    	}
	    		
	    	// set the other values
	    	for (int m = 1; m < nf - 1; m++) {
	    		double x = 0.5 - Math.cos(-m/(nf-1.0)*Math.PI) / 2.0;
	    		double logp = c + r * Math.log(x) + (n-r) * Math.log(1-x);
	    		double p = Math.exp(logp);
	    		f[m] = p;
	    	}
	    	fToA();
//	    	// slow discrete cosine transform
//	    	for (int i = 0; i < N; i++) {
//	    		double sum = 0;
//	    		for (int j = 0; j < n; j++) {
//	    			sum += f[j] * Math.cos((2*j+1)*i*Math.PI/(2*N));
//	    		}
//	    		if (i == 0) {
//	    			sum *= 1.0/n;//Math.sqrt(n);
//	    		} else {
//	    			sum *= 2.0/n;//Math.sqrt(n);
//	    		}
//	    		a[i] = sum;
//	    	}
//	    	System.out.println(Arrays.toString(f));
//			aToF();
		}
	}
	
	// convert coefficients to function values 
	public void aToF() {
		a[0] *= 2.0;
		a[a.length - 1] *= 2.0;
		f = transformer.transform(a, TransformType.FORWARD);
		a[0] /= 2.0;
		a[a.length - 1] /= 2.0;
	}

	// convert function values to coefficients
	public void fToA() {		
		a = transformer.transform(f, TransformType.INVERSE);
		a[0] /= 2.0;
		a[a.length - 1] /= 2.0;
 	}


public void transform_to_chebyshev_coef(double [] cheby_values, int K) 
{
	Complex [] cheby_complex = new Complex[2 * (K-1)];
	Complex [] cheby_complex_transform  = new Complex[2 * (K-1)];

 	for (int k = 0; k < K; ++k)
	{
		cheby_complex[k] = new Complex(cheby_values[K-1-k]);
	}
	for (int k = 1; k < K-1; ++k)
	{
		cheby_complex[K-1 + k] = new Complex(cheby_values[k]);
	}
	
	cheby_complex_transform = FFTtransformer.transform(cheby_complex, TransformType.INVERSE);
	
	for (int k = 0; k < K; ++k)
	{
		cheby_values[k] = 2 * cheby_complex_transform[k].getReal();
	}
	cheby_values[0] = cheby_values[0] / 2 ;
}

public void transform_to_chebyshev_values(double  [] cheby_coeff, int K)
{	
	Complex [] cheby_complex = new Complex[2 * (K-1)];
	Complex [] cheby_complex_transform = new Complex[2 * (K-1)];

	for (int k = 0; k < K; ++k)
	{
		cheby_coeff[k] = cheby_coeff[k] ;
	}
		cheby_coeff[0] = cheby_coeff[0] * 2 ;
	
	for (int k = 0; k < K; ++k)
	{
		cheby_complex[k] = new Complex(cheby_coeff[k]); 
	}
	for (int k = 1; k < K-1; ++k)
	{
		cheby_complex[K-1 +k] = new Complex(cheby_coeff[K-1-k]);
	}
	
	cheby_complex_transform = FFTtransformer.transform(cheby_complex, TransformType.FORWARD);

	for (int k = 0; k < K; ++k)
	{
		cheby_coeff[k] = cheby_complex_transform[K-1 - k].getReal() / 2 ;
	}
}




	
	public void setPolyFactors(double[] a) {
		if (this.a.length != a.length) {
			throw new IllegalArgumentException("Expected dimension " + this.a.length + " but got " + a.length);
		}
		System.arraycopy(a, 0, this.a, 0, a.length);
	}

	public void setPolyValues(double[] f) {
		if (this.f.length != f.length) {
			throw new IllegalArgumentException("Expected dimension " + this.f.length + " but got " + f.length);
		}
		System.arraycopy(f, 0, this.f, 0, f.length);
	}
	
    private double logBinom(int k, int n) {
    	double f = 0;
    	for (int i = k + 1; i <= n; i++) {
    		f += Math.log(i) - Math.log(n - i + 1);
    	}
		return f;
	}

    @Override
    public String toString() {
    	return "a: " + Arrays.toString(a) + "\n" +
    		   "f: " + Arrays.toString(f) + "\n";
    }

    /** use trapezoid rule to determine contributions of each interval **/
	public static double[] getInterValContributions(int N) {
		double [] f = new double[N];
    	for (int m = 1; m < N - 1; m++) {
    		f[m] = 0.5 - Math.cos(-m/(N-1.0)*Math.PI) / 2.0;
    	}
    	f[N-1] = 1.0;
		double [] w = new double[N-1];
    	for (int m = 0; m < N - 1; m++) {
    		w[m] = f[m+1] - f[m];
    	}
    	
    	double [] delta = new double[N];
    	delta[0] = w[0]/2.0;
    	for (int m = 1; m < N - 1; m++) {
    		delta[m] = 0.5*(w[m-1]+w[m]);
    	}
    	delta[N-1] = w[N-2]/2.0;
    	
		return delta;
	}
}
