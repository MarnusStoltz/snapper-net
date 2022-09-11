package snappNetProject.core;

import java.util.ArrayList;
import java.util.Arrays;
//import java.util.Collection;
import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.collect.Multiset;

import beast.core.BEASTInterface;
import beast.core.Input;
import beast.core.util.Log;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.alignment.Taxon;
import snappNetProject.matrix.QMatrix;
import snapper.ChebyshevPolynomial;
import snapper.QMatrix_snapper;
import snapper.MatrixExponentiator;
import org.apache.commons.math3.special.Gamma;
import org.apache.commons.math3.special.Beta;
import java.lang.Math.*;

public class IntegrateRoot {
    
    public int K = 9;

    public IntegrateRoot() {
    }
    
    void solver_root(double [] coeff, int K)
{
	coeff[K-1] = -16.0 * coeff[K-1];
	coeff[K-2] = -16.0 * coeff[K-2];	
	coeff[K-3] = -16.0 * (-0.125 * coeff[K-1] + coeff[K-3]);
	coeff[K-4] = -16.0 * (-0.125 * coeff[K-2] + coeff[K-4]);
	for (int  k = K - 5; k > 2; --k) 
	{
		coeff[k] = -16.0 * (-0.125 * coeff[k+2] + 0.0625*coeff[k+4] + coeff[k]); 
	}
	coeff[2] = -8.0 * (-0.125 * coeff[4] + 0.0625*coeff[6] + coeff[2]);
/*g(x) is a K-2 Chebyshev basis approximation*/
	coeff[1] = 0;
	coeff[0] = 0;	
}

double Clenshaw_Curtis_integration(double [] coeff, int K)
{
	double integral = 0;
	for (int  k = 0; k < K; ++k) 
	{
		if (k % 2 == 0)
        	integral += coeff[k]/(1-Math.pow(k,2));
	}
	return integral;
}

void get_clobatto_points(double [] x, int K) 
{
	for (int  k = 0; k < K; ++k) 
	{
		x[k] = 0.5 - Math.cos(-k/(K-1.0)*Math.PI) / 2.0;
	}
}

    public double ComputeLikelihood() {
      double sum = 0;
     return sum;   
    }

    double integrate_at_root(double [] partial_likelihood_bottom, double theta, double u, double v)
    {	
        double [] x = new double[K]; 
        double [] temp = new double[K]; 
        double root_integral = 0;
        
        double f_0 = partial_likelihood_bottom[0];
        double f_1 = partial_likelihood_bottom[K-1];
        double beta_1 = theta * (u / ( u + v));
        double beta_2 = theta * (u / ( u + v));
        
        ChebyshevPolynomial c = new ChebyshevPolynomial(K);
        //c.f = partial_likelihood_bottom;
        c.transform_to_chebyshev_coef(partial_likelihood_bottom, K);
        //c.transform_to_chebyshev_values(temp, K);
        //c.fToA()
        //transform_to_chebyshev_coef(partial_likelihood_bottom, tree);
    
        solver_root(partial_likelihood_bottom, K); 
        
        for (int i = 0; i < K-2; ++i)
        {
             temp[i] = partial_likelihood_bottom[i+2];
        } 
        
        //ChebyshevPolynomial c2 = new ChebyshevPolynomial(K);
        //c2.a = temp;
        //c2.aToF();
        c.transform_to_chebyshev_values(temp, K);

        //transform_to_chebyshev_values_root(&partial_likelihood_bottom[2], tree->root_cheby_complex, tree->p_root_forward, tree->K-2);
        get_clobatto_points(x,K);

        for (int k = 0; k < K; ++k)
        {
            temp[k] = Math.pow(x[k], beta_1) * Math.pow((1 - x[k]), beta_2) * temp[k];
        }
        
        // c2.f = temp;
        // c2.fToA();
        c.transform_to_chebyshev_coef(temp, K);
        //transform_to_chebyshev_coef_root(partial_likelihood_bottom, tree->root_cheby_complex, tree->p_root_backward, tree->K-2);;
        
        root_integral += Math.exp(Gamma.logGamma(beta_1 + beta_2) 
                        - Gamma.logGamma(beta_1) - Gamma.logGamma(beta_2))
            *Clenshaw_Curtis_integration(temp, K);
        
        root_integral += f_0 + beta_1 / (beta_1 + beta_2) * (f_1 - f_0); 
        
        return root_integral;
    }
    public static void main(String[] args) {

        // IntegrateRoot I = new IntegrateRoot(); 
        // double [] x = new double[I.K];
        // new IntegrateRoot().get_clobatto_points(x, I.K);
        // System.out.println(x);
        // ChebyshevPolynomial c = new ChebyshevPolynomial(I.K);
        // c.init(3, 10);
        // System.out.println(c.f);

        // double [] test = {0.00299497,0.00301063,0.00305787,0.00313753,0.00325098,0.00340017,0.00358756,0.00381621,0.00408969,0.0044121,0.00478805,0.0052226,0.0057212,0.00628965,0.00693398,0.00766028,0.00847465,0.00938293,0.0103905,0.0115022,0.0127218,0.014052,0.0154939,0.0170471,0.0187088,0.0204743,0.0223362,0.0242845,0.0263065,0.0283866,0.0305065,0.0326453,0.0347797,0.0368845,0.0389328,0.0408965,0.0427473,0.0444569,0.045998,0.0473446,0.0484733,0.0493635,0.0499983,0.0503645,0.0504539,0.0502625,0.0497917,0.0490475,0.0480406,0.0467865,0.0453044,0.0436172,0.0417506,0.0397327,0.0375928,0.0353612,0.033068,0.0307427,0.0284135,0.0261068,0.0238464,0.0216536,0.0195467,0.0175408,0.015648,0.0138771,0.012234,0.0107219,0.00934131,0.00809053,0.00696597,0.00596246,0.00507354,0.00429185,0.00360936,0.00301768,0.00250833,0.00207285,0.00170309,0.00139125,0.00113,0.000912592,0.00073285,0.00058521,0.000464716,0.000367001,0.000288254,0.000225186,0.000174985,0.000135265,0.000104025,7.95983e-05,6.06083e-05,4.59283e-05,3.46424e-05,2.60128e-05,1.94487e-05,1.44813e-05,1.07406e-05,7.93708e-06,5.8455e-06,4.29181e-06,3.14241e-06,2.29536e-06,1.67334e-06,1.21805e-06,8.85759e-07,6.43857e-07,4.6813e-07,3.40692e-07,2.48385e-07,1.81573e-07,1.33224e-07,9.82215e-08,7.28576e-08,5.44498e-08,4.10628e-08,3.13023e-08,2.41653e-08,1.89313e-08,1.50827e-08,1.22484e-08,1.01622e-08,8.6341e-09,7.52897e-09,6.75195e-09,6.23808e-09,5.94539e-09,5.85036e-09};
        // System.out.println(I.integrate_at_root(test, 0.1, 1, 1));  

        // double [] test = {1, 0, 0, 0, 0};
        // ChebyshevPolynomial c = new ChebyshevPolynomial(5);
        // c.a = test;
        // c.transform_to_chebyshev_values(test, 5);
        // c.transform_to_chebyshev_coef(test, 5);
        // c.a = test;
    }

}