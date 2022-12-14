package snappNetProject.core;

import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.core.util.Log;
import beast.evolution.datatype.DataType;
import beast.evolution.substitutionmodel.EigenDecomposition;
import beast.evolution.substitutionmodel.SubstitutionModel;
import beast.evolution.tree.Node;
 
public class SnappNetSubstitutionModel extends SubstitutionModel.Base {
	public Input<RealParameter> m_pU = new Input<RealParameter>("mutationRateU", "Instantaneous rate of mutating from the 0 allele to the 1 alelle");
	public Input<RealParameter> m_pV = new Input<RealParameter>("mutationRateV", "Instantaneous rate of mutating from the 1 allele to the 0 alelle");
	public Input<RealParameter> m_pCoalescenceRate = new Input<RealParameter>("coalescenceRate", "population size parameter with one value for each node in the tree");
	
	public SnappNetSubstitutionModel() {
		frequenciesInput.setRule(Validate.OPTIONAL);
	}
	
    @Override
    public void initAndValidate() {
    	double u = m_pU.get().getValue();
    	double v  = m_pV.get().getValue();
    	if (Math.abs(2*u*v-u-v) > 1e-6) {
    		Log.warning.println("WARNING: Mutation rates are not normalised. "
    				+ "This means that the tree height may not be in units of substitution any more. "
    				+ "To ensure mutation rates are normalised, check that 2 * u * v = u + v, where v and u are the mutation rates.");
    	}
    }
    
	@Override
    public void getTransitionProbabilities(Node node, double fStartTime, double fEndTime, double fRate, double[] matrix) {}

	@Override
	public double[] getFrequencies() {return null;}

	@Override
	public EigenDecomposition getEigenDecomposition(Node node) {return null;}

	@Override
	public boolean canReturnComplexDiagonalization() {return false;}

	@Override
	public boolean canHandleDataType(DataType dataType) {return true;}

}
