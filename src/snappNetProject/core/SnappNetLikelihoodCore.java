package snappNetProject.core;
 
import java.util.List;

import beast.core.util.Log;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.alignment.Taxon;
import java.util.ArrayList;
import java.util.Arrays;
import java.lang.Exception;

public class SnappNetLikelihoodCore {

	SiteProbabilityCalculator [] m_siteProbabilityCalculator;
	int numPatterns;
	
	public SnappNetLikelihoodCore(Network speciesNetwork, SnapData data) {    	    	     
		
		numPatterns = data.getPatternCount();
        m_siteProbabilityCalculator = new SiteProbabilityCalculator[numPatterns];
        
        for(int id = 0; id < numPatterns; id++) {
        	m_siteProbabilityCalculator[id]= new SiteProbabilityCalculator(speciesNetwork);
        }
        
	}
	
	public SnappNetLikelihoodCore(Network speciesNetwork) {    	    	     
		
		numPatterns = 1;
        m_siteProbabilityCalculator = new SiteProbabilityCalculator[numPatterns];
        
        for(int id = 0; id < numPatterns; id++) {
        	m_siteProbabilityCalculator[id]= new SiteProbabilityCalculator(speciesNetwork);
        }
        
	}

	public double [] computeLogLikelihood(SnapData data, Network speciesNetwork, double u, double v, 
            Double [] coalescenceRate) throws Exception
	{
						
		int numPatterns = data.getPatternCount();
		
		//Temporarily store pattern probabilities... used for numerical checks.
        double [] patternProb = new double[numPatterns];
        List<TaxonSet> taxonSets=data.m_taxonsets.get();
		
		for(int id = 0; id < numPatterns; id++) {
            int [] dataAtThisSite = data.getPattern(id);
            int [] lineageCounts = data.getPatternLineagCounts(id);	               
            patternProb[id] = m_siteProbabilityCalculator[id].computeSiteLikelihood(dataAtThisSite, taxonSets, lineageCounts, speciesNetwork, u, v, coalescenceRate);            
            m_siteProbabilityCalculator[id]=null;            
        }
		
		return patternProb;
			
	}	



	
// 	public double [] test_computeLogLikelihood(List<TaxonSet> taxonSets, Network speciesNetwork, double u, double v, 
// 	Double [] coalescenceRate) throws Exception
// {
				
// int numPatterns = 1;

// //Temporarily store pattern probabilities... used for numerical checks.
// double [] patternProb = new double[numPatterns];
		 
// for(int id = 0; id < numPatterns; id++) {
// 	int [] dataAtThisSite = {1,2,3};
// 	int [] lineageCounts = {15,15,15};	               
// 	patternProb[id] = m_siteProbabilityCalculator[id].computeSiteLikelihood(dataAtThisSite, taxonSets, lineageCounts, speciesNetwork, u, v, coalescenceRate);            
// 	m_siteProbabilityCalculator[id]=null;            
// }

// return patternProb;
	
// }	

	// public void test()
	// {
	// 	int test = 1;
	// 	test = test + 1;
	// }

	// public static void main(String[] args) {
   
	// 	List<Taxon> Taxa_name = new ArrayList<>();
	// 	Taxon species_1 = new Taxon("taxa_1");
	// 	species_1.setID("taxa_1");

	// 	Taxon species_2 = new Taxon("taxa_2");
	// 	species_2.setID("taxa_2");
		
	// 	Taxon species_3 = new Taxon("taxa_3");
	// 	species_3.setID("taxa_3");

	// 	Taxa_name.add(species_1);
	// 	Taxa_name.add(species_2);
	// 	Taxa_name.add(species_3);
	// 	TaxonSet Simple_set = new TaxonSet(Taxa_name);
		
	// 	List<Taxon> TaxaSet_S1 = new ArrayList<>();
	// 	TaxaSet_S1.add(species_1);
	// 	TaxonSet Sspecies_1 = new TaxonSet(TaxaSet_S1);
	// 	Sspecies_1.setID("taxa_1");

	// 	List<Taxon> TaxaSet_S2 = new ArrayList<>();
	// 	TaxaSet_S2.add(species_2);
	// 	TaxonSet Sspecies_2 = new TaxonSet(TaxaSet_S2);
	// 	Sspecies_2.setID("taxa_2");
		
	// 	List<Taxon> TaxaSet_S3 = new ArrayList<>();
	// 	TaxaSet_S3.add(species_3);
	// 	TaxonSet Sspecies_3 = new TaxonSet(TaxaSet_S3);
	// 	Sspecies_3.setID("taxa_3");
		


	// 	//check to see if TaxonSet loaded correctly
	// 	System.out.println(Simple_set.getTaxonCount());
		
	// 	//Set-up a simple Caterpillar tree
	// 	Network Simple_net = new Network();
	// 	Simple_net.taxonSet = Simple_set;
	// 	Simple_net.test_initAndValidate();
		
	// 	//Add Speciation
	// 	NetworkNode speciationNode = new NetworkNode(Simple_net);
	// 	speciationNode.height = 1;

	// 	//Simple_net.addSpeciationNode(speciationNode);
	// 	//Add Reticulation 
	// 	NetworkNode reticNode = new NetworkNode(Simple_net);
	// 	reticNode.height = 0.5;
		
	// 	//Simple_net.addReticulationNode(reticNode);
	// 	Simple_net.updateRelationships();
	// 	//System.out.println(reticNode.nodeNumber);
	// 	Simple_net.addReticulationBranch(reticNode, speciationNode,1,2);

	// 	double u = 1;
	// 	double v = 1;
	// 	Double pop_size = 0.05;
	// 	Double [] Pop_size = {pop_size,pop_size,pop_size,pop_size,pop_size,pop_size,pop_size,pop_size};
		
	// 	List<TaxonSet> list_Simple_set = new ArrayList<>();
	// 	list_Simple_set.add(Simple_set);
		
	// 	List<TaxonSet> list_set = new ArrayList<>();
	// 	list_set.add(Sspecies_1);
	// 	list_set.add(Sspecies_2);
	// 	list_set.add(Sspecies_3);
		
	// 	try{
	// 	//System.out.println(new SnappNetLikelihoodCore(Simple_net).test_computeLogLikelihood(list_set, Simple_net, u, v, Pop_size));
	// 	}
	// 	catch (Exception e) {
	// 		/* This is a generic Exception handler which means it can handle
	// 		 * all the exceptions. This will execute if the exception is not
	// 		 * handled by previous catch blocks.
	// 		 */
	// 		System.out.println("Exception occurred");
	// 	 }

		//System.out.println(speciationNode.getChildren());

		//System.out.println(Simple_net.getReticulationOffset() );

		//System.out.println(Simple_net.taxonSet.getTaxonCount());
		//Simple_net.initAndValidate_noinput();
		//System.out.println(Simple_net.getNodeCount());
		//System.out.println(new SiteProbabilityCalculator().test_class_import());
		//new SiteProbabilityCalculator().leafLikelihood(data, sets, counts, u, v, Pop_size); 
	 //  }

	
}
