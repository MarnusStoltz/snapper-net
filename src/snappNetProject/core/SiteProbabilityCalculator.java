//analogue  of SiteProbabilityCalulator.java of Snapp,
//but it is dedicated to SnappNet 

//author CE Rabier

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
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.TrapezoidIntegrator;


public class SiteProbabilityCalculator {
	int N = 9;
	FMatrixAugmented[] tableFMatAugmented;
	ArrayList<NetworkNode>  networkNodesReadyToBeTreated;
	ArrayList<NetworkNode>  networkNodesNotReadyToBeTreated;
	
	public SiteProbabilityCalculator() {
	}

	public SiteProbabilityCalculator(Network speciesNetwork) {    	    	     

		tableFMatAugmented = new FMatrixAugmented[speciesNetwork.getBranchCount()];
      	networkNodesReadyToBeTreated = new ArrayList <NetworkNode>();
		networkNodesNotReadyToBeTreated = new ArrayList <NetworkNode>();
		
		//
		// We will see if we need this later
        for (int j = 0; j < speciesNetwork.getBranchCount(); j++) {         	
        	tableFMatAugmented[j]= new FMatrixAugmented();
        }
		//		
				
		//Initialization of List of Nodes Ready to Be treated and
        //and also the list of nodes not Ready to Be treated       
        //Initialization of lists of Nodes Ready to Be treated 
        // i.e. only the leaves are ready to be treated !!!  
		final NetworkNode[] networkLeaves=speciesNetwork.getLeafNodes();       
        for (int j = 0; j < networkLeaves.length; j++) { 
        networkNodesReadyToBeTreated.add(networkLeaves[j]);
        }    
		
        //Initialization of  List of Nodes Not Ready to Be treated 
        // i.e. all the internal nodes are Not ready to be treated !!!     
        final NetworkNode[] networkInternalNodes=speciesNetwork.getInternalNodes(); 
        for (int j = 0; j < networkInternalNodes.length; j++) { 
            networkNodesNotReadyToBeTreated.add(networkInternalNodes[j]);
            }         
        	
}
	
	

	public void printListNodes(java.util.ListIterator<NetworkNode> listIterator, boolean ready) { //throws SAXException, JAXBException {
		//printing list of nodes just to check if our algorithm is fine !!!		
		//other methods will use this method by calling this method in the following way
		//  this.printListNodes(networkNodesReadyToBeTreated.listIterator(),true);
		//  this.printListNodes(networkNodesNotReadyToBeTreated.listIterator(),false) ;  									
		
		if (ready==true) {
			Log.debug.println("AWESOME !!! Here is my list of nodes ready !!!\n");}
    			else {
    				Log.debug.println("BOUHHHHHH !!!  Here is my list of nodes not ready  !!!\n");}            

         NetworkNode myNode;
         while(listIterator.hasNext()) {
         	myNode = listIterator.next() ;
         	Log.debug.println(myNode.getLabel() +" ");
         }
         Log.debug.println("\n");
        
    	
    }
	

	
	
	public double computeSiteLikelihood(int [] dataAtThisSite, List<TaxonSet> taxonSets, int [] lineageCounts, Network speciesNetwork, double u, double v, 
            Double [] coalescenceRate) throws Exception {
								
		leafLikelihood(dataAtThisSite, taxonSets, lineageCounts, u, v, coalescenceRate); //handle leaves branches and go to Top on those branches
		UpdateListNodesReadyOrNot(speciesNetwork); //update the lists
		
		NetworkNode temp_node = new NetworkNode();
        //System.out.println("here");
		while (!networkNodesReadyToBeTreated.isEmpty()){
			   
        		NetworkNode nodeReady = networkNodesReadyToBeTreated.listIterator().next() ; 
				//Reticulation
        		if (nodeReady.isReticulation()) {
					//System.out.println("Reticulation");
        		    //System.out.println(nodeReady.childBranchNumbers.get(0));
        			//case (3*), handle reticulation node        		 
        			reticulateLikelihood(nodeReady, u, v, coalescenceRate);       		
					updateReticulate(nodeReady); //update the nodes to be treated or not, and also FMatrix 
					//System.out.println("left branch");
					//System.out.println(tableFMatAugmented[7].left_branch_number);
					//System.out.println(tableFMatAugmented[6].left_branch_number);
				//Speciation	       		
        		}if (nodeReady.getChildCount()>1 && speciesNetwork.reticulationNodeCount==0) {
        			FMatrixAugmented FMatChild1=tableFMatAugmented[nodeReady.childBranchNumbers.get(0)]; 
        			FMatrixAugmented FMatChild2=tableFMatAugmented[nodeReady.childBranchNumbers.get(1)]; 
    	         
        													 
					internalLikelihoodTwoDifferentChildren(nodeReady, FMatChild1, FMatChild2,
									u, v, coalescenceRate);
							
					updateInternalLikTwoDiffChildren(nodeReady); 
							
							
        					
        			}
				if (nodeReady.getChildCount()>1 && speciesNetwork.reticulationNodeCount>0) {
					
        			//System.out.println("left branch");
					//System.out.println(tableFMatAugmented[7].left_branch_number);
					//System.out.println(tableFMatAugmented[6].left_branch_number);
					tableFMatAugmented[6].isLeft = true;
					tableFMatAugmented[7].isLeft = false;
					FMatrixAugmented FMatChild1=tableFMatAugmented[nodeReady.childBranchNumbers.get(0)]; 
        			FMatrixAugmented FMatChild2=tableFMatAugmented[nodeReady.childBranchNumbers.get(1)]; 
					
					
					//Assign FMatChild1 to be the reticulation node
					if (FMatChild2.left_branch_number>0){
						FMatChild1=tableFMatAugmented[nodeReady.childBranchNumbers.get(1)]; 
        			    FMatChild2=tableFMatAugmented[nodeReady.childBranchNumbers.get(0)];
					}
        			//System.out.println("Child branches");
					//System.out.println(nodeReady.childBranchNumbers.get(0));
					//System.out.println(nodeReady.childBranchNumbers.get(1));
					
					//System.out.println("left branch");
					//System.out.println(FMatChild1.left_branch_number);
					//System.out.println(FMatChild2.left_branch_number);
					speciation_bottom_bottom(nodeReady, FMatChild1, FMatChild2);
					
					//System.out.println(FMatChild1.left_branch_number);
					if (FMatChild1.retic_vist_count>1){
						tableFMatAugmented[nodeReady.gammaBranchNumber].Partial_likelihood_r = FMatChild1.Partial_likelihood_r;
						speciation_top_top(nodeReady, temp_node, u, v, coalescenceRate);
					}

					temp_node = nodeReady;

					updateInternalLikTwoDiffChildren(nodeReady); 
							
							
        					
        			}
		}
        		
        //need to handle the root        
        int branchRoot=speciesNetwork.getRoot().gammaBranchNumber;
        FMatrixAugmented rootFMatrix=tableFMatAugmented[branchRoot];
        
    	double likelihoodSite=0;
    	try {
			
			FMatrixAugmented FMatChild1=tableFMatAugmented[speciesNetwork.getRoot().childBranchNumbers.get(0)]; 
			FMatrixAugmented FMatChild2=tableFMatAugmented[speciesNetwork.getRoot().childBranchNumbers.get(1)]; 
			//System.out.println("Root");
			//System.out.println(speciesNetwork.getRoot().childBranchNumbers.get(0));
			//System.out.println(speciesNetwork.getRoot().childBranchNumbers.get(1));
			//
		//	likelihoodSite=doRootLikelihood(rootFMatrix, u, v, coalescenceRate[branchRoot], false);

		//SNAPPER-----------------------------------------------------------------------------------------------------------------
			likelihoodSite=doRootLikelihood_SNAPPER(speciesNetwork, FMatChild1, FMatChild2, u, v, coalescenceRate[branchRoot], false);
		//------------------------------------------------------------------------------------------------------------------------

		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
        
		return likelihoodSite;
		
	}
	
	
	public void leafLikelihood(int [] dataAtThisSite, List<TaxonSet> taxonSets, int [] lineageCounts, double u, double v, Double [] coalescenceRate) {
		//compute likelihood for leaves branches, and go to Top of those branches
				
		//handle leaves since leaves are ready to be treated
		//we are at bottom of these branches

		//Example on how to compute Partial likelihood using SNAPPER classes
		// ChebyshevPolynomial c = new ChebyshevPolynomial(N);
		// c.init(8, 16);
		
		// QMatrix_snapper Q = new QMatrix_snapper(N);
		// Q.setQ(1, 1, 0.01);

		// GoToTopOfBranch(0.001, Q.Q, c.a);
		

		for (int j = 0; j < networkNodesReadyToBeTreated.size(); j++) {        
        	        	
        		NetworkNode myNodeToTreat=networkNodesReadyToBeTreated.get(j);               	
        		//look for taxonset with same label as NetworkLeaves[j]
        		for (int i = 0; i < taxonSets.size(); i++) {         		
        			if (taxonSets.get(i).getID().equals(myNodeToTreat.getLabel())){ 	
        			        			        			        			
        				tableFMatAugmented[myNodeToTreat.gammaBranchNumber]= new FMatrixAugmented(lineageCounts[i],dataAtThisSite[i], N);  
        			
        			}// end if
        		
        	}
        	       	
        	tableFMatAugmented[myNodeToTreat.gammaBranchNumber].addBranchNumbersAndLocations(myNodeToTreat.gammaBranchNumber,"B");    	
       
       }                 
        
		 //Go to Top Of those branches
        //System.out.println("Let us go to TOP for leaves !!!\n");
		// let s move to the Top on those branches          
        //So let us update FMatrixAugmented for the leaves by going at the top of their branches
        for (int j = 0; j < networkNodesReadyToBeTreated.size(); j++) {         	       	
        	try {         					
        		NetworkNode myNodeToTreat=networkNodesReadyToBeTreated.get(j);
				double heightBranch=myNodeToTreat.getParentByBranch(myNodeToTreat.gammaBranchNumber).getHeight();        						
				
				tableFMatAugmented[myNodeToTreat.gammaBranchNumber].goToTopLeaf(u, v, coalescenceRate[myNodeToTreat.gammaBranchNumber], heightBranch);
				
				//SNAPPER METHOD ---------------
				QMatrix_snapper Q = new QMatrix_snapper(N);
				Q.setQ(u, v, coalescenceRate[myNodeToTreat.gammaBranchNumber]);
				//System.out.println(myNodeToTreat.gammaBranchNumber);
				tableFMatAugmented[myNodeToTreat.gammaBranchNumber].GoToTopOfLeafBranch_SNAPPER(heightBranch, Q.Q, tableFMatAugmented[myNodeToTreat.gammaBranchNumber].Partial_likelihood_c);
				// -----------------------------
				// for (int n = 0; n < N; n++){
				// 	System.out.printf(" %.3f", tableFMatAugmented[myNodeToTreat.gammaBranchNumber].Partial_likelihood_f[n] );
				// }
				// System.out.println(";");
        	} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			};    	
        }          

    
		
	}

		
	
	public void UpdateListNodesReadyOrNot(Network speciesNetwork) {		
		//remove all the elements (i.e. the leaves) from the list of nodes ready to be treated
        networkNodesReadyToBeTreated.clear();
        
        //fill the list of nodes ready with internal nodes whose children are leaves,
        //since at this time leaves have been treated	               
        int countChildLeaves;
    		int count;        
    		Multiset<NetworkNode> children;
    		final NetworkNode[] networkInternalNodes=speciesNetwork.getInternalNodes(); 
    	
        for (int j = 0; j < networkInternalNodes.length; j++) {                 	       	        	
        	//Multiset<NetworkNode> children=NetworkSpeciationNodes[j].getChildren(); 
        		children=networkInternalNodes[j].getChildren(); 
        		countChildLeaves=0;
        		count=0;
        		for (NetworkNode n: children) {        		
        			count = n.isLeaf() ? 1 : 0;       		
        			countChildLeaves=countChildLeaves + count;        		
        		}        		
        	       	
        		if ( countChildLeaves==children.size() ) {
        			//the speciation node is ready to be treated
        			networkNodesReadyToBeTreated.add(networkInternalNodes[j]);       		
        			//remove this node from the list NetworkNodesNotReadyToBeTreated
        			networkNodesNotReadyToBeTreated.remove(networkInternalNodes[j]);
        		} 
        	        	
        } 	        	
               
	//this.printListNodes(networkNodesReadyToBeTreated.listIterator(),true);
	//this.printListNodes(networkNodesNotReadyToBeTreated.listIterator(),false) ;  
	       
}
public double Evaluate_chebyshev_function( double x ,double [] c){
	double c0;
	double c1;
	double x2;
	if (c.length == 1){
	c0 = c[0];
	c1 = 0;
	}
	else if (c.length == 2){
	c0 = c[0];
	c1 = c[1];
	}
	else {
	x2 = 2*x;
	c0 = c[c.length-2];
	c1 = c[c.length-1];
	for (int i = 3; i < c.length + 1; i++)
	{
		double tmp = c0;
		c0 = c[c.length-i] - c1;
		c1 = tmp + c1*x2;
	}
	}
return c0 + c1*x;
}

public void get_clobatto_points(double [] x, int K) 
    {
        for (int  k = 0; k < K; ++k) 
        {
            x[k] = Math.cos(-k/(K-1.0)*Math.PI) ;
        }
    }
		
		
	
	public void reticulateLikelihood(NetworkNode nodeReady, double u, double v, Double [] coalescenceRate) {
		//compute likelihood at the Top of the two branches above retic node 
		//handle case 3*	
				     		
		 int belowReticBranchNumber=0;        		 
		 int lBranchNumber = nodeReady.gammaBranchNumber;            	 
		 Multiset<NetworkNode> reticChild = nodeReady.getChildren();          	 
		 for (NetworkNode n: reticChild) {        		
			 //use for because of Multiset<NetworkNode> 
			 belowReticBranchNumber = n.gammaBranchNumber;      		        		
		 }      
		 double [] coeff = tableFMatAugmented[belowReticBranchNumber].Partial_likelihood_c;  
		 double [][] A = tableFMatAugmented[belowReticBranchNumber].evaluate_joint_partial_likelihood_at_reticulation(coeff, nodeReady.inheritProb);

		 //start with right branch
		 int rBranchNumber=lBranchNumber+1;
		  tableFMatAugmented[rBranchNumber]=new FMatrixAugmented(tableFMatAugmented[belowReticBranchNumber], lBranchNumber, belowReticBranchNumber,
		 		 nodeReady.inheritProb); // version that fill extra dimension
		 
		tableFMatAugmented[rBranchNumber].isLeft = false;
		tableFMatAugmented[lBranchNumber].isLeft = true;

		tableFMatAugmented[rBranchNumber].left_branch_number = lBranchNumber;
		tableFMatAugmented[lBranchNumber].left_branch_number = lBranchNumber;
		//System.out.println("Parent Reticulation branches");
		//System.out.println(lBranchNumber);  
		//System.out.println(rBranchNumber);  
		//System.out.println("--------------------");
		 //SNAPPER--------------------------------------------------------------
		 QMatrix_snapper Qr = new QMatrix_snapper(N);
		 Qr.setQ(u, v, coalescenceRate[rBranchNumber]);
		 QMatrix_snapper Ql = new QMatrix_snapper(N);
		 Ql.setQ(u, v, coalescenceRate[lBranchNumber]);
		 
		 //----------------------------------------------------------------------
		 
		 // Let us go to TOP on the right side    	
		 double heightBranch_r=nodeReady.getParentByBranch(rBranchNumber).getHeight() - nodeReady.getHeight();     		  
		 //tableFMatAugmented[rBranchNumber].goToBotTopRetic(u, v, coalescenceRate[rBranchNumber], heightBranch);
		 
		 // Let us go to TOP on the left side            	
		  double heightBranch_l=nodeReady.getParentByBranch(lBranchNumber).getHeight() - nodeReady.getHeight();           	       	
		 //tableFMatAugmented[rBranchNumber].goToTopTopRetic(u, v, coalescenceRate[lBranchNumber], heightBranch);
		 
		 A = tableFMatAugmented[belowReticBranchNumber].go_to_top_of_reticulation(heightBranch_l, heightBranch_r, Ql.Q, Qr.Q, A);

		//  for (int i=0; i < A.length; i++){
		// 	for (int j=0; j < A.length; j++) {
		// 		System.out.println(A[i][j]);    
		// 	}
		// }
		tableFMatAugmented[rBranchNumber].isLeft = false;
		tableFMatAugmented[lBranchNumber].isLeft = true;
		 tableFMatAugmented[lBranchNumber].Partial_likelihood_r = A;
		 tableFMatAugmented[rBranchNumber].Partial_likelihood_r = A;
		 //----------------------------------------------------------------------      	
	}
		
	
		
	public void updateReticulate(NetworkNode nodeReady){
		//let us update the Fmatrices and the list of Nodes
		//it is related to case 3*
		
		  int lBranchNumber=nodeReady.gammaBranchNumber;
		  int rBranchNumber=lBranchNumber + 1;
				
		//  for (int i : tableFMatAugmented[rBranchNumber].branchNumbers) {   	
			 
		// 	boolean temp_flag = tableFMatAugmented[i].isLeft; 	

		// 	 if (i!=rBranchNumber) {
		// 		 tableFMatAugmented[i]=tableFMatAugmented[rBranchNumber];
		// 		 tableFMatAugmented[i].isLeft = temp_flag; 
		// 	 }		  		        		  		
		//  }
			
		 		   	
 		//let us update the list of nodes ready to be treated 		 
		 updateReticulateLists(nodeReady, rBranchNumber);
		 updateReticulateLists(nodeReady, lBranchNumber);
		 //need to remove the retic node 
		 networkNodesReadyToBeTreated.remove(nodeReady);
		 //this.printListNodes(networkNodesReadyToBeTreated.listIterator(),true);
		 //this.printListNodes(networkNodesNotReadyToBeTreated.listIterator(),false) ;  
		
	}
	
	 
	public void updateReticulateLists(NetworkNode nodeReady, int branchNumber){
		//related to case 3*, it is used by the method just above
		
		//let us update the list of nodes ready to be treated 		 	    
 		NetworkNode parentNode=nodeReady.getParentByBranch(branchNumber);   				
 		if (parentNode.isReticulation()) {  			
 			networkNodesReadyToBeTreated.add(parentNode);
 			networkNodesNotReadyToBeTreated.remove(parentNode);
 		}else{        			
     		int edgeNumber;    		
     		if (parentNode.childBranchNumbers.get(0)==branchNumber) {
     			edgeNumber=parentNode.childBranchNumbers.get(1);		
     		}else {
     			edgeNumber=parentNode.childBranchNumbers.get(0);    			
     		}     		    
     		if (!tableFMatAugmented[edgeNumber].branchNumbers.isEmpty()) {
     			networkNodesReadyToBeTreated.add(parentNode);
     			networkNodesNotReadyToBeTreated.remove(parentNode);  			
     		} 			   			 			
 		}
					
	}
	
	public void speciation_bottom_bottom(NetworkNode nodeReady, FMatrixAugmented FMatChild1, FMatrixAugmented FMatChild2){
		
		//System.out.println( "*");
		int branchNumber = nodeReady.gammaBranchNumber;
		int branchFirstChild=nodeReady.childBranchNumbers.get(0);
        int branchSecondChild=nodeReady.childBranchNumbers.get(1);
		
		ChebyshevPolynomial c = new ChebyshevPolynomial(N);
		double [][] A = new double [N][N];
		//System.out.println( "**");
		
		A = FMatChild1.Partial_likelihood_r;
		
		//System.out.println( "***");
		tableFMatAugmented[branchNumber]=new FMatrixAugmented();
				//System.out.println( "****");
		
				//System.out.println( FMatChild1.isLeft);
				//System.out.println( FMatChild2.isLeft);
	if (FMatChild1.isLeft){
		//System.out.println( "da_fuck1");
		for (int j=0; j < A.length; j++) {
			
			for (int i=0; i < A.length; i++) {
				c.a[i] = A[i][j];
			}
			c.aToF();
			double [] f = c.f;
			
			for (int i=0; i < A.length; i++) {
				A[i][j] = f[i] * FMatChild2.Partial_likelihood_f[i];
			}
		}
		//System.out.println( "da_fuck2");
	}
	else {
		//System.out.println( "da_fuck3");
		for (int i=0; i < N; i++) {
			//System.out.println( "*");
			for (int j=0; j < N; j++) {
				//System.out.println( "**");	
				c.a[j] = A[i][j];
			}
			//System.out.println( "***");
			c.aToF();
			//System.out.println( "****");
			double [] f = c.f;
			//System.out.println( "*****");
			for (int j=0; j < N; j++) {
				//System.out.println( "******");
				//System.out.println( f[j]);
				//System.out.println( FMatChild2.Partial_likelihood_f[j]);
				
				//System.out.println( "******");
				A[i][j] = f[j] * FMatChild2.Partial_likelihood_f[j];
			}
			
		}
		//System.out.println( "da_fuck4"); 
	}
	//System.out.println( "*****");
	FMatChild1.Partial_likelihood_r = A;
	
	//System.out.println( "******");
	tableFMatAugmented[FMatChild1.left_branch_number].retic_vist_count = tableFMatAugmented[FMatChild1.left_branch_number].retic_vist_count + 1;	
	
	//System.out.println( "*******");//System.out.println( tableFMatAugmented[branch1].Partial_likelihood_r[i][j]);   

	// for (int i=0; i < A.length; i++){
    //     for (int j=0; j < A.length; j++) {
    //         System.out.println( A[i][j]);   
    //     }
    // }

	}
	
	public void speciation_top_top(NetworkNode node1, NetworkNode node2,
	double u, double v, Double [] coalescenceRate){

        int branch1=node1.gammaBranchNumber;
        int branch2=node2.gammaBranchNumber;
		//System.out.println(branch1);
		//System.out.println(branch2);

		double [][] A = new double [N][N];
		QMatrix_snapper Q1 = new QMatrix_snapper(N);
		Q1.setQ(u, v, coalescenceRate[branch1]);
	    QMatrix_snapper Q2 = new QMatrix_snapper(N);
		Q2.setQ(u, v, coalescenceRate[branch2]);
		 
		//----------------------------------------------------------------------
		double heightBranch_1=node1.getParentByBranch(branch1).getHeight() - node1.getHeight();     		  	 	
		double heightBranch_2=node2.getParentByBranch(branch2).getHeight() - node2.getHeight();           	       	
		 
		
		// for (int i=0; i < A.length; i++){
		// 	for (int j=0; j < A.length; j++) {
		// 		System.out.println( tableFMatAugmented[branch1].Partial_likelihood_r[i][j]);   
		// 	}
		// }
		
		A = tableFMatAugmented[branch1].go_to_top_of_reticulation(heightBranch_1, heightBranch_2, Q1.Q, Q2.Q,  tableFMatAugmented[branch1].Partial_likelihood_r);

		tableFMatAugmented[branch1].Partial_likelihood_r = A;
		tableFMatAugmented[branch2].Partial_likelihood_r = A;

		// for (int i=0; i < A.length; i++){
		// 	for (int j=0; j < A.length; j++) {
		// 		System.out.println( A[i][j]);   
		// 	}
		// }
	}


	public void internalLikelihoodTwoDifferentChildren(NetworkNode nodeReady, FMatrixAugmented FMatChild1, FMatrixAugmented FMatChild2,
			double u, double v, Double [] coalescenceRate) {
		//handle case 2*
		double [] f_bottom = new double[N];		
        int branchNumber = nodeReady.gammaBranchNumber;
        int branchFirstChild=nodeReady.childBranchNumbers.get(0);
        int branchSecondChild=nodeReady.childBranchNumbers.get(1);
            
        tableFMatAugmented[branchNumber]=new FMatrixAugmented(tableFMatAugmented[branchFirstChild], 
        		tableFMatAugmented[branchSecondChild], branchFirstChild, branchSecondChild, branchNumber);
       
        double heightBranch=nodeReady.getParentByBranch(branchNumber).getHeight()-nodeReady.getHeight();        		
		tableFMatAugmented[branchNumber].goToTopInternal(u, v, coalescenceRate[branchNumber], heightBranch); 
		
		//SNAPPER--------------------------------------------------------------------------------
		//multiply partial likelihoods together at the bottom of the branch
		for (int i = 0; i < N; i++) { 
			f_bottom[i] = FMatChild1.Partial_likelihood_f[i]*FMatChild2.Partial_likelihood_f[i];
		}
		QMatrix_snapper Q = new QMatrix_snapper(N);
		Q.setQ(u, v, coalescenceRate[branchNumber]);
		tableFMatAugmented[branchNumber].GoToTopOfInternalBranch_SNAPPER(heightBranch , Q.Q, f_bottom);

		//----------------------------------------------------------------------------------------
	}
	
	public void internalLikelihoodwithReticulationChild(NetworkNode nodeReady, FMatrixAugmented FMatChild1, FMatrixAugmented FMatChild2,
			double u, double v, Double [] coalescenceRate) {
		//handle case 2*
		double [] f_bottom = new double[N];		
        int branchNumber = nodeReady.gammaBranchNumber;
        int branchFirstChild=nodeReady.childBranchNumbers.get(0);
        int branchSecondChild=nodeReady.childBranchNumbers.get(1);
            
        tableFMatAugmented[branchNumber]=new FMatrixAugmented(tableFMatAugmented[branchFirstChild], 
        		tableFMatAugmented[branchSecondChild], branchFirstChild, branchSecondChild, branchNumber);
       
        double heightBranch=nodeReady.getParentByBranch(branchNumber).getHeight()-nodeReady.getHeight();        		
		tableFMatAugmented[branchNumber].goToTopInternal(u, v, coalescenceRate[branchNumber], heightBranch); 
		
		//SNAPPER--------------------------------------------------------------------------------
		//multiply partial likelihoods together at the bottom of the branch
		for (int i = 0; i < N; i++) { 
			f_bottom[i] = FMatChild1.Partial_likelihood_f[i]*FMatChild2.Partial_likelihood_f[i];
		}
		QMatrix_snapper Q = new QMatrix_snapper(N);
		Q.setQ(u, v, coalescenceRate[branchNumber]);
		tableFMatAugmented[branchNumber].GoToTopOfInternalBranch_SNAPPER(heightBranch , Q.Q, f_bottom);
		//----------------------------------------------------------------------------------------
	}

	
	public void updateInternalLikTwoDiffChildren(NetworkNode nodeReady) {
		//related to case 2*		
		int branchNumber = nodeReady.gammaBranchNumber;		 
		// for (int i : tableFMatAugmented[branchNumber].branchNumbers) {   				
		//	boolean temp_flag = tableFMatAugmented[i].isLeft; 
		//	if (i!=branchNumber) {
		//		tableFMatAugmented[i]=tableFMatAugmented[branchNumber];
		//		tableFMatAugmented[i].isLeft = temp_flag;
		//	}	  				  		
		// }
		
		
		//check if parent is ready        				
		NetworkNode parentNode=nodeReady.getParentByBranch(branchNumber);
		if (!parentNode.isOrigin()) {
			if (parentNode.isReticulation()) {  			 
				networkNodesReadyToBeTreated.add(parentNode);
				networkNodesNotReadyToBeTreated.remove(parentNode);
			}else {
				int edgeNumber; 
				if (parentNode.childBranchNumbers.get(0)==branchNumber) {
					edgeNumber=parentNode.childBranchNumbers.get(1);		
				}else {
					edgeNumber=parentNode.childBranchNumbers.get(0);    			
				}  
			
				if (!tableFMatAugmented[edgeNumber].branchNumbers.isEmpty()) {
					networkNodesReadyToBeTreated.add(parentNode);
					networkNodesNotReadyToBeTreated.remove(parentNode);  			
				} 		
			
			}
		}

		
	
		 networkNodesReadyToBeTreated.remove(nodeReady);
		 //this.printListNodes(networkNodesReadyToBeTreated.listIterator(),true);
		 //this.printListNodes(networkNodesNotReadyToBeTreated.listIterator(),false) ;  		
	}

	
//////////////////////////////////////////////	
		
	public void internalLikelihoodTwins(NetworkNode nodeReady, FMatrixAugmented FMatChild1, double u, double v, Double [] coalescenceRate) {
		// handle case (4*)
		double [] f_bottom = new double[N];	

        ArrayList <Integer>  branchAboveDescendingLeaves = new ArrayList <Integer>();
        nodeReady.getLeafBranchNumber(branchAboveDescendingLeaves);
         
        int nMax=0; //will refer to the max number of lineages that go along this edge
        for (int i=0; i<branchAboveDescendingLeaves.size(); i++) {
        	 	 nMax += tableFMatAugmented[branchAboveDescendingLeaves.get(i)].m_n_MultiDim.get(0);
        }

        int branchNumber=nodeReady.gammaBranchNumber;
        int branchFirstChild=nodeReady.childBranchNumbers.get(0);
        int branchSecondChild=nodeReady.childBranchNumbers.get(1);       
        tableFMatAugmented[branchNumber]=new FMatrixAugmented(FMatChild1, branchFirstChild, branchSecondChild, branchNumber, nMax);
         	    	
        //need to go at the top of the branch
    	 	double heightBranch=nodeReady.getParentByBranch(branchNumber).getHeight()-nodeReady.getHeight();        	    	  
			 tableFMatAugmented[branchNumber].goToTopInternal(u, v, coalescenceRate[branchNumber], heightBranch); 	
			 
		//SNAPPER--------------------------------------------------------------------------------
		//multiply partial likelihoods together at the bottom of the branch
		for (int i = 0; i < N; i++) { 
			f_bottom[i] = FMatChild1.Partial_likelihood_f[i]*FMatChild1.Partial_likelihood_f[i];
		}
		QMatrix_snapper Q = new QMatrix_snapper(N);
		Q.setQ(u, v, coalescenceRate[branchNumber]);
		tableFMatAugmented[branchNumber].GoToTopOfInternalBranch_SNAPPER(heightBranch , Q.Q, f_bottom);
		//----------------------------------------------------------------------------------------
	}
		
		
		
	public void updateInternalLikTwins(NetworkNode nodeReady) {
		//related to case 4*		
		int branchNumber = nodeReady.gammaBranchNumber;
		for (int i : tableFMatAugmented[branchNumber].branchNumbers) {   		
		
			boolean temp_flag = tableFMatAugmented[i].isLeft; 
			if (i!=branchNumber) {
				tableFMatAugmented[i]=tableFMatAugmented[branchNumber];
				tableFMatAugmented[i].isLeft = temp_flag;
			}		  		
		  		
		}
		
		//check if parent is ready        				
		NetworkNode parentNode=nodeReady.getParentByBranch(branchNumber);
		if (!parentNode.isOrigin()) { 
			if (parentNode.isReticulation()) {  			
				Log.debug.println(" c est une reticulation !! \n");
				networkNodesReadyToBeTreated.add(parentNode);
				networkNodesNotReadyToBeTreated.remove(parentNode);
			}else {
				int edgeNumber;
				Log.debug.println(" ce n'est pas une reticulation !! \n");
				if (parentNode.childBranchNumbers.get(0)==branchNumber) {
					edgeNumber=parentNode.childBranchNumbers.get(1);		
				}else {
					edgeNumber=parentNode.childBranchNumbers.get(0);    			
				}  
					
				if (!tableFMatAugmented[edgeNumber].branchNumbers.isEmpty()) {
					networkNodesReadyToBeTreated.add(parentNode);
					networkNodesNotReadyToBeTreated.remove(parentNode);  			
				} 		
					
			}
		}
		
		 networkNodesReadyToBeTreated.remove(nodeReady);
		 //this.printListNodes(networkNodesReadyToBeTreated.listIterator(),true);
		 //this.printListNodes(networkNodesNotReadyToBeTreated.listIterator(),false) ;  
	}	
		
		 

	
	/**
	 * David Bryant 's code below
    Determines a non-zero right e-vector for the matrix Q, defined by u,v,coalescenceRate and N.
    The e-vector is normalised so that the entries for each n sum to 1.0

    //TODO: incorporate into code for abstract matrix
    */
   double [][]  findRootProbabilities(int N, double u, double v, double coalescenceRate, boolean dprint) throws Exception {
	   double [][] x;
       QMatrix Qt = new QMatrix(N,u,v,coalescenceRate);
       double [] xcol;
       xcol = Qt.findOrthogonalVector(dprint);
      
       if (dprint) {
           Log.debug.println("xcol = " +Arrays.toString(xcol));
       }
       
       int index = 1;
       x = new double[N+1][];
       for(int n=1;n<=N;n++) {
           x[n] = new double[n+1];
           double rowsum = 0.0;
           for(int r=0;r<=n;r++) {
               double xcol_index = Math.max(xcol[index], 0.0);
               rowsum += xcol_index;
               x[n][r] = xcol_index;
               index++;
           }
           for(int r=0;r<=n;r++)
               x[n][r] = x[n][r] / rowsum;
       }
       return x;
   } // findRootProbabilities

   
double doRootLikelihood(FMatrixAugmented rootFMatrix, double u, double v, double gamma, boolean dprint) throws Exception
   {	   
       
       int N = rootFMatrix.getSizeMultidDim().get(0);		
       double[][] conditional = findRootProbabilities(N, u, v, gamma, dprint);
		
       double sum = 0.0;
       for(int n=1;n<=N;n++) {
           for(int r=0;r<=n;r++) {
        	   	double term =  conditional[n][r] * 
        			   rootFMatrix.getF()[rootFMatrix.getLocationMultidDim(Arrays.asList(n),Arrays.asList(r))];       	   
        	   	sum += term;
        	   //CE we can leave this flag or not, I remove it for the moment
        	   	// System.out.println("Voici la valeur de rootFMatrix!!!" + rootFMatrix.getF()[rootFMatrix.getLocationMultidDim(Arrays.asList(n),Arrays.asList(r))]+ "\n");
            //   if (sum<0.0) {
             //      System.out.println("Numerical problems"+ sum +"\n")                   
              // }
           }
       }
       return sum;
   } // doRootLikelihood



//SNAPPER-----------------------------------------------------------------------------------------------------------------------------
   double doRootLikelihood_SNAPPER(Network speciesNetwork,FMatrixAugmented FMatrix1, FMatrixAugmented FMatrix2, double u, double v, double gamma, boolean dprint) throws Exception
   {	   
	   double sum = 0.0;
	   double [] f_bottom = new double[N];
	   
	  // TEST PRINT - BEGIN
	//    System.out.println("1--");	
	//    for (int n = 0; n < N; n++){
	// 				System.out.printf(" %.3f", FMatrix1.Partial_likelihood_f[n] );
	// 			}
	// 			System.out.println(";");	 
	// 	System.out.println("2--");	
	// System.out.println("2--");	
	// for (int n = 0; n < N; n++){
	// 			 System.out.printf(" %.3f", FMatrix2.Partial_likelihood_f[n] );
	// 		 }
	// 		 System.out.println(";");
		// TEST PRINT - END
		// if (FMatrix2.Partial_likelihood_f == null)
		// {
		if (speciesNetwork.reticulationNodeCount>0){
		    for (int i = 0; i < N; i++) { 
				f_bottom[i] = FMatrix1.Partial_likelihood_r[i][i];
				// for (int j=0; j < N; j++){
				// 	for (int k=0; k < N; k++) {
				// 		System.out.println(FMatrix1.Partial_likelihood_r[j][k]);    
				// 	}
				// }
			};
		} else{
			for (int i = 0; i < N; i++) { 
				f_bottom[i] = FMatrix1.Partial_likelihood_f[i]*FMatrix2.Partial_likelihood_f[i] ;
			};
		}
			
		// }		
	 	// else {
	//    for (int i = 0; i < N; i++) { 
	// 		f_bottom[i] = 1;
	// 	}
	//}
		IntegrateRoot I = new IntegrateRoot(); 
		sum = I.integrate_at_root(f_bottom, gamma, u, v);
		//Initialize root integration
		//System.out.println(String.format("%6.3e",sum));
       return sum;
   } // doRootLikelihood
//--------------------------------------------------------------------------------------------------------------------------------------------------
   public boolean test_class_import()
   {
	   try {
		   Class.forName("com.google.common.collect.Multiset");
		   return true;
	   } catch(Exception e) {
		   return false;
	   }
		   
   }
   
	}
	
	
	
	
	

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
