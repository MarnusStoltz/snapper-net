
/*
 * File FMatrixAugmented.java CE Rabier (thanks to V Berry), inspired by FMatrix.java from snapp ( Remco Bouckaert, David Bryant)
 * ADAPTED by MArnus Stoltz for SnappNetter
 * 
 */
package snappNetProject.core;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import beast.core.util.Log;
import snappNetProject.core.Num;
import snapper.ChebyshevPolynomial;
import snapper.MatrixExponentiator;
import org.apache.commons.math3.special.Gamma;
import org.apache.commons.math3.special.Beta;
/**FMatrixAugmented contains a  Fmatrix and extra key elements , see the paper for the different rules 
 * **/ 
/** The FMatrix will contain likelihood multiplied by the lineage probabilities Pr(Ry | n,r ) x Pr(n)  . Be careful it is multidimensional
 *  since we handle several branches simultaneously contrary to snapp
 *  **/


public class FMatrixAugmented {
		
	public List<Integer> branchNumbers;
	public List<String> branchLocations; //either T for Top or B for Bottom
	public List<Integer>  m_n_MultiDim;
	
	public List<Boolean> isBranchAboveRetic;		// either TRUE or FALSE, depending if the branch is located above 
	// a reticulation node
	int m_n;
	public double [] F;
	//SNAPPER METHOD -------------------------------------------------------------
	//Specify the number of Chebyshev coefficients
	int N = 9;
	//Stores Partial likelihoods for SNAPPER method
	//Chebyshev coefficients
	public double [] Partial_likelihood_c;
	//Function values
	public double [] Partial_likelihood_f;

	public double [][] Partial_likelihood_r;

	public boolean isLeft;

	public boolean isReticulateParentBranch;
	
	public int retic_vist_count = 0;

	public int left_branch_number = -1;
	
		//---------------------------------------------------------------------------
	
	public int getSize() {
		return m_n;
	}
	
	public List<Integer> getSizeMultidDim() {
		return m_n_MultiDim;
	}
	
	
	public int getIndex(int i) {		 
		int index = branchNumbers.indexOf(i);
		return index;		
	}
	
	public double [] getF() {						
		return F;
	}
	
	public double get(int n, int r) {return F[n*(n+1)/2-1+r];}
	public void set(int n, int r, double f) {F[n*(n+1)/2-1+r] = f;} //Needed to make this public to handle dominant markers. Could move that calculation here?


	public void addBranchNumbersAndLocations(int i, String loc){
			branchNumbers.add(i);
			branchLocations.add(loc);
			}
	
	
	
	public double [] asVectorCopy() {
		double [] copy = new double[F.length];
		System.arraycopy(F,0,copy,0,F.length);
		return copy;
	}
	public double [] asVectorCopyBase1() {
		double [] copy = new double[F.length + 1];
		System.arraycopy(F,0,copy,1,F.length);
		return copy;
		}
	
	
	public FMatrixAugmented() {
		m_n = 0;
		branchNumbers= new ArrayList<Integer>(); //be careful it is useful !
		// do not remove it !!!
	}

	public FMatrixAugmented(int n) {
		resize(n);
	}
	
	public FMatrixAugmented(FMatrixAugmented other) {
		assign(other);
	} // c'tor
	
	
		
	// CE we will see if we need that constructor
	/** construct a top-of-branch likelihood for leaf branch**/
	public FMatrixAugmented(int n, double []_F, List<Integer> _branchNumbers,
			List<String> _branchLocations) {
		
		
		Log.debug.println("Je suis ds le constructeur de FMatrixAugmented pour lexponentiar unidim \n");
		Log.debug.println("voici le F a copier" + _F +"\n");
		m_n_MultiDim = Arrays.asList(n);
		F = _F;
		branchNumbers = _branchNumbers;
		branchLocations = _branchLocations;
		branchLocations.set(0, "T");
	
	} // c'tor
	
	
	
///////////////////////////////////////////////// Case 3 Stars (Reticulation Node) 
	
		
	// Constructor for case 3*	(the right branch of retic node is located at the end of lists, the left branch is the last but one)  **/
	public FMatrixAugmented(FMatrixAugmented fBelowRetic, int leftBranchNumber, int belowReticBranchNumber, double inheritProb) {
		
		// In the lists, this.branchNumbers , this.branchLocations, this.isBranchAboveRetic, the last element will refer to the right branch above
		//the retic node, whereas the last but one will refer to the left branch
		
		m_n_MultiDim = new ArrayList<Integer>(fBelowRetic.m_n_MultiDim);		
		int indListBelowReticBranchNumber=fBelowRetic.getIndex(belowReticBranchNumber);
		int m_n_BelowRet=fBelowRetic.getSizeMultidDim().get(indListBelowReticBranchNumber);		
		
		m_n_MultiDim.remove(indListBelowReticBranchNumber);
		m_n_MultiDim.add(m_n_BelowRet);
		m_n_MultiDim.add(m_n_BelowRet);
		// the dimension for retic is handled at the end of the list 
		
						
		branchNumbers= new ArrayList<Integer>(fBelowRetic.branchNumbers);
		branchLocations= new ArrayList<String>(fBelowRetic.branchLocations);
		isBranchAboveRetic= new ArrayList<Boolean>(fBelowRetic.isBranchAboveRetic);
				
		branchNumbers.remove(indListBelowReticBranchNumber);
		branchLocations.remove(indListBelowReticBranchNumber);
		isBranchAboveRetic.remove(indListBelowReticBranchNumber);
		
					
		branchNumbers.add(leftBranchNumber);
		branchNumbers.add(leftBranchNumber+1);					
		branchLocations.add("B");
		branchLocations.add("B");								
		isBranchAboveRetic.add(true);
		isBranchAboveRetic.add(true);
		
	
		resizeMultiDim();//initialize F
				
		//fill in F						
		//be careful, n_listTop and r_listTop refers to the top below the retic node
		//that is why we use the location indListBelowReticBranchNumber for those list			 
		List<Integer> n_listTop = new ArrayList<Integer>(fBelowRetic.m_n_MultiDim);
		List<Integer> r_listTop = new ArrayList<Integer>(fBelowRetic.m_n_MultiDim);		
		
		recurseThreeStars(n_listTop, r_listTop, 0,
				fBelowRetic, indListBelowReticBranchNumber, m_n_BelowRet, inheritProb);	
		
	}
/// end constructor case 3*
	
	
		public void recurseThreeStars(List<Integer> n_listTop, List<Integer> r_listTop, int i, 
				FMatrixAugmented fBelowRetic, int indListBelowReticBranchNumber, int m_n_BelowRet, double inheritProb) {
						//it gives  a n_listTop and a r_listTop !!!
			
			//we do not set the last dimension for
			//  n_listBot and r_lisBot since the last dimension will be initialized later in fillF
			// so we use m_n_MultiDim.size()-1 as a terminal condition  !
	 
			if (i== (fBelowRetic.m_n_MultiDim.size()) ) {
								
				//we need to fill F			   
				fillFThreeStars(n_listTop, r_listTop, fBelowRetic,  indListBelowReticBranchNumber, m_n_BelowRet, inheritProb);
				
			}else if (i== indListBelowReticBranchNumber) {
				//do nothing !!!
				//we do not want to set this dimension 
				recurseThreeStars(n_listTop, r_listTop, i+1, fBelowRetic, indListBelowReticBranchNumber, m_n_BelowRet, inheritProb);							
			}else {
					int minNumberLineages= fBelowRetic.isBranchAboveRetic.get(i)==true ? 0 :1;
					// to know if at least 1 lineage go up, or if it is possible to have 0 lineage 			
					for (int j = minNumberLineages; j <= fBelowRetic.m_n_MultiDim.get(i); j++) {				
						n_listTop.set(i, j);				
						for (int k = 0; k <= j; k++) {					
							r_listTop.set(i, k);
							recurseThreeStars(n_listTop, r_listTop, i+1, fBelowRetic, indListBelowReticBranchNumber, m_n_BelowRet, inheritProb);								
						}						
					}
			
				}	
					
		}
		
	
		public void fillFThreeStars(List<Integer> n_listTop, List<Integer> r_listTop,  
				FMatrixAugmented fBelowRetic, int indListBelowReticBranchNumber, int m_n_BelowRet, double  inheritProb) {
								
			
			List<Integer> n_listBot = new ArrayList<Integer>(n_listTop);
			n_listBot.remove(indListBelowReticBranchNumber);
			n_listBot.add(0);
			n_listBot.add(0);		
			//at the end of lists we have the right branch above the retic
			//and the last but one element refers to the left branch
			//at this time these two entries are filled with 0 but will be updated later
			
			List<Integer> r_listBot = new ArrayList<Integer>(r_listTop);
			r_listBot.remove(indListBelowReticBranchNumber);
			r_listBot.add(0);
			r_listBot.add(0);
			//same comment as above but for r_listBot
			
				
			for (int s = 1; s <= m_n_BelowRet; s++) {
				//s lineages at the top of the edge below retic node
												
				n_listTop.set(indListBelowReticBranchNumber, s);			
				
				for (int j = 0; j <= s; j++) {
					// j among lineages go left, s-j go right						
							
					n_listBot.set(n_listBot.size()-2, j);  //for left branch
					n_listBot.set(n_listBot.size()-1, s-j);	 // for right branch				
				
						for (int rTopBelow = 0; rTopBelow <= s; rTopBelow++) {
														
							r_listTop.set(indListBelowReticBranchNumber, rTopBelow);
									
							//for (int rLeft = 0; rLeft<= Math.min(j,rTopBelow); rLeft++) { 													
							for (int rLeft = Math.max(0,rTopBelow-s+j); rLeft<= Math.min(j,rTopBelow); rLeft++) { 
								
								r_listBot.set(n_listBot.size()-2, rLeft);
								r_listBot.set(n_listBot.size()-1, rTopBelow-rLeft);																				
								Num C=new Num();
								double coeff= Math.pow(inheritProb,j)*Math.pow(1-inheritProb,s-j) * C.combination(s,j);							
								F[getLocationMultidDim(n_listBot, r_listBot)] =	coeff * fBelowRetic.getF()[fBelowRetic.getLocationMultidDim(n_listTop, r_listTop)];	
								
								
							}		
						}	
				}	
			
			}
					
			
		}
	
////////	/////////////////////////// END case 3* ( i.e. retic node)
	
		
		
	

//////////////////////////////////////////////////////////////////// Case Two stars

		//Constructor for case 2* 
		public  FMatrixAugmented(FMatrixAugmented fFirstChild, FMatrixAugmented fSecondChild, int branchFirstChild, 
				int branchSecondChild, int branchNumber) {
			
			
			//To begin let us initialize	
			m_n_MultiDim = new ArrayList<Integer>(fFirstChild.m_n_MultiDim);	
			branchNumbers= new ArrayList<Integer>(fFirstChild.branchNumbers);
			branchLocations= new ArrayList<String>(fFirstChild.branchLocations);
			isBranchAboveRetic= new ArrayList<Boolean>(fFirstChild.isBranchAboveRetic);		
			
		    int locationFirstChild=fFirstChild.getIndex(branchFirstChild);		
			m_n_MultiDim.remove(locationFirstChild);
			branchNumbers.remove(locationFirstChild);
			branchLocations.remove(locationFirstChild);
			isBranchAboveRetic.remove(locationFirstChild);
					
			int locationSecondChild=fSecondChild.getIndex(branchSecondChild);
			
			m_n_MultiDim.addAll(fSecondChild.m_n_MultiDim);
			m_n_MultiDim.remove(fFirstChild.m_n_MultiDim.size()-1+locationSecondChild);
			branchNumbers.addAll(fSecondChild.branchNumbers);
			branchNumbers.remove(fFirstChild.m_n_MultiDim.size()-1+locationSecondChild);
			branchLocations.addAll(fSecondChild.branchLocations);
			branchLocations.remove(fFirstChild.m_n_MultiDim.size()-1+locationSecondChild);
			isBranchAboveRetic.addAll(fSecondChild.isBranchAboveRetic);
			isBranchAboveRetic.remove(fFirstChild.m_n_MultiDim.size()-1+locationSecondChild);
		
			m_n_MultiDim.add(fFirstChild.m_n_MultiDim.get(locationFirstChild) + fSecondChild.m_n_MultiDim.get(locationSecondChild));	
			branchNumbers.add(branchNumber);
			branchLocations.add("B");
				
			if ( (fFirstChild.branchNumbers.size()>1) &   (fSecondChild.branchNumbers.size()>1) ) {
				isBranchAboveRetic.add(true);	
				//there could be none lineage on that branch 
			} else { isBranchAboveRetic.add(false);
			//at least one lineage is going up
			}			
			
			//Note: our node of interest is at the end of the list !!!
			//end of initialization
			
			//////////////////////////////
			resizeMultiDim();//initialize F
				
			
			//be careful, n_listBot and r_listBot refers to the Bottom above the internal node
			//we will use the locations at the end for those list		
			List<Integer> n_listBot = new ArrayList<Integer>(m_n_MultiDim);
			List<Integer> r_listBot = new ArrayList<Integer>(m_n_MultiDim);			
			
			
			recurseTwoStars(n_listBot, r_listBot, 0,
					fFirstChild, fSecondChild, locationFirstChild, locationSecondChild);
		}
				
	//end constructor for case 2*
		
	

	public void recurseTwoStars(List<Integer> n_listBot, List<Integer> r_listBot, int i, 
			FMatrixAugmented fFirstChild, FMatrixAugmented fSecondChild, int locationFirstChild, int locationSecondChild) {
					//it gives  a n_listBot and a r_listBot !!!
		
		//we do not set the last dimension for
		//  n_listBot and r_lisBot since the last dimension will be initialized later in fillF
		// so we use m_n_MultiDim.size()-1 as a terminal condition  !
 
			if (i== (m_n_MultiDim.size()-1) ) {
			 
			    fillFTwoStars(fFirstChild, fSecondChild, n_listBot, r_listBot, locationFirstChild, locationSecondChild);
			
			} else {
		
				int minNumberLineages= isBranchAboveRetic.get(i)==true ? 0 :1;
				// to know if at least 1 lineage go up, or if it is possible to have 0 lineage 			
				for (int j = minNumberLineages; j <= m_n_MultiDim.get(i); j++) {				
					n_listBot.set(i, j);				
					for (int k = 0; k <= j; k++) {					
						r_listBot.set(i, k);			
						recurseTwoStars(n_listBot, r_listBot, i+1, fFirstChild, fSecondChild,
							locationFirstChild, locationSecondChild);					
					}						
				}
		
			}
	
		
	}
	
	
	
	
	public void fillFTwoStars(FMatrixAugmented fFirstChild, FMatrixAugmented fSecondChild, 
			List<Integer> n_listBot, List<Integer> r_listBot, 
			int locationFirstChild, int locationSecondChild) {
		
		
		Num C=new Num(); //for computing combinations
		double fFirstChildValue= 0;
		double fSecondChildValue=0;
		double coeff=0;
		
		List<Integer> n_listTopFirstChild = new ArrayList<Integer>(fFirstChild.m_n_MultiDim);
		List<Integer> r_listTopFirstChild = new ArrayList<Integer>(fFirstChild.m_n_MultiDim);				
		List<Integer> n_listTopSecondChild = new ArrayList<Integer>(fSecondChild.m_n_MultiDim);
		List<Integer> r_listTopSecondChild = new ArrayList<Integer>(fSecondChild.m_n_MultiDim);				
		
		int minNumberLineagesFirstChild=  fFirstChild.isBranchAboveRetic.get(locationFirstChild)==true ? 0 :1;
		int minNumberLineagesSecondChild=  fSecondChild.isBranchAboveRetic.get(locationSecondChild)==true ? 0 :1;
		int maxNumberLineagesFirstChild=fFirstChild.m_n_MultiDim.get(locationFirstChild);
		int maxNumberLineagesSecondChild=fSecondChild.m_n_MultiDim.get(locationSecondChild);						
		int minNumberLineages= isBranchAboveRetic.get(isBranchAboveRetic.size()-1)==true ? 0:1;
		
		
		//initialize n_listTopFirstChild, r_listTopFirstChild, n_listTopSecondChild, r_listTopSecondChild
		// with nBot and rBot
			
		// all the branches contained in n_listTopFirstChild are contained in n_Bot
		//excepted the branch corresponding to FirstChild
		
		for (int i=0; i<fFirstChild.branchNumbers.size(); i++) {			
			int branchNb =fFirstChild.branchNumbers.get(i);			
			
			if (branchNb !=  fFirstChild.branchNumbers.get(locationFirstChild)) {
			int idBranchInBot = this.branchNumbers.indexOf(branchNb); 
			n_listTopFirstChild.set(i,n_listBot.get(idBranchInBot)); 	
			r_listTopFirstChild.set(i,r_listBot.get(idBranchInBot)); 	}		
		
		}
		
		// all the branches contained in n_listTopSecondChild are contained in n_Bot
		//excepted the branch corresponding to SecondChild
		for (int i=0; i<fSecondChild.branchNumbers.size(); i++) {			
			int branchNb =fSecondChild.branchNumbers.get(i);				
			
			if (branchNb !=  fSecondChild.branchNumbers.get(locationSecondChild)) {
				int idBranchInBot = this.branchNumbers.indexOf(branchNb); 
				n_listTopSecondChild.set(i,n_listBot.get(idBranchInBot)); 	
				r_listTopSecondChild.set(i,r_listBot.get(idBranchInBot)); 	}
			
		}
		
		
		
		for (int s = minNumberLineages; s <= m_n_MultiDim.get(m_n_MultiDim.size()-1); s++) {
			//s lineages at the bottom of the edge of interest												
			n_listBot.set(n_listBot.size()-1, s);
			
			for (int k = 0; k <= s; k++) {
				//k red lineages among the s 
				
				r_listBot.set(r_listBot.size()-1, k);
						
				//for (int j = 0; j <= s; j++) {
				for (int j = Math.max(minNumberLineagesFirstChild, s-maxNumberLineagesSecondChild); j <= Math.min(s-minNumberLineagesSecondChild, maxNumberLineagesFirstChild); j++) {
					// j among s lineages come from first child
									
					n_listTopFirstChild.set(locationFirstChild, j); 
					n_listTopSecondChild.set(locationSecondChild, s-j);
										
					//for (int m = Math.max(0,k-fSecondChild.m_n_MultiDim.get(locationSecondChild)); m <= Math.min(j,k); m++) {
					for (int m =Math.max(0, Math.max(-s+j+k,k-fSecondChild.m_n_MultiDim.get(locationSecondChild))); m <= Math.min(j,k); m++) {
						// m red lineages among the j lineages of the first child
						
						r_listTopFirstChild.set(locationFirstChild, m);
						r_listTopSecondChild.set(locationSecondChild, k-m);
							
						fFirstChildValue= fFirstChild.getF()[fFirstChild.getLocationMultidDim(n_listTopFirstChild, r_listTopFirstChild)];
						fSecondChildValue= fSecondChild.getF()[fSecondChild.getLocationMultidDim(n_listTopSecondChild, r_listTopSecondChild)];
						coeff= C.combination(k,m) * C.combination(s-k,j-m) / C.combination(s,j);		
				
						F[getLocationMultidDim(n_listBot, r_listBot)] += fFirstChildValue * fSecondChildValue * coeff;
								
 			
			       }
				}
			}
		}		
		
	}	
	//end fillFTwoStars
	
///////////////////////////////////////////////////////////////////////  End case Two Stars
	

////////////////////////////////////////////////////////////  Case 4* 
	
	//constructor for case 4*  
		public  FMatrixAugmented(FMatrixAugmented fFirstChild, int branchFirstChild, 
				int branchSecondChild, int branchNumber, int nMax) {		
			
			// in this constructor, no need to have two FMatChild since they are equal !!!
			
			//To begin let us initialize	
			m_n_MultiDim = new ArrayList<Integer>(fFirstChild.m_n_MultiDim);	
			branchNumbers= new ArrayList<Integer>(fFirstChild.branchNumbers);
			branchLocations= new ArrayList<String>(fFirstChild.branchLocations);
			isBranchAboveRetic= new ArrayList<Boolean>(fFirstChild.isBranchAboveRetic);		
		
			//where are the first and second child in fFirstChild ? 
			int locationFirstChild=fFirstChild.getIndex(branchFirstChild);
			int locationSecondChild=fFirstChild.getIndex(branchSecondChild);				
			int locFirstToRemove = locationFirstChild>locationSecondChild == true ? locationFirstChild :locationSecondChild;
			int locSecondToRemove = locationFirstChild>locationSecondChild == true ? locationSecondChild :locationFirstChild;
			
			int nbMax=nMax; //max number of lineages that go along this edge
		
			m_n_MultiDim.remove(locFirstToRemove);
			m_n_MultiDim.remove(locSecondToRemove);			
			m_n_MultiDim.add(nbMax);
			//we set the branch of interest at the end
			
			branchNumbers.remove(locFirstToRemove);
			branchNumbers.remove(locSecondToRemove);
			branchNumbers.add(branchNumber);
			
			branchLocations.remove(locFirstToRemove);
			branchLocations.remove(locSecondToRemove);
			branchLocations.add("B");
			
			isBranchAboveRetic.remove(locFirstToRemove);
			isBranchAboveRetic.remove(locSecondToRemove);
						
			if ( fFirstChild.branchNumbers.size() > 2 ) {
				isBranchAboveRetic.add(true); //since it is possible to have no lineages along the edge
			}else {
				// fFirstChild.branchNumbers.size() is equal to 2
				//at least one lineage is going up
				isBranchAboveRetic.add(false); 
			}
								
			resizeMultiDim();//initialise F
											
			//be careful, n_listBot and r_listBot refers to the Bottom above the internal node
			//we will use the locations at the end for those list		
			List<Integer> n_listBot = new ArrayList<Integer>(m_n_MultiDim);
			List<Integer> r_listBot = new ArrayList<Integer>(m_n_MultiDim);					
			recurseFourStars(n_listBot, r_listBot, 0, fFirstChild, locationFirstChild, locationSecondChild);		
		
		} //end constructor 4*
			
	
	
	public void recurseFourStars(List<Integer> n_listBot, List<Integer> r_listBot, int i, 
			FMatrixAugmented fFirstChild,  int locationFirstChild, int locationSecondChild) {
					//it gives  a n_listBot and a r_listBot !!!
		
		//we do not set the last dimension for
		//  n_listBot and r_lisBot since the last dimension will be initialized later in fillFFourStars
		// so we use m_n_MultiDim.size()-1 as a terminal condition  !
 
			if (i== (m_n_MultiDim.size()-1) ) {
			 		
				fillFFourStars(fFirstChild,  n_listBot, r_listBot, locationFirstChild, locationSecondChild);
			
			} else {
		
				int minNumberLineages= isBranchAboveRetic.get(i)==true ? 0 :1;
				// to know if at least 1 lineage go up, or if it is possible to have 0 lineage 			
				for (int j = minNumberLineages; j <= m_n_MultiDim.get(i); j++) {				
					n_listBot.set(i, j);				
					for (int k = 0; k <= j; k++) {					
						r_listBot.set(i, k);			
						recurseFourStars(n_listBot, r_listBot, i+1, fFirstChild, 
							locationFirstChild, locationSecondChild);					
					}						
				}
		
			}
	
		
	}
	
	
		public void fillFFourStars(FMatrixAugmented fFirstChild, 
				List<Integer> n_listBot, List<Integer> r_listBot, 
				int locationFirstChild, int locationSecondChild) {
			
			//locationFirstChild is the location of First Child in the Fmatrix fFirstChild
			//locationSecondChild is the location of Second Child in the Fmatrix fFirstChild
			
			Num C=new Num(); //for computing combinations			 
			double fFirstAndSecondChildValue=0;
			double coeff=0;
		
			List<Integer> n_listTopFirstAndSecondChild = new ArrayList<Integer>(fFirstChild.m_n_MultiDim);
			List<Integer> r_listTopFirstAndSecondChild = new ArrayList<Integer>(fFirstChild.m_n_MultiDim);		
			//n_listTopFirstAndSecondChild will be the n_list for First Child  and Second Child
			//r_listTopFirstAndSecondChild will be the r_list for First Child  and Second Child
						
			int minNumberLineagesFirstChild=  fFirstChild.isBranchAboveRetic.get(locationFirstChild)==true ? 0 :1;
			int minNumberLineagesSecondChild=  fFirstChild.isBranchAboveRetic.get(locationSecondChild)==true ? 0 :1;
			int maxNumberLineagesFirstChild=fFirstChild.m_n_MultiDim.get(locationFirstChild);
			int maxNumberLineagesSecondChild=fFirstChild.m_n_MultiDim.get(locationSecondChild);					
			
			int minNumberLineages= isBranchAboveRetic.get(isBranchAboveRetic.size()-1)==true ? 0:1;
									
			//initialize n_listTopFirstAndSecondChild, r_listTopFirstAndSecondChild 
			// with nBot and rBot   
					
			// all the branches contained in n_listTopFirstAndSecondChild are contained in n_Bot
			//excepted the branch corresponding to First Child and Second Child
			
			for (int i=0; i<fFirstChild.branchNumbers.size(); i++) {			
				int branchNb =fFirstChild.branchNumbers.get(i);			
				
				if ( (branchNb !=  fFirstChild.branchNumbers.get(locationFirstChild)) &&  (branchNb !=  fFirstChild.branchNumbers.get(locationSecondChild)) )
				{
				int idBranchInBot = this.branchNumbers.indexOf(branchNb); 
				n_listTopFirstAndSecondChild.set(i,n_listBot.get(idBranchInBot)); 	
				r_listTopFirstAndSecondChild.set(i,r_listBot.get(idBranchInBot)); 	}					
			}
		
			
			for (int s = minNumberLineages; s <= m_n_MultiDim.get(m_n_MultiDim.size()-1); s++) {
				//s lineages at the bottom of the edge of interest
													
				n_listBot.set(n_listBot.size()-1, s);				
				for (int k = 0; k <= s; k++) {
					//k red lineages among the s 					
					r_listBot.set(r_listBot.size()-1, k);
												
					for (int j= Math.max(minNumberLineagesFirstChild, s-maxNumberLineagesSecondChild) ; j<=Math.min(s-minNumberLineagesSecondChild, maxNumberLineagesFirstChild); j++) {						
						// j among s lineages come from first child
						
						n_listTopFirstAndSecondChild.set(locationFirstChild, j);
						n_listTopFirstAndSecondChild.set(locationSecondChild, s-j);
																	
						for (int m =Math.max(0, Math.max(-s+j+k,k-fFirstChild.m_n_MultiDim.get(locationSecondChild))); m <= Math.min(j,k); m++) {
							// m red lineages among the j lineages of the first child
														
							r_listTopFirstAndSecondChild.set(locationFirstChild, m);
							r_listTopFirstAndSecondChild.set(locationSecondChild, k-m);					 									
							fFirstAndSecondChildValue = fFirstChild.getF()[fFirstChild.getLocationMultidDim(n_listTopFirstAndSecondChild, r_listTopFirstAndSecondChild)];					 
							coeff= C.combination(k,m) * C.combination(s-k,j-m) / C.combination(s,j);							
							F[getLocationMultidDim(n_listBot, r_listBot)] += fFirstAndSecondChildValue * coeff;		
	 			
				       }
					}
				}
			}		
			

			
		}//end fillFFourSTars
		
/////////////////////////////////////////////////////////   End case 4stars
		
	
	
//////////////////////////////////////////////////////////////////////////////////////////
		

	void goToTopInternal(double u, double v, double gamma, double height ) {
		//to go at the top of the  branch located at the end of the list of branchNumbers
		
		try {
						
			if (isBranchAboveRetic.get(branchNumbers.size()-1)) { 
			// it is possible that no lineages are going up			
				F=MatrixExponentiatorMultiDim.expQTtxTopInternalNoLineagePossible(m_n_MultiDim.get(m_n_MultiDim.size()-1),u, v, gamma, height, this);
			}else {			
				//at least one lineage is going up
				F=MatrixExponentiatorMultiDim.expQTtxTopInternalAtLeastOneLineage(m_n_MultiDim.get(m_n_MultiDim.size()-1),u, v, gamma, height, this);			
			}
				
			branchLocations.set(branchLocations.size()-1, "T");
			
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	
////////////////////////////////////////////////////////////
		
	
	void goToBotTopRetic(double u, double v, double gamma, double height ) {
		//to go at the top of the right  branch, above the reticulation node
		//that is to say go from BotBot to BotTop
		
		// at this time, the right branch of the retic node is at the end of my lists
		try {
			
			F=MatrixExponentiatorMultiDim.expQTtxBotTopRetic(m_n_MultiDim.get(m_n_MultiDim.size()-1),u, v, gamma, height, this);						
			branchLocations.set(branchLocations.size()-1, "T");
			
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	
	
	
	void goToTopTopRetic(double u, double v, double gamma, double height ) {
		//to go at the top of the left  branch, above the reticulation node, assuming that we are already
		//at the top of right branch
		//that is to say go from BotTop to TopTop
			
		try {
						
		F=MatrixExponentiatorMultiDim.expQTtxTopTopRetic(m_n_MultiDim.get(m_n_MultiDim.size()-2),u, v, gamma, height, this);						
		branchLocations.set(branchLocations.size()-2, "T");
		// at this time, the right branch of the retic node is at the end of my list
		//i.e. left branch is the last but one
			
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
		
	void goToTopLeaf( double u, double v, double gamma, double height ) {
				//to go at the top of a leaf branch
		try {
			//F = MatrixExponentiatorMultiDim.expQTtxTopLeaf(10 ,u, v, gamma, height, this);
			F = MatrixExponentiatorMultiDim.expQTtxTopLeaf(m_n_MultiDim.get(0),u, v, gamma, height, this);						
			branchLocations.set(0, "T");			
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}

	//SNAPPER METHOD -------------------------------------------------------------
	//a: Chebyshev coefficients at bottom of branch
	void GoToTopOfLeafBranch_SNAPPER(double time, double[] matrix, double[] a) {


		MatrixExponentiator e = new MatrixExponentiator();
		e.expCF(time, matrix, a);
		ChebyshevPolynomial c = new ChebyshevPolynomial(N);
	    //assign Chebychev coefficients to the instance c (of class ChebyshevPolynomial)
		c.a = a;
		c.aToF();
		 //assign function values to the instance c (of class ChebyshevPolynomial)
		//System.arraycopy(c.f, 0, a, 0, N);

		Partial_likelihood_c = c.a;
		Partial_likelihood_f = c.f;

	}
	//----------------------------------------------------------------------------------

	public static double [][] evaluate_joint_partial_likelihood_at_reticulation(double [] coeff, double gamma){
		double [] x = new double [coeff.length]; 
		double [][] A = new double [coeff.length][coeff.length];
		get_clobatto_points(x, coeff.length); 
		for (int i=0; i < coeff.length; i++){
			for (int j=0; j < coeff.length; j++) {
				A[i][j] =  Evaluate_chebyshev_function( gamma*x[i] + (1-gamma)*x[j], coeff);   
			}
		}
		return A;
	}

	public static double [][] go_to_top_of_reticulation(double t_left, double t_right, double [] Q_l, double [] Q_r, double [][] A){

		ChebyshevPolynomial c = new ChebyshevPolynomial(A.length);
		MatrixExponentiator e = new MatrixExponentiator();
	
		for (int j=0; j < A.length; j++) {
	
			for (int i=0; i < A.length; i++) {
				c.f[i] = A[i][j];
			}
			c.fToA();
			double [] a = c.a;
			e.expCF(t_left, Q_l, a);
			for (int i=0; i < A.length; i++) {
				A[i][j] = a[i];
			}
		}
	
		for (int i=0; i < A.length; i++) {
	
			for (int j=0; j < A.length; j++) {
				c.f[j] = A[i][j];
			}
			c.fToA();
			double [] a = c.a;
			e.expCF(t_right, Q_r, a);
			for (int j=0; j < A.length; j++) {
				A[i][j] = a[j];
			}
	
		}
	
	
		return A;
	}
	


	//SNAPPER METHOD -------------------------------------------------------------
	//a: Chebyshev coefficients at bottom of branch
	void GoToTopOfReticilateBranch_SNAPPER(double time, double[] matrix, double[] f, double prob) {

		ChebyshevPolynomial c = new ChebyshevPolynomial(N);
		c.f = f;
		c.fToA();
		double [] a = c.a;

		MatrixExponentiator e = new MatrixExponentiator();
		e.expCF(time, matrix, a);
		
	    //assign Chebychev coefficients to the instance c (of class ChebyshevPolynomial)
		c.a = a;
		c.aToF();
		//assign function values to the instance c (of class ChebyshevPolynomial)
		//System.arraycopy(c.f, 0, a, 0, N);

		Partial_likelihood_f = c.f;
		Partial_likelihood_c = c.a;

	}
	//----------------------------------------------------------------------------------

	//SNAPPER METHOD -------------------------------------------------------------
	//a: Chebyshev coefficients at bottom of branch
	void GoToTopOfInternalBranch_SNAPPER(double time, double[] matrix, double[] f) {
		//Go from function values to Cehbyshev coefficients
		ChebyshevPolynomial c = new ChebyshevPolynomial(N);
		c.f = f;
		c.fToA();
		double [] a = c.a;

		MatrixExponentiator e = new MatrixExponentiator();
		e.expCF(time, matrix, a);
		
	    //assign Chebychev coefficients to the instance c (of class ChebyshevPolynomial)
		c.a = a;
		c.aToF();
		//assign function values to the instance c (of class ChebyshevPolynomial)
		//System.arraycopy(c.f, 0, a, 0, N);

		Partial_likelihood_f = c.f;
		Partial_likelihood_c = c.a;
	}
	//----------------------------------------------------------------------------------



	/** construct a leaf likelihood **/
	public FMatrixAugmented(int n, int nReds, int N) {
		
		m_n_MultiDim=new ArrayList<Integer>(Arrays.asList(n));
		branchNumbers = new ArrayList<Integer>();
		branchLocations = new ArrayList<String>();
		isBranchAboveRetic = new ArrayList<Boolean>(Arrays.asList(false));
		
		resizeMultiDim();
		if (n > 0) {
			setMultiDim(Arrays.asList(n), Arrays.asList(nReds), 1.0);
			//specific to leaf likelihood because of the 1
		}
//SNAPPER METHOD -----------------------------
		computeSNAPPERpartials(n, nReds, N);
//--------------------------------------------	
		
	} // c'tor

//SNAPPER METHOD -------------------------------------
	void computeSNAPPERpartials(int n, int r, int N){
		ChebyshevPolynomial c = new ChebyshevPolynomial(N);
		c.init(r, n);
		Partial_likelihood_c = c.a;
		Partial_likelihood_f = c.f;
	}
// -----------------------------------------------------------



void resize(int n) {
		if (F != null && getSize() == n) {
			// no need to resize, just init to zero
			Arrays.fill(F, 0);
		}
		m_n = n;
		F = new double[(n+1)*(n+2)/2-1];
	} // resize
	
		
	
	public double getMultiDim(List<Integer> n_List, List<Integer> r_List) 
	{		 
		return F[getLocationMultidDim(n_List, r_List)];
	}
	

	public void setMultiDim(List<Integer> n_List, List<Integer> r_List, double f)
	{
		F[getLocationMultidDim(n_List, r_List)]=f;		
	} 	
	

		void resizeMultiDim() {
			
			int product = 1; 
			int ind=0;
		    for (int i : m_n_MultiDim) {
		     
		    		if (!isBranchAboveRetic.get(ind))
		    			{		    		         	
		    				product = product * ((i+1)*(i+2)/2-1);		    	 
		    			}else {
		    				//deal with and edge above retic node 
		    				product = product * ((i+1)*(i+2)/2);	
		    			}
		    		ind ++; 
		     }
			     
		     
			F = new double[product];			
		} // end resizeMultiDim
		
		 
		public int getLocationMultidDim(List<Integer> n_List, List<Integer> r_List) 
		{		
			int n_Lineage=0;
			int r_Lineage=0;
			int location=0;
			
			for (int j = 0; j < n_List.size(); j++) {				
				n_Lineage=n_List.get(j);	 
				r_Lineage = r_List.get(j);			
			
				if (!isBranchAboveRetic.get(j)) {		    	 				
					location = location + (n_Lineage*(n_Lineage+1)/2 + r_Lineage -1) * calc(m_n_MultiDim,j+1);
				} else {
					location = location + (n_Lineage*(n_Lineage+1)/2 + r_Lineage) * calc(m_n_MultiDim,j+1);					
				}
			
			}
			 
			return location;
		}
			
		
		

		public int  calc(List<Integer> n_MultidDim, int ind) 
		{
			
			int prod=1;
			int n_Lineage=0;
			for (int j = ind; j < n_MultidDim.size(); j++) {	
				n_Lineage=n_MultidDim.get(j);
			
				if (!isBranchAboveRetic.get(j)) {
					prod=prod*((n_Lineage+1)*(n_Lineage+2)/2 -1) ;			
				} else {
					prod=prod*((n_Lineage+1)*(n_Lineage+2)/2) ;	
				}
			}		
			return prod;
					
		}
		
		

		@SuppressWarnings("null")
		public FMatrixAugmented getClone() {		
						
			FMatrixAugmented myClone=new FMatrixAugmented();			
			myClone.branchLocations=new ArrayList<String>(this.branchLocations);
			myClone.branchNumbers= new ArrayList<Integer>(this.branchNumbers);
			myClone.isBranchAboveRetic= new ArrayList<Boolean>(this.isBranchAboveRetic);
			myClone.m_n_MultiDim= new ArrayList<Integer>(this.m_n_MultiDim); 
			myClone.F=new double[this.F.length];
			System.arraycopy(this.F, 0, myClone.F, 0, F.length);
		
			return myClone;
			
			
		}
		
		
		public boolean compare (FMatrixAugmented x) {
						
			return this.branchNumbers.containsAll(x.branchNumbers) && this.branchNumbers.size() == x.branchNumbers.size() ? true :false;						
			
		}
		
		
	
	
	void rawresize(int n) {
		if (F != null && getSize() == n) {
			return;
		}
		m_n = n;
		F = new double[(n+1)*(n+2)/2-1];
	} // rawresize
	

	public void assign(FMatrixAugmented other) {
		int n = other.getSize();
		rawresize(n);
		if (getSize() != other.getSize()) {
			System.err.println("diff in length " + getSize() +"!="+ other.getSize());
		}
		System.arraycopy(other.F, 0, F , 0 , F.length);
	} // assign
	
	public String toString() {
		int n = getSize();
		StringBuffer buf = new StringBuffer(); 
		for(int i=1;i<=n;i++) {
			for(int r=0;r<=i;r++) {
				buf.append(get(i,r)+" ");
			}
			buf.append(';');
		} 
		return buf.toString();
	} // toString


	public static void get_clobatto_points(double [] x, int K) 
    {
        for (int  k = 0; k < K; ++k) 
        {
            x[k] = Math.cos(-k/(K-1.0)*Math.PI) ;
        }
    }

    public static double Evaluate_chebyshev_function( double x ,double [] c){
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

   public static void main(String[] args) {
   
	    int [] data = {1,2,3}; 
		 
	    int [] counts = {5,5,5};
	    double u = 1; 
	    double v = 1; 
	    double Pop_size = 0.05;
		double branch_height = 0.1;
	 // System.out.println(new SiteProbabilityCalculator().test_class_import());
	  //  new FMatrixAugmented().goToTopLeaf(u, v, Pop_size, branch_height); 
   }			


	
} // class FMatrixAugmented
