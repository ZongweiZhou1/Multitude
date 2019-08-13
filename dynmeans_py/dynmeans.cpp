#ifndef __DYNMEANS_IMPL_HPP
#include "dynmeans.hpp"
using namespace std;

DynMeans::DynMeans(double lambda, double Q, double tau, bool verbose){
	this->verbose = verbose;
	this->lambda = lambda;
	this->Q = Q;
	this->tau = tau;
	this->ages.clear();
	this->oldprms.clear();
	this->oldprmlbls.clear();
	this->observations.clear();
	this->weights.clear();
	this->nextLbl = 0;

	this->prms.clear();
	this->rprms.clear();
	this->cnts.clear();
	this->lbls.clear();
	this->finalParams.clear();
	this->finalCnts.clear();
	this->finalLabels.clear();
	this->mars.clear();
	this->finalMars.clear();
	this->mrds.clear();
	this->finalMrds.clear();
}


DynMeans::~DynMeans(){}


void DynMeans::reset(){
	this->ages.clear();
	this->oldprms.clear();
	this->oldprmlbls.clear();
	this->observations.clear();
	this->weights.clear();
	this->nextLbl = 0;

	this->prms.clear();
	this->rprms.clear();
	this->cnts.clear();
	this->lbls.clear();
	this->finalParams.clear();
	this->finalCnts.clear();
	this->finalLabels.clear();
	this->mars.clear();
	this->finalMars.clear();
	this->mrds.clear();
	this->finalMrds.clear();

}

// tools to pass parameters between cpp and python
void DynMeans::set_tmpVariables(int observation_num){
    this->prms.clear();
    this->lbls.clear();
    this->cnts.clear();

    for(int i=0; i < this->oldprms.size(); i++){
        this->prms.push_back(this->oldprms[i]); //placeholders for updated parameters if the old ones get instantiated
        this->cnts.push_back(0); //start with count 0
    }
    //Initialization: no label on anything
    for(int i=0; i<observation_num; i++){
        this->lbls.push_back(-1);  // start with no labels on anything
    }

}

/* **********************************************************
* assgnOrdering:        random order of input points
* lbls:                 label assigned to each point
* cnts:                 points number of each group
* prms:                 feature of each center point
* **********************************************************/
void DynMeans::assignObservations(std::vector<int> assgnOrdering){
    std::vector<int>& lbls = this->lbls;
    std::vector<int>& cnts = this->cnts;
    std::vector<V2d>& prms = this->prms;

    for (int i=0; i<assgnOrdering.size(); i++){
        // get the observation idx from the random ordering
        int idx = assgnOrdering[i];
        // store the old lbl for possibly deleting clusers later
        int oldlbl = lbls[idx];
        //calculate the distances to all the parameters
        int minind = 0;
        double mindistsq = std::numeric_limits<double>::max();
        for(int j=0; j < prms.size(); j++){
            double tmpdistsq = (prms[j] - this->observations[idx]).squaredNorm();
            if(cnts[j] == 0){//the only way cnts can get to 0 is if it's an old parameters
                double gamma = 1.0/(1.0/ this->weights[j] + this->ages[j]*this->tau);
                tmpdistsq = gamma/(1.0+gamma)*tmpdistsq + this->ages[j]*this->Q;
            }
            if(tmpdistsq < mindistsq){
                minind = j;
                mindistsq = tmpdistsq;
            }
        }

        // if the minimum distance is still greater than lambda+startup cost, start a new cluster
        if (mindistsq > this->lambda)
        {
            prms.push_back(this->observations[idx]);
            lbls[idx] = prms.size()-1;
            cnts.push_back(1);
        }else{
            if(cnts[minind] == 0){ //if we just instantiated an old cluster
                //update its parameter to the current timestep so that upcoming assignment
                //are valid
                double gamma = 1.0/(1.0/this->weights[minind] + this->ages[minind]*this->tau);
                prms[minind] = (this->oldprms[minind]*gamma + this->observations[idx])/(gamma + 1);
            }
            lbls[idx] = minind;
            cnts[minind]++;
        }

        // if obs was previously assigned to something, decrease the count of the cluster it was assigned to
        // we do cluster deletion *after* assignment to prevent corner cases with monotinicity
        if (oldlbl != -1){
            cnts[oldlbl]--;
            //if this cluster now has no observations, but was a new one (no age recording for it yet)
            //remove it and shift labels downwards
            if(cnts[oldlbl] ==0 && oldlbl >= this->oldprms.size()){
                prms.erase(prms.begin() + oldlbl);
                cnts.erase(cnts.begin() + oldlbl);
                for(int j=0; j<lbls.size(); j++){
                    if(lbls[j] > oldlbl) lbls[j]--;
                }

            }else if (cnts[oldlbl] == 0){//it was an old parameter, reset it to the oldprm
                prms[oldlbl] = this->oldprms[oldlbl];
            }
        }
    }
    return;
 }

/* *********************************************************
* prms returned in this function will be used to match with
  targets from cv
* *********************************************************/
std::vector<double> DynMeans::first_updateParams(){
    for(int i=0; i < this->prms.size(); i++){
        if (this->cnts[i] > 0){
            std::vector<V2d> obsInClus = this->getObsInCluster(i, this->lbls);
            V2d tmpvec = obsInClus[0];
            for(int j=1; j < obsInClus.size(); j++)
                tmpvec = tmpvec + obsInClus[j];
            tmpvec = tmpvec / obsInClus.size();
            if(i < this->oldprms.size()){//updating an old param
                double gamma = 1.0/(1.0/this->weights[i] + this->ages[i]*this->tau);
				this->prms[i] = (this->oldprms[i]*gamma + tmpvec*this->cnts[i])/(gamma + this->cnts[i]);
            }else{
                this->prms[i] = tmpvec;
            }
        }
    }

    std::cout<<"current cluster number: "<<this->prms.size()<<endl;
    std::vector<double> prms;
    for(int i=0; i < this->prms.size(); i++){
        prms.push_back(this->prms[i](0));
        prms.push_back(this->prms[i](1));
    }
    return prms;
}

/* **********************************************************
* setting parameters given labels
* mars:         the matching relationship between current
                group and targets obtained from cv
* mrds:         the matching degree
* return:       the objective of current assignment
* **********************************************************/
double DynMeans::setParameters(std::vector<int> mars, std::vector<double> mrds){

	if (mars.size() != mrds.size() || mars.size() != this->prms.size()){
            cout<<"Error: mrds has invalid size with mars." <<endl;
            int ddd;
            cin>>ddd;
        }
    double objective = 0;

    std::vector<int>& lbls = this->lbls;
    std::vector<int>& cnts = this->cnts;
    std::vector<V2d>& prms = this->prms;
    std::vector<V2d>& rprms = this->rprms;
    this->mars = mars;
    this->mrds = mrds;

    for (int i = 0; i < prms.size(); i++){
		if (cnts[i] > 0){
		    std::vector<V2d> obsInClus = this->getObsInCluster(i, this->lbls);
			//add cost for new clusters - lambda
			// or add cost for old clusters - Q
			if (i < this->oldprms.size()){
				objective += this->Q*this->ages[i];
			} else {
				objective += this->lambda;
			}
			// update the prms using the matching result
			if(i < this->oldprms.size()){//updating an old param
			    double gamma = 1.0/(1.0/this->weights[i] + this->ages[i]*this->tau);
			    prms[i] = (prms[i] * (gamma + cnts[i]) + cnts[i]*mrds[i]*rprms[mars[i]])/(gamma+cnts[i]*(1+mrds[i]));
			    double tmpsqdist = (prms[i] - this->oldprms[i]).squaredNorm();
			    objective += gamma*tmpsqdist;
			}else{//just updating a new param
			    prms[i] = (prms[i] + mrds[i]*rprms[mars[i]])/(1+mrds[i]);
			    //no lag cost for new params
			}
			//get cost for prms[i]
			for (int j = 0; j < obsInClus.size(); j++){
				objective += (prms[i] - obsInClus[j]).squaredNorm();
			}
			// the targets from cv
//			objective += (prms[i] - rprms[mars[i]]).squaredNorm()*cnts[i]*mrds[i];
		}
	}
	return objective;
}

/* **********************************************************
* save the best result
* **********************************************************/
void DynMeans::set_finalPrms(){
    this->finalParams = this->prms;
    this->finalLabels = this->lbls;
    this->finalCnts = this->cnts;
    this->finalMars = this->mars;
    this->finalMrds = this->mrds;
}

/* **********************************************************
* This function is used when sampling parameters - it returns
  a vector of the observations in the next cluster, along with
  the index of the cluster.
  If this is the last parameter to be sampled, the function
  returns true; otherwise, false.
* idx:          the label of requested cluster
* lbls:         the labels corresponding to points
*************************************************************/
std::vector<V2d> DynMeans::getObsInCluster(int idx, std::vector<int> lbls){
    std::vector<V2d> obsIncluster;
    obsIncluster.reserve(lbls.size());
    for(int i=0; i < lbls.size(); i++){
        if(lbls[i] == idx)
            obsIncluster.push_back(this->observations[i]);
    }
    return obsIncluster;
}

/* *********************************************************
* update state with the cluster results and the cv targets,
  only the new cluster who assigned to a targets from cv will
  be reserved.
* mars:         the matching relationship between current
                group and targets obtained from cv
* mrds:         the matching degree
* return:
* **********************************************************/
std::vector<int> DynMeans::updateState(){

    std::vector<int>& lbls = this->finalLabels;
    std::vector<int>& cnts = this->finalCnts;
    std::vector<V2d>& prms = this->finalParams;
    std::vector<V2d>& rprms= this->rprms;
    std::vector<int>& mars = this->finalMars;
    std::vector<double>& mrds = this->finalMrds;

    if (mars.size() != mrds.size() || mars.size() != prms.size()){
            cout<<"Error: mrds has invalid size with mars." <<endl;
            int ddd;
            cin>>ddd;
        }
    this->oldprms = prms;
    int old_cluster_num = this->weights.size();
    std::vector<int> outLbls;  //stores the label output using oldprmlbls
    // update the weights/ages
    for(int i=0; i < prms.size(); i++){
        if (i<this->weights.size() && cnts[i]>0){
            //this is an instantiated cluster from a previous time; set age to 0 and update weights
            this->weights[i] = 1.0/(1.0/this->weights[i] + this->ages[i]*this->tau) + cnts[i]*(1+mrds[i]);
            this->ages[i] = 0;
        }else if(i >= this->weights.size()){
            //new cluster
            //push back a 0 for the age, and set the weight to the number of observation
            this->ages.push_back(0);
            this->weights.push_back(cnts[i]*(1+mrds[i]));
            this->oldprmlbls.push_back(this->nextLbl++);
        }
        this->ages[i]++;
    }
    for(int i=0; i<lbls.size(); i++){
        if (lbls[i] < old_cluster_num || mrds[lbls[i]] > 0 )
            outLbls.push_back(this->oldprmlbls[lbls[i]]);
        else
            outLbls.push_back(-1); // the corresponding points is deemed as a noisy
    }
    // now check whether the new cluster is valid, i.e., whether it assigned to a targets from cv
    for(int i=old_cluster_num-1; i<this->oldprms.size(); i++){
        if(mrds[i]==0){
        this->oldprms.erase(this->oldprms.begin()+i);
        this->oldprmlbls.erase(this->oldprmlbls.begin()+i);
        this->weights.erase(this->weights.begin()+i);
        this->ages.erase(this->ages.begin()+i);
        mrds.erase(mrds.begin()+i);
        i--;
        }
    }

    // check to see which clusters are permanently dead
    for(int i=0; i< this->oldprms.size(); i++){
        if(this->ages[i]*this->Q > this->lambda ){
            this->oldprms.erase(this->oldprms.begin()+i);
            this->oldprmlbls.erase(this->oldprmlbls.begin()+i);
            this->weights.erase(this->weights.begin()+i);
			this->ages.erase(this->ages.begin()+i);
			i--;
        }
    }
    return outLbls;
}


void DynMeans::set_data(std::vector<double>observations, std::vector<double> prms){
    // assign data to observations
    for(int i=0; i<observations.size(); i=i+2)
    {
        V2d newData(observations[i], observations[i+1]);
        this->observations.push_back(newData);
    }
    // assign data to clusterData from newobservations
    for(int j=0; j<prms.size(); j=j+2)
    {
        V2d newData(prms[j], prms[j+1]);
        this->rprms.push_back(newData);
    }

}

/* *******************************************************8
* debug function
* ********************************************************/
void DynMeans::pin_debug(int c){
    int numinst = 0;
    if(c==1){
        std::vector<int>& cnts = this->cnts;
        std::vector<V2d>& prms = this->prms;
    }else if (c==2){
        std::vector<int>& cnts = this->finalCnts;
        std::vector<V2d>& prms = this->finalParams;
    }else {
        std::cout<<"libdynmeans: ERROR: the case in `pin_debug` is invalid. [1|2]"<<endl;
        return;
    }
    for(int i=0; i<cnts.size(); i++)
        if(cnts[i]>0)
            numinst++;
    int numnew = prms.size() - this->oldprms.size();
    int numoldinst = numinst - numnew;
    int numolduninst = cnts.size() - numinst;
    std::cout<<"Old Uninst: "<< numolduninst << " old Inst: "<< numoldinst << " New: "<< numnew << endl;
}

#define __DYNMEANS_IMPL_HPP
#endif /* __DYNMEANS_IMPL_HPP */


