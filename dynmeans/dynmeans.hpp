#ifndef __DYNMEANS_HPP
#include<vector>
#include<iostream>
#include<algorithm>
#include<boost/static_assert.hpp>
#include<boost/function.hpp>
#include<boost/bind.hpp>
#include<sys/time.h>
#include <ctime>
#include <Eigen/Dense>
using namespace std;

typedef Eigen::Vector2d V2d;


class DynMeans{
	public:
		DynMeans(double lambda, double Q, double tau, bool verbose = false);
		~DynMeans();

		//initialize a new step and cluster

		std::vector<int> cluster_wrapper(std::vector<double> newobservations, int nRestarts);

		//reset DDP chain
		void reset();
	private:
		double lambda, Q, tau;
		bool verbose;
		std::vector<V2d> observations;
		//during each step, constants which are information about the past steps
		//once each step is complete, these get updated
		int nextLbl;
		std::vector<V2d> oldprms;
		std::vector<int> oldprmlbls;
		std::vector<double> weights;
		std::vector<int> ages;

		void cluster(std::vector<V2d>& newobservations, int nRestarts, std::vector<int>& finalLabels, std::vector<V2d>& finalParams, double& finalObj, double& tTaken);

		//tools to help with kmeans
		std::vector<V2d> getObsInCluster(int idx, std::vector<int> lbls);
		void assignObservations(std::vector<int> assgnOrdering, std::vector<int>& lbls, std::vector<int>& cnts, std::vector<V2d>& prms);
		double setParameters(std::vector<int>& lbls, std::vector<int>& cnts, std::vector<V2d>& prms);
		std::vector<int> updateState(std::vector<int> lbls, std::vector<int> cnts, std::vector<V2d> prms);
};

#define __DYNMEANS_HPP
#endif /* __DYNMEANS_HPP */
