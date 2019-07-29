#ifndef __MAXMATCHING_IMPL_CPP

#include "maxmatching.hpp"

//THIS FUNCTION FINDS THE BEST CORRESPONDANCE BETWEEN LABELS1 AND LABELS2
//NOT *ALL* LABELS IN EITHER LABELS1 OR LABELS2 IS GUARANTEED TO BE USED
//THE MAP WILL NOT CONTAIN A KEY FOR ANY UNUSED LABEL IN LABELS1
//THE MAP WILL NOT MAP ANY KEY TO AN UNUSED LABEL IN LABELS2

//FURTHER, IF ANY LABEL IN LABELS1 OR LABELS2 IS < 0, THOSE ELEMENTS
//ARE REMOVED FROM THE PROCEEDINGS AT THE START (IT IS ASSUMED THAT NEGATIVE LABELS
//STAND FOR *UNKNOWN* / *UNLABELLED* DATA THAT IS "WRONG" BY DEFAULT)
map<int, int>
getMaxMatching(vector<int> labels1, vector<int> labels2){
	if (labels1.size() != labels2.size() || labels1.size() == 0){
		cout << "Error: labels have invalid size for getMaxMatching." << endl;
		cout << "Labels1 size: " << labels1.size() << " Labels2 size: " << labels2.size() << endl;
		int ddd;
		cin >> ddd;
	}
	//prune data with negative labels
	for (int i = 0; i < labels1.size(); i++){
		if(labels1[i] < 0 || labels2[i] < 0){
			labels1.erase(labels1.begin()+i);
			labels2.erase(labels2.begin()+i);
			i--;
		}
	}
	if (labels1.empty() || labels2.empty()){
		return map<int,int>();
	}
	/*cout << "Labels 1:"<< endl;
	for (int i =0 ; i < labels1.size(); i++){
		cout << labels1[i] << " ";
	}
	cout << endl;
	cout << "Labels 2:"<< endl;
	for (int i =0 ; i < labels2.size(); i++){
		cout << labels2[i] << " ";
	}
	cout << endl;*/
	//first create ordered sets of labels1 and labels2
	//this pulls out the unique labels
	set<int> l1set, l2set;
	for (int i = 0; i < labels1.size(); i++){
		l1set.insert(labels1[i]);
	}
	for (int i = 0; i < labels2.size(); i++){
		l2set.insert(labels2[i]);
	}
	/*
	cout << "Unique Labels 1" << endl;
	for (set<int>::iterator it = l1set.begin(); it != l1set.end(); it++){
		cout << *it << " ";
	}
	cout << endl;
	cout << "Unique Labels 2" << endl;
	for (set<int>::iterator it = l2set.begin(); it != l2set.end(); it++){
		cout << *it << " ";
	}
	cout << endl;*/
	//now create maps from the label indices in l1 and l2 to variable indices in the LP
	//and simultaneously create the weight map
	map<pair<int, int>, int> varMap;
	map<int, pair<int, int> > invvarMap;
	map<int, double> weightMap;
	int k = 1;
	for (set<int>::iterator it1 = l1set.begin(); it1 != l1set.end(); it1++){
	for (set<int>::iterator it2 = l2set.begin(); it2 != l2set.end(); it2++){
			//variable maps
			varMap[pair<int, int>(*it1, *it2)] =  k;
			invvarMap[k] = pair<int,int>(*it1, *it2);
			//weight map
			//increment the weight of this variable (graph edge) based on the number
			//of matching observations between the two
			double weight = 0;
			for (int i = 0; i < labels1.size(); i++){
				if (labels1[i] == *it1 && labels2[i] == *it2){
					weight++;
				}
			}
			weightMap[k] = weight;

			//cout << "Var map (" << *it1 << "," << *it2 << ") = " << k << ", weight map (" << k << ") = " << weight << endl;
			k++;
	}
	}

	//construct the linear program
	lprec *lp;
	lp = make_lp(0, varMap.size()); //there are l1Map.size()*l2Map.size() variables to optimize over
	if (lp == NULL){
		cout << "Error: Couldn't create linear program." << endl;
		int ddd;
		cin >> ddd;
	}

	//add constraints
	set_add_rowmode(lp, true);
	int* varnos = new int[varMap.size()];
	double* varweights = new double[varMap.size()];
	//constrant type 1: the sum of outgoing edges from each of the A vertices = 1
	for (set<int>::iterator it1 = l1set.begin(); it1 != l1set.end(); it1++){
		int j = 0;
		for (set<int>::iterator it2 = l2set.begin(); it2 != l2set.end(); it2++){
			varweights[j] = 1;
			varnos[j] = varMap[pair<int,int>(*it1, *it2)];
			j++;
		}

		if(!add_constraintex(lp, j, varweights, varnos, LE, 1)){
			cout << "Error adding constraint in LP" << endl;
			cout << "During sum over edges connecting to L1 label no. " << *it1 << endl;
			int ddd;
			cin >> ddd;
		}
	}
	//constraint type 2: the sum of incoming edges to each fo the B vertices = 1
	for (set<int>::iterator it2 = l2set.begin(); it2 != l2set.end(); it2++){
		int j = 0;
		for (set<int>::iterator it1 = l1set.begin(); it1 != l1set.end(); it1++){
			varweights[j] = 1;
			varnos[j] = varMap[pair<int,int>(*it1, *it2)];
			j++;
		}

		if(!add_constraintex(lp, j, varweights, varnos, LE, 1)){
			cout << "Error adding constraint in LP" << endl;
			cout << "During sum over edges connecting to L2 label no. " << *it2 << endl;
			int ddd;
			cin >> ddd;
		}
	}
	//add objective and set to maximization mode
	set_add_rowmode(lp, false);
	int j = 0;
	for (map<pair<int, int>, int>::iterator it = varMap.begin(); it != varMap.end(); it++){
		varnos[j] = it->second;
		varweights[j] = weightMap[it->second];
		j++;
	}
	if(!set_obj_fnex(lp, j, varweights, varnos)){
			cout << "Error setting objective in LP" << endl;
			int ddd;
			cin >> ddd;
	}
	set_maxim(lp);

	//write_LP(lp, stdout);//print the LP to stdout

	set_verbose(lp, 3);//only display important messages (warnings and errors)
	//solve
	if (solve(lp) != OPTIMAL){
		cout << "Error solving maximization problem." << endl;
		cout << "Could not find maximum." << endl;
		int ddd;
		cin >> ddd;
	}

	//double objective = get_objective(lp);//unused
	get_variables(lp, varweights);

	/*cout << "raw output: " << endl;
	for (int i = 0; i <invvarMap.size(); i++){
		cout << varweights[i] << endl;
	}*/
	//create the output
	map<int, int> retMap;
	set<int> usedL1Labels, usedL2Labels;
	for (int i = 1; i < invvarMap.size()+1; i++){ //we have to use weird indexing here due to lp_solve using 1-indexing like an idiot
		if (fabs(varweights[i-1]-1.0) < 1e-6){
			//ensure no duplicate L1 labels
			if(usedL1Labels.find(invvarMap[i].first) != usedL1Labels.end()){
				cout << "Error: Tried to map a label in L1 to multiple labels in L2." << endl;
				cout << "L1 label: " << invvarMap[i].first
						  << " L2 label1: " << retMap[invvarMap[i].first]
						  << " replacing with L2 label2: " << invvarMap[i].second << endl;
				int ddd;
				cin >> ddd;
			}
			//ensure no duplicate L2 labels
			if(usedL2Labels.find(invvarMap[i].second) != usedL2Labels.end()){
				cout << "Error: Tried to map a label in L2 to multiple labels in L1." << endl;
				int l1label1 = -1;
				for (map<int, int>::iterator it = retMap.begin(); it != retMap.end(); it++){
					if (it->second == invvarMap[i].second){
						l1label1 = it->first;
						break;
					}
				}
				cout << "L2 label: " << invvarMap[i].second
						  << " L1 label1: " << l1label1
						  << " replacing with L1 label2: " << invvarMap[i].first << endl;
				int ddd;
				cin >> ddd;
			}
			//cout<<invvarMap[i].first<<endl;
			retMap[invvarMap[i].first] = invvarMap[i].second;
			usedL1Labels.insert(invvarMap[i].first);
			usedL2Labels.insert(invvarMap[i].second);
		}
	    else
	    {
	        retMap[invvarMap[i].first] = -1e2;
	        //cout<<"varweights[i-1]-1.0  "<< varweights[i-1]-1.0 <<endl;
	    }

	}

	//clean up
	delete_lp(lp);
	delete[] varnos;
	delete[] varweights;
	return retMap;
}


map<int, int> getMaxMatchingConsistentWithOldMatching(vector<int> labels1, vector<int> labels2, map<int, int> oldmatchings){
	if (labels1.size() != labels2.size() || labels1.size() == 0){
		cout << "Error: labels have invalid size for getMaxMatching." << endl;
		cout << "Labels1 size: " << labels1.size() << " Labels2 size: " << labels2.size() << endl;
		int ddd;
		cin >> ddd;
	}
	//prune data with negative labels
	for (int i = 0; i < labels1.size(); i++){
		if(labels1[i] < 0 || labels2[i] < 0){
			labels1.erase(labels1.begin()+i);
			labels2.erase(labels2.begin()+i);
			i--;
		}
	}
	//prune data with labels already in the oldmatching
	for (map<int, int>::iterator it = oldmatchings.begin(); it != oldmatchings.end(); ++it){
		int l1 = it->first, l2 = it->second;
		for (int i = 0; i < labels1.size(); i++){
			if(labels1[i] == l1 || labels2[i] == l2){
				labels1.erase(labels1.begin()+i);
				labels2.erase(labels2.begin()+i);
				i--;
			}
		}
	}
	if (labels1.empty() || labels2.empty()){
		return oldmatchings;
	}
	/*cout << "Labels 1:"<< endl;
	for (int i =0 ; i < labels1.size(); i++){
		cout << labels1[i] << " ";
	}
	cout << endl;
	cout << "Labels 2:"<< endl;
	for (int i =0 ; i < labels2.size(); i++){
		cout << labels2[i] << " ";
	}
	cout << endl;*/
	//first create ordered sets of labels1 and labels2
	//this pulls out the unique labels
	set<int> l1set, l2set;
	for (int i = 0; i < labels1.size(); i++){
		l1set.insert(labels1[i]);
	}
	for (int i = 0; i < labels2.size(); i++){
		l2set.insert(labels2[i]);
	}
	/*cout << "Unique Labels 1" << endl;
	for (set<int>::iterator it = l1set.begin(); it != l1set.end(); it++){
		cout << *it << " ";
	}
	cout << endl;
	cout << "Unique Labels 2" << endl;
	for (set<int>::iterator it = l2set.begin(); it != l2set.end(); it++){
		cout << *it << " ";
	}
	cout << endl;*/
	//now create maps from the label indices in l1 and l2 to variable indices in the LP
	//and simultaneously create the weight map
	map<pair<int, int>, int> varMap;
	map<int, pair<int, int> > invvarMap;
	map<int, double> weightMap;
	int k = 1;
	for (set<int>::iterator it1 = l1set.begin(); it1 != l1set.end(); it1++){
	for (set<int>::iterator it2 = l2set.begin(); it2 != l2set.end(); it2++){
			//variable maps
			varMap[pair<int, int>(*it1, *it2)] =  k;
			invvarMap[k] = pair<int,int>(*it1, *it2);
			//weight map
			//increment the weight of this variable (graph edge) based on the number
			//of matching observations between the two
			double weight = 0;
			for (int i = 0; i < labels1.size(); i++){
				if (labels1[i] == *it1 && labels2[i] == *it2){
					weight++;
				}
			}
			weightMap[k] = weight;

			//cout << "Var map (" << *it1 << "," << *it2 << ") = " << k << ", weight map (" << k << ") = " << weight << endl;
			k++;
	}
	}

	//construct the linear program
	lprec *lp;
	lp = make_lp(0, varMap.size()); //there are l1Map.size()*l2Map.size() variables to optimize over
	if (lp == NULL){
		cout << "Error: Couldn't create linear program." << endl;
		int ddd;
		cin >> ddd;
	}

	//add constraints
	set_add_rowmode(lp, true);
	int* varnos = new int[varMap.size()];
	double* varweights = new double[varMap.size()];
	//constrant type 1: the sum of outgoing edges from each of the A vertices = 1
	for (set<int>::iterator it1 = l1set.begin(); it1 != l1set.end(); it1++){
		int j = 0;
		for (set<int>::iterator it2 = l2set.begin(); it2 != l2set.end(); it2++){
			varweights[j] = 1;
			varnos[j] = varMap[pair<int,int>(*it1, *it2)];
			j++;
		}

		if(!add_constraintex(lp, j, varweights, varnos, LE, 1)){
			cout << "Error adding constraint in LP" << endl;
			cout << "During sum over edges connecting to L1 label no. " << *it1 << endl;
			int ddd;
			cin >> ddd;
		}
	}
	//constraint type 2: the sum of incoming edges to each fo the B vertices = 1
	for (set<int>::iterator it2 = l2set.begin(); it2 != l2set.end(); it2++){
		int j = 0;
		for (set<int>::iterator it1 = l1set.begin(); it1 != l1set.end(); it1++){
			varweights[j] = 1;
			varnos[j] = varMap[pair<int,int>(*it1, *it2)];
			j++;
		}

		if(!add_constraintex(lp, j, varweights, varnos, LE, 1)){
			cout << "Error adding constraint in LP" << endl;
			cout << "During sum over edges connecting to L2 label no. " << *it2 << endl;
			int ddd;
			cin >> ddd;
		}
	}
	//add objective and set to maximization mode
	set_add_rowmode(lp, false);
	int j = 0;
	for (map<pair<int, int>, int>::iterator it = varMap.begin(); it != varMap.end(); it++){
		varnos[j] = it->second;
		varweights[j] = weightMap[it->second];
		j++;
	}
	if(!set_obj_fnex(lp, j, varweights, varnos)){
			cout << "Error setting objective in LP" << endl;
			int ddd;
			cin >> ddd;
	}
	set_maxim(lp);

	//write_LP(lp, stdout);//print the LP to stdout

	set_verbose(lp, 3);//only display important messages (warnings and errors)
	//solve
	if (solve(lp) != OPTIMAL){
		cout << "Error solving maximization problem." << endl;
		cout << "Could not find maximum." << endl;
		int ddd;
		cin >> ddd;
	}

	//double objective = get_objective(lp);//unused
	get_variables(lp, varweights);

	/*cout << "raw output: " << endl;
	for (int i = 0; i <invvarMap.size(); i++){
		cout << varweights[i] << endl;
	}*/
	//create the output
	map<int, int> resMap;
	set<int> usedL1Labels, usedL2Labels;
	for (int i = 1; i < invvarMap.size()+1; i++){ //we have to use weird indexing here due to lp_solve using 1-indexing like an idiot
		if (fabs(varweights[i-1]-1.0) < 1e-6){
			//ensure no duplicate L1 labels
			if(usedL1Labels.find(invvarMap[i].first) != usedL1Labels.end()){
				cout << "Error: Tried to map a label in L1 to multiple labels in L2." << endl;
				cout << "L1 label: " << invvarMap[i].first
						  << " L2 label1: " << resMap[invvarMap[i].first]
						  << " replacing with L2 label2: " << invvarMap[i].second << endl;
				int ddd;
				cin >> ddd;
			}
			//ensure no duplicate L2 labels
			if(usedL2Labels.find(invvarMap[i].second) != usedL2Labels.end()){
				cout << "Error: Tried to map a label in L2 to multiple labels in L1." << endl;
				int l1label1 = -1;
				for (map<int, int>::iterator it = resMap.begin(); it != resMap.end(); it++){
					if (it->second == invvarMap[i].second){
						l1label1 = it->first;
						break;
					}
				}
				cout << "L2 label: " << invvarMap[i].second
						  << " L1 label1: " << l1label1
						  << " replacing with L1 label2: " << invvarMap[i].first << endl;
				int ddd;
				cin >> ddd;
			}
			resMap[invvarMap[i].first] = invvarMap[i].second;
			usedL1Labels.insert(invvarMap[i].first);
			usedL2Labels.insert(invvarMap[i].second);
		}
	}

	//update the oldMatchings
	for (map<int, int>::iterator it = resMap.begin(); it != resMap.end(); ++it){
		oldmatchings[it->first] = it->second;
	}

	//clean up
	delete_lp(lp);
	delete[] varnos;
	delete[] varweights;

	//output the updated oldmatchings
	return oldmatchings;
}

map<int, int>
getWeightedMaxMatching(vector<int> labels1, vector<int> labels2, vector<double> weights){
	if (labels1.size() != labels2.size() || labels1.size() == 0 || weights.size() != labels2.size()){
		cout << "Error: labels/weights have invalid size for getMaxMatching." << endl;
		cout << "Labels1 size: " << labels1.size() << " Labels2 size: " << labels2.size() << " Weights size: " << weights.size() << endl;
		int ddd;
		cin >> ddd;
	}
	//prune data with negative labels
	for (int i = 0; i < labels1.size(); i++){
		if(labels1[i] < 0 || labels2[i] < 0){
			labels1.erase(labels1.begin()+i);
			labels2.erase(labels2.begin()+i);
			i--;
		}
	}
	if (labels1.empty() || labels2.empty()){
		return map<int,int>();
	}
	/*cout << "Labels 1:"<< endl;
	for (int i =0 ; i < labels1.size(); i++){
		cout << labels1[i] << " ";
	}
	cout << endl;
	cout << "Labels 2:"<< endl;
	for (int i =0 ; i < labels2.size(); i++){
		cout << labels2[i] << " ";
	}
	cout << endl;*/
	//first create ordered sets of labels1 and labels2
	//this pulls out the unique labels
	set<int> l1set, l2set;
	for (int i = 0; i < labels1.size(); i++){
		l1set.insert(labels1[i]);
	}
	for (int i = 0; i < labels2.size(); i++){
		l2set.insert(labels2[i]);
	}
	/*cout << "Unique Labels 1" << endl;
	for (set<int>::iterator it = l1set.begin(); it != l1set.end(); it++){
		cout << *it << " ";
	}
	cout << endl;
	cout << "Unique Labels 2" << endl;
	for (set<int>::iterator it = l2set.begin(); it != l2set.end(); it++){
		cout << *it << " ";
	}
	cout << endl;*/
	//now create maps from the label indices in l1 and l2 to variable indices in the LP
	//and simultaneously create the weight map
	map<pair<int, int>, int> varMap;
	map<int, pair<int, int> > invvarMap;
	map<int, double> weightMap;
	int k = 1;
	for (set<int>::iterator it1 = l1set.begin(); it1 != l1set.end(); it1++){
	for (set<int>::iterator it2 = l2set.begin(); it2 != l2set.end(); it2++){
			//variable maps
			varMap[pair<int, int>(*it1, *it2)] =  k;
			invvarMap[k] = pair<int,int>(*it1, *it2);
			//weight map
			//increment the weight of this variable (graph edge) based on the number
			//of matching observations between the two
			double weight = 0;
			for (int i = 0; i < labels1.size(); i++){
				if (labels1[i] == *it1 && labels2[i] == *it2){
					weight+= weights[i]; //weighted matching
				}
			}
			weightMap[k] = weight;

			//cout << "Var map (" << *it1 << "," << *it2 << ") = " << k << ", weight map (" << k << ") = " << weight << endl;
			k++;
	}
	}

	//construct the linear program
	lprec *lp;
	lp = make_lp(0, varMap.size()); //there are l1Map.size()*l2Map.size() variables to optimize over
	if (lp == NULL){
		cout << "Error: Couldn't create linear program." << endl;
		int ddd;
		cin >> ddd;
	}

	//add constraints
	set_add_rowmode(lp, true);
	int* varnos = new int[varMap.size()];
	double* varweights = new double[varMap.size()];
	//constrant type 1: the sum of outgoing edges from each of the A vertices = 1
	for (set<int>::iterator it1 = l1set.begin(); it1 != l1set.end(); it1++){
		int j = 0;
		for (set<int>::iterator it2 = l2set.begin(); it2 != l2set.end(); it2++){
			varweights[j] = 1;
			varnos[j] = varMap[pair<int,int>(*it1, *it2)];
			j++;
		}

		if(!add_constraintex(lp, j, varweights, varnos, LE, 1)){
			cout << "Error adding constraint in LP" << endl;
			cout << "During sum over edges connecting to L1 label no. " << *it1 << endl;
			int ddd;
			cin >> ddd;
		}
	}
	//constraint type 2: the sum of incoming edges to each fo the B vertices = 1
	for (set<int>::iterator it2 = l2set.begin(); it2 != l2set.end(); it2++){
		int j = 0;
		for (set<int>::iterator it1 = l1set.begin(); it1 != l1set.end(); it1++){
			varweights[j] = 1;
			varnos[j] = varMap[pair<int,int>(*it1, *it2)];
			j++;
		}

		if(!add_constraintex(lp, j, varweights, varnos, LE, 1)){
			cout << "Error adding constraint in LP" << endl;
			cout << "During sum over edges connecting to L2 label no. " << *it2 << endl;
			int ddd;
			cin >> ddd;
		}
	}
	//add objective and set to maximization mode
	set_add_rowmode(lp, false);
	int j = 0;
	for (map<pair<int, int>, int>::iterator it = varMap.begin(); it != varMap.end(); it++){
		varnos[j] = it->second;
		varweights[j] = weightMap[it->second];
		j++;
	}
	if(!set_obj_fnex(lp, j, varweights, varnos)){
			cout << "Error setting objective in LP" << endl;
			int ddd;
			cin >> ddd;
	}
	set_maxim(lp);

	//write_LP(lp, stdout);//print the LP to stdout

	set_verbose(lp, 3);//only display important messages (warnings and errors)
	//solve
	if (solve(lp) != OPTIMAL){
		cout << "Error solving maximization problem." << endl;
		cout << "Could not find maximum." << endl;
		int ddd;
		cin >> ddd;
	}

	//double objective = get_objective(lp);//unused
	get_variables(lp, varweights);

	/*cout << "raw output: " << endl;
	for (int i = 0; i <invvarMap.size(); i++){
		cout << varweights[i] << endl;
	}*/
	//create the output
	map<int, int> retMap;
	set<int> usedL1Labels, usedL2Labels;
	for (int i = 1; i < invvarMap.size()+1; i++){ //we have to use weird indexing here due to lp_solve using 1-indexing like an idiot
		if (fabs(varweights[i-1]-1.0) < 1e-6){
			//ensure no duplicate L1 labels
			if(usedL1Labels.find(invvarMap[i].first) != usedL1Labels.end()){
				cout << "Error: Tried to map a label in L1 to multiple labels in L2." << endl;
				cout << "L1 label: " << invvarMap[i].first
						  << " L2 label1: " << retMap[invvarMap[i].first]
						  << " replacing with L2 label2: " << invvarMap[i].second << endl;
				int ddd;
				cin >> ddd;
			}
			//ensure no duplicate L2 labels
			if(usedL2Labels.find(invvarMap[i].second) != usedL2Labels.end()){
				cout << "Error: Tried to map a label in L2 to multiple labels in L1." << endl;
				int l1label1 = -1;
				for (map<int, int>::iterator it = retMap.begin(); it != retMap.end(); it++){
					if (it->second == invvarMap[i].second){
						l1label1 = it->first;
						break;
					}
				}
				cout << "L2 label: " << invvarMap[i].second
						  << " L1 label1: " << l1label1
						  << " replacing with L1 label2: " << invvarMap[i].first << endl;
				int ddd;
				cin >> ddd;
			}
			retMap[invvarMap[i].first] = invvarMap[i].second;
			usedL1Labels.insert(invvarMap[i].first);
			usedL2Labels.insert(invvarMap[i].second);
		}
	}

	//clean up
	delete_lp(lp);
	delete[] varnos;
	delete[] varweights;

	return retMap;
}


#define __MAXMATCHING_IMPL_CPP
#endif
