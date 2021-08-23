#include<iostream>
using namespace std;

string chainsClust(int nbElts, double *MS, string *seqsNames, int loopMinSize=2){
	double minSimil=1.0, maxSimil=0.0, miniMaxSimil=1.0, localMax, localSMax, threshold, *maxSimils, *clustsSimil;
	int *links, *clusterID, curClust=0, clusteredElts=0;
	string clustersString="";
	bool allEqualSimils=true;

	if(loopMinSize<2){//Invalid
		loopMinSize=2;
		cout<<"Cannot form a loop with less than 2 elements. The minimum value of 2 will be used.\n";
	}

	clusterID = new int[nbElts];
		
	//Check if all sequences are identical (all similarities are equal)
	for(int i=0;i<nbElts;i++){
		for(int j=i+1;j<nbElts;j++){
			if(MS[i*nbElts+j]!=MS[2]){
				allEqualSimils=false;
				break;
			}
		}
		clusterID[i]=-1; //initialize all the clusters flags to -1 (unclustered)	
	}

	if(allEqualSimils){//all sequences are identical or equidistant and belong to the same cluster		
		for(int i=0;i<nbElts;i++){
			clusterID[i]=0;
		}
		clusteredElts=nbElts;
	}else{//Link the each element to the most similar but not identical ones
		if(loopMinSize>=nbElts){//Impossible to find chains of size higher than the number of elements
			loopMinSize=nbElts;
		}
		maxSimils = new double[(loopMinSize-1)*nbElts];
		links = new int[(loopMinSize-1)*nbElts];

		for(int i=0;i<nbElts;i++){//Initialise the first row of maxSimils[] and links[]
			maxSimils[i]=0.0;
			for(int j=0;j<nbElts;j++){
				if(MS[i*nbElts+j]!=1 && MS[i*nbElts+j]>maxSimils[i]){
					maxSimils[i]=MS[i*nbElts+j];
					links[i]=j;
				}
			}
		}
		for(int r=1; r<(loopMinSize-1); r++){
			for(int i=0;i<nbElts;i++){
				maxSimils[r*nbElts+i]=0.0;
				links[r*nbElts+i]=-1; //Initialize to -1 for the case where all the sequences are identical except 'r' only
				for(int j=0;j<nbElts;j++){
					if(MS[i*nbElts+j]<maxSimils[(r-1)*nbElts+i]  && MS[i*nbElts+j]>maxSimils[r*nbElts+i]){
						maxSimils[r*nbElts+i]=MS[i*nbElts+j];
						links[r*nbElts+i]=j;
					}
				}
			}
		}

		//Satisfy the loopMinSize for the first row of links[]
		bool satisfied=false;
		while(!satisfied){
			satisfied=true;
			for(int i=0; i<nbElts && satisfied; i++){
				int nextLink=links[i];
				for(int j=1; j<loopMinSize && satisfied; j++){
					if(nextLink==i && links[1*nbElts+i]!=-1){
						links[i]=links[1*nbElts+i];
						for(int r=2; r<(loopMinSize-1); r++){
							links[(r-1)*nbElts+i]=links[r*nbElts+i];
						}
						links[(loopMinSize-2)*nbElts+i]=-1;
						satisfied=false;
					}
					nextLink=links[nextLink];
				}
			}
		}

		//Cluster the elements based on the loops and chains formed by the previous links
		int curClustID;
		while(clusteredElts<nbElts){
			for(int i=0; i<nbElts; i++){
				if(clusterID[i]==-1){
					int j = links[i];
					if(clusterID[j]!=-1){
						clusterID[i]=clusterID[j];
						clusteredElts++;
					} else {
						int loopCounter=nbElts-clusteredElts;
						while(j!=i && clusterID[j]==-1 && loopCounter>0){
							j=links[j];
							loopCounter--;
						}
						if(j==i||loopCounter==0){
							curClustID = curClust;
							curClust++;
						} else {
							curClustID = clusterID[j];
						}
						clusterID[i] = curClustID;
						clusteredElts++;
						j = links[i];
						while(clusterID[j]==-1){
							clusterID[j] = curClustID;
							clusteredElts++;
							j=links[j];
						}
					}
				}
			}
		}

		if(maxSimils){delete [] maxSimils;} 
		if(links){delete [] links;} 
	}
	//cout<<curClust<<" clusters generated.\n";
	for(int i=0; i<curClust; i++){
		clustersString+="[";
		bool fisrtClustrElt = true;
		for(int j=0; j<nbElts; j++){
			if(clusterID[j]==i){
				if(fisrtClustrElt){
					clustersString+=seqsNames[j];
					fisrtClustrElt=false;
				} else {
					clustersString+=", ";
					clustersString+=seqsNames[j];
				}
			}
		}
		clustersString+="]\n\n";
	}

	delete [] clusterID;
	return clustersString;
}