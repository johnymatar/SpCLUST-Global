#include <iostream>
#include <cstring>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <math.h>
//#include <cstring>
//#include <conio.h>
using namespace std;

int countMotifs(int nbElts, int nd1, int nd2, double *mat, string motif){
	int count=0;
	//Check for motifs using the links existance
	if(motif=="M4" || motif=="m4"){
		for(int i=0;i<nbElts;i++){
			if((mat[i*nbElts+nd1]!=0) && (mat[i*nbElts+nd2]!=0))
				count++;
		}
	}else{//M13 motif
		for(int i=0;i<nbElts;i++){
			if(mat[i*nbElts+nd1]!=0)
				count++;
			if(mat[i*nbElts+nd2]!=0)
				count++;
		}
	}
	return count;
}

void calcWeights(int nbElts, double *w, bool useWeights, string motif){

	//Fill the upper part of the weights matrix
	for(int i=0;i<nbElts;i++){
		for(int j=i;j<nbElts;j++){
			if(j==i){
				w[i*nbElts+j]=0.0;
			} else {
				if(w[i*nbElts+j]!=0)
					if(useWeights){
						w[i*nbElts+j]=w[i*nbElts+j]*countMotifs(nbElts,i,j,w,motif);
						//cout<<w[i*nbElts+j]<<" ";
					}else{
						w[i*nbElts+j]=countMotifs(nbElts,i,j,w,motif);
						//cout<<w[i*nbElts+j]<<" ";
					}
			}
		}
	}
	//Fill the lower part of the weights matrix
	for(int i=0;i<nbElts;i++)
		for(int j=i;j<nbElts;j++){
			w[j*nbElts+i]=w[i*nbElts+j];
			//cout<<w[j*nbElts+i]<<"\t";
		}

}

void dispMat(int nbElts, double *mat){
	for(int i=0;i<nbElts;i++){
		for(int j=0;j<nbElts;j++){
		std::cout<<setw(16);
		std::cout << std::fixed;
		std::cout << std::setprecision(12);
		std::cout <<mat[i*nbElts+j];
		}
	std::cout<<endl;
	}
	cout<<endl<<endl;
}

void matLap(int nbElts, double *mat){
	double deg;
	for(int i=0;i<nbElts;i++){
		deg=0.0;
		for(int j=0;j<nbElts;j++){
			deg+=mat[i*nbElts+j];
		}
		for(int j=0;j<nbElts;j++){
			if(i==j){
				mat[i*nbElts+j]=sqrt(deg)*(deg-mat[i*nbElts+j])*sqrt(deg);
			} else {
				mat[i*nbElts+j]=sqrt(deg)*(0-mat[i*nbElts+j])*sqrt(deg);
			}
		}
	}
}

void motifs(int nbElts, double *MS, bool useWeights = true, string motif = "M4"){
	
	//Phase 1: Eliminate the links based on the threshold
	
	//Performed in the main function

	//Phase 2: Calculate the motifs-based weights
	calcWeights(nbElts, MS, useWeights, motif);

	//dispMat(nbElts,MS);

	//Phase 3: Calculate the Laplacian of the weights matrix
	matLap(nbElts,MS);
	//dispMat(nbElts,MS);
}