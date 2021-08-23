#include <cstring>
#include <algorithm>
#include "GaussianMixture.h"
#include "eigen_jacobi.h"
#include "motifs.h"
#include "edlib.h"
#include "chains.h"
#include "dbscan.h"
#include "hdbscan.h"
//From Main
#include<cstdlib>
#include<fstream>
#include<sstream>
#include<iomanip>
#include<cmath>
#include<vector>
#include<mpi.h>
#include<time.h>
using namespace std;

//Validation flags
bool profilingOn=true; //Set to true to print the timestamp at the starting and ending of each computation phase
bool internalValidationOn=true; //Set to true to print the internal validation indices of the resulting clustering
bool dispSimilStat=true; //Set to true to displays the minimum, maximum, and average similarity between the sequences 

//---------------------------------------------------Global Variables (for spclust.cpp)--------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------

bool installed=false;
int Aligned_seq_length; //Alligned sequence length
char progPath[255]=""; //Path of the program directory
string mDist; //Scoring matrix used for distance calculation
string alignMode; //Alignment mode
string fListe; //Sequences' input filename
double gapOpenI, gapExtendI; //Parameters used for distance calculation
string EVMethod="delta+"; //For the choice of the number of features: delta, log, or delta+
float EVThreshold=0.01;
string *nc; //Array that will hold the input sequences' names
string *rawSeqs; //Array that will hold the unaligned input sequences
double *MatSimil; //MatSimil will hold the similarity matrix
int nbSequences; //The number of sequences in the input fasta file initialized in similarity() function
char *buffChar1,*buffChar2; //Buffer for sequences distance calculation
int *dCount; //Nunber of distances to be calculated by each process
bool motifsParallelCalcInProgress = false; //Used to avoid MPI deadlocks while paralelly calculating the best proportion for motifs clustering
//MPI communication buffers and variables
int myid;
char *aliSeqs; //Aligned sequences buffer
int numprocs; //Nunber of processes
MPI_Status status;
double *buffDouble; //Calculated distances or similarities buffer
int taskID; //The task ID for the slave
int progEnd=0; //Exit flag
float thresholdProp=1.0; //Threshold for motifs clustering or chains clustering or DBSCAN epsilon or HDBSCAN minimum elts per cluster (can be altered by argument -tp)
string thresholdingTech="pclosest"; //Thresholding technique used to nullify some similarities in the matrix
float techParam=-1.0; //Parameter used as minimum loop size for CHAINS or epsilon for DBSCAN or minimum elements per cluster for HDBSCAN
bool outputAlig=false; //flag to output the pairewise aligned sequences
string outputPath=""; //The output path of the aligned sequences. Will be filled by truncating the output clustering file path


//--------------------------------------------------------Independent Functions----------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------

void printTimestamp(string msg){
	if(profilingOn){
		// declaring argument of time() 
		time_t my_time = time(NULL); 
  
		// ctime() used to give the present time 
		cout<<msg<<" "<<ctime(&my_time)<<flush; 
	}
}

double **alloc_2d_double(int rows, int cols) {
    /*double *data = (double *)malloc(rows*cols*sizeof(double));
    double **array= (double **)malloc(rows*sizeof(double*));
    for (int i=0; i<rows; i++)
        array[i] = &(data[cols*i]);*/

	double **array = (double **)malloc(rows*sizeof(double*));
	array[0] = (double*)malloc(cols*rows*sizeof(double));
	for(int i=1; i<rows; i++)
		array[i]=array[0]+cols*i;

    return array;
}

char tolowerchar(char in){
	if(in>='A' && in <= 'Z')
		in = in + 'a' - 'A';
	return in;
}

string tolowercstr(char *in){
	string lowres="";
	for(int i=0; i<strlen(in); i++)
		lowres+=tolowerchar(in[i]);
	return lowres;
}

string tolowerstr(string in){
	string lowres="";
	for(int i=0; i<in.length(); i++)
		lowres+=tolowerchar(in[i]);
	return lowres;
}

double distanceS(string mat, char *ali1, char *ali2, double gapOpen = -10, double gapExtend = -0.5) {
	/*
	This function returns the distance between two aligned sequences. 
	This function is a subfunction of "matriceDistances".
	*/
	//Validating the choice of distance matrix
	/*if(mat!="EDNAFULL" && mat!="BLOSUM62" && mat!="PAM250"){
		cout<<"Invalid distance matrix. Available matrices are: EDNAFULL, BLOSUM62, and PAM250.\n";
		return 999;
	}*/
	
	int l1 = strlen(ali1), l2 = strlen(ali2), k, scoreMax, dist;
	double score = 0.0;
	
	//Nucleotides and proteins comparison scores matrices
	int EDNAFULL[225][3] = { { 'A', 'A', 5 },{ 'A', 'C', -4 },{ 'A', 'B', -4 },{ 'A', 'D', -1 },{ 'A', 'G', -4 },{ 'A', 'H', -1 },{ 'A', 'K', -4 },{ 'A', 'M', 1 },{ 'A', 'N', -2 },{ 'A', 'S', -4 },{ 'A', 'R', 1 },{ 'A', 'T', -4 },{ 'A', 'W', 1 },{ 'A', 'V', -1 },{ 'A', 'Y', -4 },{ 'C', 'A', -4 },{ 'C', 'C', 5 },{ 'C', 'B', -1 },{ 'C', 'D', -4 },{ 'C', 'G', -4 },{ 'C', 'H', -1 },{ 'C', 'K', -4 },{ 'C', 'M', 1 },{ 'C', 'N', -2 },{ 'C', 'S', 1 },{ 'C', 'R', -4 },{ 'C', 'T', -4 },{ 'C', 'W', -4 },{ 'C', 'V', -1 },{ 'C', 'Y', 1 },{ 'B', 'A', -4 },{ 'B', 'C', -1 },{ 'B', 'B', -1 },{ 'B', 'D', -2 },{ 'B', 'G', -1 },{ 'B', 'H', -2 },{ 'B', 'K', -1 },{ 'B', 'M', -3 },{ 'B', 'N', -1 },{ 'B', 'S', -1 },{ 'B', 'R', -3 },{ 'B', 'T', -1 },{ 'B', 'W', -3 },{ 'B', 'V', -2 },{ 'B', 'Y', -1 },{ 'D', 'A', -1 },{ 'D', 'C', -4 },{ 'D', 'B', -2 },{ 'D', 'D', -1 },{ 'D', 'G', -1 },{ 'D', 'H', -2 },{ 'D', 'K', -1 },{ 'D', 'M', -3 },{ 'D', 'N', -1 },{ 'D', 'S', -3 },{ 'D', 'R', -1 },{ 'D', 'T', -1 },{ 'D', 'W', -1 },{ 'D', 'V', -2 },{ 'D', 'Y', -3 },{ 'G', 'A', -4 },{ 'G', 'C', -4 },{ 'G', 'B', -1 },{ 'G', 'D', -1 },{ 'G', 'G', 5 },{ 'G', 'H', -4 },{ 'G', 'K', 1 },{ 'G', 'M', -4 },{ 'G', 'N', -2 },{ 'G', 'S', 1 },{ 'G', 'R', 1 },{ 'G', 'T', -4 },{ 'G', 'W', -4 },{ 'G', 'V', -1 },{ 'G', 'Y', -4 },{ 'H', 'A', -1 },{ 'H', 'C', -1 },{ 'H', 'B', -2 },{ 'H', 'D', -2 },{ 'H', 'G', -4 },{ 'H', 'H', -1 },{ 'H', 'K', -3 },{ 'H', 'M', -1 },{ 'H', 'N', -1 },{ 'H', 'S', -3 },{ 'H', 'R', -3 },{ 'H', 'T', -1 },{ 'H', 'W', -1 },{ 'H', 'V', -2 },{ 'H', 'Y', -1 },{ 'K', 'A', -4 },{ 'K', 'C', -4 },{ 'K', 'B', -1 },{ 'K', 'D', -1 },{ 'K', 'G', 1 },{ 'K', 'H', -3 },{ 'K', 'K', -1 },{ 'K', 'M', -4 },{ 'K', 'N', -1 },{ 'K', 'S', -2 },{ 'K', 'R', -2 },{ 'K', 'T', 1 },{ 'K', 'W', -2 },{ 'K', 'V', -3 },{ 'K', 'Y', -2 },{ 'M', 'A', 1 },{ 'M', 'C', 1 },{ 'M', 'B', -3 },{ 'M', 'D', -3 },{ 'M', 'G', -4 },{ 'M', 'H', -1 },{ 'M', 'K', -4 },{ 'M', 'M', -1 },{ 'M', 'N', -1 },{ 'M', 'S', -2 },{ 'M', 'R', -2 },{ 'M', 'T', -4 },{ 'M', 'W', -2 },{ 'M', 'V', -1 },{ 'M', 'Y', -2 },{ 'N', 'A', -2 },{ 'N', 'C', -2 },{ 'N', 'B', -1 },{ 'N', 'D', -1 },{ 'N', 'G', -2 },{ 'N', 'H', -1 },{ 'N', 'K', -1 },{ 'N', 'M', -1 },{ 'N', 'N', -1 },{ 'N', 'S', -1 },{ 'N', 'R', -1 },{ 'N', 'T', -2 },{ 'N', 'W', -1 },{ 'N', 'V', -1 },{ 'N', 'Y', -1 },{ 'S', 'A', -4 },{ 'S', 'C', 1 },{ 'S', 'B', -1 },{ 'S', 'D', -3 },{ 'S', 'G', 1 },{ 'S', 'H', -3 },{ 'S', 'K', -2 },{ 'S', 'M', -2 },{ 'S', 'N', -1 },{ 'S', 'S', -1 },{ 'S', 'R', -2 },{ 'S', 'T', -4 },{ 'S', 'W', -4 },{ 'S', 'V', -1 },{ 'S', 'Y', -2 },{ 'R', 'A', 1 },{ 'R', 'C', -4 },{ 'R', 'B', -3 },{ 'R', 'D', -1 },{ 'R', 'G', 1 },{ 'R', 'H', -3 },{ 'R', 'K', -2 },{ 'R', 'M', -2 },{ 'R', 'N', -1 },{ 'R', 'S', -2 },{ 'R', 'R', -1 },{ 'R', 'T', -4 },{ 'R', 'W', -2 },{ 'R', 'V', -1 },{ 'R', 'Y', -4 },{ 'T', 'A', -4 },{ 'T', 'C', -4 },{ 'T', 'B', -1 },{ 'T', 'D', -1 },{ 'T', 'G', -4 },{ 'T', 'H', -1 },{ 'T', 'K', 1 },{ 'T', 'M', -4 },{ 'T', 'N', -2 },{ 'T', 'S', -4 },{ 'T', 'R', -4 },{ 'T', 'T', 5 },{ 'T', 'W', 1 },{ 'T', 'V', -4 },{ 'T', 'Y', 1 },{ 'W', 'A', 1 },{ 'W', 'C', -4 },{ 'W', 'B', -3 },{ 'W', 'D', -1 },{ 'W', 'G', -4 },{ 'W', 'H', -1 },{ 'W', 'K', -2 },{ 'W', 'M', -2 },{ 'W', 'N', -1 },{ 'W', 'S', -4 },{ 'W', 'R', -2 },{ 'W', 'T', 1 },{ 'W', 'W', -1 },{ 'W', 'V', -3 },{ 'W', 'Y', -2 },{ 'V', 'A', -1 },{ 'V', 'C', -1 },{ 'V', 'B', -2 },{ 'V', 'D', -2 },{ 'V', 'G', -1 },{ 'V', 'H', -2 },{ 'V', 'K', -3 },{ 'V', 'M', -1 },{ 'V', 'N', -1 },{ 'V', 'S', -1 },{ 'V', 'R', -1 },{ 'V', 'T', -4 },{ 'V', 'W', -3 },{ 'V', 'V', -1 },{ 'V', 'Y', -3 },{ 'Y', 'A', -4 },{ 'Y', 'C', 1 },{ 'Y', 'B', -1 },{ 'Y', 'D', -3 },{ 'Y', 'G', -4 },{ 'Y', 'H', -1 },{ 'Y', 'K', -2 },{ 'Y', 'M', -2 },{ 'Y', 'N', -1 },{ 'Y', 'S', -2 },{ 'Y', 'R', -4 },{ 'Y', 'T', 1 },{ 'Y', 'W', -2 },{ 'Y', 'V', -3 },{ 'Y', 'Y', -1 } };
	int PAM250[529][3]={{'A','A',2},{'R','A',-2},{'N','A',0},{'D','A',0},{'C','A',-2},{'Q','A',0},{'E','A',0},{'G','A',1},{'H','A',-1},{'I','A',-1},{'L','A',-2},{'K','A',-1},{'M','A',-1},{'F','A',-3},{'P','A',1},{'S','A',1},{'T','A',1},{'W','A',-6},{'Y','A',-3},{'V','A',0},{'B','A',0},{'Z','A',0},{'X','A',0},{'A','R',-2},{'R','R',6},{'N','R',0},{'D','R',-1},{'C','R',-4},{'Q','R',1},{'E','R',-1},{'G','R',-3},{'H','R',2},{'I','R',-2},{'L','R',-3},{'K','R',3},{'M','R',0},{'F','R',-4},{'P','R',0},{'S','R',0},{'T','R',-1},{'W','R',2},{'Y','R',-4},{'V','R',-2},{'B','R',-1},{'Z','R',0},{'X','R',-1},{'A','N',0},{'R','N',0},{'N','N',2},{'D','N',2},{'C','N',-4},{'Q','N',1},{'E','N',1},{'G','N',0},{'H','N',2},{'I','N',-2},{'L','N',-3},{'K','N',1},{'M','N',-2},{'F','N',-3},{'P','N',0},{'S','N',1},{'T','N',0},{'W','N',-4},{'Y','N',-2},{'V','N',-2},{'B','N',2},{'Z','N',1},{'X','N',0},{'A','D',0},{'R','D',-1},{'N','D',2},{'D','D',4},{'C','D',-5},{'Q','D',2},{'E','D',3},{'G','D',1},{'H','D',1},{'I','D',-2},{'L','D',-4},{'K','D',0},{'M','D',-3},{'F','D',-6},{'P','D',-1},{'S','D',0},{'T','D',0},{'W','D',-7},{'Y','D',-4},{'V','D',-2},{'B','D',3},{'Z','D',3},{'X','D',-1},{'A','C',-2},{'R','C',-4},{'N','C',-4},{'D','C',-5},{'C','C',12},{'Q','C',-5},{'E','C',-5},{'G','C',-3},{'H','C',-3},{'I','C',-2},{'L','C',-6},{'K','C',-5},{'M','C',-5},{'F','C',-4},{'P','C',-3},{'S','C',0},{'T','C',-2},{'W','C',-8},{'Y','C',0},{'V','C',-2},{'B','C',-4},{'Z','C',-5},{'X','C',-3},{'A','Q',0},{'R','Q',1},{'N','Q',1},{'D','Q',2},{'C','Q',-5},{'Q','Q',4},{'E','Q',2},{'G','Q',-1},{'H','Q',3},{'I','Q',-2},{'L','Q',-2},{'K','Q',1},{'M','Q',-1},{'F','Q',-5},{'P','Q',0},{'S','Q',-1},{'T','Q',-1},{'W','Q',-5},{'Y','Q',-4},{'V','Q',-2},{'B','Q',1},{'Z','Q',3},{'X','Q',-1},{'A','E',0},{'R','E',-1},{'N','E',1},{'D','E',3},{'C','E',-5},{'Q','E',2},{'E','E',4},{'G','E',0},{'H','E',1},{'I','E',-2},{'L','E',-3},{'K','E',0},{'M','E',-2},{'F','E',-5},{'P','E',-1},{'S','E',0},{'T','E',0},{'W','E',-7},{'Y','E',-4},{'V','E',-2},{'B','E',3},{'Z','E',3},{'X','E',-1},{'A','G',1},{'R','G',-3},{'N','G',0},{'D','G',1},{'C','G',-3},{'Q','G',-1},{'E','G',0},{'G','G',5},{'H','G',-2},{'I','G',-3},{'L','G',-4},{'K','G',-2},{'M','G',-3},{'F','G',-5},{'P','G',0},{'S','G',1},{'T','G',0},{'W','G',-7},{'Y','G',-5},{'V','G',-1},{'B','G',0},{'Z','G',0},{'X','G',-1},{'A','H',-1},{'R','H',2},{'N','H',2},{'D','H',1},{'C','H',-3},{'Q','H',3},{'E','H',1},{'G','H',-2},{'H','H',6},{'I','H',-2},{'L','H',-2},{'K','H',0},{'M','H',-2},{'F','H',-2},{'P','H',0},{'S','H',-1},{'T','H',-1},{'W','H',-3},{'Y','H',0},{'V','H',-2},{'B','H',1},{'Z','H',2},{'X','H',-1},{'A','I',-1},{'R','I',-2},{'N','I',-2},{'D','I',-2},{'C','I',-2},{'Q','I',-2},{'E','I',-2},{'G','I',-3},{'H','I',-2},{'I','I',5},{'L','I',2},{'K','I',-2},{'M','I',2},{'F','I',1},{'P','I',-2},{'S','I',-1},{'T','I',0},{'W','I',-5},{'Y','I',-1},{'V','I',4},{'B','I',-2},{'Z','I',-2},{'X','I',-1},{'A','L',-2},{'R','L',-3},{'N','L',-3},{'D','L',-4},{'C','L',-6},{'Q','L',-2},{'E','L',-3},{'G','L',-4},{'H','L',-2},{'I','L',2},{'L','L',6},{'K','L',-3},{'M','L',4},{'F','L',2},{'P','L',-3},{'S','L',-3},{'T','L',-2},{'W','L',-2},{'Y','L',-1},{'V','L',2},{'B','L',-3},{'Z','L',-3},{'X','L',-1},{'A','K',-1},{'R','K',3},{'N','K',1},{'D','K',0},{'C','K',-5},{'Q','K',1},{'E','K',0},{'G','K',-2},{'H','K',0},{'I','K',-2},{'L','K',-3},{'K','K',5},{'M','K',0},{'F','K',-5},{'P','K',-1},{'S','K',0},{'T','K',0},{'W','K',-3},{'Y','K',-4},{'V','K',-2},{'B','K',1},{'Z','K',0},{'X','K',-1},{'A','M',-1},{'R','M',0},{'N','M',-2},{'D','M',-3},{'C','M',-5},{'Q','M',-1},{'E','M',-2},{'G','M',-3},{'H','M',-2},{'I','M',2},{'L','M',4},{'K','M',0},{'M','M',6},{'F','M',0},{'P','M',-2},{'S','M',-2},{'T','M',-1},{'W','M',-4},{'Y','M',-2},{'V','M',2},{'B','M',-2},{'Z','M',-2},{'X','M',-1},{'A','F',-3},{'R','F',-4},{'N','F',-3},{'D','F',-6},{'C','F',-4},{'Q','F',-5},{'E','F',-5},{'G','F',-5},{'H','F',-2},{'I','F',1},{'L','F',2},{'K','F',-5},{'M','F',0},{'F','F',9},{'P','F',-5},{'S','F',-3},{'T','F',-3},{'W','F',0},{'Y','F',7},{'V','F',-1},{'B','F',-4},{'Z','F',-5},{'X','F',-2},{'A','P',1},{'R','P',0},{'N','P',0},{'D','P',-1},{'C','P',-3},{'Q','P',0},{'E','P',-1},{'G','P',0},{'H','P',0},{'I','P',-2},{'L','P',-3},{'K','P',-1},{'M','P',-2},{'F','P',-5},{'P','P',6},{'S','P',1},{'T','P',0},{'W','P',-6},{'Y','P',-5},{'V','P',-1},{'B','P',-1},{'Z','P',0},{'X','P',-1},{'A','S',1},{'R','S',0},{'N','S',1},{'D','S',0},{'C','S',0},{'Q','S',-1},{'E','S',0},{'G','S',1},{'H','S',-1},{'I','S',-1},{'L','S',-3},{'K','S',0},{'M','S',-2},{'F','S',-3},{'P','S',1},{'S','S',2},{'T','S',1},{'W','S',-2},{'Y','S',-3},{'V','S',-1},{'B','S',0},{'Z','S',0},{'X','S',0},{'A','T',1},{'R','T',-1},{'N','T',0},{'D','T',0},{'C','T',-2},{'Q','T',-1},{'E','T',0},{'G','T',0},{'H','T',-1},{'I','T',0},{'L','T',-2},{'K','T',0},{'M','T',-1},{'F','T',-3},{'P','T',0},{'S','T',1},{'T','T',3},{'W','T',-5},{'Y','T',-3},{'V','T',0},{'B','T',0},{'Z','T',-1},{'X','T',0},{'A','W',-6},{'R','W',2},{'N','W',-4},{'D','W',-7},{'C','W',-8},{'Q','W',-5},{'E','W',-7},{'G','W',-7},{'H','W',-3},{'I','W',-5},{'L','W',-2},{'K','W',-3},{'M','W',-4},{'F','W',0},{'P','W',-6},{'S','W',-2},{'T','W',-5},{'W','W',17},{'Y','W',0},{'V','W',-6},{'B','W',-5},{'Z','W',-6},{'X','W',-4},{'A','Y',-3},{'R','Y',-4},{'N','Y',-2},{'D','Y',-4},{'C','Y',0},{'Q','Y',-4},{'E','Y',-4},{'G','Y',-5},{'H','Y',0},{'I','Y',-1},{'L','Y',-1},{'K','Y',-4},{'M','Y',-2},{'F','Y',7},{'P','Y',-5},{'S','Y',-3},{'T','Y',-3},{'W','Y',0},{'Y','Y',10},{'V','Y',-2},{'B','Y',-3},{'Z','Y',-4},{'X','Y',-2},{'A','V',0},{'R','V',-2},{'N','V',-2},{'D','V',-2},{'C','V',-2},{'Q','V',-2},{'E','V',-2},{'G','V',-1},{'H','V',-2},{'I','V',4},{'L','V',2},{'K','V',-2},{'M','V',2},{'F','V',-1},{'P','V',-1},{'S','V',-1},{'T','V',0},{'W','V',-6},{'Y','V',-2},{'V','V',4},{'B','V',-2},{'Z','V',-2},{'X','V',-1},{'A','B',0},{'R','B',-1},{'N','B',2},{'D','B',3},{'C','B',-4},{'Q','B',1},{'E','B',3},{'G','B',0},{'H','B',1},{'I','B',-2},{'L','B',-3},{'K','B',1},{'M','B',-2},{'F','B',-4},{'P','B',-1},{'S','B',0},{'T','B',0},{'W','B',-5},{'Y','B',-3},{'V','B',-2},{'B','B',3},{'Z','B',2},{'X','B',-1},{'A','Z',0},{'R','Z',0},{'N','Z',1},{'D','Z',3},{'C','Z',-5},{'Q','Z',3},{'E','Z',3},{'G','Z',0},{'H','Z',2},{'I','Z',-2},{'L','Z',-3},{'K','Z',0},{'M','Z',-2},{'F','Z',-5},{'P','Z',0},{'S','Z',0},{'T','Z',-1},{'W','Z',-6},{'Y','Z',-4},{'V','Z',-2},{'B','Z',2},{'Z','Z',3},{'X','Z',-1},{'A','X',0},{'R','X',-1},{'N','X',0},{'D','X',-1},{'C','X',-3},{'Q','X',-1},{'E','X',-1},{'G','X',-1},{'H','X',-1},{'I','X',-1},{'L','X',-1},{'K','X',-1},{'M','X',-1},{'F','X',-2},{'P','X',-1},{'S','X',0},{'T','X',0},{'W','X',-4},{'Y','X',-2},{'V','X',-1},{'B','X',-1},{'Z','X',-1},{'X','X',-1}};
	int BLOSUM62[400][3]={{'C','C',9},{'S','C',-1},{'T','C',-1},{'P','C',-3},{'A','C',0},{'G','C',-3},{'N','C',-3},{'D','C',-3},{'E','C',-4},{'Q','C',-3},{'H','C',-3},{'R','C',-3},{'K','C',-3},{'M','C',-1},{'I','C',-1},{'L','C',-1},{'V','C',-1},{'F','C',-2},{'Y','C',-2},{'W','C',-2},{'C','S',-1},{'S','S',4},{'T','S',1},{'P','S',-1},{'A','S',1},{'G','S',0},{'N','S',1},{'D','S',0},{'E','S',0},{'Q','S',0},{'H','S',-1},{'R','S',-1},{'K','S',0},{'M','S',-1},{'I','S',-2},{'L','S',-2},{'V','S',-2},{'F','S',-2},{'Y','S',-2},{'W','S',-3},{'C','T',-1},{'S','T',1},{'T','T',4},{'P','T',1},{'A','T',-1},{'G','T',1},{'N','T',0},{'D','T',1},{'E','T',0},{'Q','T',0},{'H','T',0},{'R','T',-1},{'K','T',0},{'M','T',-1},{'I','T',-2},{'L','T',-2},{'V','T',-2},{'F','T',-2},{'Y','T',-2},{'W','T',-3},{'C','P',-3},{'S','P',-1},{'T','P',1},{'P','P',7},{'A','P',-1},{'G','P',-2},{'N','P',-1},{'D','P',-1},{'E','P',-1},{'Q','P',-1},{'H','P',-2},{'R','P',-2},{'K','P',-1},{'M','P',-2},{'I','P',-3},{'L','P',-3},{'V','P',-2},{'F','P',-4},{'Y','P',-3},{'W','P',-4},{'C','A',0},{'S','A',1},{'T','A',-1},{'P','A',-1},{'A','A',4},{'G','A',0},{'N','A',-1},{'D','A',-2},{'E','A',-1},{'Q','A',-1},{'H','A',-2},{'R','A',-1},{'K','A',-1},{'M','A',-1},{'I','A',-1},{'L','A',-1},{'V','A',-2},{'F','A',-2},{'Y','A',-2},{'W','A',-3},{'C','G',-3},{'S','G',0},{'T','G',1},{'P','G',-2},{'A','G',0},{'G','G',6},{'N','G',-2},{'D','G',-1},{'E','G',-2},{'Q','G',-2},{'H','G',-2},{'R','G',-2},{'K','G',-2},{'M','G',-3},{'I','G',-4},{'L','G',-4},{'V','G',0},{'F','G',-3},{'Y','G',-3},{'W','G',-2},{'C','N',-3},{'S','N',1},{'T','N',0},{'P','N',-2},{'A','N',-2},{'G','N',0},{'N','N',6},{'D','N',1},{'E','N',0},{'Q','N',0},{'H','N',-1},{'R','N',0},{'K','N',0},{'M','N',-2},{'I','N',-3},{'L','N',-3},{'V','N',-3},{'F','N',-3},{'Y','N',-2},{'W','N',-4},{'C','D',-3},{'S','D',0},{'T','D',1},{'P','D',-1},{'A','D',-2},{'G','D',-1},{'N','D',1},{'D','D',6},{'E','D',2},{'Q','D',0},{'H','D',-1},{'R','D',-2},{'K','D',-1},{'M','D',-3},{'I','D',-3},{'L','D',-4},{'V','D',-3},{'F','D',-3},{'Y','D',-3},{'W','D',-4},{'C','E',-4},{'S','E',0},{'T','E',0},{'P','E',-1},{'A','E',-1},{'G','E',-2},{'N','E',0},{'D','E',2},{'E','E',5},{'Q','E',2},{'H','E',0},{'R','E',0},{'K','E',1},{'M','E',-2},{'I','E',-3},{'L','E',-3},{'V','E',-3},{'F','E',-3},{'Y','E',-2},{'W','E',-3},{'C','Q',-3},{'S','Q',0},{'T','Q',0},{'P','Q',-1},{'A','Q',-1},{'G','Q',-2},{'N','Q',0},{'D','Q',0},{'E','Q',2},{'Q','Q',5},{'H','Q',0},{'R','Q',1},{'K','Q',1},{'M','Q',0},{'I','Q',-3},{'L','Q',-2},{'V','Q',-2},{'F','Q',-3},{'Y','Q',-1},{'W','Q',-2},{'C','H',-3},{'S','H',-1},{'T','H',0},{'P','H',-2},{'A','H',-2},{'G','H',-2},{'N','H',1},{'D','H',1},{'E','H',0},{'Q','H',0},{'H','H',8},{'R','H',0},{'K','H',-1},{'M','H',-2},{'I','H',-3},{'L','H',-3},{'V','H',-2},{'F','H',-1},{'Y','H',2},{'W','H',-2},{'C','R',-3},{'S','R',-1},{'T','R',-1},{'P','R',-2},{'A','R',-1},{'G','R',-2},{'N','R',0},{'D','R',-2},{'E','R',0},{'Q','R',1},{'H','R',0},{'R','R',5},{'K','R',2},{'M','R',-1},{'I','R',-3},{'L','R',-2},{'V','R',-3},{'F','R',-3},{'Y','R',-2},{'W','R',-3},{'C','K',-3},{'S','K',0},{'T','K',0},{'P','K',-1},{'A','K',-1},{'G','K',-2},{'N','K',0},{'D','K',-1},{'E','K',1},{'Q','K',1},{'H','K',-1},{'R','K',2},{'K','K',5},{'M','K',-1},{'I','K',-3},{'L','K',-2},{'V','K',-3},{'F','K',-3},{'Y','K',-2},{'W','K',-3},{'C','M',-1},{'S','M',-1},{'T','M',-1},{'P','M',-2},{'A','M',-1},{'G','M',-3},{'N','M',-2},{'D','M',-3},{'E','M',-2},{'Q','M',0},{'H','M',-2},{'R','M',-1},{'K','M',-1},{'M','M',5},{'I','M',1},{'L','M',2},{'V','M',-2},{'F','M',0},{'Y','M',-1},{'W','M',-1},{'C','I',-1},{'S','I',-2},{'T','I',-2},{'P','I',-3},{'A','I',-1},{'G','I',-4},{'N','I',-3},{'D','I',-3},{'E','I',-3},{'Q','I',-3},{'H','I',-3},{'R','I',-3},{'K','I',-3},{'M','I',1},{'I','I',4},{'L','I',2},{'V','I',1},{'F','I',0},{'Y','I',-1},{'W','I',-3},{'C','L',-1},{'S','L',-2},{'T','L',-2},{'P','L',-3},{'A','L',-1},{'G','L',-4},{'N','L',-3},{'D','L',-4},{'E','L',-3},{'Q','L',-2},{'H','L',-3},{'R','L',-2},{'K','L',-2},{'M','L',2},{'I','L',2},{'L','L',4},{'V','L',3},{'F','L',0},{'Y','L',-1},{'W','L',-2},{'C','V',-1},{'S','V',-2},{'T','V',-2},{'P','V',-2},{'A','V',0},{'G','V',-3},{'N','V',-3},{'D','V',-3},{'E','V',-2},{'Q','V',-2},{'H','V',-3},{'R','V',-3},{'K','V',-2},{'M','V',1},{'I','V',3},{'L','V',1},{'V','V',4},{'F','V',-1},{'Y','V',-1},{'W','V',-3},{'C','F',-2},{'S','F',-2},{'T','F',-2},{'P','F',-4},{'A','F',-2},{'G','F',-3},{'N','F',-3},{'D','F',-3},{'E','F',-3},{'Q','F',-3},{'H','F',-1},{'R','F',-3},{'K','F',-3},{'M','F',0},{'I','F',0},{'L','F',0},{'V','F',-1},{'F','F',6},{'Y','F',3},{'W','F',1},{'C','Y',-2},{'S','Y',-2},{'T','Y',-2},{'P','Y',-3},{'A','Y',-2},{'G','Y',-3},{'N','Y',-2},{'D','Y',-3},{'E','Y',-2},{'Q','Y',-1},{'H','Y',2},{'R','Y',-2},{'K','Y',-2},{'M','Y',-1},{'I','Y',-1},{'L','Y',-1},{'V','Y',-1},{'F','Y',3},{'Y','Y',7},{'W','Y',2},{'C','W',-2},{'S','W',-3},{'T','W',-3},{'P','W',-4},{'A','W',-3},{'G','W',-2},{'N','W',-4},{'D','W',-4},{'E','W',-3},{'Q','W',-2},{'H','W',-2},{'R','W',-3},{'K','W',-3},{'M','W',-1},{'I','W',-3},{'L','W',-2},{'V','W',-3},{'F','W',1},{'Y','W',2},{'W','W',11}};
	
	k = 0;
	while (k<l1) {
		if (ali1[k] == '-' && ali2[k] == '-') {
			l1--;
			l2--;
			for (int i = k; i<l1; i++)
				ali1[i] = ali1[i + 1];
			ali1[l1] = '\0';
			for (int i = k; i<l2; i++)
				ali2[i] = ali2[i + 1];
			ali2[l2] = '\0';
		}
		else { k++; }
	}

	//Calculate the maximum score based on the selected scoring matrix
	if(tolowerstr(mat)=="ednafull"){
		scoreMax = 5 * l1;
	} else if(tolowerstr(mat)=="pam250"){
		scoreMax = 17 * l1;
	} else { //blosum62
		scoreMax = 11 * l1;
	}

	//scoreMax = 5 * l1;
	for (k = 160; k>0; k--) {
		int l = 0;
		while (l<l1 - k + 1) {
			bool m1 = true, m2 = true;
			for (int i = l; i<l + k; i++) {
				if (ali1[i] != '-') {
					m1 = false;
					break;
				}
			}
			for (int i = l; i<l + k; i++) {
				if (ali2[i] != '-') {
					m2 = false;
					break;
				}
			}
			if (m1 || m2) {
				score += gapOpen + gapExtend*k;
				for (int i = l; i<l1 - k; i++) {
					ali1[i] = ali1[i + k];
				}
				l1 = l1 - k;
				ali1[l1] = '\0';
				for (int i = l; i<l2 - k; i++) {
					ali2[i] = ali2[i + k];
				}
				l2 = l2 - k;
				ali2[l2] = '\0';
			}
			else { l++; }
		}
	}
	
	if(tolowerstr(mat)=="ednafull"){
		for (k = 0; k<l1; k++) {
			dist=0;
			for (int i = 0; i<225; i++)
				if (EDNAFULL[i][0] == ali1[k] && EDNAFULL[i][1] == ali2[k]) {
					dist = EDNAFULL[i][2];
					break;
				}
			score+=dist;
		}
	}else if(tolowerstr(mat)=="pam250"){
		for (k = 0; k<l1; k++) {
			dist=-8;
			for (int i = 0; i<529; i++)
				if (PAM250[i][0] == ali1[k] && PAM250[i][1] == ali2[k]) {
					dist = PAM250[i][2];
					break;
				}
			score+=dist;
		}
	}else{ //mat=="BLOSUM62"
		for (k = 0; k<l1; k++) {
			dist=-4;
			for (int i = 0; i<400; i++)
				if (BLOSUM62[i][0] == ali1[k] && BLOSUM62[i][1] == ali2[k]) {
					dist = BLOSUM62[i][2];
					break;
				}
			score+=dist;
		}
	}
	
	return double(scoreMax - score) / scoreMax;
}

std::string exec(const char* cmd) {
	/*
	This function executes an external shell command and retrieves its output in a string.
	Input:
	-cmd : The command to execute.
	Output:
	-result : A string holding the command output.
	*/
	char buffer[128];
	std::string result = "";
#if !(defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)) //Not under Windows OS
	FILE* pipe = popen(cmd, "r");
#else
	FILE* pipe = _popen(cmd, "r");
#endif
	if (!pipe) throw std::runtime_error("popen() failed!");
	try {
		while (!feof(pipe)) {
			if (fgets(buffer, 128, pipe) != NULL)
				result += buffer;
		}
	}
	catch (...) {
#if !(defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)) //Not under Windows OS
		pclose(pipe);
#else			
		_pclose(pipe);
#endif
		throw;
	}
#if !(defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)) //Not under Windows OS
	pclose(pipe);
#else			
	_pclose(pipe);
#endif
	return result;
}

string getStringFromDouble(double num){
	ostringstream streamObj;
	streamObj << std::fixed;
	streamObj << std::setprecision(18);
	streamObj << num;
	string strObj = streamObj.str();
	return strObj;
}

void toAffinityMat(double *mat, int n, string type="rwnl"){
	if(tolowerstr(type)=="rwnl"){ //Random Walk Normalized Laplacian
		
		double deg;
		for(int l=0; l<n; l++){
			deg=0.0;
			for(int m=0; m<n; m++)
				deg+=mat[l*n+m];
			for(int m=0; m<n; m++)
				if(l==m)
					mat[l*n+m]=(deg-mat[l*n+m])/deg;
				else
					mat[l*n+m]=(0.0-mat[l*n+m])/deg;
		}
		
	} else if(tolowerstr(type)=="ul"){ //Unnormalized laplacian
		double deg;
		for(int l=0; l<n; l++){
			deg=0.0;
			for(int m=0; m<n; m++)
				deg+=mat[l*n+m];
			for(int m=0; m<n; m++)
				if(l==m)
					mat[l*n+m]=deg-mat[l*n+m];
				else
					mat[l*n+m]=0.0-mat[l*n+m];
		}
	} else if(tolowerstr(type)=="mod"){ //Modularity
		double K=0.0,*k;
		k=new double [n];
		for(int l=0; l<n; l++){ //Calculating the degrees vector k and the total degree K
			k[l]=0.0;
			for(int m=0; m<n; m++){
				k[l]+=mat[l*n+m];
				K+=mat[l*n+m];
			}
		}
		for(int l=0; l<n; l++){//M=(A-kkT/K)/K
			for(int m=0; m<n; m++){
				mat[l*n+m]=(mat[l*n+m]-(k[l]*k[m])/K)/K;
			}
		}
		delete k;
	} else if(tolowerstr(type)=="bh"){ //Bethe Hessian
		//r= (sqrt(c)) where c is the average degree of the graph
		//In our case, all the nodes are linked so the degree of each node is the number of nodes - 1
		//therefore r = sqrt(c) = sqrt(n-1) and r^2 = n-1
		double r=sqrt(n-1), *k;
		k=new double [n]; // 1 + SUM i!=j ((Wij^2)/(r^2 - Wij^2)) this constitutes the first part of the equation where i=j. k=all the nodes except i because all the nodes are neighbors to i
		for(int l=0; l<n; l++){ //Calculating the degrees vector k and the total degree K
			k[l]=1.0;
			for(int m=0; m<n; m++){
				if(m!=l){
					k[l]+=(pow(mat[l*n+m],2))/((n-1)-(pow(mat[l*n+m],2)));
				}
			}
		}
		for(int l=0; l<n; l++){//H(r)=(r^2-1)I+D-rA
			for(int m=0; m<n; m++){
				mat[l*n+m]=-(r*mat[l*n+m])/((n-1)-(pow(mat[l*n+m],2)));
				if(l==m){
					mat[l*n+m]=mat[l*n+m]+k[l];
				}
			}
		}
		delete[] k;
	} else if(tolowerstr(type)=="mot4"){ //M4 motifs clustering
		motifs(n,mat, false, "M4");
	} else if(tolowerstr(type)=="wmot4"){ //Weighted M4 motifs clustering
		motifs(n,mat, true, "M4");
	} else if(tolowerstr(type)=="mot13"){ //M13 motifs clustering
		motifs(n,mat, false, "M13");
	} else if(tolowerstr(type)=="wmot13"){ //Weighted M13 motifs clustering
		motifs(n,mat, true, "M13");
	}else { //Non-Backtracking
	}
}

bool is_number(const string& s)
{
    string::const_iterator it = s.begin();
    while (it != s.end() && isdigit(*it)) ++it;
    return !s.empty() && it == s.end();
}

int featuresCalc(double *affMat, double **&featuresMat, int& nbFeatures, int& nbClusters, int nbSamples, string EVMethod="delta", float EVCutOff=0.01, string affMatType="rwnl"){
	
	//Check the compatibility between the EVMethod and affMatType
	if((tolowerstr(affMatType)=="mod" || tolowerstr(affMatType)=="mot4" || tolowerstr(affMatType)=="mot13" || tolowerstr(affMatType)=="wmot4" || tolowerstr(affMatType)=="wmot13" || tolowerstr(affMatType)=="fmot4" || tolowerstr(affMatType)=="fmot13" || tolowerstr(affMatType)=="fwmot4" || tolowerstr(affMatType)=="fwmot13") && (EVMethod=="delta" || EVMethod=="delta+")){
		EVMethod="log";
	}

	//Calculate the Eigenvalues and Eigenvectors
	double *eigVals, *eigVects; //a holds the matrix, d holds the eigenvalues, v holds the eigenvectors
	int it_max=100; //maximum number of iterations
	int it_num; //used number of iterations
	int rot_num; //used number of rotations
		
	//Allocate space for eigenvalues
	if (!(eigVals = new double[nbSamples])){
		cout << "Allocation for 'eigVals' failed. \n";
		return 1;
	} 
	//Allocate space for eigenvectors
	if (!(eigVects = new double[nbSamples*nbSamples])){
		cout << "Allocation for 'eigVects' failed. \n";
		delete []eigVals;
		return 1;
	}
		
	jacobi_eigenvalue ( nbSamples, affMat, it_max, eigVects, eigVals, it_num, rot_num );
	//Print the eingenvalues
	/*for(int i=0; i<nbSamples; i++)
		cout<<eigVals[i]<<"  ";*/

	//Calculating the number of eigenvectors to use based on the chosen method and EVCutOff
	if(EVMethod=="delta"){
		nbFeatures=1;
		while (eigVals[nbFeatures+1] > eigVals[nbFeatures] + EVCutOff){
            nbFeatures+=1;
            if (nbFeatures == (nbSamples-2))
               break;
		}
	} else if(EVMethod=="log"){
		nbFeatures=int(log(nbSamples));
	} else if(EVMethod=="delta+"){
		nbFeatures=1;
		while (eigVals[nbFeatures+1] > eigVals[nbFeatures] + EVCutOff){
            nbFeatures+=1;
            if (nbFeatures == (nbSamples-2))
               break;
		}
		int minFeat = 2;
		if(nbFeatures<minFeat && nbSamples>=minFeat){
			nbFeatures=minFeat;
		}
	} else if(EVMethod=="bh"){//Bethe Hessian: count the number of negative eigenvectors as nbClusters and use their eigenvalues
		nbFeatures=1;
		while (eigVals[nbFeatures]<0){
			nbFeatures+=1;
            if (nbFeatures == (nbSamples-1))
               break;
		}
		nbClusters=nbFeatures;
	}
	
	//Allocating and initializing the features matrix
	featuresMat = alloc_2d_double(nbSamples,nbFeatures);
	if(featuresMat==NULL){
		cout << "Allocation for 'featuresMat' failed. \n";
		return 1;
	}

	for(int i=0; i<nbSamples; i++){
		for(int j=0; j<nbFeatures; j++){
			if(tolowerstr(affMatType)!="mod"){ //If NOT using the modularity matrix then take the eigenvectors of the smallest eigenvalues
				featuresMat[i][j]=eigVects[j*nbSamples+i];
			} else { //If using the modularity matrix then take the eigenvectors of the highest eigenvalues
				featuresMat[i][j]=eigVects[(nbSamples-1-j)*nbSamples+i];
			}
		}
	}

	return 0;
}

int getBestNbClusters(double **data, int number_data, int number_features, int rndseed, string criterion = "BIC", string type_covariance="full", int max_iterations=1000){
	int *p = new int[number_data];
	int nbClusters;
	GaussianMixture *gmm;
	double crit=INFTY;
	if (tolowerstr(criterion)=="bic"){
       for(int j = 1;j<number_data+1;j++){
			gmm = new GaussianMixture(j,number_data, number_features, type_covariance,max_iterations,0.01,0.001,rndseed);
			(*gmm).train(data,p,number_data,number_features);
			/*cout<<"For nb = "<<j<<endl;
			cout<<"ICL: "<<(*gmm).icl()<<endl;
			cout<<"BIC: "<<(*gmm).bic()<<endl<<endl;*/
			if(crit>(*gmm).bic()){
				crit=(*gmm).bic();
				nbClusters = j;
			}else{
				delete gmm;
				break;
			}
			delete gmm;
		}
	}else if (tolowerstr(criterion)=="aic"){
        for(int j = 1;j<number_data+1;j++){
			gmm = new GaussianMixture(j,number_data, number_features, type_covariance,max_iterations,0.01,0.001,rndseed);
			(*gmm).train(data,p,number_data,number_features);
			if(crit>(*gmm).aic()){
				crit=(*gmm).aic();
				nbClusters = j;
			}else{
				delete gmm;
				break;
			} 
			delete gmm;
		}    
	}else if (tolowerstr(criterion)=="icl"){
        for(int j = 1;j<number_data+1;j++){
			gmm = new GaussianMixture(j,number_data, number_features, type_covariance,max_iterations,0.01,0.001,rndseed);
			(*gmm).train(data,p,number_data,number_features);
			if(crit>(*gmm).icl()){
				crit=(*gmm).icl();
				nbClusters = j;
			}else{
				delete gmm;
				break;
			} 
			delete gmm;
		}    
	}
	delete[] p;
	return nbClusters;
}

string getClusteringString(int number_data, int nbClusters, string *refs, int *p){
	/*for (int i = 0; i < nbClusters; i++){
		clustersString +="[";
		for (int j = 0; j < number_data; j++)
			if(p[j]==i)
				clustersString+="'" + refs[j] + "', ";
		clustersString = clustersString.substr(0, clustersString.size()-2);
		clustersString +="]\r\n";
	}*/
	vector<vector<string>> clusters(nbClusters);
	string clustersString="";

	for(int i = 0;i<number_data;i++){
		clusters[int(p[i])].push_back(refs[i]);
	}

	for (int i = 0; i < nbClusters; i++){
		clustersString +="[";
		for (int j = 0; j < clusters[i].size(); j++)
			clustersString += "'" + clusters[i][j] + "', ";
		clustersString = clustersString.substr(0, clustersString.size()-2);
		clustersString +="]\r\n\r\n";
	}
	return clustersString;
}

int countReadShuffleSeqs(string data, string *&seqsNames, string *&seqsData);
int readSetGlobalParamsAP(int argc, char* argv[]){
	for(int i=1; i<argc; i+=2){//Read the output file path to be used to output the alignments. The output file with its path will be stored later
		if(strcmp(tolowercstr(argv[i]).c_str(),"-out")==0){
			if(i+1>=argc){
				cout<<"Error: No data specified for argument "<<argv[i]<<".\n\n";
				return 1;
			} else{
				outputPath=argv[i+1];
				auto pos = outputPath.rfind('/');
				if (pos != std::string::npos) {
					outputPath.erase(pos);
				}else{
					pos = outputPath.rfind('\\');
					if (pos != std::string::npos) {
						outputPath.erase(pos);
					}
				}
			}
		}
	}
	if(outputPath==""){
		outputPath+="aliSeqs.txt";
	}else{
		outputPath+="/aliSeqs.txt";
	}
	
	
	for(int i=1; i<argc; i+=2){
		if(strcmp(tolowercstr(argv[i]).c_str(),"-outputalig")==0){
			if(i+1>=argc){
				if(myid==0)
					cout<<"Error: No data specified for argument "<<argv[i]<<".\n\n";
				return 1;
			} else{
				outputAlig=(tolowerstr(argv[i+1])=="true");
			}
		}
	}
	
	for(int i=1; i<argc; i+=2){
		if(strcmp(tolowercstr(argv[i]).c_str(),"-tt")==0){
			if(i+1>=argc){
				if(myid==0)
					cout<<"Error: No data specified for argument "<<argv[i]<<".\n\n";
				return 1;
			} else{
				thresholdingTech=argv[i+1];
			}
		}
	}

	string ths;
	for(int i=1; i<argc; i+=2){
		if(strcmp(tolowercstr(argv[i]).c_str(),"-th")==0){
			if(i+1>=argc){
				if(myid==0)
					cout<<"Error: No data specified for argument "<<argv[i]<<".\n\n";
				return 1;
			} else{
				ths=argv[i+1];
				try{
					thresholdProp=stof(ths);
				}catch(exception e){
					if(myid==0)
						cout<<"Error: Non numeric value specified for argument "<<argv[i]<<".\n\n";
					return 1;
				}
			}
		}
	}

	string tps;
	for(int i=1; i<argc; i+=2){
		if(strcmp(tolowercstr(argv[i]).c_str(),"-tp")==0){
			if(i+1>=argc){
				if(myid==0)
					cout<<"Error: No data specified for argument "<<argv[i]<<".\n\n";
				return 1;
			} else{
				tps=argv[i+1];
				try{
					techParam=stof(tps);
				}catch(exception e){
					if(myid==0)
						cout<<"Error: Non numeric value specified for argument "<<argv[i]<<".\n\n";
					return 1;
				}
			}
		}
	}

	mDist="none";
	for(int i=1; i<argc; i+=2)
		if(strcmp(tolowercstr(argv[i]).c_str(),"-mdist")==0){
			if(i+1>=argc){
				if(myid==0)
					cout<<"Error: No data specified for argument "<<argv[i]<<".\n\n";
				return 1;
			} else
				mDist=argv[i+1];
		}
	gapOpenI=-10;
	for(int i=1; i<argc; i+=2)
		if(strcmp(tolowercstr(argv[i]).c_str(),"-gapopen")==0)
			if(i+1>=argc){
				if(myid==0)
					cout<<"Error: No data specified for argument "<<argv[i]<<".\n\n";
				return 1;
			} else
				gapOpenI=atof(argv[i+1]);
	
	
	gapExtendI=-0.5;
	for(int i=1; i<argc; i+=2)
		if(strcmp(tolowercstr(argv[i]).c_str(),"-gapextend")==0)
			if(i+1>=argc){
				if(myid==0)
					cout<<"Error: No data specified for argument "<<argv[i]<<".\n\n";
				return 1;
			} else
				gapExtendI=atof(argv[i+1]);
			
	alignMode = "none";
	for(int i=1; i<argc; i+=2)
		if(strcmp(tolowercstr(argv[i]).c_str(),"-alignmode")==0)
			if(i+1>=argc){
				if(myid==0)
					cout<<"Error: No data specified for argument "<<argv[i]<<".\n\n";
				return 1;
			} else
				alignMode=argv[i+1];
		
	if(!installed){
	//Setting the running path of our executables
	#if !(defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)) //Not under Windows OS
		realpath(argv[0], progPath);
	#else			
		_fullpath(progPath, argv[0], sizeof(progPath));
	#endif
		progPath[strlen(progPath)-7]='\0';
	}

	fListe = string(progPath);
	fListe += "sequences.fasta";
	for(int i=1; i<argc; i+=2)
		if(strcmp(tolowercstr(argv[i]).c_str(),"-in")==0)
			if(i+1>=argc){
				if(myid==0)
					cout<<"Error: No data specified for argument "<<argv[i]<<".\n\n";
				return 1;
			} else
				fListe=argv[i+1];

	if(argc!=1){//Validate the extension and existance of the input file then read it if alignMode is none
		if (fListe.substr(fListe.find_last_of(".") + 1) != "fasta" && fListe.substr(fListe.find_last_of(".") + 1) != "dat"){
			if(myid==0)
				cout << "Error: Invalid input filename! Please input a fasta file containing the genomes using the argument -in.\n";
			return 1;
		} else {	
			std::ifstream infile(fListe.c_str());
			if(!infile.good()){
				if(myid==0)
					cout<<"Error: The input fasta file is missing or not accessible.\n\n";
				return 1;
			} else if(alignMode=="none"){ //Count the input sequences and fill NC and rawSeqs[]
				
				// Read the data from input file
				ifstream ifs(fListe.c_str());
				string indata( (std::istreambuf_iterator<char>(ifs) ),
								   (std::istreambuf_iterator<char>()    ) );
				
				//Count the input sequences, fill NC and rawSeqs, and validate the number of sequences
				nbSequences=countReadShuffleSeqs(indata,nc,rawSeqs);
				if(nbSequences==-1){
					if(myid==0)
						cout << "Allocation for 'nc' or 'rawSeqs' failed (probably no enough memory). \n";
					return 1;
				} else if(nbSequences<3){
					if(myid==0)
						cout << "The input file either contains less than 3 sequences or contains ill-formatted data. \n";
					return 1;
				}
			}
		}
	}
	return 0;
}

int readSetParamsMP(int argc, char* argv[], string &fGroupes, string &matType, string &cTech, string &ccCriterion, string &nbcCriterion, string &nbRunsStr, string &neStopStr){

	//Setting the default values for the arguments
	fGroupes = string(progPath);
	fGroupes += "Clustering.txt";
	ccCriterion = "bestBIC";
	nbcCriterion = "BIC";
	nbRunsStr = "500";
	neStopStr = "50";
	matType = "RWNL";
	cTech = "GMM";

	//Reading the input arguments specific to the main process only
	for(int i=1; i<argc; i+=2)
		if(strcmp(tolowercstr(argv[i]).c_str(),"-out")==0)
			if(i+1>=argc){
				cout<<"Error: No data specified for argument "<<argv[i]<<".\n\n";
				return 1;
			} else
				fGroupes=argv[i+1];

	for(int i=1; i<argc; i+=2)
		if(strcmp(tolowercstr(argv[i]).c_str(),"-cccriterion")==0)
			if(i+1>=argc){
				cout<<"Error: No data specified for argument "<<argv[i]<<".\n\n";
				return 1;
			} else
				ccCriterion=argv[i+1];

	for(int i=1; i<argc; i+=2)
		if(strcmp(tolowercstr(argv[i]).c_str(),"-nbccriterion")==0)
			if(i+1>=argc){
				cout<<"Error: No data specified for argument "<<argv[i]<<".\n\n";
				return 1;
			} else
				nbcCriterion=argv[i+1];

	for(int i=1; i<argc; i+=2)
		if(strcmp(tolowercstr(argv[i]).c_str(),"-nbruns")==0)
			if(i+1>=argc){
				cout<<"Error: No data specified for argument "<<argv[i]<<".\n\n";
				return 1;
			} else
				nbRunsStr=argv[i+1];

	for(int i=1; i<argc; i+=2)
		if(strcmp(tolowercstr(argv[i]).c_str(),"-nestop")==0)
			if(i+1>=argc){
				cout<<"Error: No data specified for argument "<<argv[i]<<".\n\n";
				return 1;
			} else
				neStopStr=argv[i+1];

	for(int i=1; i<argc; i+=2)
		if(strcmp(tolowercstr(argv[i]).c_str(),"-mattype")==0)
			if(i+1>=argc){
				cout<<"Error: No data specified for argument "<<argv[i]<<".\n\n";
				return 1;
			} else
				matType=argv[i+1];

	for(int i=1; i<argc; i+=2)
		if(strcmp(tolowercstr(argv[i]).c_str(),"-ctech")==0)
			if(i+1>=argc){
				cout<<"Error: No data specified for argument "<<argv[i]<<".\n\n";
				return 1;
			} else
				cTech=argv[i+1];


	return 0;
}

int validateInputParams(string mDist, string alignMode, string fGroupes, string ccCriterion, string nbcCriterion, string nbRunsStr, int &nbRuns, string neStopStr, int &neStop, string matType, string cTech, string thresholdingTech){

	if(tolowerstr(mDist)!="ednafull" && tolowerstr(mDist)!="blosum62" && tolowerstr(mDist)!="pam250" && tolowerstr(mDist)!="none"){
		cout<<"Error: Invalid distance matrix. Available matrices are: none, EDNAFULL, BLOSUM62, and PAM250.\n";
		return 1;
	}
		
	if(tolowerstr(alignMode)!="fast" && tolowerstr(alignMode)!="moderate" && tolowerstr(alignMode)!="maxprecision" && tolowerstr(alignMode)!="done" && tolowerstr(alignMode)!="none"){
		cout<<"Error: Invalid alignmnet mode. Available modes are: none, fast, moderate, maxPrecision, and done.\n";
		return 1;
	} else if(tolowerstr(alignMode)!="none"  && tolowerstr(mDist)=="none"){
		cout<<"Error: The scoring matrix (mDist) can be 'none' only when the alignMode is 'none'.\n";
		return 1;
	}
		
	/*// Removed restriction on the output file extension
	if (fGroupes.substr(fGroupes.find_last_of(".") + 1) != "txt" && fGroupes.substr(fGroupes.find_last_of(".") + 1) != "dat"){
		cout << "Error: Invalid output filename! Please name a txt file for the results using the argument -out.\n";
		return 1;
	}
	*/
		
	if(tolowerstr(ccCriterion)!="fast" && tolowerstr(ccCriterion)!="bestbic" && tolowerstr(ccCriterion)!="bestaic" && tolowerstr(ccCriterion)!="besticl" && tolowerstr(ccCriterion)!="mostfreq" && tolowerstr(ccCriterion)!="besticl-r" && tolowerstr(ccCriterion)!="bestbic-r"){
		cout<<"Error: Invalid clustering choice criterion. Available criteria are: bestBIC, bestAIC, bestICL, mostFreq, fast, bestICL-R, and bestBIC-R.\n";
		return 1;
	}
		
	if(tolowerstr(nbcCriterion)!="bic" && tolowerstr(nbcCriterion)!="aic" && tolowerstr(nbcCriterion)!="icl"){
		cout<<"Error: Invalid number of clusters choice criterion. Available criteria are: BIC, AIC, and ICL.\n";
		return 1;
	}
		
	if(!is_number(nbRunsStr)){
		cout<<"Error: The input number of runs must be a positive integer.\n";
		return 1;
	}
	nbRuns = atoi(nbRunsStr.c_str());
		
	if(!is_number(neStopStr)){
		cout<<"Error: The input no-emprovement stop parameter must be a positive integer.\n";
		return 1;
	}
	neStop = atoi(neStopStr.c_str());
		
	if(tolowerstr(matType)!="ul" && tolowerstr(matType)!="rwnl" && tolowerstr(matType)!="mod" && tolowerstr(matType)!="bh" && tolowerstr(matType)!="mot4" && tolowerstr(matType)!="wmot4" && tolowerstr(matType)!="fmot4" && tolowerstr(matType)!="fwmot4" && tolowerstr(matType)!="mot13" && tolowerstr(matType)!="wmot13" && tolowerstr(matType)!="fmot13" && tolowerstr(matType)!="fwmot13"){
		cout<<"Error: Invalid affinity matrix choice. Available choices are: LU, RWNL, MOD, BH, MOT4, WMOT4, FMOT4, FWMOT4, MOT13, WMOT13, FMOT13, and FWMOT13.\n";
		return 1;
	}
		
	if(tolowerstr(cTech)!="gmm" && tolowerstr(cTech)!="chains" && tolowerstr(cTech)!="dbscan" && tolowerstr(cTech)!="hdbscan"){
		cout<<"Error: Invalid clustering technique choice. Available choices are: GMM, DBSCAN, HDBSCAN and CHAINS.\n";
		return 1;
	}
		
	if(tolowerstr(thresholdingTech)!="keptprop" && tolowerstr(thresholdingTech)!="propavg" && tolowerstr(thresholdingTech)!="avg" && tolowerstr(thresholdingTech)!="pclosest" && tolowerstr(thresholdingTech)!="deltas"){
		cout<<"Error: Invalid clustering technique choice. Available choices are: pClosest, AVG, PropAVG, and deltas.\n";
		return 1;
	}

	return 0;
}

void orderLabels(int *p, int nbc, int nbd){
	int *cl=new int[nbc];
	for(int i=0; i<nbc; i++)
		cl[i]=-1;
	for(int i=0; i<nbd;i++){
		for(int j=0; j<nbc; j++){
			if(cl[j]==p[i]){
				break;
			} else if(cl[j]==-1){
				cl[j]=p[i];
				break;
			}
		}
	}
	for(int i=0; i<nbd;i++){
		for(int j=0; j<nbc; j++){
			if(cl[j]==p[i]){
				p[i]=j;
				break;
			}
		}
	}
	delete[] cl;
}

double silhouette(double **coords, int *clusterNbs, int nbPts, int dim){
	double sumSi=0.0, Si, ai, bi, dij, sdij;
	int nbClstrs=0, Ci, Ck;

	for(int i=0; i<nbPts; i++)
		if(clusterNbs[i]>nbClstrs)
			nbClstrs=clusterNbs[i];

	for(int i=0; i<nbPts; i++){
		Ci=0;
		sdij=0.0;
		for(int j=0; j<nbPts; j++)
			if(clusterNbs[j]==clusterNbs[i])
				Ci++;
		if(Ci!=1){
			for(int j=0; j<nbPts; j++){
				if(i!=j && clusterNbs[j]==clusterNbs[i]){
					dij=0.0;
					for(int l=0; l<dim; l++){
						dij=dij+pow((coords[i][l]-coords[j][l]),2);
					}
					dij=sqrt(dij);
					sdij=sdij+dij;
				}
			}
			ai=(1.0/double(Ci-1))*sdij;
			
			bi=INFTY;
			for(int k=0; k<=nbClstrs; k++){
				if(k!=clusterNbs[i]){
					Ck=0;
					sdij=0.0;
					for(int j=0; j<nbPts; j++)
						if(clusterNbs[j]==k)
							Ck++;
					if(Ck>0){
						for(int j=0; j<nbPts; j++){
							if(clusterNbs[j]==k){
								dij=0.0;
								for(int l=0; l<dim; l++){
									dij=dij+pow((coords[i][l]-coords[j][l]),2);
								}
								dij=sqrt(dij);
								sdij=sdij+dij;
							}
						}
						if(bi>((1.0/double(Ck))*sdij))
							bi=(1.0/double(Ck))*sdij;
					}
				}
			}
			if(ai<bi)
				Si=1.0-(ai/bi);
			else
				Si=(bi/ai)-1.0;
			sumSi=sumSi+Si;
		}
	}
	return sumSi/double(nbPts);
}

int getMedoidIndex(int clstrNb, double **coords, int *clusterNbs, int nbPts, int dim){ //Returns the medoid index of cluster clstrNb
	double sdyxi = INFTY;
	int mIndex = -1;
	double dij, sdij;
	for(int i=0; i<nbPts; i++){
		if(clusterNbs[i]==clstrNb){
			sdij=0.0;
			for(int j=0; j<nbPts; j++){
				if(i!=j && clusterNbs[j]==clusterNbs[i]){
					dij=0.0;
					for(int l=0; l<dim; l++){
						dij=dij+pow((coords[i][l]-coords[j][l]),2);
					}
					dij=sqrt(dij);
					sdij=sdij+dij;
				}
			}
			if(sdij<sdyxi){
				sdyxi=sdij;
				mIndex=i;
			}
		}
	}
	return mIndex;
}

int getDsCentroidIndex(double **coords, int nbPts, int dim){ //Returns the centoid index of the data set
	double sdyxi = INFTY;
	int cIndex = -1;
	double dij, sdij;
	for(int i=0; i<nbPts; i++){
		sdij=0.0;
		for(int j=0; j<nbPts; j++){
			if(i!=j){
				dij=0.0;
				for(int l=0; l<dim; l++){
					dij=dij+pow((coords[i][l]-coords[j][l]),2);
				}
				dij=sqrt(dij);
				sdij=sdij+dij;
			}
		}
		if(sdij<sdyxi){
			sdyxi=sdij;
			cIndex=i;
		}
	}
	return cIndex;
}

double daviesbouldin(double **coords, int *clusterNbs, int nbPts, int dim){
	double *S, **R, DB, dmij, sdmij;
	int nbClstrs=0, realNbClstrs=0, clstrSize;

	for(int i=0; i<nbPts; i++)
		if(clusterNbs[i]>nbClstrs)
			nbClstrs=clusterNbs[i];

	S=new double[nbClstrs+1]; //nbClstrs+1 because nbClstrs hold the highest cluser number and the numbers could be starting from 0.
	R=new double*[nbClstrs+1];
	for(int i=0; i<=nbClstrs; i++)
		R[i]=new double[nbClstrs+1];
	for(int i=0; i<=nbClstrs; i++){
		int mi = getMedoidIndex(i,coords,clusterNbs,nbPts,dim);
		clstrSize=0;
		sdmij=0.0;
		if(mi!=(-1)){//Cluster i not empty;
			realNbClstrs++;//The number of non-empty clusters to take into consideration
			for(int j=0; j<nbPts; j++){
				if(clusterNbs[j]==i){
					clstrSize++;
					dmij=0.0;
					for(int l=0; l<dim; l++){
						dmij=dmij+pow((coords[mi][l]-coords[j][l]),2);
					}
					dmij=sqrt(dmij);
					sdmij=sdmij+dmij;
				}
			}
			S[i]=sdmij/double(clstrSize);
		} else {S[i]=-1;}//Cluster i is empty
	}
	for(int i=0; i<=nbClstrs; i++){
		for(int j=i; j<=nbClstrs; j++){
			if(j==i){
				R[i][j]=-INFTY; //Will never be max because we have max i!=j in the DB sum
			} else{
				int mi = getMedoidIndex(i,coords,clusterNbs,nbPts,dim);
				int mj = getMedoidIndex(j,coords,clusterNbs,nbPts,dim);
				if(mi!=-1 && mj!=-1){//Both clusters are not empty
					double dmimj=0.0;
					for(int l=0; l<dim; l++){
						dmimj=dmimj+pow((coords[mi][l]-coords[mj][l]),2);
					}
					dmimj=sqrt(dmimj);
					R[i][j]=(S[i]+S[j])/dmimj;
					R[j][i]=R[i][j];
				} else {R[i][j]=-INFTY; R[j][i]=-INFTY;} //Will never be taken into consideration when chosing max
			}
		}
	}
	double sMaxRij=0.0;
	double maxRij;
	for(int i=0; i<=nbClstrs; i++){
		maxRij=-INFTY;
		for(int j=0; j<=nbClstrs; j++)
			if(R[i][j]>maxRij)
				maxRij=R[i][j];
		if(maxRij!=-INFTY) //Omit empty clusters
			sMaxRij = sMaxRij + maxRij;
	}
	DB=sMaxRij/double(realNbClstrs);

	for(int i=0; i<=nbClstrs; i++)
		delete [] R[i];
	delete [] R;
	delete [] S;
	return DB;
}

double calinskiharabasz(double **coords, int *clusterNbs, int nbPts, int dim){
	double CH, sqrdmij, ssqrdc, SSw=0.0, tss=0.0, SSb;
	int nbClstrs=0, realNbClstrs=0, ciMedInd, dsCInd;

	for(int i=0; i<nbPts; i++)
		if(clusterNbs[i]>nbClstrs)
			nbClstrs=clusterNbs[i];

	for(int i=0; i<=nbClstrs; i++){
		ssqrdc=0.0;
		ciMedInd=getMedoidIndex(i,coords,clusterNbs,nbPts,dim);
		if(ciMedInd!=-1){ //Cluster exists (not empty)
			realNbClstrs++;
			for(int j=0; j<nbPts; j++){
				sqrdmij=0.0;
				if(clusterNbs[j]==i){
					for(int l=0; l<dim; l++){
						sqrdmij=sqrdmij+pow((coords[ciMedInd][l]-coords[j][l]),2);
					}
				}
				ssqrdc = ssqrdc + sqrdmij;
			}
		}
		SSw = SSw + ssqrdc;
	}

	dsCInd=getDsCentroidIndex(coords,nbPts,dim);
	if(dsCInd!=-1){ //Data set not empty
		for(int j=0; j<nbPts; j++){
			sqrdmij=0.0;
			for(int l=0; l<dim; l++){
				sqrdmij=sqrdmij+pow((coords[dsCInd][l]-coords[j][l]),2);
			}
			tss = tss + sqrdmij;
		}
	}
	SSb=tss-SSw;

	CH=(SSb/SSw)*(double(nbPts-realNbClstrs)/double(realNbClstrs-1));

	return CH;
}

double dunn(double **coords, int *clusterNbs, int nbPts, int dim){
	double D, minSep = INFTY, maxDmtr = -INFTY, dij;
	int nbClstrs=0;

	for(int i=0; i<nbPts; i++)
		if(clusterNbs[i]>nbClstrs)
			nbClstrs=clusterNbs[i];

	for(int i=0; i<=nbPts; i++){
		for(int j=i+1; j<nbPts; j++){
			if(clusterNbs[i]!=clusterNbs[j]){
				dij=0.0;
				for(int l=0; l<dim; l++){
					dij=dij+pow((coords[i][l]-coords[j][l]),2);
				}
				dij = sqrt(dij);
				if(minSep>dij)
					minSep=dij;
			}
		}
	}

	for(int i=0; i<=nbPts; i++){
		for(int j=i+1; j<nbPts; j++){
			if(clusterNbs[i]==clusterNbs[j]){
				dij=0.0;
				for(int l=0; l<dim; l++){
					dij=dij+pow((coords[i][l]-coords[j][l]),2);
				}
				dij = sqrt(dij);
				if(maxDmtr<dij)
					maxDmtr=dij;
			}
		}
	}

	D=minSep/maxDmtr;

	return D;
}

void printInternalValidationIndices(double **coords, int *clusterNbs, int nbPts, int dim){	
	if(internalValidationOn){
		cout<<"The mean Silhouette is: "<<silhouette(coords,clusterNbs,nbPts,dim)<<endl;
		cout<<"Davies-Bouldin index is: "<<daviesbouldin(coords,clusterNbs,nbPts,dim)<<endl;
		cout<<"Calinski-Harabasz index is: "<<calinskiharabasz(coords,clusterNbs,nbPts,dim)<<endl;
		cout<<"Dunn index is: "<<dunn(coords,clusterNbs,nbPts,dim)<<endl;
	}
}

int dummyAliLength;
double dummyDist;
double calcSimil(string seq1, string seq2, int coord1, int coord2, bool align=false, string mat="none", int &aliLength=dummyAliLength, double &dist=dummyDist, string aliPath=outputPath){
	int LevDist,seq1GapsSize=0,seq2GapsSize=0,EDNAFULLMax=5,PAM250Max=17,BLOSUM62Max=11;
	double simil,mismatchScore=0.0;
	int EDNAFULL[12][3] = { { 'A', 'W',1 },{ 'A', 'R', 1 },{ 'A', 'M', 1 },{ 'T', 'W', 1 },{ 'T', 'Y', 1 },{ 'T', 'K', 1 },{ 'G', 'S', 1 },{ 'G', 'R', 1 },{ 'G', 'K', 1 },{ 'C', 'S', 1 },{ 'C', 'Y', 1 },{ 'C', 'M', 1 }};
	int PAM250[47][3] = { { 'S', 'T',1 },{ 'S', 'P', 1 },{ 'S', 'A', 1 },{ 'S', 'G', 1 },{ 'S', 'N', 1 },{ 'T', 'A', 1 },{ 'P', 'A', 1 },{ 'A', 'G', 1 },{ 'A', 'D', 1 },{ 'G', 'D', 2 },{ 'N', 'D', 2 },{ 'N', 'E', 1 },{ 'N', 'Q', 1 },{ 'N', 'H', 2 },{ 'N', 'K', 1 },{ 'N', 'B', 2 },{ 'N', 'Z', 1 },{ 'D', 'E', 3 },{ 'D', 'Q', 2 },{ 'D', 'H', 1 },{ 'D', 'B', 3 },{ 'D', 'Z', 3 },{ 'E', 'Q', 2 },{ 'E', 'H', 1 },{ 'E', 'B', 2 },{ 'E', 'Z', 3 },{ 'Q', 'H', 3 },{ 'Q', 'R', 1 },{ 'Q', 'K', 1 },{ 'Q', 'B', 1 },{ 'Q', 'Z', 3 },{ 'H', 'R', 2 },{ 'H', 'B', 1 },{ 'H', 'Z', 2 },{ 'R', 'K', 3 },{ 'R', 'W', 2 },{ 'K', 'B', 1 },{ 'M', 'I', 2 },{ 'M', 'L', 4 },{ 'M', 'V', 2 },{ 'I', 'L', 2 },{ 'I', 'V', 4 },{ 'I', 'F', 1 },{ 'L', 'V', 2 },{ 'L', 'F', 2 },{ 'F', 'Y', 7 },{ 'B', 'Z', 2 }};
	int BLOSUM62[21][3] = { { 'A', 'S',1 },{ 'R', 'Q', 1 },{ 'R', 'K', 2 },{ 'N', 'D', 1 },{ 'N', 'H', 1 },{ 'N', 'S', 1 },{ 'D', 'E', 2 },{ 'Q', 'E', 2 },{ 'Q', 'K', 1 },{ 'E', 'K', 1 },{ 'H', 'Y', 2 },{ 'I', 'L', 2 },{ 'I', 'M', 1 },{ 'I', 'V', 3 },{ 'L', 'M', 2 },{ 'L', 'V', 1 },{ 'M', 'V', 1 },{ 'F', 'W', 1 },{ 'F', 'Y', 3 },{ 'S', 'T', 1 },{ 'W', 'Y', 2 }};
	EdlibAlignResult result = edlibAlign(seq1.c_str(), seq1.length(), seq2.c_str(), seq2.length(), edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0));
    if (result.status == EDLIB_STATUS_OK) {		
		for(int i=0; i<result.alignmentLength; i++){
			//cout<<int(result.alignment[i]);
			switch(int(result.alignment[i])){
			case 1:
				if(align){seq2.insert(i,1,'-');} else {seq2GapsSize++;}
				break;
			case 2:
				if(align){seq1.insert(i,1,'-');} else {seq1GapsSize++;}
				break;
			case 3:
				if(tolowerstr(mat)=="ednafull"){
					for (int j = 0; j<12; j++)
						if ((EDNAFULL[j][0] == toupper(seq1[i-seq1GapsSize]) && EDNAFULL[j][1] == toupper(seq2[i-seq2GapsSize]))||(EDNAFULL[j][0] == toupper(seq2[i-seq2GapsSize]) && EDNAFULL[j][1] == toupper(seq1[i-seq1GapsSize]))) {
							mismatchScore += double(EDNAFULL[j][2])/EDNAFULLMax;
							break;
						}
				} else if(tolowerstr(mat)=="pam250"){
					for (int j = 0; j<47; j++)
						if ((PAM250[j][0] == toupper(seq1[i-seq1GapsSize]) && PAM250[j][1] == toupper(seq2[i-seq2GapsSize]))||(PAM250[j][0] == toupper(seq2[i-seq2GapsSize]) && PAM250[j][1] == toupper(seq1[i-seq1GapsSize]))) {
							mismatchScore += double(PAM250[j][2])/PAM250Max;
							break;
						}
				} else if(tolowerstr(mat)=="blosum62"){
					for (int j = 0; j<21; j++)
						if ((BLOSUM62[j][0] == toupper(seq1[i-seq1GapsSize]) && BLOSUM62[j][1] == toupper(seq2[i-seq2GapsSize]))||(BLOSUM62[j][0] == toupper(seq2[i-seq2GapsSize]) && BLOSUM62[j][1] == toupper(seq1[i-seq1GapsSize]))) {
							mismatchScore += double(BLOSUM62[j][2])/BLOSUM62Max;
							break;
						}
				} //else: mat=none
			default:
				break;
			}
		}
		LevDist=result.editDistance;
		dist=LevDist-mismatchScore;
		aliLength=result.alignmentLength;
		/*int seqLength;
		if(seq1.length()<seq2.length()){
			seqLength=seq1.length();
		}else{
			seqLength=seq2.length();
		}*/
		simil=1-(dist/aliLength);
	} else { return -1; }
	edlibFreeAlignResult(result);

	if(align){
		string infos;
		infos = ">" + nc[coord1] + "\n" + seq1 + "\n>" + nc[coord2] + "\n" + seq2 + "\nThe Levenstein distance between the above sequences is: " + to_string(LevDist) + "\nTheir biological distance is: " + to_string(dist) + "\nTheir similarity index is: " + to_string(simil) + "\n\n";
		std::ofstream outfile;
		outfile.open(aliPath, std::ios_base::app); //Open file in appendmode 
		outfile << infos;
		outfile.close();
	}
	return simil;
}

int countReadShuffleSeqs(string data, string *&seqsNames, string *&seqsData){
	//Count the sequences
	int nbSeqs=0;
	for (int i = 0; i<data.size(); i++)
		if ((data[i] == '>' && data[i+1] != '>') || (data[i] == ';' && data[i+1] != ';')) //Taking into consideration the error of duplication by checking liste[i+1]
			nbSeqs++;

	//Read the seqeunces in the arrays
	if (!(seqsData = new string[nbSeqs])){
		return -1;
	}
	if (!(seqsNames = new string[nbSeqs])){
		return -1;
	}
	istringstream inFile(data);
	string line;
	int i=0;
	while(getline(inFile, line)){
		line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
		if(line=="")
			continue;
        if(line.at(0)=='>'||line.at(0)==';'){
			line.erase(0,1);
			seqsNames[i]=line.c_str();
			seqsData[i]="";
			i++;
		} else {
			seqsData[i-1]+=line;
		}
	}
	
	//Shuffle the sequences order to minimize the ordering interference on GMM
	int *arrIndex = new int[nbSeqs];
	for(int i=0;i<nbSeqs;i++)
		arrIndex[i]=i;
	customRND myrnd;
	random_shuffle(arrIndex,arrIndex+nbSeqs-1,[&](int i){return (myrnd.rnd()%i);});
	for(int i=0; i<nbSeqs; i++){
		line=seqsNames[i];
		seqsNames[i]=seqsNames[arrIndex[i]];
		seqsNames[arrIndex[i]]=line;
				
		line=seqsData[i];
		seqsData[i]=seqsData[arrIndex[i]];
		seqsData[arrIndex[i]]=line;
	}

	return nbSeqs;
}

void slaveSimilCalc(int n, int *&dCount, double *&buffDouble, string *&rawSeqs, string mDist, int myid, bool outputAlignment=outputAlig){
	//Determine the number of similarities to be calculated by each slave
	int nCalcs=((n*n)+n)/2; //Number of calculations
	for(int i=0; i<numprocs; i++){
		if(nCalcs%numprocs>i)
			dCount[i]=ceil(float(nCalcs)/float(numprocs));
		else
			dCount[i]=floor(float(nCalcs)/float(numprocs));
	}
	
	//Allocate the similarities vector of this slave
	buffDouble=new double[dCount[myid]];
	
	//Fill this slave's similarities vector
	int sInd=0; //Starting calculation item's index
	for(int i=0; i<myid; i++)
		sInd+=dCount[i];
	int cInd=0; //Calculation item's index
	for (int i = 0; i<n; i++)
		for (int j = i; j<n; j++,cInd++){
			if(cInd>=sInd && cInd<sInd+dCount[myid]){
				buffDouble[cInd-sInd]=calcSimil(rawSeqs[i],rawSeqs[j],i,j,outputAlignment,mDist);
			}
		}
}

int similarityNoAlig(double *&MatSimil, int numprocs, int nbSequences, string *&rawSeqs, string mDist, int &taskID, int &progEnd, double *&buffDouble, int *&dCount, int myid, bool outputAlig){
	//Allocate MatSimil
	if (!(MatSimil = new double[nbSequences*nbSequences])){
		return 1;
	}
	
	if(numprocs == 1){
		int aliLength;
		int maxAliSeqSize=0;
		int avgAliSeqSize=0;
		int nbAliPairs = 0;
		for (int i = 0; i<nbSequences; i++){
			for (int j = i; j<nbSequences; j++){
				MatSimil[i*nbSequences+j] = calcSimil(rawSeqs[i],rawSeqs[j],i,j,outputAlig,mDist,aliLength);
				if(outputAlig){
					if(maxAliSeqSize<aliLength){
						maxAliSeqSize=aliLength;
					}
					avgAliSeqSize+=aliLength;
					nbAliPairs++;
				}
			}
		}
		if(outputAlig){
			avgAliSeqSize=avgAliSeqSize/nbAliPairs;
			cout<<"The maximum aligned sequence size is "<<maxAliSeqSize<<" and the average aligned sequences size is "<<avgAliSeqSize<<endl;
		}
	} else {
		//Send the required instructions to slaves
		taskID=5;
		for(int i=1; i<numprocs; i++){
			MPI_Send(&progEnd,1,MPI_INT,i,99,MPI_COMM_WORLD);
			MPI_Send(&taskID,1,MPI_INT,i,98,MPI_COMM_WORLD);
		}

		//Do the required tasks from this process
		slaveSimilCalc(nbSequences, dCount, buffDouble, rawSeqs, mDist, myid);

		//Collect the results from slaves and fill the upper part of the distance matrix
		int curProc=0; //Current process' results
		int toInd=dCount[curProc];
		int curInd=0, vectInd=0;
		for (int i = 0; i<nbSequences; i++){
			for (int j = i; j<nbSequences; j++,curInd++,vectInd++){
				if(curInd==toInd){ //Get and continue with the next process' results
					curProc++;
					toInd+=dCount[curProc];
					vectInd=0;
					MPI_Recv(buffDouble,dCount[curProc],MPI_DOUBLE,curProc,curProc*100+17,MPI_COMM_WORLD,&status);
				}
				MatSimil[i*nbSequences+j] = buffDouble[vectInd];
			}
		}
		delete buffDouble; //Free memory from dynamic allocations
	}
	//Fill the lower part of the symmetric similarity matrix
	for (int i = 0; i<nbSequences; i++) 
		for (int j = 0; j<i; j++)
			MatSimil[i*nbSequences+j] = MatSimil[j*nbSequences+i];

	delete [] rawSeqs; //Free up the allocated memory

	return 0;
}

int zeroNSS(double *sMat, int dim, float thresholdProportion, string technique){
	int counter=0;//counter for the zeroed elements
	if(tolowerstr(technique)=="keptprop"){//Elimination based on a proportion of the lower-value elements
		int nbKeptElts=(((dim*dim)-dim)/2)*thresholdProportion;
		for(int i=0;i<dim;i++){
			for(int j=i+1;j<dim;j++){
				int supCount=0;
				for(int k=0;k<dim&&supCount<nbKeptElts;k++){
					for(int l=k+1;l<dim&&supCount<nbKeptElts;l++){
						if(sMat[i*dim+j]<sMat[k*dim+l]){
							supCount++;
						}
						if(supCount==nbKeptElts){
							sMat[i*dim+j]=0.0;
							counter++;
							sMat[j*dim+i]=0.0;
							counter++;
						}
					}
				}
			}
		}
	}else if((tolowerstr(technique)=="avg")||(tolowerstr(technique)=="propavg")){//Elimination of similarities under the average
		double threshold;
		//Calculate the average of the similarities
		double avgSimil=0.0;
		for(int i=0;i<dim;i++){
			for(int j=i+1;j<dim;j++){
				avgSimil+=sMat[i*dim+j];
			}		
		}
		avgSimil=avgSimil/(((dim*dim)-dim)/2);
		
		if(tolowerstr(technique)=="avg"){
			threshold=avgSimil;
		}else{
			threshold=avgSimil*thresholdProportion;
		}

		//Set the simlarities under the threshold to zero
		for(int i=0;i<dim;i++){
			for(int j=0;j<dim;j++){
				if(sMat[i*dim+j]<=threshold){
					counter++;
					sMat[i*dim+j]=0.0;
				}
			}
		}
	}else if(tolowerstr(technique)=="pclosest"){//Keep only the p portion of the closest similarities for each element
		int nbKeptEltsPR=dim*thresholdProportion;
		for(int i=0;i<dim;i++){
			for(int j=i;j<dim;j++){
				int supCount=0;
				for(int k=0;k<dim&&supCount<nbKeptEltsPR;k++){
					if(sMat[i*dim+j]<sMat[i*dim+k]){
						supCount++;
					}
					if(supCount==nbKeptEltsPR){
						sMat[i*dim+j]=0.0;
						counter++;
						sMat[j*dim+i]=0.0;
						counter++;
					}
				}
			}
		}		
	}else{//Elimination based on similarities distances (deltas) * proportion
		double *deltas, localMax;//, *maxSimils, localSMax;
		//maxSimils = new double[dim];
		deltas = new double[dim];
		//Find the highest similarity (not identical) for each elements and the distance between it and the second highest one
		for(int i=0;i<dim;i++){
			localMax=0.0;
			//localSMax=0.0;
			for(int j=0;j<dim;j++){
				if(sMat[i*dim+j]!=1 && sMat[i*dim+j]>localMax){
					//localSMax=localMax;
					localMax=sMat[i*dim+j];
				}
			}
			/*if(localMax=0.0){
				maxSimils[i]=1;
			}else{
				maxSimils[i]=localMax;
			}
			if(localSMax==0.0){
				deltas[i]=0.0;
			} else {
				deltas[i]=localMax-localSMax;
			}*/
			deltas[i]=1-localMax;
		}
	
		//Set the similarities under the the |maximum - delta * prop| to zero
		for(int i=0;i<dim;i++){
			for(int j=i+1;j<dim;j++){			
				if(sMat[i*dim+j]<1-thresholdProportion*deltas[i] || sMat[i*dim+j]<1-thresholdProportion*deltas[j]){//if(sMat[i*dim+j]<maxSimils[i]-thresholdProportion*deltas[i] || sMat[i*dim+j]<maxSimils[j]-thresholdProportion*deltas[j]){
					sMat[i*dim+j]=0.0;
					counter++;
					sMat[j*dim+i]=0.0;
					counter++;
				}
			}
		}

		delete [] deltas; //delete [] maxSimils; 
	}
	return counter;
}


//---------------------------------------Interdependant Functions (from spclust.cpp)-----------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------

void sendDataToWorkers(){ //Send the required data to workors processes
	taskID=1;
	for(int i=1; i<numprocs; i++){
		MPI_Send(&progEnd,1,MPI_INT,i,99,MPI_COMM_WORLD);
		MPI_Send(&taskID,1,MPI_INT,i,98,MPI_COMM_WORLD);
		MPI_Send(&nbSequences,1,MPI_INT,i,97,MPI_COMM_WORLD);
		MPI_Send(&Aligned_seq_length,1,MPI_INT,i,96,MPI_COMM_WORLD);
		MPI_Send(aliSeqs,nbSequences*(Aligned_seq_length+1)+1,MPI_CHAR,i,95,MPI_COMM_WORLD);
	}
}

void slaveDistCalc(int n=nbSequences){
	//Determine the number of distances to be calculated by each slave
	int nCalcs=((n*n)+n)/2; //Number of calculations
	for(int i=0; i<numprocs; i++){
		if(nCalcs%numprocs>i)
			dCount[i]=ceil(float(nCalcs)/float(numprocs));
		else
			dCount[i]=floor(float(nCalcs)/float(numprocs));
	}
	
	//Allocate the distance vector of this slave
	buffDouble=new double[dCount[myid]];
	
	//Fill this slave's distances vector
	int sInd=0; //Starting calculation item's index
	for(int i=0; i<myid; i++)
		sInd+=dCount[i];
	int cInd=0; //Calculation item's index
	for (int i = 0; i<n; i++)
		for (int j = i; j<n; j++,cInd++){
			if(cInd>=sInd && cInd<sInd+dCount[myid]){
				strncpy(buffChar1,aliSeqs+(i*(Aligned_seq_length+1)),Aligned_seq_length); //Copy the i th sequence
				buffChar1[Aligned_seq_length]='\0';
				strncpy(buffChar2,aliSeqs+(j*(Aligned_seq_length+1)),Aligned_seq_length); //Copy the j th sequence
				buffChar2[Aligned_seq_length]='\0';
				buffDouble[cInd-sInd]=distanceS(mDist,buffChar1,buffChar2,gapOpenI,gapExtendI);
			}
		}
}

void matriceDistances(string **dicoMuscle, double *matDist, int n = nbSequences) {
	/*
	This function is a subfunction of "similarity".
	Input:
	-dicoMuscle : A dictionary of sequences in which the keys are the names of the sequences.
	-nc : The names of the sequences.
	Output:
	-Matrix : A matrix of distance.
	*/
	
	Aligned_seq_length=dicoMuscle[0][1].length();
	buffChar1=new char[Aligned_seq_length+1];
	buffChar2=new char[Aligned_seq_length+1];
	if(numprocs == 1){
		for (int i = 0; i<n; i++)
			for (int j = i; j<n; j++){
				strcpy(buffChar1,dicoMuscle[i][1].c_str());
				strcpy(buffChar2,dicoMuscle[j][1].c_str());
				matDist[i*n+j] = distanceS(mDist,buffChar1,buffChar2,gapOpenI,gapExtendI);
			}
	} else {
		// Vectorize the aligned sequences
		aliSeqs=new char[nbSequences*(Aligned_seq_length+1)+1];
		strcpy(aliSeqs,dicoMuscle[0][1].c_str());
		strcat(aliSeqs,"+");
		for(int i=1; i<nbSequences; i++){
			strcat(aliSeqs,dicoMuscle[i][1].c_str());
			strcat(aliSeqs,"+");
		}
		strcat(aliSeqs,"\0");

		//Send the required data to slaves
		sendDataToWorkers();

		//Do the required calculations from this slave
		slaveDistCalc();

		//Collect the results from slaves and fill the lower part of the distance matrix
		int curProc=0; //Current process' results
		int toInd=dCount[curProc];
		int curInd=0, vectInd=0;
		for (int i = 0; i<n; i++)
			for (int j = i; j<n; j++,curInd++,vectInd++){
				if(curInd==toInd){ //Get and continue with the next process' results
					curProc++;
					toInd+=dCount[curProc];
					vectInd=0;
					MPI_Recv(buffDouble,dCount[curProc],MPI_DOUBLE,curProc,curProc*100+11,MPI_COMM_WORLD,&status);
				}
				matDist[i*n+j] = buffDouble[vectInd];
			}
	}

	//Fill the upper part of the symmetric distance matrix
	for (int i = 0; i<n; i++) 
		for (int j = 0; j<i; j++)
			matDist[i*n+j] = matDist[j*n+i];

	delete aliSeqs; delete buffDouble; delete buffChar1; delete buffChar2; //Free memory from dynamic allocations
}

void listSequences(string seqs){
	/*ifstream inFile(fListe.c_str());
	if(!inFile)
	{
		cout<<"Couldn't open input fasta file"<<endl;
		exit(1);
	}*/
	istringstream inFile(seqs);
	string line;
	int i=0;
	while(getline(inFile, line)){
		if(line=="")
			continue;
        if(line.at(0)=='>'||line.at(0)==';'){
			line.erase(0,1);
			#if !(defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)) //Not under Windows OS
				if(line[strlen(line.c_str())-1]=='\r')
					line[strlen(line.c_str())-1]='\0';
			#endif
			nc[i]=line.c_str();
			i++;
		}
	}
}

int similarity(string fListe, string alignMode="maxPrecision") {
	/*
	This function returns the similarity matrix.
	Input:
	-fListe : Name of the fasta file holding a list of sequences with their names.
	Output:
	-MatSimil : The similarity matrix (global variable).
	-nc : Sequential number with the names of the sequences (global variable).
	Output:
	-0 = success or 1 = alignment failed
	*/

	char cmd[300] = "";
	if(!installed){
	#if (defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)) //Under Windows OS
		strcat(cmd,"\"");
	#endif
		strcat(cmd,"\"");
		strcat(cmd,progPath);
		strcat(cmd,"muscle\" ");
	}else
		strcat(cmd,"muscle ");
	if(tolowerstr(alignMode)=="maxprecision")
		strcat(cmd,"-quiet -in ");
	else if(tolowerstr(alignMode)=="moderate")
		strcat(cmd,"-quiet -maxiters 4 -in ");
	else if(tolowerstr(alignMode)=="fast")
		strcat(cmd,"-quiet -maxiters 2 -in ");

	string liste="";
	
	if(tolowerstr(alignMode)!="done"){
		//Perform the alignment
		strcat(cmd,"\"");
		strcat(cmd, fListe.c_str());
		strcat(cmd,"\"");
		//cout<<cmd<<endl;
		liste = exec(cmd);
		printTimestamp("Sequences' alignment ended at: ");
		// Save the result
		/*std::ofstream out("aligned.txt");
		out << liste;
		out.close();*/
	} else {
		// Read the aligned sequences fron input file
		ifstream ifs(fListe.c_str());
		string inliste( (std::istreambuf_iterator<char>(ifs) ),
						   (std::istreambuf_iterator<char>()    ) );
		liste=inliste;
	}
	if(liste==""){
		return 1;
	}
	int index;
	string **dicoMuscle,tempGen;
	double max = 0.0;
	nbSequences=0;
	for (int i = 0; i<liste.size(); i++)
		if ((liste[i] == '>' && liste[i+1] != '>') || (liste[i] == ';' && liste[i+1] != ';')) //Taking into consideration the error of duplication by checking liste[i+1]
			nbSequences++;
	//Allocate space for pointers to rows of a matrix
	if (!(MatSimil = new double[nbSequences*nbSequences])){
		cout << "Allocation for 'MatSimil' failed. \n";
		return 3;
	}
	if (!(dicoMuscle = new string*[nbSequences])){
		cout << "Allocation for 'dicoMuscle' failed. \n";
		return 3;
	}
	if (!(nc = new string[nbSequences])){
		cout << "Allocation for 'nc' failed. \n";
		return 3;
	}

	//initialize the list of sequences references as in input file order
	listSequences(liste);
	
	//Shuffle the sequences order to minimize the ordering interference on GMM
	customRND myrnd;
	random_shuffle(nc,nc+nbSequences-1,[&](int i){return (myrnd.rnd()%i);});
	
	for (int i = 0; i<nbSequences; i++) {
		dicoMuscle[i] = new string[2];
		dicoMuscle[i][0]="";//Initialize to empty string in order to detect if the current row is used in the case of genes with duplicate name
		dicoMuscle[i][1]="";
	}
	for (int i = 0; i<liste.size();) {
		if (liste[i] == '>'||liste[i]==';') {
			tempGen="";
			index=0;
			i++;
			do {
				tempGen += liste[i];
				i++;
			} while (liste[i] != '\n'&&liste[i]!='\r');
			while(tempGen.compare(nc[index])!=0 || dicoMuscle[index][0]!=""){//strcmp(nc[index].c_str(),tempGen.c_str())!=0){
				index++;
			}
			dicoMuscle[index][0]=tempGen;
		}
		else {
			do {
				if (liste[i] != '\n' && liste[i] != '\r') {
					dicoMuscle[index][1] += liste[i];
				}
				i++;
			} while (liste[i] != '>'&&liste[i] != ';'&&i<liste.size());
		}
	}

	//Validate that the aligned sequences all have the same size
	int aligseqsize = dicoMuscle[0][1].length();
	
	for (int i = 1; i<nbSequences; i++){
		if(dicoMuscle[i][1].length()!=aligseqsize){
			//Free memory from dynamically allocated matrixes
			for (int i = 0; i<nbSequences; i++) {
				delete[] dicoMuscle[i];
			}
			delete[] dicoMuscle;
			delete[] nc;
			return 2;
		}
	}
	
	matriceDistances(dicoMuscle, MatSimil, nbSequences); //Calculates the distance matrix
	
	for (int i = 0; i<nbSequences; i++)
		for (int j = 0; j<nbSequences; j++)
			if (MatSimil[i*nbSequences+j]>max)
				max = MatSimil[i*nbSequences+j];
	for (int i = 0; i<nbSequences; i++)
		for (int j = 0; j<nbSequences; j++)
			MatSimil[i*nbSequences+j] /= max;
	for (int i = 0; i<nbSequences; i++)
		for (int j = 0; j<nbSequences; j++)
			MatSimil[i*nbSequences+j] = 1 - MatSimil[i*nbSequences+j];

	//Free memory from dynamically allocated matrixes
	for (int i = 0; i<nbSequences; i++) {
		delete[] dicoMuscle[i];
	}
	delete[] dicoMuscle;
	return 0;
}

void killWorkers(){ //Send a message to the worker processes that the work is done
	progEnd=1;
	for(int i=1; i<numprocs; i++)
		MPI_Send(&progEnd,1,MPI_INT,i,99,MPI_COMM_WORLD);
}

string GMM_Clustering(double **data, string *refs, int number_data, int number_features, string ccCriterion = "bestBIC", int nbRuns = 500, int neStop = 50, string nbcCriterion = "BIC", string type_covariance="full", int max_iterations=1000, int nbClusters=-1, double *ch=NULL){ //ch returns Calinski-Harabasz index by ref

	int *p = (int*)malloc(number_data*sizeof(int));
	int sentNbClusters;
	sentNbClusters = nbClusters; //Preserve the sent number to be sent for the workers processes
	bool calculateBestNbClusters=(nbClusters==-1); //if nbClusters != -1 then it has not been calculated and provided by other methods. It will be calculated here based on the best BIC.
	GaussianMixture *gmm;
	string clustersString="";

	if(tolowerstr(ccCriterion)=="fast"){
		
		int rndseed=320; // Fixed seed chosen randomly in order to avoid getting different clustering at each run
		if(calculateBestNbClusters)
			nbClusters = getBestNbClusters(data, number_data, number_features, rndseed, nbcCriterion, type_covariance, max_iterations);
		gmm = new GaussianMixture(nbClusters, number_data, number_features, type_covariance, max_iterations, 0.01, 0.001, rndseed);    
		(*gmm).train(data,p,number_data,number_features);
		clustersString = getClusteringString(number_data, nbClusters, refs, p);
		delete gmm;
		
	} else if(tolowerstr(ccCriterion)=="bestbic" || tolowerstr(ccCriterion)=="bestaic" || tolowerstr(ccCriterion)=="besticl"){

		int neTries=0;
		double crit, bestCrit=INFTY;
		int bestSeed;
		//int totalTries=0;
		for(int i=0; i<nbRuns && neTries<neStop; i++){ // "i" will be used as the random seed for each loop
			
			//Set the task ID to be sent to the workers processes if they exist and distribute the data to them
			if(!motifsParallelCalcInProgress){
				taskID=2;
				int ccCriterionCode, nbcCriterionCode;
				if(tolowerstr(ccCriterion)=="bestbic") //bestBIC
					ccCriterionCode=1;
				else if(tolowerstr(ccCriterion)=="bestaic") //bestAIC
					ccCriterionCode=2;
				else //bestICL
					ccCriterionCode=3;
				
				if(tolowerstr(nbcCriterion)=="bic") //BIC
					nbcCriterionCode=1;
				else if(tolowerstr(nbcCriterion)=="aic") //AIC
					nbcCriterionCode=2;
				else //ICL
					nbcCriterionCode=3;

				for(int j=1; j<numprocs; j++){
					MPI_Send(&progEnd,1,MPI_INT,j,99,MPI_COMM_WORLD);
					MPI_Send(&taskID,1,MPI_INT,j,98,MPI_COMM_WORLD);
					MPI_Send(&number_data,1,MPI_INT,j,97,MPI_COMM_WORLD);
					MPI_Send(&sentNbClusters,1,MPI_INT,j,96,MPI_COMM_WORLD);
					MPI_Send(&number_features,1,MPI_INT,j,95,MPI_COMM_WORLD);
					MPI_Send(&(data[0][0]),number_data*number_features,MPI_DOUBLE,j,94,MPI_COMM_WORLD);
					MPI_Send(&i,1,MPI_INT,j,93,MPI_COMM_WORLD);
					MPI_Send(&max_iterations,1,MPI_INT,j,92,MPI_COMM_WORLD);
					MPI_Send(&ccCriterionCode,1,MPI_INT,j,91,MPI_COMM_WORLD);
					MPI_Send(&nbcCriterionCode,1,MPI_INT,j,90,MPI_COMM_WORLD);
				}
			}
			
			//Master as worker
			if(calculateBestNbClusters)
				nbClusters = getBestNbClusters(data, number_data, number_features, i, nbcCriterion, type_covariance, max_iterations);
			gmm = new GaussianMixture(nbClusters, number_data, number_features, type_covariance, max_iterations, 0.01, 0.001, i);    
			(*gmm).train(data,p,number_data,number_features);
			
			if(tolowerstr(nbcCriterion)=="bic") //BIC
				crit=(*gmm).bic();
			else if(tolowerstr(nbcCriterion)=="aic") //AIC
				crit=(*gmm).aic();
			else //ICL
				crit=(*gmm).icl();

			if(bestCrit>crit){
				bestCrit=crit;
				neTries=1;
				bestSeed=i;
			} else {
				neTries++;
			}
			delete gmm;
			
			//Receive the data from the workers if they exist
			if(!motifsParallelCalcInProgress){
				for(int j=1; j<numprocs; j++){
					i++;
					double receivedCrit;
					MPI_Recv(&receivedCrit, 1,MPI_DOUBLE,j,j*100+12,MPI_COMM_WORLD,&status);
					if(bestCrit>receivedCrit){
						bestCrit=receivedCrit;
						neTries=1;
						bestSeed=i;
					} else {
						neTries++;
					}
				}
			}
		}
		
		if(calculateBestNbClusters)
			nbClusters = getBestNbClusters(data, number_data, number_features, bestSeed, nbcCriterion, type_covariance, max_iterations);
		gmm = new GaussianMixture(nbClusters, number_data, number_features, type_covariance, max_iterations, 0.01, 0.001, bestSeed);    
		(*gmm).train(data,p,number_data,number_features);
		clustersString = getClusteringString(number_data, nbClusters, refs, p);
		delete gmm;

	} else if(tolowerstr(ccCriterion)=="mostfreq"){ // The choice is mostFreq

		string *cStrings = new string[nbRuns];
		int *occurences = new int[nbRuns];
		double crit, *bestCrit = new double[nbRuns];
		for(int i=0;i<nbRuns;i++){
			cStrings[i]="";
		}
		int highestFreq=0, bestCluteringIndex, nbcCriterionCode;

		for(int i=0; i<nbRuns; i++){ // i will be used as the random seed for each loop
			
			//Set the task ID to be sent to the workers processes if they exist and distribute the data to them
			if(!motifsParallelCalcInProgress){
				taskID=3;
				
				if(tolowerstr(nbcCriterion)=="bic") //BIC
					nbcCriterionCode=1;
				else if(tolowerstr(nbcCriterion)=="aic") //AIC
					nbcCriterionCode=2;
				else //ICL
					nbcCriterionCode=3;

				for(int j=1; j<numprocs; j++){
					MPI_Send(&progEnd,1,MPI_INT,j,99,MPI_COMM_WORLD);
					MPI_Send(&taskID,1,MPI_INT,j,98,MPI_COMM_WORLD);
					MPI_Send(&number_data,1,MPI_INT,j,97,MPI_COMM_WORLD);
					MPI_Send(&sentNbClusters,1,MPI_INT,j,96,MPI_COMM_WORLD);
					MPI_Send(&number_features,1,MPI_INT,j,95,MPI_COMM_WORLD);
					MPI_Send(&(data[0][0]),number_data*number_features,MPI_DOUBLE,j,94,MPI_COMM_WORLD);
					MPI_Send(&i,1,MPI_INT,j,93,MPI_COMM_WORLD);
					MPI_Send(&max_iterations,1,MPI_INT,j,92,MPI_COMM_WORLD);
					MPI_Send(&nbcCriterionCode,1,MPI_INT,j,91,MPI_COMM_WORLD);
				}
			}
			
			//Master as worker
			if(calculateBestNbClusters)
				nbClusters = getBestNbClusters(data, number_data, number_features, i, nbcCriterion, type_covariance, max_iterations);
			gmm = new GaussianMixture(nbClusters, number_data, number_features, type_covariance, max_iterations, 0.01, 0.001, i);    
			(*gmm).train(data,p,number_data,number_features);
			orderLabels(p,nbClusters,number_data);
			clustersString = getClusteringString(number_data, nbClusters, refs, p);

			for(int j=0; j<nbRuns; j++){					
				if(tolowerstr(nbcCriterion)=="bic") //BIC
					crit=(*gmm).bic();
				else if(tolowerstr(nbcCriterion)=="aic") //AIC
					crit=(*gmm).aic();
				else //ICL
					crit=(*gmm).icl();

				if(clustersString==cStrings[j]){
					occurences[j]++;

					if(bestCrit[j]>crit)
						bestCrit[j]=crit;
					break;
				} else if(cStrings[j]==""){
					cStrings[j]=clustersString;
					occurences[j]=1;
					bestCrit[j]=crit;
					break;
				}
			}
			delete gmm;
			
			//Receive the data from the workers if they exist
			if(!motifsParallelCalcInProgress){
				for(int k=1; k<numprocs; k++){
					i++;
					double receivedCrit;
					MPI_Recv(&receivedCrit, 1,MPI_DOUBLE,k,k*100+13,MPI_COMM_WORLD,&status);				
					MPI_Recv(&nbClusters, 1,MPI_INT,k,k*100+14,MPI_COMM_WORLD,&status);			
					MPI_Recv(p, number_data, MPI_INT, k, k*100+15, MPI_COMM_WORLD, &status);

					orderLabels(p,nbClusters,number_data);
					clustersString = getClusteringString(number_data, nbClusters, refs, p);

					for(int j=0; j<nbRuns; j++){
						if(clustersString==cStrings[j]){
							occurences[j]++;
							if(bestCrit[j]>receivedCrit)
								bestCrit[j]=receivedCrit;
							break;
						} else if(cStrings[j]==""){
							cStrings[j]=clustersString;
							occurences[j]=1;
							bestCrit[j]=receivedCrit;
							break;
						}
					}
				}
			}
		}

		for(int i=0; i<nbRuns; i++){ // Getting the index of the best clustering: the most frequently occurent one, and if many we take the one with the best BIC
			if(cStrings[i]==""){//cout<<"stopped at "<<i<<endl;
				break;
			}else if(occurences[i]>highestFreq){
				highestFreq=occurences[i];
				bestCluteringIndex=i;
			} else if(occurences[i]==highestFreq){
				if(bestCrit[i]<bestCrit[bestCluteringIndex])
					bestCluteringIndex=i;
			}
		}//cout<<"\nThe best index is "<<bestCluteringIndex<<" and the highest freq is "<<highestFreq<<endl;

		clustersString=cStrings[bestCluteringIndex];
		delete[] cStrings; delete[] occurences; delete[] bestCrit;

	} else { //The choice is bestICL-R or bestBIC-R

		//Output the features data to text files
		ofstream features;
		ifstream clustVect;
		string tempElt;
		char fname[20]="featuresMat";
		strcat(fname,(to_string(myid)).c_str());
		strcat(fname,".txt");
		
		features.open(fname);
		features.precision(17);
		for(int i=0; i<number_data; i++){
			for(int j=0; j<number_features; j++){
				features<<data[i][j]<<" ";
			}
			features<<endl;
		}
		features.close();

		//Perform the clustering using the R script
		//exec("\"C:/Program Files/R/R-4.0.3/bin/Rscript.exe\" GMM_ICL.R");
		char cmd[300] = "Rscript \"";
		strcat(cmd,progPath);
		if(tolowerstr(ccCriterion)=="besticl-r"){ //The choice is bestICL-R
			strcat(cmd,"GMM_ICL.R\" ");
		} else { //The choice is bestBIC-R
			strcat(cmd,"GMM_BIC.R\" ");
		}
		strcat(cmd,(to_string(myid)).c_str());
		
		string res = exec(cmd);
		
		if(res!=""){ //Get sure that the R scrpit returned a result		
			istringstream isres(res);
		
			//Read the labels vector
			//clustVect.open("clustVect.txt");
			string elt;
			for(int i=0; i<number_data; i++){
				isres>>elt;
				p[i]=stoi(elt)-1; //The labels should start from 0 (not from 1 like output from mclust)
				if(p[i]>nbClusters){
					nbClusters=p[i];
				}
				//cout<<p[i]<<" ";
			}
			nbClusters++; //The highest label + 1 because the labels start from 0
			clustersString = getClusteringString(number_data, nbClusters, refs, p);
		} else {
			clustersString = "";
		}
		//Cleanup
		while(remove(fname) != 0);
	}
	
	if(ch!=NULL){//a ch different than null is sent by ref
		if(clustersString == ""){//Clustering failed
			*ch=-9999999999;
		} else {
			*ch=calinskiharabasz(data,p,number_data,number_features);
		}
	}
	if(clustersString != ""){//Prevent printing nonesense internal validation indices if the clustering failed
		printInternalValidationIndices(data,p,number_data,number_features);
	}

	free(p);

	return clustersString;	
}

int calcBestThreshold(float &thresholdProp, string matType, string ccCriterion, int nbRuns, int neStop, string nbcCriterion, int &errres){

	//Preserve the state of "internalValidationOn" and avoid printing during iterations
	bool internalValidationOnState=internalValidationOn;
	internalValidationOn=false;

	double bestCHIndex=0.1, CHIndex=0.1;

	for(float prop=0.05; prop<=0.8; prop+=0.05){//Iterate to find the best proportion
				
		//Allocate a new matrix to preserve the original values of MatSimil
		double *MatSimilCopy = (double*)malloc(nbSequences*nbSequences*sizeof(double)); 
		if (MatSimilCopy == NULL){
			return 1;
		}

		//Initialize MatSimilCopy
		for(int i=0; i<nbSequences*nbSequences; i++){
			MatSimilCopy[i]=MatSimil[i];
		}
				
		//Lock the usage of the workers processes
		motifsParallelCalcInProgress=true;
		//Send data to workers if available
		taskID=4;
		for(int j=1; j<numprocs; j++){
			MPI_Send(&progEnd,1,MPI_INT,j,99,MPI_COMM_WORLD);
			MPI_Send(&taskID,1,MPI_INT,j,98,MPI_COMM_WORLD);
			MPI_Send(&nbSequences,1,MPI_INT,j,97,MPI_COMM_WORLD);
			MPI_Send(MatSimilCopy,nbSequences*nbSequences,MPI_DOUBLE,j,96,MPI_COMM_WORLD);
			int matTypeCode, ccCriterionCode, nbcCriterionCode;
			if(tolowerstr(matType)=="mot4"){ //M4 motifs
				matTypeCode=1;
			}else if(tolowerstr(matType)=="wmot4"){ //Weighted M4 motifs
				matTypeCode=2;
			}else if(tolowerstr(matType)=="mot13"){ //M13 motifs
				matTypeCode=3;
			}else if(tolowerstr(matType)=="wmot13"){ //Weighted M13 motifs
				matTypeCode=4;
			}
			MPI_Send(&matTypeCode,1,MPI_INT,j,95,MPI_COMM_WORLD);
			MPI_Send(&prop,1,MPI_FLOAT,j,94,MPI_COMM_WORLD);
			if(tolowerstr(ccCriterion)=="bestbic")
				ccCriterionCode=1;
			else if(tolowerstr(ccCriterion)=="bestaic")
				ccCriterionCode=2;
			else if(tolowerstr(ccCriterion)=="besticl")
				ccCriterionCode=3;
			else if(tolowerstr(ccCriterion)=="fast")
				ccCriterionCode=4;
			else if(tolowerstr(ccCriterion)=="mostfreq")
				ccCriterionCode=5;
			else if(tolowerstr(ccCriterion)=="besticl-r")
				ccCriterionCode=6;
			else if(tolowerstr(ccCriterion)=="bestbic-r")
				ccCriterionCode=7;
			MPI_Send(&ccCriterionCode,1,MPI_INT,j,93,MPI_COMM_WORLD);
			MPI_Send(&nbRuns,1,MPI_INT,j,92,MPI_COMM_WORLD);
			MPI_Send(&neStop,1,MPI_INT,j,91,MPI_COMM_WORLD);
			if(tolowerstr(nbcCriterion)=="bic")
				nbcCriterionCode=1;
			else if(tolowerstr(nbcCriterion)=="aic")
				nbcCriterionCode=2;
			else
				nbcCriterionCode=3;
			MPI_Send(&nbcCriterionCode,1,MPI_INT,j,90,MPI_COMM_WORLD);
		}

		//Master as worker

		zeroNSS(MatSimilCopy, nbSequences, prop, "pclosest"); //Apply the threshold annulation on matsimil

		toAffinityMat(MatSimilCopy, nbSequences,matType);
				
		int chNbEigenVectsCopy, nbClustersCopy=-1; //if nbClusters is provided as -1 to GMM_Clustering() then it will be calculated within GMM_Clustering() based on the best BIC.
				
		double **vecPropTMot; //Will be allocated and initialized in featuresCalc()	
				
		errres=featuresCalc(MatSimilCopy,vecPropTMot,chNbEigenVectsCopy,nbClustersCopy,nbSequences,EVMethod,EVThreshold,matType);
		//Check for errors
		if(errres==1){
			return 0;
		}
				
		//Free the memory used by MatSimil to leave enough memory for the GMM
		free(MatSimilCopy);

		if(!(chNbEigenVectsCopy==1 && (tolowerstr(ccCriterion)=="bestbic-r" || tolowerstr(ccCriterion)=="besticl-r"))){ //mclust cannot cluster a single feature matrix
			GMM_Clustering(vecPropTMot, nc, nbSequences, chNbEigenVectsCopy, ccCriterion, nbRuns, neStop, nbcCriterion, "full", 1000, nbClustersCopy, &CHIndex); //Got the Calinsky-Harabasz index for this proportion in CHIndex
									
			if(CHIndex>9999999999 || CHIndex<-9999999999)//If the index is infinity then the clustering is a single cluster, then bad
				CHIndex = -1.0;

			if(CHIndex>bestCHIndex){//Found a better CH => set the new proportion
				bestCHIndex=CHIndex;
				thresholdProp=prop;
			}
		}
				
		//Receive the data from the workers if they exist
		for(int k=1; k<numprocs; k++){
			prop+=0.05;
			MPI_Recv(&CHIndex, 1,MPI_DOUBLE, k, k*100+16, MPI_COMM_WORLD, &status);		
										
			if(CHIndex>bestCHIndex){//Found a better CH => set the new proportion
				bestCHIndex=CHIndex;
				thresholdProp=prop;
			}
		}
		//unlock the usage of the workers processes
		motifsParallelCalcInProgress=false;
				
				
		//Free danamically allocated memory
		free(vecPropTMot[0]); free(vecPropTMot);
	}

	//Set back "internalValidationOn" to its initial state
	internalValidationOn=internalValidationOnState;

	return 0;

}

void printRefsMatSimilToFile(){
	std::ofstream out1("matSimil.txt");
	std::ofstream out2("refs.txt");
	
	for(int i=0; i<nbSequences; i++){
		for(int j=0; j<nbSequences; j++){
			out1 << getStringFromDouble(MatSimil[i*nbSequences+j]);
			out1 << " ";
		}
		out1 << "\r\n";

		out2 << nc[i];
		out2 << ",";
	}
	out1.close();
	out2.close();
}