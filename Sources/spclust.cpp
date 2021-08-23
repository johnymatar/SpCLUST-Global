//#include<amconf.h>
#include<iostream>
#include"spclustFunctions.h"
using namespace std;


int main(int argc, char* argv[]) {
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD,&numprocs); 
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	
	dCount=new int[numprocs];
	//Setting default value and reading the input arguments needed by all processes
	if(readSetGlobalParamsAP(argc, argv)==1){//An invalid input is detected
		MPI_Finalize();
		return 0;
	}

	if(myid==0){ //This is the main process
		string fGroupes;
		string groupes;
		string matType;
		string cTech;
		string ccCriterion;
		string nbcCriterion;
		string nbRunsStr;
		int nbRuns;
		string neStopStr;
		int neStop;
		int errres; //Used later to capture the error number from several functions

		//Displaying usage message in case called without arguments
		if(argc==1){
			//cout<<"You are running "<<PACKAGE_STRING<<endl;
			cout<<"SpClust performs biological sequences clustering using GMM.\n\n"
			<<"usage: mpispclust -in [input fasta file] -out [output clustering file] -cTech [clustering technique] -alignMode [alignment mode] -mdist [scoring matrix] -ccCriterion [clustering choice criterion] -nbcCriterion [number of clusters choice criterion] -nbRuns [maximum number of GMM runs] -neStop [nb of runs with no emprovement to stop] -matType [affinity matrix type] -tt [thresholding technique] -th [threshold used with the selected thresholding technique] -tp [clustering technique parameter used as epsilon for DBSCAN or minimum elements per cluster for HDBSCAN or minimum loop elements for CHAINS] -outputAlig [true or not]\n"
			<<"   or: mpiexec -n [number of slave processes] spclust -in [input fasta file] -out [output clustering file] -alignMode [alignment mode] -mdist [scoring matrix] -ccCriterion [clustering choice criterion] -nbcCriterion [number of clusters choice criterion] -nbRuns [maximum number of GMM runs] -neStop [nb of runs with no emprovement to stop] -matType [affinity matrix type] -tt [thresholding technique] -th [threshold used with the selected thresholding technique] -tp [clustering technique parameter used as epsilon for DBSCAN or minimum elements per cluster for HDBSCAN or minimum loop elements for CHAINS] -outputAlig [true or not]\n\n"
			<<"Available clustering techniques are: GMM, CHAINS, DBSCAN, and HDBSCAN. Defaults to 'GMM' if not specified.\n"
			<<"Available alignment modes are: none, done, fast, moderate, and maxPrecision. done considers the input sequences aligned and does not perforn any alignment. fast and moderate limit the number of iterations for the alignment to 2 and 4 respectively while maxPrecision does not set such limit (using MUSCLE). Defaults to 'none' if not specified.\n"
			<<"Available scoring matrices are: none, EDNAFULL, BLOSUM62, and PAM250. Defaults to 'none' if not specified.\n"
			<<"Available clustering choice criteria are: bestBIC, bestAIC, bestICL, mostFreq, fast, bestICL-R, and bestBIC-R. Defaults to bestBIC.\n"
			<<"   nbRuns and neStop are both used for the clustering choice criterion bestBIC to indicate the stop condition for the best BIC choice.\n"
			<<"   neStop is ignored for the clustering choice criterion mostFreq. GMM will be run nbRuns times and the most occurent clustering will be chosen.\n"
			<<"   nbRuns and neStop are ignored for the clustering choice criterion fast. GMM will be run only once.\n"
			<<"   If not specified, nbRuns defaults to 500 and neStop defaults to 50.\n"
			<<"Available number of clusters choice criteria are: BIC, AIC, and ICL. Defaults to BIC.\n"
			<<"Available affinity matrices types are UL (Unnormalized Laplacian), RWNL (Random Walk Normalized Laplacian), MOD (Modularity), BH (Bethe Hessian), MOT4 (M4 Motifs), WMOT4 (Weighted M4 Motifs), FMOT4 (Fast MOT4 with user defined threshold), FWMOT4 (Fast WMOT4 with user defined threshold), MOT13 (M13 Motifs), WMOT13 (Weighted M13 Motifs), FMOT13 (Fast MOT13 with user defined threshold), and FWMOT13 (Fast WMOT13 with user defined threshold). Defaults to RWNL if not specified.\n"
			<<"Available thresholding techniques are pClosest, AVG, PropAVG, and deltas. Defaults to pClosest if not specified.\n"
			<<"The 'th' is used as threshold for thresholding techniques. It defaults to 1.0 if not specified.\n"
			<<"The 'tp' is used as epsilon for DBSCAN or minimum elements per cluster for HDBSCAN or minimum loop elements for CHAINS. If not specified, it defaults to 0.1 for DBSCAN, 5 for HDBSCAN and 2 for CHAINS.\n"
			<<"Only in the case where the 'alignMode' is set to none, the 'outputAlig' enables the output of the pairwise aligments generated when set to 'true'. It is either ignored or considered false otherwise.\n\n"
			/*<<"Note: parameters are case sensitive, e.g. using -alignmode instead of -alignMode will cause this parameter to be disregarded, and using blosum62 instead of BLOSUM62 will be mentioned as an error.\n\n"*/;
			killWorkers();
			MPI_Finalize();
			return 0;
		}

		//Checking if our program is installed
	#if !(defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)) //Not under Windows OS
		char chkcmd[20];
		string spclustres, /*spclustGMMres, */muscleres;
		strcpy(chkcmd,"whereis spclust");
		spclustres = exec(chkcmd);
		/*strcpy(chkcmd,"whereis spclustGMM");
		spclustGMMres = exec(chkcmd);*/
		strcpy(chkcmd,"whereis muscle");
		muscleres = exec(chkcmd);
		if(spclustres.length()>9 && /*spclustGMMres.length()>12 && */muscleres.length()>9){
			installed=true;
			//Setting the working directory
			strcpy(progPath,exec("pwd").c_str());
			progPath[strlen(progPath)-1]='/';
		}
	#endif

		//Setting default value and reading the input arguments needed by the master process only
		if(readSetParamsMP(argc, argv, fGroupes, matType, cTech, ccCriterion, nbcCriterion, nbRunsStr, neStopStr)==1){//An invalid input is detected
			killWorkers();
			MPI_Finalize();
			return 0;
		}

		//Validate the input arguments
		if(validateInputParams(mDist, alignMode, fGroupes, ccCriterion, nbcCriterion, nbRunsStr, nbRuns, neStopStr, neStop, matType, cTech, thresholdingTech)==1){//An invalid input is detected
			killWorkers();
			MPI_Finalize();
			return 0;
		}
	
		if(tolowerstr(alignMode)!="none"){//Calculate the similarities in the original GCLUST's python method

			//Check if the required modules are present if needed
			char modsPath[300]="";//For storing the paths of the required files in order to check their existance		
		#if !(defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)) //Not under Windows OS
			//Check muscle
			strcpy(modsPath,progPath);
			strcat(modsPath,"muscle");
			std::ifstream infile2(modsPath);
			if(!infile2.good() && !installed){
				cout<<"Error: muscle executable is missing or not accessible.\n\n";
				killWorkers();
				MPI_Finalize();
				return 0;
			}

		#else
			//Check muscle
			strcpy(modsPath,progPath);
			strcat(modsPath,"muscle.exe");
			std::ifstream infile2(modsPath);
			if(!infile2.good()){
				cout<<"Error: muscle executable is missing or not accessible.\n\n";
				killWorkers();
				MPI_Finalize();
				return 0;
			}
		#endif

			//Calculate the similarity matrix
			printTimestamp("Alignment and similarity matrix calculation started at: ");
			errres = similarity(fListe,alignMode); //Capture the return value to check for errors

			//Check for alignment errors
			if(errres==1){
				killWorkers();
				cout<<"Sequences alignment failed. Job aborted!\n";
				MPI_Finalize();
				return 0;
			} else if(errres==2){
				killWorkers();
				cout<<"Input aligned sequences must have the same size. Job aborted!\n";
				MPI_Finalize();
				return 0;
			} else if(errres==3){
				killWorkers();
				cout<<"Memory allocation error. Job aborted!\n";
				MPI_Finalize();
				return 0;
			}
			printTimestamp("Alignment and similarity matrix calculation ended at: ");

		} else { //Calculate the similarities using edlib library
			
			printTimestamp("Similarity matrix calculation started at: ");

			ifstream f(outputPath.c_str());
			if(outputAlig && f.good()){//Delete the output alignment files if an alignment output is requested
				f.close();
				while( remove( outputPath.c_str() ) != 0 );
			}

			if(similarityNoAlig(MatSimil, numprocs, nbSequences, rawSeqs, mDist, taskID, progEnd, buffDouble, dCount, myid, outputAlig)==1){
				killWorkers();
				cout<<"Memory allocation error. Job aborted!\n";
				MPI_Finalize();
				return 0;
			}

			printTimestamp("Similarity matrix calculation ended at: ");

		}
		//Display the similarity report
		if(dispSimilStat){
			double minSimil=1.0, maxSimil=0.0, avgSimil=0.0;
			vector<int> duplicates;		
			for (int i = 0; i<nbSequences; i++){
				for (int j = i+1; j<nbSequences; j++){
					if(MatSimil[i*nbSequences+j]>maxSimil)
						maxSimil=MatSimil[i*nbSequences+j];
					if(MatSimil[i*nbSequences+j]<minSimil){
						minSimil=MatSimil[i*nbSequences+j];
						//cout<<nc[j]<<" and "<<nc[i]<<" scored a new minimum of "<<minSimil<<"\n";
					}
					//Print a list of identical sequences
					if(MatSimil[i*nbSequences+j]==1){
						if (std::find(duplicates.begin(), duplicates.end(),j)==duplicates.end()){
							duplicates.push_back(j);
							cout<<nc[j]<<" is identical to "<<nc[i]<<"\n";
						}
					}
					avgSimil+=MatSimil[i*nbSequences+j];
					//cout<<MatSimil[i*nbSequences+j]<<"\t";
				}
				//cout<<endl<<endl;
			}
			if(duplicates.size()>0){
				cout<<duplicates.size()<<" sequences indetical to others were detected!\n";
			}
			avgSimil=avgSimil/(((nbSequences*nbSequences)-nbSequences)/2);
			cout<<"The maximum similarity is: "<<maxSimil<<endl;
			cout<<"The minimum similarity is: "<<minSimil<<endl;
			cout<<"The average similarity is: "<<avgSimil<<endl;
		}

		//printRefsMatSimilToFile();

		if(tolowerstr(cTech)=="chains"){
			if(techParam==-1.0){
				techParam=2;
			}else if(techParam<2){
				techParam=2;
				cout<<"A minimum of 2 elements is required to form a loop. 'tp' is set to 2.\n";
			}
			printTimestamp("Clustering started at: ");
			string clustersString = chainsClust(nbSequences,MatSimil,nc,techParam);
			ofstream outfile;
			outfile.open(fGroupes.c_str());
			outfile<<clustersString;
			outfile.close();
			printTimestamp("Clustering ended at: ");
		} else {

			if((tolowerstr(cTech)=="gmm")&&(techParam!=-1.0)){
				cout<<"The 'tp' parameter is ignored for the GMM 'cTech'.\n";
			}else{
				if((tolowerstr(cTech)=="dbscan")&&(techParam==-1.0)){
					techParam=0.1;
				}else if((tolowerstr(cTech)=="hdbscan")&&(techParam==-1.0)){
					techParam=5;
				}else if((tolowerstr(cTech)=="dbscan")&&(techParam<=0)){
					cout<<"The 'tp' parameter must be a positive float for the chosen 'cTech'. The default value of 0.1 will be used.\n";
					techParam=0.1;
				}else if((tolowerstr(cTech)=="hdbscan")&&(techParam<=0)){
					cout<<"The 'tp' parameter must be a positive integer for the chosen 'cTech'. The default value of 5 will be used.\n";
					techParam=5;
				}
			}
		
			if(tolowerstr(matType)=="mot4" || tolowerstr(matType)=="wmot4" || tolowerstr(matType)=="mot13" || tolowerstr(matType)=="wmot13"){//Iterate to find the threshold for motifs clustering
			
				if(tolowerstr(thresholdingTech)!="pclosest"){
					cout<<"The thresholding technique has been set to pClosest.\n";
					thresholdingTech="pClosest";
				}
				if(thresholdProp!=1){
					cout<<"The input threshold is ignored. This threshold is automatically calculated.\n";
				}
				printTimestamp("Best threshold calculation started at: ");

				if(calcBestThreshold(thresholdProp, matType, ccCriterion, nbRuns, neStop, nbcCriterion, errres)==1){
					killWorkers();
					cout<<"Allocation for 'MatSimilCopy' failed. \n";
					MPI_Finalize();
					return 0;
				}			
				if(errres==1){
					killWorkers();
					cout<<"Memory allocation error. Job aborted!\n";
					MPI_Finalize();
					return 0;
				}
			
				printTimestamp("Best threshold calculation ended at: ");
			} else if(tolowerstr(matType)=="fmot4"){//Use the user input threshold and annulation technique
				matType="mot4";
			} else if(tolowerstr(matType)=="fwmot4"){//Use the user input threshold and annulation technique
				matType="wmot4";
			} else if(tolowerstr(matType)=="fmot13"){//Use the user input threshold and annulation technique
				matType="mot13";
			} else if(tolowerstr(matType)=="fwmot13"){//Use the user input threshold and annulation technique
				matType="wmot13";
			}

			//cout<<endl<<"Threshold: "<<thresholdProp<<endl;

			//Validate the threshold and eliminate the similarities under a certain threshold
			if(thresholdProp<=0){
				cout<<"The 'th' parameter must be a positive float. The default value of 1.0 will be used.\n";
				thresholdProp=1.0;
			}

			if(tolowerstr(thresholdingTech)!="keptprop" && tolowerstr(thresholdingTech)!="pclosest" || thresholdProp!=1){ //Nothing to do if this condition is not satisfied
				printTimestamp("Thresholding started at: ");
				if(dispSimilStat){
					cout<<zeroNSS(MatSimil, nbSequences, thresholdProp, thresholdingTech)<<" elements were zeroed.\n";
				} else{
					zeroNSS(MatSimil, nbSequences, thresholdProp, thresholdingTech);
				}
				printTimestamp("Thresholding ended at: ");
			}
			
			printTimestamp("Affinity matrix calculation started at: ");

			//Calculate the affinity matrix
			toAffinityMat(MatSimil, nbSequences,matType);
		
			printTimestamp("Affinity matrix calculation ended at: ");
		
			printTimestamp("Features matrix calculation started at: ");
		
			int chNbEigenVects, nbClusters=-1; //if nbClusters is provided as -1 to GMM_Clustering() then it will be calculated within GMM_Clustering() based on the provided criterion.
		
			double **vecPropT; //Will be allocated and initialized in featuresCalc()

			errres=featuresCalc(MatSimil,vecPropT,chNbEigenVects,nbClusters,nbSequences,EVMethod,EVThreshold,matType);
			//Check for errors
			if(errres==1){
				killWorkers();
				cout<<"Memory allocation error. Job aborted!\n";
				MPI_Finalize();
				return 0;
			}

			//Free the memory used by MatSimil to leave enough memory for the GMM
			if(MatSimil!=NULL)
				delete [] MatSimil;

			//Output the features matrix to a file
			/*string featuresString = to_string(nbSequences) + "\n";
			for(int i=0; i<nbSequences; i++){
				for(int j=0; j<chNbEigenVects; j++){
					if(j!=0){
						featuresString+=",";
					}
					featuresString+=to_string(vecPropT[i][j]);
				}
				featuresString+="\n";
			}
			ofstream outFeatures;
			outFeatures.open("featuresMatrix.txt");
			outFeatures<<featuresString;
			outFeatures.close();*/

			printTimestamp("Features calculation ended at: ");

			if(!(chNbEigenVects==1 && (tolowerstr(ccCriterion)=="bestbic-r" || tolowerstr(ccCriterion)=="besticl-r"))){ //mclust cannot cluster a single feature matrix
				printTimestamp("Clustering started at: ");
				string clustersString;
				
				if(tolowerstr(cTech)=="dbscan"){
					clustersString=dbscanClust(nbSequences,chNbEigenVects,vecPropT,nc,techParam);
				} else if(tolowerstr(cTech)=="hdbscan"){
					clustersString=hdbscanClust(vecPropT,nc,chNbEigenVects,nbSequences,1,techParam,"Euclidean");
				} else {//cTech is GMM
					clustersString= GMM_Clustering(vecPropT, nc, nbSequences, chNbEigenVects, ccCriterion, nbRuns, neStop, nbcCriterion, "full", 1000, nbClusters);
				}

				if((tolowerstr(ccCriterion)=="bestbic-r" || tolowerstr(ccCriterion)=="besticl-r") && clustersString == ""){//mclust failed to produce a clustering
					cout<<"Clustering failed! Mclust failed to cluster this set.\nPlease try a different ccCriterion.";
				} else {
					printTimestamp("Clustering ended at: ");

					ofstream outfile;
					outfile.open(fGroupes.c_str());
					outfile<<clustersString;
					outfile.close();
				}
			} else {
				cout<<"Clustering failed! Mclust cannot perform the clustering using a single feature.\nPlease switch to a ccCriterion other than bestBIC-R and bestICL-R, or try using a different matType.\n";
			}

			//Free danamically allocated memory
			free(vecPropT[0]); free(vecPropT);
		}

		//Close the workers processes
		killWorkers();

		//cout<<"done";system("pause");
		

		//Free danamically allocated memory
		if(nc!=NULL)
			delete[] nc;

	} else { //This is a calculation worker process
		//The slave does not print any profiling or validations
		profilingOn=false; internalValidationOn=false;

		//Check if there is no work to do
		MPI_Recv(&progEnd, 1, MPI_INT, 0, 99, MPI_COMM_WORLD, &status);
		//Variables used by the slave
		int number_data, *p, sentNbClusters, nbClusters, number_features, seed, max_iterations;
		int matTypeCode, ccCriterionCode, nbRuns, neStop, nbcCriterionCode;
		int chNbEigenVectsCopy, nbClustersCopy=-1; //if nbClusters is provided as -1 to GMM_Clustering() then it will be calculated within GMM_Clustering() based on the best BIC.
		float prop;
		double **data, myCrit, *MatSimilCopy, **vecPropTMot, CHIndex=0.1;
		bool calculateBestNbClusters;
		string matType, ccCriterion, nbcCriterion, clString="";
		GaussianMixture *gmm;

		while(progEnd!=1){
			MPI_Recv(&taskID, 1, MPI_INT, 0, 98, MPI_COMM_WORLD, &status); //Receive the aligned sequence size
			switch(taskID){
				case 1:
					MPI_Recv(&nbSequences, 1, MPI_INT, 0, 97, MPI_COMM_WORLD, &status); //Receive the aligned sequence size
					MPI_Recv(&Aligned_seq_length, 1, MPI_INT, 0, 96, MPI_COMM_WORLD, &status); //Receive the aligned sequence size
					aliSeqs=new char[nbSequences*(Aligned_seq_length+1)+1]; //Allocate aliSeqs to receive the aligned sequences
					MPI_Recv(aliSeqs, nbSequences*(Aligned_seq_length+1)+1, MPI_CHAR, 0, 95, MPI_COMM_WORLD, &status); //Receive the aligned sequences
		
					buffChar1=new char[Aligned_seq_length+1];
					buffChar2=new char[Aligned_seq_length+1];
		
					//Do the required calculations from this slave
					slaveDistCalc();

					//Send the results to the master
					MPI_Send(buffDouble,dCount[myid],MPI_DOUBLE,0,myid*100+11,MPI_COMM_WORLD);

					delete [] aliSeqs; delete [] buffDouble; delete [] buffChar1; delete [] buffChar2; //Free dynamic allocations
					break;
				case 2:
					MPI_Recv(&number_data, 1, MPI_INT, 0, 97, MPI_COMM_WORLD, &status); //Receive the number of data			
					p = (int*)malloc(number_data*sizeof(int));

					MPI_Recv(&sentNbClusters, 1, MPI_INT, 0, 96, MPI_COMM_WORLD, &status); //Receive the user sent number of clusters
					calculateBestNbClusters=(sentNbClusters==-1); //if nbClusters != -1 then it has not been calculated and provided by other methods. It will be calculated here based on the best BIC.
					
					nbClusters=sentNbClusters;

					MPI_Recv(&number_features, 1, MPI_INT, 0, 95, MPI_COMM_WORLD, &status); //Receive the number of features
					
					data=alloc_2d_double(number_data,number_features);
					MPI_Recv(&(data[0][0]), number_data*number_features, MPI_DOUBLE, 0, 94, MPI_COMM_WORLD, &status); //Receive the data
					
					MPI_Recv(&seed, 1, MPI_INT, 0, 93, MPI_COMM_WORLD, &status); //Receive the round seed
					seed+=myid; //Calculate my seed

					MPI_Recv(&max_iterations, 1, MPI_INT, 0, 92, MPI_COMM_WORLD, &status); //Receive the seed

					MPI_Recv(&ccCriterionCode, 1, MPI_INT, 0, 91, MPI_COMM_WORLD, &status); //Receive the clustering criterion code

					MPI_Recv(&nbcCriterionCode, 1, MPI_INT, 0, 90, MPI_COMM_WORLD, &status); //Receive the nb clusters criterion code
										
					if(nbcCriterionCode==1)
						nbcCriterion="bic";
					else if(ccCriterionCode==2)
						nbcCriterion="aic";
					else //bestICL
						nbcCriterion="icl";

					if(calculateBestNbClusters)
						nbClusters = getBestNbClusters(data, number_data, number_features, seed, nbcCriterion, "full", max_iterations);
					
					gmm = new GaussianMixture(nbClusters, number_data, number_features, "full", max_iterations, 0.01, 0.001, seed);    
					(*gmm).train(data,p,number_data,number_features);
					
					if(ccCriterionCode==1)
						myCrit=(*gmm).bic();
					else if(ccCriterionCode==2)
						myCrit=(*gmm).aic();
					else //bestICL
						myCrit=(*gmm).icl();
					
					//Send the results to the master
					MPI_Send(&myCrit, 1, MPI_DOUBLE, 0, myid*100+12, MPI_COMM_WORLD);

					//Free dynamically allocated data
					delete gmm;
					free(data[0]); free(data);
					free(p);
					break;
				case 3:
					MPI_Recv(&number_data, 1, MPI_INT, 0, 97, MPI_COMM_WORLD, &status); //Receive the number of data			
					p = (int*)malloc(number_data*sizeof(int));

					MPI_Recv(&sentNbClusters, 1, MPI_INT, 0, 96, MPI_COMM_WORLD, &status); //Receive the user sent number of clusters
					calculateBestNbClusters=(sentNbClusters==-1); //if nbClusters != -1 then it has not been calculated and provided by other methods. It will be calculated here based on the best BIC.
					
					nbClusters=sentNbClusters;

					MPI_Recv(&number_features, 1, MPI_INT, 0, 95, MPI_COMM_WORLD, &status); //Receive the number of features
					
					data=alloc_2d_double(number_data,number_features);
					MPI_Recv(&(data[0][0]), number_data*number_features, MPI_DOUBLE, 0, 94, MPI_COMM_WORLD, &status); //Receive the data
					
					MPI_Recv(&seed, 1, MPI_INT, 0, 93, MPI_COMM_WORLD, &status); //Receive the round seed
					seed+=myid; //Calculate my seed

					MPI_Recv(&max_iterations, 1, MPI_INT, 0, 92, MPI_COMM_WORLD, &status); //Receive the seed
					
					MPI_Recv(&nbcCriterionCode, 1, MPI_INT, 0, 91, MPI_COMM_WORLD, &status); //Receive the nb clusters criterion code
										
					if(nbcCriterionCode==1)
						nbcCriterion="bic";
					else if(nbcCriterionCode==2)
						nbcCriterion="aic";
					else //bestICL
						nbcCriterion="icl";


					if(calculateBestNbClusters)
						nbClusters = getBestNbClusters(data, number_data, number_features, seed, nbcCriterion, "full", max_iterations);
					
					gmm = new GaussianMixture(nbClusters, number_data, number_features, "full", max_iterations, 0.01, 0.001, seed);    
					(*gmm).train(data,p,number_data,number_features);
					
					if(nbcCriterionCode==1)
						myCrit=(*gmm).bic();
					else if(nbcCriterionCode==2)
						myCrit=(*gmm).aic();
					else //bestICL
						myCrit=(*gmm).icl();
					
					//Send the results to the master
					MPI_Send(&myCrit, 1, MPI_DOUBLE, 0, myid*100+13, MPI_COMM_WORLD);
					MPI_Send(&nbClusters, 1, MPI_INT, 0, myid*100+14, MPI_COMM_WORLD);
					MPI_Send(p, number_data, MPI_INT, 0, myid*100+15, MPI_COMM_WORLD);

					//Free dynamically allocated data
					delete gmm;
					free(data[0]); free(data);
					free(p);
					break;
				case 4:
					//Lock the usage of the workers processes
					motifsParallelCalcInProgress=true;

					MPI_Recv(&number_data, 1, MPI_INT, 0, 97, MPI_COMM_WORLD, &status);
					MatSimilCopy = (double*)malloc(number_data*number_data*sizeof(double));
					MPI_Recv(MatSimilCopy, number_data*number_data, MPI_DOUBLE, 0, 96, MPI_COMM_WORLD, &status);
					MPI_Recv(&matTypeCode, 1, MPI_INT, 0, 95, MPI_COMM_WORLD, &status);
					MPI_Recv(&prop, 1, MPI_FLOAT, 0, 94, MPI_COMM_WORLD, &status);
					prop=prop+(myid*0.05);
					MPI_Recv(&ccCriterionCode, 1, MPI_INT, 0, 93, MPI_COMM_WORLD, &status);
					MPI_Recv(&nbRuns, 1, MPI_INT, 0, 92, MPI_COMM_WORLD, &status);
					MPI_Recv(&neStop, 1, MPI_INT, 0, 91, MPI_COMM_WORLD, &status);
					MPI_Recv(&nbcCriterionCode, 1, MPI_INT, 0, 90, MPI_COMM_WORLD, &status);
					nbSequences=number_data;

					if(matTypeCode==1){ //M4 motifs
						matType="mot4";
					}else if(matTypeCode==2){//Weighted M4 motifs
						matType="wmot4";
					}else if(matTypeCode==3){//M13 motifs
						matType="mot13";
					}else if(matTypeCode==4){//Weighted M13 motifs
						matType="wmot13";
					}

					if(ccCriterionCode==1)
						ccCriterion="bestbic";
					else if(ccCriterionCode==2)
						ccCriterion="bestaic";
					else if(ccCriterionCode==3)
						ccCriterion="besticl";
					else if(ccCriterionCode==4)
						ccCriterion="fast";
					else if(ccCriterionCode==5)
						ccCriterion="mostfreq";
					else if(ccCriterionCode==6)
						ccCriterion="besticl-r";
					else if(ccCriterionCode==7)
						ccCriterion="bestbic-r";

					if(nbcCriterionCode==1)
						nbcCriterion="bic";
					else if(nbcCriterionCode==2)
						nbcCriterion="aic";
					else
						nbcCriterion="icl";

					zeroNSS(MatSimilCopy, number_data, prop, "pclosest"); //Apply the threshold annulation on matsimil
					
					toAffinityMat(MatSimilCopy, number_data,matType);
				
					featuresCalc(MatSimilCopy,vecPropTMot,chNbEigenVectsCopy,nbClustersCopy,number_data,EVMethod,EVThreshold,matType);
		
					//Free the memory used by MatSimil to leave enough memory for the GMM
					free(MatSimilCopy);
					
					//Fill NC with dummy data to make GMM_Clustering() work
					nc=new string[number_data];
					for(int i=0; i<number_data; i++){
						nc[i]=i;
					}

					if(!(chNbEigenVectsCopy==1 && (tolowerstr(ccCriterion)=="bestbic-r" || tolowerstr(ccCriterion)=="besticl-r"))){ //mclust cannot cluster a single feature matrix
						clString=GMM_Clustering(vecPropTMot, nc, nbSequences, chNbEigenVectsCopy, ccCriterion, nbRuns, neStop, nbcCriterion, "full", 1000, nbClustersCopy, &CHIndex); //Got the Calinsky-Harabasz index for this proportion in CHIndex
					
						if(CHIndex>9999999999 || CHIndex<-9999999999) //If the index is infinity then the clustering is a single cluster, then bad
							CHIndex=-1.0;
					}else{CHIndex=-1.0;};
					//Send the results to the master
					MPI_Send(&CHIndex, 1, MPI_DOUBLE, 0, myid*100+16, MPI_COMM_WORLD);
					
					//Free danamically allocated memory
					delete [] nc;
					free(vecPropTMot[0]); free(vecPropTMot);
					
					//Unlock the usage of the workers processes
					motifsParallelCalcInProgress=false;
					break;
				case 5:					
					//Do the required calculations from this slave
					slaveSimilCalc(nbSequences, dCount, buffDouble, rawSeqs, mDist, myid);

					//Send the results to the master
					MPI_Send(buffDouble,dCount[myid],MPI_DOUBLE,0,myid*100+17,MPI_COMM_WORLD);

					delete [] buffDouble; delete [] rawSeqs; delete [] nc;//Free dynamic allocations
					break;
				default:
					cout<<"Invalid task ID received by worker process "<<myid<<".\n";
					break;
			}
			MPI_Recv(&progEnd, 1, MPI_INT, 0, 99, MPI_COMM_WORLD, &status); //Check if there is no work to do

		}
	}
	delete []dCount;
	MPI_Finalize();
	return 0;
}
