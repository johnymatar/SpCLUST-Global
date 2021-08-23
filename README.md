# SpCLUST-Global
SpCLUST-Global is a package embedding multiple clustering techniques for clustering biological sequences, including those presenting a high level of divergence.
SpCLUST-Global includes the option of using Edgar, R.C.'s MUSCLE module (www.drive5.com) for sequences alignment.

# Prerequisite
SpCLUST-Global uses MPI for parallel computation and the executable building for the installation package. Below are some basic instructions for installing MPI on your system.
- For Linux users:
  •	install mpich
  •	install openmpi
  •	install openmpi-devel
  •	echo "export PATH=$PATH:/usr/lib64/openmpi/bin" >> ~/.bashrc
- For Windows users:
  •	Download MS-MPI SDK and Redist installers from Microsoft's website: https://msdn.microsoft.com/en-us/library/bb524831.aspx
  •	Install the downloaded packages

# Usage instructions
- Download the executables from the "Executables" directory 
  or download the source files from the "Sources" directory and build the executables on your own machine by using "mpic++" or "M.S. Visual Studio"
- Call the executables with the desired arguments
- For serial computation use "spclust" with the desired arguments
- For parallel computation use "mpispclust" with the desired arguments
- To use the graphical interface, install mono (run "apt install mono-complete" as a sudoer) and then call "guispclust"

# Current version features
- Cross-platform: tested on Linux and Windows. Although not tested, the source files should also compile and run on MAC OS
- Embeds the traditional EM-GMM, MOTIFS GMM, DBSCAN, HDBSCAN, and CHAINS clustering techniques
- Parallel computation for the distance matrix using MPI, thus enabling its use on computation clusters.
- usage: mpispclust -in [input fasta file] -out [output clustering file] -cTech [clustering technique] -alignMode [alignment mode] -mdist [scoring matrix] -ccCriterion [clustering choice criterion] -nbcCriterion [number of clusters choice criterion] -nbRuns [maximum number of GMM runs] -neStop [nb of runs with no emprovement to stop] -matType [affinity matrix type] -tt [thresholding technique] -th [threshold used with the selected thresholding technique] -tp [clustering technique parameter used as epsilon for DBSCAN or minimum elements per cluster for HDBSCAN or minimum loop elements for CHAINS] -outputAlig [true or not]
     or: mpiexec -n [number of slave processes] spclust -in [input fasta file] -out [output clustering file] -alignMode [alignment mode] -mdist [scoring matrix] -ccCriterion [clustering choice criterion] -nbcCriterion [number of clusters choice criterion] -nbRuns [maximum number of GMM runs] -neStop [nb of runs with no emprovement to stop] -matType [affinity matrix type] -tt [thresholding technique] -th [threshold used with the selected thresholding technique] -tp [clustering technique parameter used as epsilon for DBSCAN or minimum elements per cluster for HDBSCAN or minimum loop elements for CHAINS] -outputAlig [true or not]

		Available clustering techniques are: GMM, CHAINS, DBSCAN, and HDBSCAN. Defaults to 'GMM' if not specified.
		Available alignment modes are: none, done, fast, moderate, and maxPrecision. done considers the input sequences aligned and does not perforn any alignment. fast and moderate limit the number of iterations for the alignment to 2 and 4 respectively while maxPrecision does not set such limit (using MUSCLE). Defaults to 'none' if not specified.
		Available scoring matrices are: none, EDNAFULL, BLOSUM62, and PAM250. Defaults to 'none' if not specified.
		Available clustering choice criteria are: bestBIC, bestAIC, bestICL, mostFreq, fast, bestICL-R, and bestBIC-R. Defaults to bestBIC.
		   nbRuns and neStop are both used for the clustering choice criterion bestBIC to indicate the stop condition for the best BIC choice.
		   neStop is ignored for the clustering choice criterion mostFreq. GMM will be run nbRuns times and the most occurent clustering will be chosen.
		   nbRuns and neStop are ignored for the clustering choice criterion fast. GMM will be run only once.
		   If not specified, nbRuns defaults to 500 and neStop defaults to 50.
		Available number of clusters choice criteria are: BIC, AIC, and ICL. Defaults to BIC.
		Available affinity matrices types are UL (Unnormalized Laplacian), RWNL (Random Walk Normalized Laplacian), MOD (Modularity), BH (Bethe Hessian), MOT4 (M4 Motifs), WMOT4 (Weighted M4 Motifs), FMOT4 (Fast MOT4 with user defined threshold), FWMOT4 (Fast WMOT4 with user defined threshold), MOT13 (M13 Motifs), WMOT13 (Weighted M13 Motifs), FMOT13 (Fast MOT13 with user defined threshold), and FWMOT13 (Fast WMOT13 with user defined threshold). Defaults to RWNL if not specified.
		Available thresholding techniques are pClosest, AVG, PropAVG, and deltas. Defaults to pClosest if not specified.
		The 'th' is used as threshold for thresholding techniques. It defaults to 1.0 if not specified.
		The 'tp' is used as epsilon for DBSCAN or minimum elements per cluster for HDBSCAN or minimum loop elements for CHAINS. If not specified, it defaults to 0.1 for DBSCAN, 5 for HDBSCAN and 2 for CHAINS.
		Only in the case where the 'alignMode' is set to none, the 'outputAlig' enables the output of the pairwise aligments generated when set to 'true'. It is either ignored or considered false otherwise.