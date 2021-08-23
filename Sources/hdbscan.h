#pragma once
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <set>
#include <map>
#include <cstdint>
#include <tuple> 
#include <cmath>
#include <limits>
#include <algorithm>
#include <list>
#include <stdexcept>
#include<cstdio>
using namespace std;

class outlierScore
{
private:
	double coreDistance;
public:
	double score;
	int id;
	/// <summary>
	/// Creates a new OutlierScore for a given point.
	/// </summary>
	/// <param name="score">The outlier score of the point</param>
	/// <param name="coreDistance">The point's core distance</param>
	/// <param name="id">The id (index) of the point</param>
	outlierScore(double score, double coreDistance, int id);
	outlierScore();
	/// <summary>
	/// Method Overridden to compare two objects.
	/// </summary>
	bool operator<(const outlierScore& other) const;


};

class hdbscanResult
{
public:
	vector <int> labels;
	vector <outlierScore> outliersScores;
	vector <double> membershipProbabilities;
	bool hasInfiniteStability;
	hdbscanResult();
	hdbscanResult(vector<int> pLables, vector<outlierScore> pOutlierScores, vector <double> pmembershipProbabilities, bool pHsInfiniteStability);
};

class HDBSCANCluster
{
private:
	int _id;
	double _birthLevel;
	double _deathLevel;
	int _numPoints;
	double _propagatedStability;
	int _numConstraintsSatisfied;
	int _propagatedNumConstraintsSatisfied;
	std::set<int> _virtualChildCluster;
	static int counter;

public:
	std::vector<HDBSCANCluster*> PropagatedDescendants;
	double PropagatedLowestChildDeathLevel;
	HDBSCANCluster* Parent;
	double Stability;
	bool HasChildren;
	int Label;
	int HierarchyPosition;   //First level where points with this HDBSCANCluster's label appear

	HDBSCANCluster();

	HDBSCANCluster(int label, HDBSCANCluster *parent, double birthLevel, int numPoints);
	bool operator==(const HDBSCANCluster& other) const;
	void detachPoints(int numPoints, double level);
	void propagate();
	void addPointsToVirtualChildCluster(std::set<int> points);
	
	bool virtualChildClusterConstraintsPoint(int point);

	void addVirtualChildConstraintsSatisfied(int numConstraints);
	

	void addConstraintsSatisfied(int numConstraints);


	void releaseVirtualChildCluster();

	int getClusterId();

};

enum hdbscanConstraintType{mustLink, cannotLink};
/// <summary>
/// A clustering constraint (either a must-link or cannot-link constraint between two points).
/// </summary>
class hdbscanConstraint
{
private : 
	hdbscanConstraintType _constraintType;
	int _pointA;
	int _pointB;
/// <summary>
/// Creates a new constraint.
/// </summary>
/// <param name="pointA">The first point involved in the constraint</param>
/// <param name="pointB">The second point involved in the constraint</param>
/// <param name="type">The constraint type</param>
public:
	hdbscanConstraint(int pointA, int pointB, hdbscanConstraintType type);

	int getPointA();

	int getPointB();

	hdbscanConstraintType getConstraintType();

};

class undirectedGraph
{
private:
	int _numVertices;
	std::vector<int> _verticesA;
	std::vector<int> _verticesB;
	std::vector<double> _edgeWeights;
	std::vector<std::vector<int>> _edges;

public:
	undirectedGraph(int numVertices, std::vector<int> verticesA, std::vector<int> verticesB, std::vector<double> edgeWeights)
	{
		_numVertices = numVertices;
		_verticesA = verticesA;
		_verticesB = verticesB;
		_edgeWeights = edgeWeights;
		_edges.resize(numVertices);
		int _edgesLength = _edges.size();
		int _edgeWeightsLength = _edgeWeights.size();
		for (int i = 0; i < _edgeWeightsLength; i++)
		{
			_edges[_verticesA[i]].push_back(_verticesB[i]);

			if (_verticesA[i] != _verticesB[i])
				_edges[_verticesB[i]].push_back(_verticesA[i]);
		}

	}

	void quicksortByEdgeWeight();
	int getNumVertices();

	int getNumEdges();

	int getFirstVertexAtIndex(int index);
	int getSecondVertexAtIndex(int index);

	double getEdgeWeightAtIndex(int index);
	std::vector<int> &getEdgeListForVertex(int vertex);
private:
	int selectPivotIndex(int startIndex, int endIndex);

	int partition(int startIndex, int endIndex, int pivotIndex);
	void swapEdges(int indexOne, int indexTwo);

};

class IDistanceCalculator
{
	/// <summary>
	/// Computes the distance between two points.
	/// Note that larger values indicate that the two points are farther apart.
	/// </summary>
	/// <param name="attributesOne">The attributes of the first point</param>
	/// <param name="attributesTwo">The attributes of the second point</param>
	/// <returns>A double for the distance between the two points</returns>
public:
	virtual double computeDistance(std::vector<double> attributesOne, std::vector<double> attributesTwo)=0;
};

class EuclideanDistance : IDistanceCalculator
{
public:
	double computeDistance(std::vector<double> attributesOne, std::vector<double> attributesTwo);

};

class ManhattanDistance : IDistanceCalculator
{
public:
	double computeDistance(std::vector<double> attributesOne, std::vector<double> attributesTwo);
};

class hdbscanAlgorithm
{
public:
	/// <summary>
	/// Calculates the core distances for each point in the data set, given some value for k.
	/// </summary>
	/// <param name="distances">A vector of vectors where index [i][j] indicates the jth attribute of data point i</param>
	/// <param name="k">Each point's core distance will be it's distance to the kth nearest neighbor</param>
	/// <returns> An array of core distances</returns>
	static std::vector<double> calculateCoreDistances(std::vector<std::vector<double>> distances, int k);
	
	static undirectedGraph constructMst(std::vector<std::vector<double>> distances, std::vector<double> coreDistances, bool selfEdges);
	

	/// <summary>
	/// Propagates constraint satisfaction, stability, and lowest child death level from each child
	/// HDBSCANCluster to each parent HDBSCANCluster in the tree.  This method must be called before calling
	/// findProminentClusters() or calculateOutlierScores().
	/// </summary>
	/// <param name="clusters">A list of Clusters forming a HDBSCANCluster tree</param>
	/// <returns>true if there are any clusters with infinite stability, false otherwise</returns>


	static void computeHierarchyAndClusterTree(undirectedGraph *mst, int minClusterSize, std::vector<hdbscanConstraint> constraints, std::vector<std::vector<int>> &hierarchy, std::vector<double> &pointNoiseLevels, std::vector<int> &pointLastClusters, std::vector<HDBSCANCluster*> &clusters);
	
	static std::vector<int> findProminentClusters(std::vector<HDBSCANCluster*> &clusters, std::vector<std::vector<int>> &hierarchy, int numPoints);

	static std::vector<double> findMembershipScore(std::vector<int> clusterids, std::vector<double> coreDistances);
	
	static bool propagateTree(std::vector<HDBSCANCluster*> &sclusters);
	
	/// <summary>
	/// Produces the outlier score for each point in the data set, and returns a sorted list of outlier
	/// scores.  propagateTree() must be called before calling this method.
	/// </summary>
	/// <param name="clusters">A list of Clusters forming a HDBSCANCluster tree which has already been propagated</param>
	/// <param name="pointNoiseLevels">A double[] with the levels at which each point became noise</param>
	/// <param name="pointLastClusters">An int[] with the last label each point had before becoming noise</param>
	/// <param name="coreDistances">An array of core distances for each data point</param>
	/// <returns>An List of OutlierScores, sorted in descending order</returns>
	static std::vector<outlierScore> calculateOutlierScores(
		std::vector<HDBSCANCluster*> &clusters,
		std::vector<double> &pointNoiseLevels,
		std::vector<int> &pointLastClusters,
		std::vector<double> coreDistances);
	
	/// <summary>
	/// Removes the set of points from their parent HDBSCANCluster, and creates a new HDBSCANCluster, provided the
	/// clusterId is not 0 (noise).
	/// </summary>
	/// <param name="points">The set of points to be in the new HDBSCANCluster</param>
	/// <param name="clusterLabels">An array of HDBSCANCluster labels, which will be modified</param>
	/// <param name="parentCluster">The parent HDBSCANCluster of the new HDBSCANCluster being created</param>
	/// <param name="clusterLabel">The label of the new HDBSCANCluster </param>
	/// <param name="edgeWeight">The edge weight at which to remove the points from their previous HDBSCANCluster</param>
	/// <returns>The new HDBSCANCluster, or null if the clusterId was 0</returns>
	static HDBSCANCluster* createNewCluster(
		std::set<int>& points,
		std::vector<int> &clusterLabels,
		HDBSCANCluster *parentCluster,
		int clusterLabel,
		double edgeWeight);
	
	/// <summary>
	/// Calculates the number of constraints satisfied by the new clusters and virtual children of the
	/// parents of the new clusters.
	/// </summary>
	/// <param name="newClusterLabels">Labels of new clusters</param>
	/// <param name="clusters">An List of clusters</param>
	/// <param name="constraints">An List of constraints</param>
	/// <param name="clusterLabels">An array of current HDBSCANCluster labels for points</param>
	static void calculateNumConstraintsSatisfied(
		std::set<int>& newClusterLabels,
		std::vector<HDBSCANCluster*>& clusters,
		std::vector<hdbscanConstraint>& constraints,
		std::vector<int>& clusterLabels);
	
};

class bitSet
{
private:
	std::vector<bool> _bits;
public:
	bool get(int pos);
	
	void set(int pos);

	void ensure(int pos);
};

class hdbscanParameters
{
public:

	/// <summary>
	/// Parameters to be Passed to the HDBSCAN Algorithm
	/// </summary>
	/// <param name="distances">The attributes of the first point</param>
	/// <param name="dataset">The attributes of the second point</param>
	/// <param name="dataset">The attributes of the second point</param>
	/// <param name="distanceFunction">Defines the type of distance measure to use : Euclidean, Manhattan ,..</param>
	/// <param name="minPoints">Min Points in the HDBSCANCluster</param>
	/// <param name="minClusterSize">The minimum number of points which a HDBSCANCluster needs to be a valid HDBSCANCluster</param>
	vector< vector <double> > distances;
	vector< vector <double> > dataset;
	string distanceFunction;
	uint32_t minPoints;
	uint32_t minClusterSize;
	vector<hdbscanConstraint> constraints;
};

class hdbscanRunner
{
public:
	static hdbscanResult run(hdbscanParameters parameters);
};

class Hdbscan
{
private:

	hdbscanResult result;

public:

	vector < vector <double > > dataset;

	std::vector<int> labels_;

	std::vector<int> normalizedLabels_;

	std::vector<outlierScore>outlierScores_;

	std::vector <double> membershipProbabilities_;

	uint32_t noisyPoints_;

	uint32_t numClusters_;


	Hdbscan(){}
	Hdbscan(vector < vector <double > > dataset, int minPoints, int minClusterSize, string distanceMetric) {
		this->dataset = dataset;
		execute(minPoints, minClusterSize, distanceMetric);
	}
	Hdbscan(double **features, int nbFeatures, int nbData, int minPoints, int minClusterSize, string distanceMetric) {
		loadFeatures(features, nbFeatures, nbData);
		execute(minPoints, minClusterSize, distanceMetric);
	}	

	void loadFeatures(double **features, int nbFeatures, int nbData);
	void execute(int minPoints, int minClusterSize, string distanceMetric);
	hdbscanResult getResult();

};


void Hdbscan::execute(int minPoints, int minClusterSize, string distanceMetric) {
	//Call The Runner Class here
	hdbscanRunner runner;
	hdbscanParameters parameters;
	uint32_t noisyPoints = 0;
	set<int> numClustersSet;
	map<int, int> clustersMap;
	vector<int> normalizedLabels;

	parameters.dataset = this->dataset;
	parameters.minPoints = minPoints;
	parameters.minClusterSize = minClusterSize;
	parameters.distanceFunction = distanceMetric;
    	this->result = runner.run(parameters);
	this->labels_ = result.labels;
	this->outlierScores_ = result.outliersScores;
	for (uint32_t i = 0; i < result.labels.size(); i++) {
		if (result.labels[i] == 0) {
			noisyPoints++;
		}
		else {
			numClustersSet.insert(result.labels[i]);
		}
	}
	this->numClusters_ = numClustersSet.size();
	this->noisyPoints_ = noisyPoints;
	int iNdex = 1;
	for (auto it = numClustersSet.begin(); it != numClustersSet.end(); it++) {
		clustersMap[*it] = iNdex++;
	}
	for (int i = 0; i < labels_.size(); i++) {
		if (labels_[i] != 0)
			normalizedLabels.push_back(clustersMap[labels_[i]]);
		else if (labels_[i] == 0) {
			normalizedLabels.push_back(-1);
		}

	}
	this->normalizedLabels_ = normalizedLabels;
	this->membershipProbabilities_ = result.membershipProbabilities;
}

hdbscanResult Hdbscan::getResult(){
	return this->result;
}

void Hdbscan::loadFeatures(double **features, int nbFeatures, int nbData) {
	
	vector<vector<double> > dataset;
	for(int i=0; i<nbData; i++) {
		vector<double> row;
		for(int j=0; j<nbFeatures; j++){
			row.push_back(features[i][j]);
		}
		dataset.push_back(row);
	}	
	this->dataset=dataset;
}

hdbscanResult hdbscanRunner::run(hdbscanParameters parameters) {
	int numPoints = parameters.dataset.size() != 0 ? parameters.dataset.size() : parameters.distances.size();

	hdbscanAlgorithm algorithm;
	hdbscanResult result;
	if (parameters.distances.size() == 0) {
		std::vector<std::vector<double>> distances(numPoints);
		for (int i = 0; i < numPoints; i++) {
			distances[i].resize(numPoints);
			//distances[i]=std::vector<double>(numPoints);
			for (int j = 0; j < i; j++) {
				if (parameters.distanceFunction.length() == 0) {
					//Default to Euclidean
					EuclideanDistance EDistance;
					double distance;
					distance = EDistance.computeDistance(parameters.dataset[i], parameters.dataset[j]);
					distances[i][j] = distance;
					distances[j][i] = distance;

				}
				else if (parameters.distanceFunction == "Euclidean") {
					EuclideanDistance EDistance;
					double distance;
					distance = EDistance.computeDistance(parameters.dataset[i], parameters.dataset[j]);
					distances[i][j] = distance;
					distances[j][i] = distance;
				}
				else if (parameters.distanceFunction == "Manhattan") {
					ManhattanDistance MDistance;
					double distance;
					distance = MDistance.computeDistance(parameters.dataset[i], parameters.dataset[j]);
					distances[i][j] = distance;
					distances[j][i] = distance;
				}
			}
		}

		parameters.distances = distances;
	}

	std::vector <double> coreDistances = algorithm.calculateCoreDistances(
		parameters.distances,
		parameters.minPoints);

	undirectedGraph mst = algorithm.constructMst(
		parameters.distances,
		coreDistances,
		true);
	mst.quicksortByEdgeWeight();

	std::vector<double> pointNoiseLevels(numPoints);
	std::vector<int> pointLastClusters(numPoints);

	std::vector< std::vector <int> > hierarchy;

	std::vector<HDBSCANCluster*> clusters;
	algorithm.computeHierarchyAndClusterTree(
		&mst,
		parameters.minClusterSize,
		parameters.constraints,
		hierarchy,
		pointNoiseLevels,
		pointLastClusters,
		clusters);
	bool infiniteStability = algorithm.propagateTree(clusters);

	std::vector<int> prominentClusters = algorithm.findProminentClusters(clusters, hierarchy, numPoints);
	std::vector<double> membershipProbabilities = algorithm.findMembershipScore(prominentClusters, coreDistances);
	std::vector<outlierScore> scores = algorithm.calculateOutlierScores(
		clusters,
		pointNoiseLevels,
		pointLastClusters,
		coreDistances);

	return hdbscanResult(prominentClusters, scores, membershipProbabilities,  infiniteStability);
}

hdbscanResult::hdbscanResult() {
	;
}
hdbscanResult::hdbscanResult(vector<int> pLables, vector<outlierScore> pOutlierScores, vector<double> pmembershipProbabilities, bool pHsInfiniteStability) {
	labels = pLables;
	outliersScores = pOutlierScores;
	membershipProbabilities = pmembershipProbabilities;
	hasInfiniteStability = pHsInfiniteStability;
}

outlierScore::outlierScore() {
	;
}

outlierScore::outlierScore(double score, double coreDistance, int id) {
	outlierScore::score = score;
	outlierScore::coreDistance = coreDistance;
	outlierScore::id = id;
}

bool outlierScore::operator<(const outlierScore& other) const {
	/*
	if (score < other.score)
		return score < other.score;
	else if (coreDistance < other.coreDistance)
		return coreDistance < other.coreDistance;
	else if (id < other.id)
		return id < other.id;
	else
		return false;*/
	return std::tie(score, coreDistance, id) < std::tie(other.score, other.coreDistance, other.id);
}

double EuclideanDistance::computeDistance(std::vector<double> attributesOne, std::vector<double> attributesTwo) {
	double distance = 0;
	for (uint32_t i = 0; i < attributesOne.size() && i < attributesTwo.size(); i++) {
		distance += ((attributesOne[i] - attributesTwo[i]) * (attributesOne[i] - attributesTwo[i]));
	}

	return sqrt(distance);
}

double ManhattanDistance::computeDistance(std::vector<double> attributesOne, std::vector<double> attributesTwo) {
	double distance = 0;
	for (uint32_t i = 0; i < attributesOne.size() && i < attributesTwo.size(); i++) {
		distance += fabs(attributesOne[i] - attributesTwo[i]);
	}

	return distance;
}


std::vector<double> hdbscanAlgorithm::calculateCoreDistances(std::vector<std::vector<double>> distances, int k)
{
	int length = distances.size();

	int numNeighbors = k - 1;
	std::vector<double>coreDistances(length);
	if (k == 1)
	{
		for (int point = 0; point < length; point++)
		{
			coreDistances[point] = 0;
		}
		return coreDistances;
	}
	for (int point = 0; point < length; point++)
	{
		std::vector<double> kNNDistances(numNeighbors);  //Sorted nearest distances found so far
		for (int i = 0; i < numNeighbors; i++)
		{
			kNNDistances[i] = std::numeric_limits<double>::max();
		}

		for (int neighbor = 0; neighbor < length; neighbor++)
		{
			if (point == neighbor)
				continue;
			double distance = distances[point][neighbor];
			int neighborIndex = numNeighbors;
			//Check at which position in the nearest distances the current distance would fit:
			while (neighborIndex >= 1 && distance < kNNDistances[neighborIndex - 1])
			{
				neighborIndex--;
			}
			//Shift elements in the array to make room for the current distance:
			if (neighborIndex < numNeighbors)
			{
				for (int shiftIndex = numNeighbors - 1; shiftIndex > neighborIndex; shiftIndex--)
				{
					kNNDistances[shiftIndex] = kNNDistances[shiftIndex - 1];
				}
				kNNDistances[neighborIndex] = distance;
			}

		}
		coreDistances[point] = kNNDistances[numNeighbors - 1];
	}
	return coreDistances;
}
undirectedGraph hdbscanAlgorithm::constructMst(std::vector<std::vector<double>> distances, std::vector<double> coreDistances, bool selfEdges)
{
	int length = distances.size();
	int selfEdgeCapacity = 0;
	if (selfEdges)
		selfEdgeCapacity = length;
	bitSet attachedPoints;

	std::vector<int> nearestMRDNeighbors(length - 1 + selfEdgeCapacity);
	std::vector<double> nearestMRDDistances(length - 1 + selfEdgeCapacity);

	for (int i = 0; i < length - 1; i++)
	{
		nearestMRDDistances[i] = std::numeric_limits<double>::max();
	}

	int currentPoint = length - 1;
	int numAttachedPoints = 1;
	attachedPoints.set(length - 1);

	while (numAttachedPoints < length)
	{

		int nearestMRDPoint = -1;
		double nearestMRDDistance = std::numeric_limits<double>::max();
		for (int neighbor = 0; neighbor < length; neighbor++)
		{
			if (currentPoint == neighbor)
				continue;
			if (attachedPoints.get(neighbor) == true)
				continue;
			double distance = distances[currentPoint][neighbor];
			double mutualReachabiltiyDistance = distance;
			if (coreDistances[currentPoint] > mutualReachabiltiyDistance)
				mutualReachabiltiyDistance = coreDistances[currentPoint];

			if (coreDistances[neighbor] > mutualReachabiltiyDistance)
				mutualReachabiltiyDistance = coreDistances[neighbor];

			if (mutualReachabiltiyDistance < nearestMRDDistances[neighbor])
			{
				nearestMRDDistances[neighbor] = mutualReachabiltiyDistance;
				nearestMRDNeighbors[neighbor] = currentPoint;
			}

			if (nearestMRDDistances[neighbor] <= nearestMRDDistance)
			{
				nearestMRDDistance = nearestMRDDistances[neighbor];
				nearestMRDPoint = neighbor;
			}

		}
		attachedPoints.set(nearestMRDPoint);
		numAttachedPoints++;
		currentPoint = nearestMRDPoint;
	}
	std::vector<int> otherVertexIndices(length - 1 + selfEdgeCapacity);
	for (int i = 0; i < length - 1; i++)
	{
		otherVertexIndices[i] = i;
	}
	if (selfEdges)
	{
		for (int i = length - 1; i < length * 2 - 1; i++)
		{
			int vertex = i - (length - 1);
			nearestMRDNeighbors[i] = vertex;
			otherVertexIndices[i] = vertex;
			nearestMRDDistances[i] = coreDistances[vertex];
		}
	}
	undirectedGraph undirectedGraphObject(length, nearestMRDNeighbors, otherVertexIndices, nearestMRDDistances);
	return undirectedGraphObject;

}

void hdbscanAlgorithm::computeHierarchyAndClusterTree(undirectedGraph* mst, int minClusterSize, std::vector<hdbscanConstraint> constraints, std::vector<std::vector<int>>& hierarchy, std::vector<double>& pointNoiseLevels, std::vector<int>& pointLastClusters, std::vector<HDBSCANCluster*>& clusters)
{
	int hierarchyPosition = 0;

	//The current edge being removed from the MST:
	int currentEdgeIndex = mst->getNumEdges() - 1;
	int nextClusterLabel = 2;
	bool nextLevelSignificant = true;

	//The previous and current HDBSCANCluster numbers of each point in the data set:
	std::vector<int> previousClusterLabels(mst->getNumVertices());
	std::vector<int> currentClusterLabels(mst->getNumVertices());

	for (int i = 0; i < currentClusterLabels.size(); i++)
	{
		currentClusterLabels[i] = 1;
		previousClusterLabels[i] = 1;
	}
	//std::vector<HDBSCANCluster *> clusters;
	clusters.push_back(NULL);
	//HDBSCANCluster cluster_object(1, NULL, std::numeric_limits<double>::quiet_NaN(),  mst->getNumVertices());
	clusters.push_back(new HDBSCANCluster(1, NULL, std::numeric_limits<double>::quiet_NaN(), mst->getNumVertices()));

	std::set<int> clusterOne;
	clusterOne.insert(1);
	calculateNumConstraintsSatisfied(
		clusterOne,
		clusters,
		constraints,
		currentClusterLabels);
	std::set<int> affectedClusterLabels;
	std::set<int> affectedVertices;
	while (currentEdgeIndex >= 0)
	{
		double currentEdgeWeight = mst->getEdgeWeightAtIndex(currentEdgeIndex);
		std::vector<HDBSCANCluster*> newClusters;
		while (currentEdgeIndex >= 0 && mst->getEdgeWeightAtIndex(currentEdgeIndex) == currentEdgeWeight)
		{
			int firstVertex = mst->getFirstVertexAtIndex(currentEdgeIndex);
			int secondVertex = mst->getSecondVertexAtIndex(currentEdgeIndex);
			std::vector<int>& firstVertexEdgeList = mst->getEdgeListForVertex(firstVertex);
			std::vector<int>::iterator secondVertexInFirstEdgeList = std::find(firstVertexEdgeList.begin(), firstVertexEdgeList.end(), secondVertex);
			if (secondVertexInFirstEdgeList != mst->getEdgeListForVertex(firstVertex).end())
				mst->getEdgeListForVertex(firstVertex).erase(secondVertexInFirstEdgeList);
			std::vector<int>& secondVertexEdgeList = mst->getEdgeListForVertex(secondVertex);
			std::vector<int>::iterator firstVertexInSecondEdgeList = std::find(secondVertexEdgeList.begin(), secondVertexEdgeList.end(), firstVertex);
			if (firstVertexInSecondEdgeList != mst->getEdgeListForVertex(secondVertex).end())
				mst->getEdgeListForVertex(secondVertex).erase(firstVertexInSecondEdgeList);

			if (currentClusterLabels[firstVertex] == 0)
			{
				currentEdgeIndex--;
				continue;
			}
			affectedVertices.insert(firstVertex);
			affectedVertices.insert(secondVertex);
			affectedClusterLabels.insert(currentClusterLabels[firstVertex]);
			currentEdgeIndex--;
		}
		if (!affectedClusterLabels.size())
			continue;
		while (affectedClusterLabels.size())
		{
			int examinedClusterLabel = *prev(affectedClusterLabels.end());
			affectedClusterLabels.erase(prev(affectedClusterLabels.end()));
			std::set<int> examinedVertices;
			//std::set<int>::iterator affectedIt;
			for (auto affectedIt = affectedVertices.begin(); affectedIt != affectedVertices.end();)
			{
				int vertex = *affectedIt;
				if (currentClusterLabels[vertex] == examinedClusterLabel)
				{
					examinedVertices.insert(vertex);
					affectedIt = affectedVertices.erase(affectedIt);

				}
				else {
					++affectedIt;
				}
			}
			std::set<int> firstChildCluster;
			std::list<int> unexploredFirstChildClusterPoints;
			int numChildClusters = 0;
			while (examinedVertices.size())
			{

				std::set<int> constructingSubCluster;
				int iters = 0;
				std::list<int> unexploredSubClusterPoints;
				bool anyEdges = false;
				bool incrementedChildCount = false;
				int rootVertex = *prev(examinedVertices.end());
				constructingSubCluster.insert(rootVertex);
				unexploredSubClusterPoints.push_back(rootVertex);
				examinedVertices.erase(prev(examinedVertices.end()));
				while (unexploredSubClusterPoints.size())
				{
					int vertexToExplore = *unexploredSubClusterPoints.begin();
					unexploredSubClusterPoints.erase(unexploredSubClusterPoints.begin());
					std::vector<int>& vertexToExploreEdgeList = mst->getEdgeListForVertex(vertexToExplore);
					for (std::vector<int>::iterator it = vertexToExploreEdgeList.begin(); it != vertexToExploreEdgeList.end();)
					{
						int neighbor = *it;
						anyEdges = true;
						if (std::find(constructingSubCluster.begin(), constructingSubCluster.end(), neighbor) == constructingSubCluster.end())
						{
							constructingSubCluster.insert(neighbor);
							unexploredSubClusterPoints.push_back(neighbor);
							if (std::find(examinedVertices.begin(), examinedVertices.end(), neighbor) != examinedVertices.end())
								examinedVertices.erase(std::find(examinedVertices.begin(), examinedVertices.end(), neighbor));

						}
						else {
							++it;
						}
					}
					if (!incrementedChildCount && constructingSubCluster.size() >= minClusterSize && anyEdges)
					{
						incrementedChildCount = true;
						numChildClusters++;

						//If this is the first valid child HDBSCANCluster, stop exploring it:
						if (firstChildCluster.size() == 0)
						{
							firstChildCluster = constructingSubCluster;
							unexploredFirstChildClusterPoints = unexploredSubClusterPoints;
							break;
						}
					}

				}
				//If there could be a split, and this child HDBSCANCluster is valid:
				if (numChildClusters >= 2 && constructingSubCluster.size() >= minClusterSize && anyEdges)
				{
					//Check this child HDBSCANCluster is not equal to the unexplored first child HDBSCANCluster:
					int firstChildClusterMember = *prev(firstChildCluster.end());
					if (std::find(constructingSubCluster.begin(), constructingSubCluster.end(), firstChildClusterMember) != constructingSubCluster.end())
						numChildClusters--;
					//Otherwise, c a new HDBSCANCluster:
					else
					{
						HDBSCANCluster* newCluster = createNewCluster(constructingSubCluster, currentClusterLabels,
							clusters[examinedClusterLabel], nextClusterLabel, currentEdgeWeight);
						newClusters.push_back(newCluster);
						clusters.push_back(newCluster);
						nextClusterLabel++;
					}
				}
				else if (constructingSubCluster.size() < minClusterSize || !anyEdges)
				{
					createNewCluster(constructingSubCluster, currentClusterLabels,
						clusters[examinedClusterLabel], 0, currentEdgeWeight);

					for (std::set<int>::iterator it = constructingSubCluster.begin(); it != constructingSubCluster.end(); it++)
					{
						int point = *it;
						pointNoiseLevels[point] = currentEdgeWeight;
						pointLastClusters[point] = examinedClusterLabel;
					}
				}
			}
			if (numChildClusters >= 2 && currentClusterLabels[*firstChildCluster.begin()] == examinedClusterLabel)
			{
				while (unexploredFirstChildClusterPoints.size())
				{
					int vertexToExplore = *unexploredFirstChildClusterPoints.begin();
					unexploredFirstChildClusterPoints.pop_front();
					for (std::vector<int>::iterator it = mst->getEdgeListForVertex(vertexToExplore).begin(); it != mst->getEdgeListForVertex(vertexToExplore).end(); it++)
					{
						int neighbor = *it;
						if (std::find(firstChildCluster.begin(), firstChildCluster.end(), neighbor) == firstChildCluster.end())
						{
							firstChildCluster.insert(neighbor);
							unexploredFirstChildClusterPoints.push_back(neighbor);
						}
					}
				}
				HDBSCANCluster* newCluster = createNewCluster(firstChildCluster, currentClusterLabels,
					clusters[examinedClusterLabel], nextClusterLabel, currentEdgeWeight);
				newClusters.push_back(newCluster);
				clusters.push_back(newCluster);
				nextClusterLabel++;
			}
		}
		if (nextLevelSignificant || newClusters.size())
		{
			std::vector<int> lineContents(previousClusterLabels.size());
			for (int i = 0; i < previousClusterLabels.size(); i++)
				lineContents[i] = previousClusterLabels[i];
			hierarchy.push_back(lineContents);
			hierarchyPosition++;
		}
		std::set<int> newClusterLabels;
		for (std::vector<HDBSCANCluster*>::iterator it = newClusters.begin(); it != newClusters.end(); it++)
		{
			HDBSCANCluster* newCluster = *it;
			newCluster->HierarchyPosition = hierarchyPosition;
			newClusterLabels.insert(newCluster->Label);
		}
		if (newClusterLabels.size())
			calculateNumConstraintsSatisfied(newClusterLabels, clusters, constraints, currentClusterLabels);

		for (int i = 0; i < previousClusterLabels.size(); i++)
		{
			previousClusterLabels[i] = currentClusterLabels[i];
		}
		if (!newClusters.size())
			nextLevelSignificant = false;
		else
			nextLevelSignificant = true;
	}

	{
		std::vector<int> lineContents(previousClusterLabels.size() + 1);
		for (int i = 0; i < previousClusterLabels.size(); i++)
			lineContents[i] = 0;
		hierarchy.push_back(lineContents);
	}
}
std::vector<int> hdbscanAlgorithm::findProminentClusters(std::vector<HDBSCANCluster*>& clusters, std::vector<std::vector<int>>& hierarchy, int numPoints)
{
	//Take the list of propagated clusters from the root HDBSCANCluster:
	std::vector<HDBSCANCluster*> solution = clusters[1]->PropagatedDescendants;
	std::vector<int> flatPartitioning(numPoints);

	//Store all the hierarchy positions at which to find the birth points for the flat clustering:
	std::map<int, std::vector<int>> significantHierarchyPositions;

	std::vector<HDBSCANCluster*>::iterator it = solution.begin();
	while (it != solution.end())
	{
		int hierarchyPosition = (*it)->HierarchyPosition;
		if (significantHierarchyPositions.count(hierarchyPosition) > 0)
			significantHierarchyPositions[hierarchyPosition].push_back((*it)->Label);
		else
			significantHierarchyPositions[hierarchyPosition].push_back((*it)->Label);
		it++;
	}

	//Go through the hierarchy file, setting labels for the flat clustering:
	while (significantHierarchyPositions.size())
	{
		std::map<int, std::vector<int>>::iterator entry = significantHierarchyPositions.begin();
		std::vector<int> clusterList = entry->second;
		int hierarchyPosition = entry->first;
		significantHierarchyPositions.erase(entry->first);

		std::vector<int> lineContents = hierarchy[hierarchyPosition];

		for (int i = 0; i < lineContents.size(); i++)
		{
			int label = lineContents[i];
			if (std::find(clusterList.begin(), clusterList.end(), label) != clusterList.end())
				flatPartitioning[i] = label;
		}
	}
	return flatPartitioning;
}
std::vector<double> hdbscanAlgorithm::findMembershipScore(std::vector<int> clusterids, std::vector<double> coreDistances)
{
	
	int length = clusterids.size();
	std::vector<double> prob(length, std::numeric_limits<double>::max());
	int i=0;
	
	while(i<length)
	{
		if(prob[i]==std::numeric_limits<double>::max())
		{
			
			int clusterno = clusterids[i];
			std::vector<int>::iterator iter = clusterids.begin()+i;
			std::vector<int> indices;
			while ((iter = std::find(iter, clusterids.end(), clusterno)) != clusterids.end())
			{
				
				indices.push_back(distance(clusterids.begin(), iter));
				iter++;
				if(iter==clusterids.end())
					break;

			}
			if(clusterno==0)
			{
				for(int j=0; j<indices.size();j++)
				{
					prob[indices[j]] = 0;
				}
				i++;
				continue;
			}
			std::vector<double> tempCoreDistances(indices.size());
			for(int j=0; j<indices.size();j++)
			{
				tempCoreDistances[j] = coreDistances[j];
			}
			double maxCoreDistance = *max_element(tempCoreDistances.begin(), tempCoreDistances.end());
			for(int j=0; j<tempCoreDistances.size();j++)
			{
				prob[indices[j]] = ( maxCoreDistance - tempCoreDistances[j] ) / maxCoreDistance;
			}

		}
		
		i++;
	}
	return prob;
	
}

bool hdbscanAlgorithm::propagateTree(std::vector<HDBSCANCluster*>& clusters)
{
	std::map<int, HDBSCANCluster*> clustersToExamine;
	bitSet addedToExaminationList;
	bool infiniteStability = false;

	//Find all leaf clusters in the HDBSCANCluster tree:
	for (HDBSCANCluster* cluster : clusters)
	{
		if (cluster != NULL && !cluster->HasChildren)
		{
			int label = cluster->Label;
			clustersToExamine.erase(label);
			clustersToExamine.insert(std::pair<int,HDBSCANCluster*>(label,cluster));
			addedToExaminationList.set(label);
		}
	}
	//Iterate through every HDBSCANCluster, propagating stability from children to parents:
	while (clustersToExamine.size())
	{
		std::map<int, HDBSCANCluster*>::iterator currentKeyValue = prev(clustersToExamine.end());
		HDBSCANCluster* currentCluster = currentKeyValue->second;
		clustersToExamine.erase(currentKeyValue->first);
		currentCluster->propagate();

		if (currentCluster->Stability == std::numeric_limits<double>::infinity())
			infiniteStability = true;

		if (currentCluster->Parent != NULL)
		{
			HDBSCANCluster* parent = currentCluster->Parent;
			int label = parent->Label;

			if (!addedToExaminationList.get(label))
			{
				clustersToExamine.erase(label);
				clustersToExamine.insert(std::pair<int,HDBSCANCluster*>(label,parent));
				addedToExaminationList.set(label);
			}
		}
	}

	return infiniteStability;
}

/// <summary>
/// Produces the outlier score for each point in the data set, and returns a sorted list of outlier
/// scores.  propagateTree() must be called before calling this method.
/// </summary>
/// <param name="clusters">A list of Clusters forming a HDBSCANCluster tree which has already been propagated</param>
/// <param name="pointNoiseLevels">A double[] with the levels at which each point became noise</param>
/// <param name="pointLastClusters">An int[] with the last label each point had before becoming noise</param>
/// <param name="coreDistances">An array of core distances for each data point</param>
/// <returns>An List of OutlierScores, sorted in descending order</returns>
std::vector<outlierScore> hdbscanAlgorithm::calculateOutlierScores(
	std::vector<HDBSCANCluster*>& clusters,
	std::vector<double>& pointNoiseLevels,
	std::vector<int>& pointLastClusters,
	std::vector<double> coreDistances)
{
	int numPoints = pointNoiseLevels.size();
	std::vector<outlierScore> outlierScores;

	//Iterate through each point, calculating its outlier score:
	for (int i = 0; i < numPoints; i++)
	{
		double epsilonMax = clusters[pointLastClusters[i]]->PropagatedLowestChildDeathLevel;
		double epsilon = pointNoiseLevels[i];
		double score = 0;

		if (epsilon != 0)
			score = 1 - (epsilonMax / epsilon);

		outlierScores.push_back(outlierScore(score, coreDistances[i], i));
	}
	//Sort the outlier scores:
	sort(outlierScores.begin(), outlierScores.end());

	return outlierScores;
}

/// <summary>
/// Removes the set of points from their parent HDBSCANCluster, and creates a new HDBSCANCluster, provided the
/// clusterId is not 0 (noise).
/// </summary>
/// <param name="points">The set of points to be in the new HDBSCANCluster</param>
/// <param name="clusterLabels">An array of HDBSCANCluster labels, which will be modified</param>
/// <param name="parentCluster">The parent HDBSCANCluster of the new HDBSCANCluster being created</param>
/// <param name="clusterLabel">The label of the new HDBSCANCluster </param>
/// <param name="edgeWeight">The edge weight at which to remove the points from their previous HDBSCANCluster</param>
/// <returns>The new HDBSCANCluster, or null if the clusterId was 0</returns>
HDBSCANCluster* hdbscanAlgorithm::createNewCluster(
	std::set<int>& points,
	std::vector<int>& clusterLabels,
	HDBSCANCluster* parentCluster,
	int clusterLabel,
	double edgeWeight)
{
	std::set<int>::iterator it = points.begin();
	while (it != points.end())
	{
		clusterLabels[*it] = clusterLabel;
		++it;
	}
	parentCluster->detachPoints(points.size(), edgeWeight);

	if (clusterLabel != 0)
	{
		return new HDBSCANCluster(clusterLabel, parentCluster, edgeWeight, points.size());
	}

	parentCluster->addPointsToVirtualChildCluster(points);
	return NULL;
}
/// <summary>
/// Calculates the number of constraints satisfied by the new clusters and virtual children of the
/// parents of the new clusters.
/// </summary>
/// <param name="newClusterLabels">Labels of new clusters</param>
/// <param name="clusters">An List of clusters</param>
/// <param name="constraints">An List of constraints</param>
/// <param name="clusterLabels">An array of current HDBSCANCluster labels for points</param>
void hdbscanAlgorithm::calculateNumConstraintsSatisfied(
	std::set<int>& newClusterLabels,
	std::vector<HDBSCANCluster*>& clusters,
	std::vector<hdbscanConstraint>& constraints,
	std::vector<int>& clusterLabels)
{

	if (constraints.size() == 0)
		return;

	std::vector<HDBSCANCluster> parents;
	std::vector<HDBSCANCluster> ::iterator it;
	for (int label : newClusterLabels)
	{
		HDBSCANCluster* parent = clusters[label]->Parent;
		if (parent != NULL && !(find(parents.begin(), parents.end(), *parent) != parents.end()))
			parents.push_back(*parent);
	}

	for (hdbscanConstraint constraint : constraints)
	{
		int labelA = clusterLabels[constraint.getPointA()];
		int labelB = clusterLabels[constraint.getPointB()];

		if (constraint.getConstraintType() == hdbscanConstraintType::mustLink && labelA == labelB)
		{
			if (find(newClusterLabels.begin(), newClusterLabels.end(), labelA) != newClusterLabels.end())
				clusters[labelA]->addConstraintsSatisfied(2);
		}
		else if (constraint.getConstraintType() == hdbscanConstraintType::cannotLink && (labelA != labelB || labelA == 0))
		{
			if (labelA != 0 && find(newClusterLabels.begin(), newClusterLabels.end(), labelA) != newClusterLabels.end())
				clusters[labelA]->addConstraintsSatisfied(1);
			if (labelB != 0 && (find(newClusterLabels.begin(), newClusterLabels.end(), labelA) != newClusterLabels.end()))
				clusters[labelB]->addConstraintsSatisfied(1);
			if (labelA == 0)
			{
				for (HDBSCANCluster parent : parents)
				{
					if (parent.virtualChildClusterConstraintsPoint(constraint.getPointA()))
					{
						parent.addVirtualChildConstraintsSatisfied(1);
						break;
					}
				}
			}
			if (labelB == 0)
			{
				for (HDBSCANCluster parent : parents)
				{
					if (parent.virtualChildClusterConstraintsPoint(constraint.getPointB()))
					{
						parent.addVirtualChildConstraintsSatisfied(1);
						break;
					}
				}
			}
		}
	}

	for (HDBSCANCluster parent : parents)
	{
		parent.releaseVirtualChildCluster();
	}
}

void undirectedGraph::quicksortByEdgeWeight()
{
	int _edgeWeightsLength = _edgeWeights.size();
	if (_edgeWeightsLength <= 1)
		return;

	std::vector<int> startIndexStack(_edgeWeightsLength / 2);
	std::vector<int> endIndexStack(_edgeWeightsLength / 2);

	startIndexStack[0] = 0;
	endIndexStack[0] = _edgeWeightsLength - 1;

	int stackTop = 0;
	while (stackTop >= 0)
	{
		int startIndex = startIndexStack[stackTop];
		int endIndex = endIndexStack[stackTop];
		stackTop--;
		int pivotIndex = selectPivotIndex(startIndex, endIndex);
		pivotIndex = partition(startIndex, endIndex, pivotIndex);
		if (pivotIndex > startIndex + 1)
		{
			startIndexStack[stackTop + 1] = startIndex;
			endIndexStack[stackTop + 1] = pivotIndex - 1;
			stackTop++;
		}
		if (pivotIndex < endIndex - 1)
		{
			startIndexStack[stackTop + 1] = pivotIndex + 1;
			endIndexStack[stackTop + 1] = endIndex;
			stackTop++;
		}

	}
}
int undirectedGraph::selectPivotIndex(int startIndex, int endIndex)
{
	if (startIndex - endIndex <= 1)
		return startIndex;

	double first = _edgeWeights[startIndex];
	double middle = _edgeWeights[startIndex + (endIndex - startIndex) / 2];
	double last = _edgeWeights[endIndex];

	if (first <= middle)
	{
		if (middle <= last)
			return startIndex + (endIndex - startIndex) / 2;

		if (last >= first)
			return endIndex;

		return startIndex;
	}

	if (first <= last)
		return startIndex;

	if (last >= middle)
		return endIndex;

	return startIndex + (endIndex - startIndex) / 2;
}

int undirectedGraph::partition(int startIndex, int endIndex, int pivotIndex)
{
	double pivotValue = _edgeWeights[pivotIndex];
	swapEdges(pivotIndex, endIndex);
	int lowIndex = startIndex;
	for (int i = startIndex; i < endIndex; i++)
	{
		if (_edgeWeights[i] < pivotValue)
		{
			swapEdges(i, lowIndex);
			lowIndex++;
		}
	}
	swapEdges(lowIndex, endIndex);
	return lowIndex;
}

void undirectedGraph::swapEdges(int indexOne, int indexTwo)
{
	if (indexOne == indexTwo)
		return;

	int tempVertexA = _verticesA[indexOne];
	int tempVertexB = _verticesB[indexOne];
	double tempEdgeDistance = _edgeWeights[indexOne];
	_verticesA[indexOne] = _verticesA[indexTwo];
	_verticesB[indexOne] = _verticesB[indexTwo];
	_edgeWeights[indexOne] = _edgeWeights[indexTwo];
	_verticesA[indexTwo] = tempVertexA;
	_verticesB[indexTwo] = tempVertexB;
	_edgeWeights[indexTwo] = tempEdgeDistance;
}

int undirectedGraph::getNumVertices()
{
	return _numVertices;
}

int undirectedGraph::getNumEdges()
{
	return _edgeWeights.size();
}

int undirectedGraph::getFirstVertexAtIndex(int index)
{
	return _verticesA[index];
}

int undirectedGraph::getSecondVertexAtIndex(int index)
{
	return _verticesB[index];
}

double undirectedGraph::getEdgeWeightAtIndex(int index)
{
	return _edgeWeights[index];
}

std::vector<int>& undirectedGraph::getEdgeListForVertex(int vertex)
{
	return _edges[vertex];
}


int HDBSCANCluster::counter = 0;
HDBSCANCluster::HDBSCANCluster()
{

	_id = ++counter;
}

HDBSCANCluster::HDBSCANCluster(int label, HDBSCANCluster* parent, double birthLevel, int numPoints) //:Label(label), Parent(parent), _birthLevel(birthLevel), _numPoints(numPoints)
{
	_id = ++counter;
	_deathLevel = 0;

	_propagatedStability = 0;
	_numConstraintsSatisfied = 0;
	_propagatedNumConstraintsSatisfied = 0;

	Parent = parent;
	Label = label;
	_birthLevel = birthLevel;
	_numPoints = numPoints;
	HierarchyPosition = 0;
	Stability = 0;
	PropagatedLowestChildDeathLevel = std::numeric_limits<double>::max();

	if (Parent != NULL)
		Parent->HasChildren = true;
	HasChildren = false;
	PropagatedDescendants.resize(0);
}
bool HDBSCANCluster ::operator==(const HDBSCANCluster& other) const {
	return (this->_id == other._id);
}
void HDBSCANCluster::detachPoints(int numPoints, double level)
{
	_numPoints -= numPoints;
	Stability += (numPoints * (1 / level - 1 / _birthLevel));

	if (_numPoints == 0)
		_deathLevel = level;
	else if (_numPoints < 0)
		throw std::invalid_argument("HDBSCANCluster cannot have less than 0 points.");
}

void HDBSCANCluster::propagate()
{
	if (Parent != NULL)
	{
		if (PropagatedLowestChildDeathLevel == std::numeric_limits<double>::max())
			PropagatedLowestChildDeathLevel = _deathLevel;
		if (PropagatedLowestChildDeathLevel < Parent->PropagatedLowestChildDeathLevel)
			Parent->PropagatedLowestChildDeathLevel = PropagatedLowestChildDeathLevel;
		if (!HasChildren)
		{
			Parent->_propagatedNumConstraintsSatisfied += _numConstraintsSatisfied;
			Parent->_propagatedStability += Stability;
			Parent->PropagatedDescendants.push_back(this);
		}
		else if (_numConstraintsSatisfied > _propagatedNumConstraintsSatisfied)
		{
			Parent->_propagatedNumConstraintsSatisfied += _numConstraintsSatisfied;
			Parent->_propagatedStability += Stability;
			Parent->PropagatedDescendants.push_back(this);
		}
		else if (_numConstraintsSatisfied < _propagatedNumConstraintsSatisfied)
		{
			Parent->_propagatedNumConstraintsSatisfied += _propagatedNumConstraintsSatisfied;
			Parent->_propagatedStability += _propagatedStability;
			Parent->PropagatedDescendants.insert(Parent->PropagatedDescendants.end(), PropagatedDescendants.begin(), PropagatedDescendants.end());
		}
		else if (_numConstraintsSatisfied == _propagatedNumConstraintsSatisfied)
		{
			//Chose the parent over descendants if there is a tie in stability:
			if (Stability >= _propagatedStability)
			{
				Parent->_propagatedNumConstraintsSatisfied += _numConstraintsSatisfied;
				Parent->_propagatedStability += Stability;
				Parent->PropagatedDescendants.push_back(this);
			}
			else
			{
				Parent->_propagatedNumConstraintsSatisfied += _propagatedNumConstraintsSatisfied;
				Parent->_propagatedStability += _propagatedStability;
				Parent->PropagatedDescendants.insert(Parent->PropagatedDescendants.end(), PropagatedDescendants.begin(), PropagatedDescendants.end());
			}
		}
	}
}
void HDBSCANCluster::addPointsToVirtualChildCluster(std::set<int> points)
{
	for (std::set<int>::iterator it = points.begin(); it != points.end(); ++it) {
		_virtualChildCluster.insert(*it);
	}
}
bool HDBSCANCluster::virtualChildClusterConstraintsPoint(int point)
{
	return (_virtualChildCluster.find(point) != _virtualChildCluster.end());
}

void HDBSCANCluster::addVirtualChildConstraintsSatisfied(int numConstraints)
{
	_propagatedNumConstraintsSatisfied += numConstraints;
}

void HDBSCANCluster::addConstraintsSatisfied(int numConstraints)
{
	_numConstraintsSatisfied += numConstraints;
}

void HDBSCANCluster::releaseVirtualChildCluster()
{
	_virtualChildCluster.clear();
}

int HDBSCANCluster::getClusterId() {
	return this->_id;
}

hdbscanConstraint::hdbscanConstraint(int pointA, int pointB, hdbscanConstraintType type) {
	_pointA = pointA;
	_pointB = pointB;
	_constraintType = type;
}

int hdbscanConstraint::getPointA() {
	return _pointA;
}

int hdbscanConstraint::getPointB() {
	return _pointB;
}

hdbscanConstraintType hdbscanConstraint::getConstraintType() {
	return _constraintType;
}

bool bitSet::get(int pos) {
		return pos < _bits.size() && _bits[pos];
}

void bitSet::set(int pos) {
		ensure(pos);
		_bits[pos] = true;
}

void bitSet::ensure(int pos) {
	if (pos >= _bits.size())
	{
		_bits.resize(pos + 64);
	}
}

string hdbscanClust(double **features, string *nc, int nbFeatures, int nbData, int minPoints, int minClusterSize, string distanceMetric){
	string clustersString="";
	Hdbscan hdbscan(features, nbFeatures, nbData, minPoints, minClusterSize, distanceMetric);
	vector<int> labels = hdbscan.getResult().labels;
	int nbClusters = hdbscan.numClusters_;
	vector<int> normLabels = hdbscan.normalizedLabels_;

	//Add the clustered elements
	for(int i=1; i<=nbClusters; i++){
		clustersString+="[";
		bool firstElt=true;
		for(int j=0; j<nbData; j++){
			if(normLabels[j]==i){
				if(firstElt){
					clustersString+="'"+nc[j]+"'";
					firstElt=false;
				} else {
					clustersString+=", '"+nc[j]+"'";
				}
			}
		}
		clustersString+="]\n\n";
	}

	//Add the noisy points each to a seperate cluster
	for(int j=0; j<nbData; j++){
		if(normLabels[j]==-1){
			clustersString+="['"+nc[j]+"']\n\n";
		}
	}

	return clustersString;
}