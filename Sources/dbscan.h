#ifndef DBSCAN_H
#define DBSCAN_H

#include <vector>
#include <cmath>

#include <stdio.h>
#include <iostream>

#define UNCLASSIFIED -1
#define CORE_POINT 1
#define BORDER_POINT 2
#define NOISE -2
#define DBSCAN_SUCCESS 0
#define FAILURE -3

using namespace std;

typedef struct DBSPoint_
{
	int nbCoords;
    double c0,c1,c2,c3,c4,c5,c6,c7,c8,c9;// point coordinates in a space of dimension nbCoords<=10
    int clusterID;  // cluster ID of the point
}DBSPoint;

class DBSCAN {
public:    
    DBSCAN(unsigned int minPts, float eps, vector<DBSPoint> points){
        m_minPoints = minPts;
        m_epsilon = eps;
        m_points = points;
        m_pointSize = points.size();
    }
    ~DBSCAN(){}

    int run();
    vector<int> calculateCluster(DBSPoint point);
    int expandCluster(DBSPoint point, int clusterID);
    inline double calculateDistance(DBSPoint pointCore, DBSPoint pointTarget);

    int getTotalPointSize() {return m_pointSize;}
    int getMinimumClusterSize() {return m_minPoints;}
    int getEpsilonSize() {return m_epsilon;}
    vector<DBSPoint> m_points;
private:
    unsigned int m_pointSize;
    unsigned int m_minPoints;
    float m_epsilon;
};

#endif // DBSCAN_H


int DBSCAN::run()
{
    int clusterID = 1;
    vector<DBSPoint>::iterator iter;
    for(iter = m_points.begin(); iter != m_points.end(); ++iter)
    {
        if ( iter->clusterID == UNCLASSIFIED )
        {
            if ( expandCluster(*iter, clusterID) != FAILURE )
            {
                clusterID += 1;
            }
        }
    }

    return 0;
}

int DBSCAN::expandCluster(DBSPoint point, int clusterID)
{    
    vector<int> clusterSeeds = calculateCluster(point);

    if ( clusterSeeds.size() < m_minPoints )
    {
        point.clusterID = NOISE;
        return FAILURE;
    }
    else
    {
        int index = 0, indexCorePoint = 0;
        vector<int>::iterator iterSeeds;
        for( iterSeeds = clusterSeeds.begin(); iterSeeds != clusterSeeds.end(); ++iterSeeds)
        {
            m_points.at(*iterSeeds).clusterID = clusterID;
			bool equalCoords = true;
			for(int i=0; i<point.nbCoords && equalCoords; i++){
				if(i==0 && m_points.at(*iterSeeds).c0 != point.c0){
					equalCoords = false;
				}else if(i==1 && m_points.at(*iterSeeds).c1 != point.c1){
					equalCoords = false;
				}else if(i==2 && m_points.at(*iterSeeds).c2 != point.c2){
					equalCoords = false;
				}else if(i==3 && m_points.at(*iterSeeds).c3 != point.c3){
					equalCoords = false;
				}else if(i==4 && m_points.at(*iterSeeds).c4 != point.c4){
					equalCoords = false;
				}else if(i==5 && m_points.at(*iterSeeds).c5 != point.c5){
					equalCoords = false;
				}else if(i==6 && m_points.at(*iterSeeds).c6 != point.c6){
					equalCoords = false;
				}else if(i==7 && m_points.at(*iterSeeds).c7 != point.c7){
					equalCoords = false;
				}else if(i==8 && m_points.at(*iterSeeds).c8 != point.c8){
					equalCoords = false;
				}else if(i==9 && m_points.at(*iterSeeds).c9 != point.c9){
					equalCoords = false;
				}
			}
            if (equalCoords)
            {
                indexCorePoint = index;
            }
            ++index;
        }
        clusterSeeds.erase(clusterSeeds.begin()+indexCorePoint);

        for( vector<int>::size_type i = 0, n = clusterSeeds.size(); i < n; ++i )
        {
            vector<int> clusterNeighors = calculateCluster(m_points.at(clusterSeeds[i]));

            if ( clusterNeighors.size() >= m_minPoints )
            {
                vector<int>::iterator iterNeighors;
                for ( iterNeighors = clusterNeighors.begin(); iterNeighors != clusterNeighors.end(); ++iterNeighors )
                {
                    if ( m_points.at(*iterNeighors).clusterID == UNCLASSIFIED || m_points.at(*iterNeighors).clusterID == NOISE )
                    {
                        if ( m_points.at(*iterNeighors).clusterID == UNCLASSIFIED )
                        {
                            clusterSeeds.push_back(*iterNeighors);
                            n = clusterSeeds.size();
                        }
                        m_points.at(*iterNeighors).clusterID = clusterID;
                    }
                }
            }
        }

        return DBSCAN_SUCCESS;
    }
}

vector<int> DBSCAN::calculateCluster(DBSPoint point)
{
    int index = 0;
    vector<DBSPoint>::iterator iter;
    vector<int> clusterIndex;
    for( iter = m_points.begin(); iter != m_points.end(); ++iter)
    {
        if ( calculateDistance(point, *iter) <= m_epsilon )
        {
            clusterIndex.push_back(index);
        }
        index++;
    }
    return clusterIndex;
}

inline double DBSCAN::calculateDistance( DBSPoint pointCore, DBSPoint pointTarget)
{
	double dist=0.0;
	for(int i=0; i<pointCore.nbCoords; i++){
		if(i==0){
			dist+=pow(pointCore.c0 - pointTarget.c0, 2);
		} else if(i==1){
			dist+=pow(pointCore.c1 - pointTarget.c1, 2);
		} else if(i==2){
			dist+=pow(pointCore.c2 - pointTarget.c2, 2);
		} else if(i==3){
			dist+=pow(pointCore.c3 - pointTarget.c3, 2);
		} else if(i==4){
			dist+=pow(pointCore.c4 - pointTarget.c4, 2);
		} else if(i==5){
			dist+=pow(pointCore.c5 - pointTarget.c5, 2);
		} else if(i==6){
			dist+=pow(pointCore.c6 - pointTarget.c6, 2);
		} else if(i==7){
			dist+=pow(pointCore.c7 - pointTarget.c7, 2);
		} else if(i==8){
			dist+=pow(pointCore.c8 - pointTarget.c8, 2);
		} else if(i==9){
			dist+=pow(pointCore.c9 - pointTarget.c9, 2);
		}
	}
    return dist;
}

string dbscanClust(int nbElts, int nbFeatures, double **features, string *labels, float eps, int minClusterNbElts = 1){
	string clustersString="";
	if(nbFeatures>10){
		cout<<"Warning: only 10 coordinates out of "<<nbFeatures<<" were used!\n";
		nbFeatures=10;
	} else if (nbFeatures == 0){
		cout<<"Error: empty coordinates table is provided!\n";
		clustersString+="[";
		for(int i=0; i<nbElts; i++){
			if(i==0){
				clustersString+="'"+labels[i]+"'";
			}else{
				clustersString+=", '"+labels[i]+"'";
			}
		}
		clustersString+="]";
		return clustersString;
	}

	vector<DBSPoint> points;

	//Initialize the points vector from the features matrix

	DBSPoint *p = (DBSPoint *)calloc(nbElts, sizeof(DBSPoint));
	
    for(int i=0; i < nbElts; i++)
    {
		p[i].nbCoords=nbFeatures;
        p[i].clusterID = UNCLASSIFIED;
		if(nbFeatures >= 1){
			p[i].c0 = features[i][0];
		}
		if(nbFeatures >= 2){
			p[i].c1 = features[i][1];
		}
		if(nbFeatures >= 3){
			p[i].c2 = features[i][2];
		}
		if(nbFeatures >= 4){
			p[i].c3 = features[i][3];
		}
		if(nbFeatures >= 5){
			p[i].c4 = features[i][4];
		}
		if(nbFeatures >= 6){
			p[i].c5 = features[i][5];
		}
		if(nbFeatures >= 7){
			p[i].c6 = features[i][6];
		}
		if(nbFeatures >= 8){
			p[i].c7 = features[i][7];
		}
		if(nbFeatures >= 9){
			p[i].c8 = features[i][8];
		}
		if(nbFeatures == 10){
			p[i].c9 = features[i][9];
		}
        points.push_back(p[i]);
	}

    free(p);

	//Perform the clustering

	DBSCAN ds(minClusterNbElts, eps*eps, points);
	ds.run();

	//Get the produced number of clusters

	int nbClusters = 0;
	for (int i=0; i<nbElts; i++)
    {
		//cout<<ds.m_points[i].c0<<"\t"<<ds.m_points[i].c1<<"\t"<<ds.m_points[i].c2<<"\t"<<ds.m_points[i].c3<<"\t"<<ds.m_points[i].clusterID<<endl;
		if(ds.m_points[i].clusterID>nbClusters){
			nbClusters=ds.m_points[i].clusterID;
		}
    }

	//Generate the clusters string

	for(int i=1; i<=nbClusters; i++){
		bool firstElt=true;
		clustersString+="[";
		for(int j=0; j<nbElts; j++){
			if(ds.m_points[j].clusterID == i){
				if(firstElt){
					clustersString+="'"+labels[j]+"'";
					firstElt=false;
				}else{
					clustersString+=", '"+labels[j]+"'";
				}
			}
		}
		clustersString+="]\n\n";
	}

	return clustersString;
}

