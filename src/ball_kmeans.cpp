//
// Created by sumit on 3/16/20.
//

#include "ball_kmeans.h"
#include "general_functions.h"
#include <cassert>
#include <cstring>
#include <vector>
#include <cmath>
#include <algorithm>

int BallKmeans::runThread(int threadId, int maxIterations) {


    // track the number of iterations the algorithm performs
    int iterations = 0;

    int startNdx = start(threadId);
    int endNdx = end(threadId);

    int dataSize = x->n;
    int dimensions = x->d;

    Dataset *clusterCentroids = new Dataset(k, dimensions);

    int *clusterMemberCount = new int[k]{0};
    std::vector<int> *clusterMembers = new std::vector<int>[k];
    int *clusterRadius = new int[k]{0};

    //holds the index of the neighbouring cluster and distance from the current centroid.
    std::vector<pair<int,double>> *neighborClusters = new std::vector<pair<int,double>>[k];


    double *stableAreaRadius = new double[k]{0};
    std::vector<int> *annulusAreaRadius = new std::vector<int>[k];


    while ((iterations < maxIterations) && (! converged)) {
        ++iterations;


        for (int i = startNdx; i < endNdx; ++i){
            int clusterIndex = assignment[i];
            clusterMembers[clusterIndex].push_back(i);

            clusterMemberCount[clusterIndex]+=1;
            for (int j=0;j<dimensions;j++){
                //may need to initialize clusterCentroids to 0
                clusterCentroids->data[clusterIndex+j]+=x->data[i+j];

            }
        }

        for (int i=0;i<k;i++){
            for (int j=0;j<dimensions;j++){
                clusterCentroids->data[i+j]/=clusterMemberCount[i];
            }
        }


        for (int iter=0;iter<k;iter++){

            for (int i=0;i<clusterMembers[iter].size();i++){
                int dataIndex = clustermembers[iter].at(i);
                double temp = 0;
                for (int j=0;j<dimensions;j++){
                    double diff = x->data[dataIndex+j]-clusterCentroids->data[iter+j];
                    temp+=diff*diff;
                }
                temp = sqrt(temp);
                if (temp>clusterRadius[iter])
                    clusterRadius[iter]=temp;
            }

        }


        // call findNeighbors function here

        for (int iter=0;iter<k;iter++){
            sort(neighborClusters[iter].begin(), neighborClusters[iter].end(), sortPair);
        }

        for (int iter=0; iter<k;iter++){
                int closestClusterIndex = neighborClusters[iter].at(0).first;
                double temp = 0;
                for (int j=0;j<dimensions;j++){
                    double diff = clusterCentroids->data[iter+j] - clusterCentroids->data[closestClusterIndex+j];
                    temp+=diff*diff;
                }
                stableAreaRadius[iter] = 0.5 * sqrt(temp);

                //find annulus area
                int numberOfAnnulusArea = neighborClusters[iter].size();

                for (int areaIdx=0;areaIdx<numberOfAnnulusArea-1;areaIdx++){
                    // need to merge stable area radius here as well
                    annulusAreaRadius[iter].push_back(0.5 * neighborClusters[iter].at(areaIdx+1).second);
                }

                //TODO fill last annulus area


                for (int j=0;j<clusterMembers[iter].size();j++){
                    int dataIdx = clusterMembers[iter].at(j);
                    //calculate distance
                    double distancePointToCentroid = 0;
                    for (int idx=0;idx<k;idx++){
                        double diff = x->data[dataIdx+idx] - clusterCentroids->data[iter+j];
                        distancePointToCentroid+=diff * diff;
                    }
                    distancePointToCentroid = sqrt(distancePointToCentroid);

                    //step 12
                    if (distancePointToCentroid>stableAreaRadius[iter]){

                        //this means it is in annulus area, not stable area

                        int annulusAreaIdx=0;
                        while (distancePointToCentroid>annulusAreaRadius[iter].at(annulusAreaIdx)){

                            double distPointToNeighbour = 0;
                            int currentNeighbourIndex = neighborClusters[iter].at(annulusAreaIdx).first;
                            for (int idx=0;idx<k;idx++){
                                double diff = x->data[dataIdx+idx] - clusterCentroids->data[currentNeighbourIndex+j];
                                distPointToNeighbour+=diff * diff;
                            }
                            distPointToNeighbour = sqrt(distPointToNeighbour);

                            if(distPointToNeighbour<distancePointToCentroid){
                                assignment[dataIdx] = currentNeighbourIndex;
                            }

                            annulusAreaIdx++;


                        }





                    }
                }
        }


        // TODO fix this block

        synchronizeAllThreads();

        if (threadId == 0) {
            int furthestMovingCenter = move_centers();
            converged = (0.0 == centerMovement[furthestMovingCenter]);
        }

        synchronizeAllThreads();
    }

    return iterations;

}


void searchNeighborClusters(Dataset *centroids, Dataset *oldCentroids, int k, int *clusterRadius, int j,
        int dimension, std::vector<pair<int,double>> *neighborClusters){

    neighborClusters[j].clear();

    for (int i=0;i<k;i++){
        if (i!=j){
            double distance_t = calculateDistance(centroids->data,centroids->data,i,j,dimension);
            double deltaI = calculateDistance(centroids->data,oldCentroids->data,i,i,dimension);
            double deltaJ = calculateDistance(centroids->data,oldCentroids->data,j,j,dimension);
            double radius = clusterRadius[i];
            double distance_t1;

            if (distance_t>=2 * radius * deltaI * deltaJ)
                distance_t1 = distance_t - deltaI - deltaJ;
            else{
                distance_t1 = calculateDistance(centroids->data,centroids->data,i,j,dimension);
                if (distance_t1<2*radius){
                    //add the index of the cluster and its distance from the queried cluster
                    neighborClusters[j].push_back(std::make_pair(i,distance_t1));

                }
            }


        }

    }


}

bool sortPair(const pair<int,double> &a,
               const pair<int,double> &b)
{
    return (a.second < b.second);
}

double BallKmeans::calculateDistance(double *centroids, double *anotherCentroids int i, int j, int dimension){
    double distance = 0;
    for (int it=0:it<dimension;it++){
        double diff = centroids[i+it] - anotherCentroids[j+it];
        distance+=diff * diff;
    }
    return sqrt(distance);
}