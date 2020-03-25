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

    //int dataSize = x->n;
    int dimensions = x->d;

    Dataset *clusterCentroids = new Dataset(k, dimensions);
    Dataset *oldCentroids = new Dataset(k, dimensions);

    int *clusterMemberCount = new int[k]{0};
    std::vector<int> *clusterMembers = new std::vector<int>[k];
    double *clusterRadius = new double[k]{0};

    //holds the index of the neighbouring cluster and distance from the current centroid.
    std::vector<std::pair<int,double>> *neighborClusters = new std::vector<std::pair<int,double>>[k];


    double *stableAreaRadius = new double[k]{0};
    std::vector<double> *annulusAreaRadius = new std::vector<double>[k];

    double **centerDistances = new double *[k];
    for (int iter = 0; iter < k; iter++) {
        centerDistances[iter] = new double[k];
    }


    while ((iterations < maxIterations) && (! converged)) {
        ++iterations;

        oldCentroids = clusterCentroids;

        //step 3-5: calculatng centroids of ball clusters
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

        // step 6 begin
        for (int iter=0;iter<k;iter++){

            for (int i=0;i<(int)(clusterMembers[iter].size());i++){
                int dataIndex = clusterMembers[iter].at(i);
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

        for (int iter=0;iter<k;iter++){
            //step 8,9,10 in the algorithm: neighbour cluster search and sort
            searchNeighborClusters(clusterCentroids,oldCentroids,k,clusterRadius,iter,
                    dimensions,neighborClusters,centerDistances);
            std::vector<std::pair<int,double>> currentNeighbors = neighborClusters[iter];
            std::sort(currentNeighbors.begin(), currentNeighbors.end(), [](const std::pair<int,double> &left, const std::pair<int,double> &right) {
                return left.second < right.second;
            });
        }

        for (int iter=0; iter<k;iter++){
            // step 11: determining stableArea points
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
                //this corresponds the actual radius of the current cluster.
                annulusAreaRadius[iter].push_back(clusterRadius[iter]);


                //step 12: assigning points to annulus areas
                for (int j=0;j<(int)(clusterMembers[iter].size());j++){
                    int dataIdx = clusterMembers[iter].at(j);
                    //calculate distance
                    double distancePointToCentroid = 0;
                    for (int idx=0;idx<k;idx++){
                        double diff = x->data[dataIdx+idx] - clusterCentroids->data[iter+j];
                        distancePointToCentroid+=diff * diff;
                    }
                    distancePointToCentroid = sqrt(distancePointToCentroid);

                    //step 12: checking if point is closer to any neighbour cluster
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


void BallKmeans::searchNeighborClusters(Dataset *centroids, Dataset *oldCentroids, int k, double *clusterRadius, int j,
        int dimension, std::vector<std::pair<int,double>> *neighborClusters, double **centerDistances){

    neighborClusters[j].clear();

    for (int i=0;i<k;i++){
        if (i!=j){
           // double distance_t = calculateDistance(centroids->data,centroids->data,i,j,dimension);
            double distance_t = centerDistances[j][i];
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
            centerDistances[j][i] = distance_t1;



        }

    }


}

/*
bool BallKmeans::sortPairs(const std::pair<int,double> &a,const std::pair<int,double> &b)
{
    return (a.second < b.second);
}
 */

double BallKmeans::calculateDistance(double *centroids, double *anotherCentroids, int i, int j, int dimension) {
    double distance = 0;
    for (int it=0;it<dimension;it++){
        double diff = centroids[i+it] - anotherCentroids[j+it];
        distance+=diff * diff;
    }
    return sqrt(distance);
}