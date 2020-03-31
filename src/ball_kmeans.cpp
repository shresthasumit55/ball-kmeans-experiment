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

    //initializing clusterCentroids with initial centers
    for (int iter = 0; iter < k; iter++) {
        for (int j = 0; j < dimensions; j++) {
            clusterCentroids->data[iter + j] = centers->data[iter + j];
        }
    }

    for (int i=0;i<k;i++){
        for (int j=i+1;j<k;j++){
            centerDistances[i][j] = calculateDistance(clusterCentroids->data,clusterCentroids->data,i,j,dimensions);
            centerDistances[j][i] = centerDistances[i][j];
        }
    }

    while ((iterations < maxIterations) && (! converged)) {
        ++iterations;

        for (int i=0;i<k;i++){
            for (int j = 0; j < dimensions; j++) {
                oldCentroids->data[i+j] = clusterCentroids->data[i+j];
            }
            clusterMembers[i].clear();
            clusterMemberCount[i] = 0;
            clusterRadius[i]=0;
        }

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
        for (int iter=0;iter<k;++iter) {

            int currentClusterSize = (int) (clusterMembers[iter].size());

            for (int i = 0; i < currentClusterSize; ++i) {
                int dataIndex = clusterMembers[iter].at(i);
                double temp = 0;
                for (int j = 0; j < dimensions; j++) {
                    double diff = x->data[dataIndex + j] - clusterCentroids->data[iter + j];
                    temp += diff * diff;
                }
                temp = sqrt(temp);
                if (temp > clusterRadius[iter]) {
                    clusterRadius[iter] = temp;
                }
            }

            //}

            //for (int iter=0;iter<k;iter++){

            //step 8,9,10 in the algorithm: neighbour cluster search and sort
            neighborClusters[iter].clear();
            searchNeighborClusters(clusterCentroids, oldCentroids, k, clusterRadius, iter,
                                   dimensions, neighborClusters, centerDistances);

            if (neighborClusters[iter].size() != 0) {

            std::vector <std::pair<int, double>> currentNeighbors = neighborClusters[iter];
            std::sort(currentNeighbors.begin(), currentNeighbors.end(),
                      [](const std::pair<int, double> &left, const std::pair<int, double> &right) {
                          return left.second < right.second;
                      });
            //}


            //for (int iter=0; iter<k;iter++){
            // step 11: determining stableArea points



            int closestClusterIndex = neighborClusters[iter].at(0).first;
            double temp = 0;

            for (int j = 0; j < dimensions; j++) {
                double diff = clusterCentroids->data[iter + j] - clusterCentroids->data[closestClusterIndex + j];
                temp += diff * diff;
            }
            stableAreaRadius[iter] = 0.5 * sqrt(temp);


                //find annulus area
            int numberOfAnnulusArea = neighborClusters[iter].size();
            annulusAreaRadius[iter].clear();

            for (int areaIdx = 0; areaIdx < (numberOfAnnulusArea-1); areaIdx++) {
                // need to merge stable area radius here as well
                annulusAreaRadius[iter].push_back(0.5 * neighborClusters[iter].at(areaIdx + 1).second);
            }

                //this corresponds the actual radius of the current cluster.
            annulusAreaRadius[iter].push_back(clusterRadius[iter]);

            //step 12: assigning points to annulus areas
            for (int j = 0; j < currentClusterSize; j++) {
                int dataIdx = clusterMembers[iter].at(j);
                //calculate distance
                double distancePointToCentroid = 0;
                for (int idx = 0; idx < dimensions; idx++) {
                    double diff = x->data[dataIdx + idx] - clusterCentroids->data[iter + idx];
                    distancePointToCentroid += diff * diff;
                }
                distancePointToCentroid = sqrt(distancePointToCentroid);

                //step 12: checking if point is closer to any neighbour cluster
                if (distancePointToCentroid > stableAreaRadius[iter]) {

                    //this means it is in annulus area, not stable area
                    int annulusAreaIdx = 0;

                    bool stopNeighborCompare = false;
                    while (!stopNeighborCompare) {

                        double distPointToNeighbour = 0;
                        int currentNeighbourIndex = neighborClusters[iter].at(annulusAreaIdx).first;

                        for (int idx = 0; idx < k; idx++) {
                            double diff = x->data[dataIdx + idx] - clusterCentroids->data[currentNeighbourIndex + j];
                            distPointToNeighbour += diff * diff;
                        }
                        distPointToNeighbour = sqrt(distPointToNeighbour);

                        if (distPointToNeighbour < distancePointToCentroid) {
                            assignment[dataIdx] = currentNeighbourIndex;
                        }

                        if (distancePointToCentroid > annulusAreaRadius[iter].at(annulusAreaIdx)) {
                            annulusAreaIdx++;
                        }else{
                            stopNeighborCompare = true;
                        }
                }

                }

            }
        }
        }



        // TODO fix this block

        synchronizeAllThreads();

        if (threadId == 0) {
            converged = true;

            int iter=0;
            while ((iter<k) && (converged)){
                double centersDistance = 0;
                for (int j = 0; j < dimensions; j++) {
                    double delta = clusterCentroids->data[iter+j] - oldCentroids->data[iter + j];
                    double delta2 = delta * delta;
                    centersDistance += delta2;
                }
                if (centersDistance != 0) {
                    converged = false;
                }
                centerMovement[iter] = sqrt(centersDistance);
                iter++;
            }

            //int furthestMovingCenter = move_centers();
            //converged = (0.0 == centerMovement[furthestMovingCenter]);
        }

        synchronizeAllThreads();
    }

    /*
     * Note: Since we are working in ball clusters, will need to assign the final centers to its corresponding
     *   data structure in the framework
     */

    for (int iter = 0; iter < k; iter++) {
        for (int j = 0; j < dimensions; j++) {
            centers->data[iter + j] = clusterCentroids->data[iter + j];
        }
    }
    return iterations;

}


void BallKmeans::searchNeighborClusters(Dataset *centroids, Dataset *oldCentroids, int k, double *clusterRadius, int j,
        int dimension, std::vector<std::pair<int,double>> *neighborClusters, double **centerDistances){

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

double BallKmeans::calculateDistance(double *centroids, double *anotherCentroids, int i, int j, int dimension) {
    double distance = 0;
    for (int it=0;it<dimension;it++){
        double diff = centroids[i+it] - anotherCentroids[j+it];
        distance+=diff * diff;
    }
    return sqrt(distance);
}