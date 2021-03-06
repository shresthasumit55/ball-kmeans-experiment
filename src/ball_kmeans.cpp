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
   /* for (int iter = 0; iter < k; iter++) {
        for (int j = 0; j < dimensions; j++) {
            clusterCentroids->data[iter*dimensions + j] = centers->data[iter*dimensions + j];
        }
    }
    */

    for (int i = startNdx; i < endNdx; ++i) {
        // look for the closest center to this example
        int closest = 0;
        double closestDist2 = std::numeric_limits<double>::max();
        for (int j = 0; j < k; ++j) {
            double d2 = pointCenterDist2(i, j);
            if (d2 < closestDist2) {
                closest = j;
                closestDist2 = d2;
            }
        }
        if (assignment[i] != closest) {
            changeAssignment(i, closest, threadId);
        }
    }





    while ((iterations < maxIterations) && (! converged)) {
        ++iterations;

        for (int iter = 0; iter < k; iter++) {
            for (int j = 0; j < dimensions; j++) {
                oldCentroids->data[iter*dimensions + j] = centers->data[iter*dimensions + j];
            }
            clusterMembers[iter].clear();
            clusterMemberCount[iter] = 0;
            clusterRadius[iter]=0;
        }

        if (threadId == 0) {
            /*
            converged = true;

            int iter=0;
            while ((iter<k) && (converged)){
                double centersDistance = 0;
                for (int j = 0; j < dimensions; j++) {
                    double delta = clusterCentroids->data[iter*dimensions+j] - oldCentroids->data[iter*dimensions + j];
                    double delta2 = delta * delta;
                    centersDistance += delta2;
                }
                if (centersDistance != 0) {
                    converged = false;
                }
                centerMovement[iter] = sqrt(centersDistance);
                iter++;
            }
             */

            int furthestMovingCenter = move_centers();
            converged = (0.0 == centerMovement[furthestMovingCenter]);
        }




        for (int iter = 0; iter < k; iter++) {
            for (int j = 0; j < dimensions; j++) {
                clusterCentroids->data[iter*dimensions + j] = centers->data[iter*dimensions + j];
            }
        }

        for (int i=0;i<k;i++){
            for (int j=i+1;j<k;j++){
                centerDistances[i][j] = calculateDistance(clusterCentroids->data,clusterCentroids->data,i,j,dimensions);
                centerDistances[j][i] = centerDistances[i][j];
            }
        }

        //step 3-5: calculatng centroids of ball clusters
        for (int i = startNdx; i < endNdx; ++i){
            int clusterIndex = assignment[i];
            clusterMembers[clusterIndex].push_back(i);

            clusterMemberCount[clusterIndex]+=1;

            /*
            for (int j=0;j<dimensions;j++){
                clusterCentroids->data[clusterIndex*dimensions+j]+=x->data[i*dimensions+j];
            }
             */
        }

/*
        for (int i=0;i<k;i++){
            for (int j=0;j<dimensions;j++){
                clusterCentroids->data[i*dimensions+j]/=clusterMemberCount[i];
            }
        }

        */

        // step 6 begin
        for (int iter=0;iter<k;++iter) {

            int currentClusterSize = (int) (clusterMembers[iter].size());

            for (int i = 0; i < currentClusterSize; ++i) {
                int dataIndex = clusterMembers[iter].at(i);
                double temp = 0;
                for (int j = 0; j < dimensions; j++) {
                    double diff = x->data[dataIndex*dimensions + j] - clusterCentroids->data[iter*dimensions + j];
                    temp += diff * diff;
                }
                temp = sqrt(temp);
                if (temp > clusterRadius[iter]) {
                    clusterRadius[iter] = temp;
                }
            }


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

            // step 11: determining stableArea points



            int closestClusterIndex = neighborClusters[iter].at(0).first;
            double temp = 0;

            for (int j = 0; j < dimensions; j++) {
                double diff = clusterCentroids->data[iter*dimensions + j] - clusterCentroids->data[closestClusterIndex*dimensions + j];
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
                    double diff = x->data[dataIdx*dimensions + idx] - clusterCentroids->data[iter*dimensions + idx];
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
                            double diff = x->data[dataIdx*dimensions + idx] - clusterCentroids->data[currentNeighbourIndex*dimensions + j];
                            distPointToNeighbour += diff * diff;
                        }
                        distPointToNeighbour = sqrt(distPointToNeighbour);

                        if (distPointToNeighbour < distancePointToCentroid) {
                            changeAssignment(dataIdx,currentNeighbourIndex,threadId);
                            //assignment[dataIdx] = currentNeighbourIndex;
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


    }

    /*
     * Note: Since we are working in ball clusters, will need to assign the final centers to its corresponding
     *   data structure in the framework
     */

    /*
    for (int iter = 0; iter < k; iter++) {
        for (int j = 0; j < dimensions; j++) {
            //(*centers)(iter,j) = (*clusterCentroids)(iter,j);
            centers->data[iter*dimensions + j] = clusterCentroids->data[iter*dimensions + j];
        }
    }
     */
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
        double diff = centroids[i*dimension+it] - anotherCentroids[j*dimension+it];
        distance+=diff * diff;
    }
    return sqrt(distance);
}