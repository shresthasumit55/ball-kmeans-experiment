//
// Created by sumit on 3/16/20.
//

#include "ball_kmeans.h"
#include "general_functions.h"
#include <cassert>
#include <cstring>
#include <vector>
#include <cmath>

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


        // find neightbour clusters from step 9

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

                //TODO last annulus area fill


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



        // loop over all examples
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

        verifyAssignment(iterations, startNdx, endNdx);

        synchronizeAllThreads();

        if (threadId == 0) {
            int furthestMovingCenter = move_centers();
            converged = (0.0 == centerMovement[furthestMovingCenter]);
        }

        synchronizeAllThreads();
    }

    return iterations;

}


void searchNeighborClusters(Dataset *centroids, Dataset *oldCentroids, int k, int *clusterRadius, int iter){



    for (int i=0;i<k;i++){
        if (i!=iter){

        }

    }


}

bool sortPair(const pair<int,double> &a,
               const pair<int,double> &b)
{
    return (a.second < b.second);
}
