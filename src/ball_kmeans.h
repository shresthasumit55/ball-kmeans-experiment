//
// Created by sumit on 3/16/20.
//

#ifndef BALL_K_MEANS_BALL_KMEANS_H
#define BALL_K_MEANS_BALL_KMEANS_H


#include "original_space_kmeans.h"

class BallKmeans : public OriginalSpaceKmeans {
public:
    virtual std::string getName() const { return "ball"; }
    virtual ~BallKmeans() { free(); }
protected:
    virtual int runThread(int threadId, int maxIterations);
private:
    double calculateDistance(double *centroids, int i, int j, int dimension);
    void searchNeighborClusters(Dataset *centroids, Dataset *oldCentroids, int k, int *clusterRadius, int j,
                                int dimension, std::vector<pair<int,double>> *neighborClusters);
    double BallKMeans::bool sortPair(const pair<int,double> &a,
                                     const pair<int,double> &b);
};

#endif //BALL_K_MEANS_BALL_KMEANS_H
