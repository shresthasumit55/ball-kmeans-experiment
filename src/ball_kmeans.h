//
// Created by sumit on 3/16/20.
//

#ifndef BALL_K_MEANS_BALL_KMEANS_H
#define BALL_K_MEANS_BALL_KMEANS_H


#include "original_space_kmeans.h"
#include <vector>

class BallKmeans : public OriginalSpaceKmeans {
public:
    virtual std::string getName() const { return "ball"; }
    virtual ~BallKmeans() { free(); }
protected:
    virtual int runThread(int threadId, int maxIterations);
private:
    double calculateDistance(double *centroids, double *anotherCentroids, int i, int j, int dimension);
    void searchNeighborClusters(Dataset *centroids, Dataset *oldCentroids, int k, double *clusterRadius, int j,
                                int dimension, std::vector<std::pair<int,double>> *neighborClusters,
                                double **centerDistances);
    //bool sortPairs(const std::pair<int,double> &a,const std::pair<int,double> &b);
};

#endif //BALL_K_MEANS_BALL_KMEANS_H
