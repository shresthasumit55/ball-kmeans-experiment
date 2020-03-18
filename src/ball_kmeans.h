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
};

#endif //BALL_K_MEANS_BALL_KMEANS_H
