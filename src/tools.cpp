#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
    
    VectorXd rmse(4);
    rmse << 0,0,0,0;
    
    if (estimations.size() == 0 || estimations.size() != ground_truth.size()){
        return rmse;
    }
    
    int n = estimations.size();
    
    //std::cout << "Size : " << estimations.size() << std::endl;
    
    //accumulate squared residuals
    for(int i=0; i < estimations.size(); ++i){
        VectorXd e = estimations[i];
        VectorXd g = ground_truth[i];

        VectorXd diff = e - g;

        VectorXd diffSqrd = diff.array() * diff.array();
        rmse += diffSqrd;
        
    }
    
    //calculate the mean
    // ... your code here
    rmse = rmse.array()/n;
    rmse = rmse.array().sqrt();
    //calculate the squared root
    // ... your code here
    
    //return the result
    return rmse;
}
