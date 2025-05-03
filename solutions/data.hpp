#ifndef DATA_HPP
#define DATA_HPP

#include "model.hpp"

class Data {
public:
  Data(const Model &model) {
    iXj_.resize(model.actuated_joints_.size());
    for (size_t i = 0; i < model.actuated_joints_.size(); i++) {
      iXj_[i] = Transform::Identity();
    }
  };

  // the position of the body frame j relative to the joint frame i
  std::vector<Transform> iXj_;

  // // links with respect to world frame
  // std::vector<Transform> oTb;

  // // joints with respect to world frame
  // std::vector<Transform> oTj;

  // // frames velocities
  // std::vector<Eigen::Vector3d> bodyLinVel_w;
  // std::vector<Eigen::Vector3d> bodyAngVel_w;

  // // link frames
  // std::vector<Transform> oTw;

  // // jacobian
  // Eigen::MatrixXd posJacobian;
  // Eigen::MatrixXd rotJacobian;

  // // crba
  // std::vector<Transform> comW; // center of mass in world frame for each link
  // std::vector<Eigen::Matrix3d> inertiaW; // inertia in world frame for each
  // link

  // // composite inertia in world frame of each subtree rooted at the link
  // std::vector<Eigen::Matrix3d> compositeInertiaW;
  // std::vector<double> compositeMassW; // mass of each subtree rooted at the
  // link std::vector<Eigen::Vector3d>
  //     compositeComW; // com of each subtree rooted at the link

  // // composite inertia in world frame of each subtree rooted at the link
  // std::vector<Eigen::MatrixXd> compositeSpatialInertia6;
  // std::vector<Eigen::Matrix<double, 6, 1>> S_W;

  // Eigen::MatrixXd massMatrix; // mass matrix
};

#endif // DATA_HPP