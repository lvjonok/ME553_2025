#ifndef DATA_HPP
#define DATA_HPP

#include "model.hpp"

class Data {
public:
  Data(const Model &model) {};

  // links with respect to world frame
  std::vector<Transform> oTb;

  // joints with respect to world frame
  std::vector<Transform> oTj;

  // frames velocities
  std::vector<Eigen::Vector3d> bodyLinVel_w;
  std::vector<Eigen::Vector3d> bodyAngVel_w;

  // link frames
  std::vector<Transform> oTw;

  // jacobian
  Eigen::MatrixXd posJacobian;
  Eigen::MatrixXd rotJacobian;
};

#endif // DATA_HPP