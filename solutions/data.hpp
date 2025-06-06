#ifndef DATA_HPP
#define DATA_HPP

#include "model.hpp"

namespace algorithms {

class Data {
public:
  Data(const Model &model) {
    iTj_.resize(model.actuated_joints_.size());
    iXj_.resize(model.actuated_joints_.size());
    iX0_.resize(model.actuated_joints_.size());
    vJ_.resize(model.actuated_joints_.size());
    for (size_t i = 0; i < model.actuated_joints_.size(); i++) {
      iTj_[i] = Transform::Identity();
      iXj_[i] = SpatialTransform::Identity();
      iX0_[i] = SpatialTransform::Identity();
      vJ_[i].setZero();
    }

    rot_WB.resize(model.actuated_joints_.size());
    jointPos_W.resize(model.actuated_joints_.size());
    jointAxis_W.resize(model.actuated_joints_.size());
    bodyLinVel_w.resize(model.actuated_joints_.size());
    bodyAngVel_w.resize(model.actuated_joints_.size());
    joint2joint_W.resize(model.actuated_joints_.size());
    S.resize(model.actuated_joints_.size());
    for (size_t i = 0; i < model.actuated_joints_.size(); i++) {
      rot_WB[i] = Eigen::Matrix3d::Identity();
      jointPos_W[i] = Eigen::Vector3d::Zero();
      jointAxis_W[i] = Eigen::Vector3d::Zero();
      bodyLinVel_w[i] = Eigen::Vector3d::Zero();
      bodyAngVel_w[i] = Eigen::Vector3d::Zero();
      joint2joint_W[i] = Eigen::Vector3d::Zero();
    }

    ST.resize(model.actuated_joints_.size());
    for (size_t i = 0; i < model.actuated_joints_.size(); i++) {
      ST[i].setZero();
    }
  };

  // the position of the body frame j relative to the joint frame i
  std::vector<Transform> iTj_;
  std::vector<SpatialTransform> iXj_;

  std::vector<SpatialTransform> iX0_; // the position of the body frame i
  // relative to the world frame

  std::vector<Motion> vJ_; // spatial velocity of the body frame j

  // EVERYTHING IN THE WORLD FRAME BELOW
  std::vector<Eigen::Matrix3d>
      rot_WB; // rotation of the body i in the world frame
  std::vector<Eigen::Vector3d>
      jointPos_W; // position of the joint i in the world frame
  std::vector<Eigen::Vector3d> jointAxis_W; // axis of the joint i in the world
  // vector between the joint i and its parent joint in the world frame
  std::vector<Eigen::Vector3d> joint2joint_W;
  // motion subspace of the body i in the world frame
  std::vector<Eigen::MatrixXd> S;

  // velocities
  std::vector<Eigen::Vector3d> bodyLinVel_w; // linear velocity of the body
  std::vector<Eigen::Vector3d> bodyAngVel_w; // angular velocity of the body

  // dynamic quantities
  std::vector<Eigen::Matrix3d>
      inertiaW; // inertia of the body i in the world frame
  std::vector<Eigen::Vector3d>
      comW; // center of mass of the body i in the world frame
  std::vector<Eigen::MatrixXd> spatialInertia6; // spatial inertia of the body i

  // composite mass, com and inertia of i-th composite body in the world frame
  std::vector<Eigen::Matrix3d> compositeInertiaW;
  std::vector<Eigen::Vector3d> compositeComW;
  std::vector<double> compositeMassW;

  // crba
  std::vector<Eigen::MatrixXd> spatialCompositeInertia6;
  Eigen::MatrixXd massMatrix; // mass matrix

  // rnea
  std::vector<Eigen::Vector3d> bodyLinAcc;
  std::vector<Eigen::Vector3d> bodyAngAcc;
  std::vector<Eigen::MatrixXd> dS;
  Eigen::MatrixXd nonlinearities;

  // aba

  // velocity and acceleration matched with pinocchio
  std::vector<Eigen::VectorXd> ov, oa_gf, oh, of;
  Eigen::VectorXd u, dv;
  std::vector<Eigen::MatrixXd> aXb, U, Dinv, Jcols;

  std::vector<Eigen::Matrix<double, 1, 6>> ST; // S^T term for each joint

  // plucker transform with [[I, -skew(r)], [0, I]]
  std::vector<Eigen::Matrix<double, 6, 6>> XT;

  // articulated inertias and bias
  std::vector<Eigen::Matrix<double, 6, 6>> Ma;
  std::vector<Eigen::Matrix<double, 6, 1>> Pa;
  std::vector<Eigen::MatrixXd> STMaXT; // S^T * Ma * X^T
  std::vector<Eigen::MatrixXd> STMa;   // S^T * Ma
  std::vector<Eigen::MatrixXd> STMaSinvSTMaXT;
  std::vector<Eigen::MatrixXd> SdotUpXdotTV;
  std::vector<double> STMaSinv;
  std::vector<Eigen::MatrixXd> Xbp, dXbpT; // plucker
  std::vector<Eigen::MatrixXd> W;          // simply twist
  std::vector<Eigen::VectorXd> Wdot;       // acc of body
  Eigen::VectorXd udot; // acceleration of the generalized coordinates

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

} // namespace algorithms

#endif // DATA_HPP