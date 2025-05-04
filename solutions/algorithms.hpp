#ifndef ME553_2025_ALGORITHMS_HPP_
#define ME553_2025_ALGORITHMS_HPP_

#include "data.hpp"
#include "model.hpp"
#include <Eigen/Core>
#include <cstddef>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <tinyxml_rai/tinyxml_rai.h>
#include <vector>

namespace algorithms {

// TODO: this is a test, what if we compute the forward kinematics completely
// bare-bone, directly utilize model and joints and use only vectors and
// matrices
inline void setState(const Model &model, Data &data, const Eigen::VectorXd &gc,
                     const Eigen::VectorXd &gv) {
  Eigen::Vector3d pPosW = Eigen::Vector3d::Zero();
  Eigen::Matrix3d pRotW = Eigen::Matrix3d::Identity();
  for (size_t i = 0; i < model.nbodies_; i++) {
    auto parentId = model.parents_[i];
    if (parentId != -1) {
      pPosW = data.jointPos_W[parentId];
      pRotW = data.rot_WB[parentId];
    }

    auto joint = model.actuated_joints_[i];
    // this is the static rotation of the joint from parent body to the child
    // (joint) frame
    Eigen::Matrix3d R_PJ = joint->jointPlacement().block<3, 3>(0, 0);
    Eigen::Vector3d P_PJ = joint->jointPlacement().block<3, 1>(0, 3);

    auto motion = joint->motion(gc, model.gc_idx_[i], joint->gc_length());
    Eigen::Matrix3d R_J = motion.block<3, 3>(0, 0);
    Eigen::Vector3d P_J = motion.block<3, 1>(0, 3);

    // perform the transformation
    data.rot_WB[i] = pRotW * R_PJ * R_J;
    data.jointPos_W[i] = pPosW + pRotW * (P_PJ + R_PJ * P_J);
  }
}

// inline void setState(const Model &model, Data &data, const Eigen::VectorXd
// &gc,
//                      const Eigen::VectorXd &gv) {
//   for (size_t i = 0; i < model.nbodies_; i++) {
//     // joint model is identity right now, we don't consider the generalized
//     // coordinate yet
//     auto joint = model.actuated_joints_[i];
//     auto body = model.bodies_[i];

//     Transform T_J = Transform::Identity();
//     SpatialTransform X_J = SpatialTransform::Identity();
//     // Eigen::MatrixXd S;
//     Motion vJ = Motion::Zero();
//     if (joint->getType() != JointType::FIXED) {
//       // std::cout << "joint " << joint->getName() << " uses gc "
//       //           << gc[model.gc_idx_[i]] << " gv " <<
//       gv[model.gv_idx_[i]]
//       //           << std::endl;

//       auto res = joint->jcalc(gc, model.gc_idx_[i], gv, model.gv_idx_[i]);
//       // T_J = joint->jcalcOld(gc, model.gc_idx_[i]);

//       X_J = std::get<0>(res);
//       vJ = std::get<1>(res);
//     } else {
//       std::cout << "Use identity for joint " << joint->getName() <<
//       std::endl;
//     }
//     Transform iTp = model.T_T_[i];
//     SpatialTransform iXp = model.X_T_[i];
//     size_t parentId = model.parents_[i];

//     // std::cout << "body " << body->getName() << " iTp:\n"
//     //           << iTp << "\n "
//     //           << "iXp:\n"
//     //           << iXp << std::endl;

//     if (parentId != -1) {
//       data.iTj_[i] = data.iTj_[parentId] * iTp * T_J;
//       data.iXj_[i] = data.iXj_[parentId] * iXp * X_J;

//       SpatialTransform Xup = iXp * X_J;

//       data.vJ_[i] = data.vJ_[parentId] + data.iXj_[i] * vJ;
//     } else {
//       data.iTj_[i] = iTp * T_J;
//       data.iXj_[i] = iXp * X_J;

//       // velocity
//       data.vJ_[i] = Motion::Zero();
//     }
//   }

//   return;
// }

inline void getBodyPose(const Model &model, Data &data, size_t bodyId,
                        Eigen::Matrix3d &R, Eigen::Vector3d &p) {
  // get the pose of the body in the world frame
  R = data.rot_WB[bodyId];
  p = data.jointPos_W[bodyId];
}

// inline void getBodyTwist(const Model &model, Data &data, size_t bodyId,
//                          Motion &twist) {
//   // get the twist of the body in the world frame
//   auto body = model.bodies_[bodyId];

//   auto iXj = data.iXj_[bodyId];
//   auto R = iXj.block<3, 3>(0, 0);

//   twist = data.vJ_[bodyId];
// }

// inline void framesForwardKinematics(const Model &model, Data &data,
//                                     const Eigen::VectorXd &gc) {
//   // compute the forward kinematics of the robot
//   data.oTb.resize(model.links_.size());
//   data.oTb[0] = Transform::Identity();

//   data.oTj.resize(model.joints_.size());
//   data.oTj[0] = Transform::Identity();

//   // first, update all the placements of the joints
//   for (size_t i = 0; i < model.joints_.size(); i++) {
//     auto joint = model.joints_[i];
//     auto parentFrame = data.oTb[joint->parent->getIndex()];
//     Transform jointFrame = parentFrame * joint->jointPlacement();

//     // if joint is not fixed, apply the motion
//     if (model.idx_qs_[i] != -1) {
//       Transform motion =
//           joint->motion(gc, model.idx_qs_[i], joint->gc_length());
//       jointFrame = jointFrame * motion;
//     }

//     data.oTb[joint->child->getIndex()] = jointFrame;

//     data.oTj[i] = jointFrame;
//   }
// }

// inline void jointAxisW(const Model &model, Data &data, size_t jointId,
//                        Eigen::Vector3d &axisW) {
//   // get the axis of the joint in the world frame
//   auto joint = data.oTj[jointId];
//   Eigen::Matrix3d R = joint.block<3, 3>(0, 0);
//   axisW = R * model.joints_[jointId]->getAxis();
// }

// inline void checkStructure(const Model &model, Data &data, size_t bodyId) {
//   // we simulate as we are doing backward pass from the bodyid to the root
//   // check that with parentDof we only get the actuated joints (i.e. not
//   fixed)

//   data.rotJacobian = Eigen::MatrixXd::Zero(3, model.nv_);
//   data.posJacobian = Eigen::MatrixXd::Zero(3, model.nv_);
//   auto joint = model.joints_[bodyId];

//   //   // skip fixed bodies
//   //   while (joint->getType() == JointType::FIXED) {
//   //     std::cout << "skipping joint " << joint->getName() << std::endl;
//   //     joint = model.joints_[joint->getParentDof()];
//   //   }

//   Eigen::Vector3d axisW;

//   auto parentDof = joint->getParentDof();
//   while (parentDof >= 0) {
//     auto columnIdx = model.idx_vs_[parentDof];

//     std::cout << "joint name " << joint->getName()
//               << " parentdof: " << joint->getParentDof()
//               << " column: " << columnIdx << std::endl;

//     jointAxisW(model, data, parentDof, axisW);

//     auto posCurrent = data.oTj[bodyId].block<3, 1>(0, 3);
//     auto posParent = data.oTj[joint->getParentDof()].block<3, 1>(0, 3);
//     auto r = posCurrent - posParent;

//     // fill if the joint is revolute
//     data.rotJacobian.block<3, 1>(0, columnIdx) = axisW;
//     // if (joint->getType() == JointType::REVOLUTE) {
//     // } else {
//     //   std::cout << "skipping joint " << joint->getName()
//     //             << " column: " << columnIdx << std::endl;
//     // }

//     // positional jacobian
//     if (joint->getType() == JointType::REVOLUTE) {
//       data.posJacobian.block<3, 1>(0, model.idx_vs_[joint->getParentDof()])
//       =
//           axisW.cross(r);
//     } else if (joint->getType() == JointType::PRISMATIC) {
//       data.posJacobian.block<3, 1>(0, model.idx_vs_[joint->getParentDof()])
//       =
//           axisW;
//     } else {
//       std::cout << "skipping joint " << joint->getName()
//                 << " column: " << columnIdx << std::endl;
//     }

//     // kinematic chain iteration
//     joint = model.joints_[parentDof];
//     parentDof = joint->getParentDof();
//   }

//   //   for (auto joint = model.joints_[bodyId]; joint->parent->getName() !=
//   //   "world";
//   //        joint = model.joints_[joint->getParentDof()]) {

//   //     std::cout << "joint name " << joint->getName() << " parentdof: "
//   //               << joint->getParentDof()
//   //               //   << " start gv: " <<
//   model.idx_vs_[joint->getParentDof()]
//   //               << std::endl;
//   //   }
// }

// inline void frameRotJacobian(const Model &model, Data &data,
//                              const Eigen::VectorXd &gc, size_t bodyId) {

//   // framesForwardKinematics should have been called
//   data.rotJacobian = Eigen::MatrixXd::Zero(3, model.nv_);
//   Eigen::Vector3d axisW;

//   auto joint = model.joints_[bodyId];

//   for (auto joint = model.joints_[bodyId]; joint->getParentDof() != 0;
//        joint = model.joints_[joint->getParentDof()]) {

//     // if this is not generalized, skip
//     if (model.idx_vs_[joint->getParentDof()] == -1) {
//       // std::cout << "skipping joint " << joint->getName() << std::endl;
//       continue;
//     }

//     auto posCurrent = data.oTj[bodyId].block<3, 1>(0, 3);
//     auto posParent = data.oTj[joint->getParentDof()].block<3, 1>(0, 3);
//     auto r = posCurrent - posParent;

//     jointAxisW(model, data, joint->getParentDof(), axisW);

//     // std::cout << "fill the data for joint name " << joint->getName()
//     //           << " parentdof: " << joint->getParentDof()
//     //           << " start gv: " << model.idx_vs_[joint->getParentDof()]
//     //           << " axisW: " << axisW.transpose() << std::endl;

//     // find the type of the generalized joint

//     if (joint->getType() == JointType::REVOLUTE || joint->getParentDof() ==
//     7) {
//       data.rotJacobian.block<3, 1>(0, model.idx_vs_[joint->getParentDof()])
//       =
//           axisW;
//     }
//   }
// }

// inline void framePosJacobian(const Model &model, Data &data,
//                              const Eigen::VectorXd &gc, size_t bodyId) {

//   // framesForwardKinematics should have been called
//   data.posJacobian = Eigen::MatrixXd::Zero(3, model.nv_);
//   Eigen::Vector3d axisW;

//   auto joint = model.joints_[bodyId];

//   for (auto joint = model.joints_[bodyId]; joint->getParentDof() != 0;
//        joint = model.joints_[joint->getParentDof()]) {

//     // if this is not generalized, skip
//     if (model.idx_vs_[joint->getParentDof()] == -1) {
//       // std::cout << "skipping joint " << joint->getName() << std::endl;
//       continue;
//     }

//     auto posCurrent = data.oTj[bodyId].block<3, 1>(0, 3);
//     auto posParent = data.oTj[joint->getParentDof()].block<3, 1>(0, 3);
//     auto r = posCurrent - posParent;

//     jointAxisW(model, data, joint->getParentDof(), axisW);

//     // std::cout << "fill the data for joint name " << joint->getName()
//     //           << " parentdof: " << joint->getParentDof()
//     //           << " start gv: " << model.idx_vs_[joint->getParentDof()]
//     //           << std::endl;

//     // find the type of the generalized joint

//     if (joint->getParentDof() == 11) {
//       // TODO: hack for the mess at the end effector

//       jointAxisW(model, data, 11, axisW);
//       data.posJacobian.block<3, 1>(0, model.idx_vs_[joint->getParentDof()])
//       =
//           axisW;
//     } else {
//       data.posJacobian.block<3, 1>(0, model.idx_vs_[joint->getParentDof()])
//       =
//           axisW.cross(r);
//     }
//   }
// }

// inline void forwardKinematics(const Model &model, Data &data,
//                               const Eigen::VectorXd &gc) {
//   Eigen::VectorXd v = Eigen::VectorXd::Zero(model.nv_);
//   algorithms::framesForwardKinematics(model, data, gc);
// }

// // CRBA algorithm

// // this function should be called after the forward kinematics
// // it computes the inertia of each link in the world frame
// inline void inertiaUpdate(const Model &model, Data &data,
//                           const Eigen::VectorXd &gc) {
//   data.inertiaW.resize(model.links_.size());
//   data.comW.resize(model.links_.size());

//   // find bodies Inertia expressed in the world frame
//   for (size_t i = 0; i < model.links_.size(); i++) {
//     auto link = model.links_[i];
//     auto inertia = link->getInertia();

//     auto worldInertia = data.oTb[i].block<3, 3>(0, 0) * inertia *
//                         data.oTb[i].block<3, 3>(0, 0).transpose();

//     data.inertiaW[i] = worldInertia;

//     // find the center of mass in the world frame
//     Transform com = link->getTransform();
//     Transform comW = data.oTb[i] * com;

//     data.comW[i] = comW; //.block<3, 1>(0, 3);
//   }

//   // now we compute the composite inertia of each body
//   data.compositeInertiaW.resize(model.links_.size());
//   data.compositeMassW.resize(model.links_.size());
//   data.compositeComW.resize(model.links_.size());
//   data.compositeSpatialInertia6.resize(model.links_.size());

//   auto skew = [](const Eigen::Vector3d &v) {
//     Eigen::Matrix3d skew;
//     skew.setZero();
//     skew(0, 1) = -v(2);
//     skew(0, 2) = v(1);
//     skew(1, 0) = v(2);
//     skew(1, 2) = -v(0);
//     skew(2, 0) = -v(1);
//     skew(2, 1) = v(0);
//     return skew;
//   };

//   for (int i = model.links_.size() - 1; i >= 1; i--) {
//     auto link = model.links_[i];
//     Eigen::Vector3d com = data.comW[i].block<3, 1>(0, 3);

//     auto mass = link->getMass();

//     data.compositeMassW[i] = mass;
//     data.compositeComW[i] = com;
//     data.compositeInertiaW[i] = data.inertiaW[i];

//     for (auto childJoint : model.getAttachedJoints(link)) {
//       Eigen::Vector3d c1 = data.compositeComW[i];
//       auto m1 = data.compositeMassW[i];
//       Eigen::Matrix3d I1 = data.compositeInertiaW[i];

//       auto childIdx = childJoint->child->getIndex();

//       Eigen::Vector3d c2 = data.compositeComW[childIdx];
//       auto m2 = data.compositeMassW[childIdx];
//       Eigen::Matrix3d I2 = data.compositeInertiaW[childIdx];

//       Eigen::Vector3d r_com_new = (m1 * c1 + m2 * c2) / (m1 + m2);
//       Eigen::Vector3d r1 = c1 - r_com_new;
//       Eigen::Vector3d r2 = c2 - r_com_new;

//       Eigen::Matrix3d I_new =
//           I1 + I2 - m1 * skew(r1) * skew(r1) - m2 * skew(r2) * skew(r2);

//       data.compositeMassW[i] = m1 + m2;
//       data.compositeComW[i] = r_com_new;
//       data.compositeInertiaW[i] = I_new;
//     }

//     // std::cout << "updated composite mass of index " << i << ": "
//     //           << data.compositeMassW[i] << std::endl;
//     // std::cout << "updated composite com of index " << i << ": "
//     //           << data.compositeComW[i].transpose() << std::endl;

//     // std::cout << "updated composite inertia of index " << i << ": " <<
//     // std::endl
//     //           << data.compositeInertiaW[i] << std::endl;
//   }

//   Eigen::Vector3d d = data.compositeComW[1];
//   Eigen::Matrix3d I_aboutOrigin =
//       data.compositeInertiaW[1] +
//       data.compositeMassW[1] *
//           (d.squaredNorm() * Eigen::Matrix3d::Identity() - d *
//           d.transpose());

//   data.compositeInertiaW[1] = I_aboutOrigin;

//   // std::cout << "updated composite inertia about origin of index " << 1
//   <<
//   ":
//   // "
//   //           << std::endl
//   //           << I_aboutOrigin << std::endl;
// }

// inline Eigen::Matrix<double, 6, 6> spatialTransform(const Transform &T) {
//   Eigen::Matrix<double, 6, 6> X;
//   X.setZero();
//   const auto R = T.block<3, 3>(0, 0);
//   const auto p = T.block<3, 1>(0, 3);
//   Eigen::Matrix3d px;
//   px << 0, -p.z(), p.y(), p.z(), 0, -p.x(), -p.y(), p.x(), 0;
//   X.topLeftCorner<3, 3>() = R;
//   X.topRightCorner<3, 3>().setZero();
//   X.bottomLeftCorner<3, 3>() = px * R;
//   X.bottomRightCorner<3, 3>() = R;
//   return X;
// }

// inline void computeMotionSubspaces(const Model &model, Data &data) {
//   data.S_W.assign(model.nv_, Eigen::Matrix<double, 6, 1>::Zero());
//   for (size_t i = 0; i < model.joints_.size(); ++i) {
//     int dof = model.idx_vs_[i];
//     if (dof < 0)
//       continue; // fixed joint → no generalized DOF

//     // 1) build S_i in joint‐local
//     Eigen::Matrix<double, 6, 1> S_loc;
//     const auto axis = model.joints_[i]->getAxis();
//     if (model.joints_[i]->getType() == JointType::REVOLUTE) {
//       S_loc.head<3>() = axis;    // ω
//       S_loc.tail<3>().setZero(); //   0
//     } else {                     // PRISMATIC
//       S_loc.head<3>().setZero();
//       S_loc.tail<3>() = axis; // v
//     }

//     // 2) transform into world
//     const auto X = spatialTransform(data.oTj[i]);
//     data.S_W[dof] = X * S_loc;
//   }
// }

// inline void buildCompositeSpatialInertia(const Model &model, Data &data) {
//   data.compositeSpatialInertia6.resize(model.links_.size());
//   auto skew = [&](const Eigen::Vector3d &v) {
//     Eigen::Matrix3d S;
//     S << 0, -v.z(), v.y(), v.z(), 0, -v.x(), -v.y(), v.x(), 0;
//     return S;
//   };

//   for (size_t i = 1; i < model.links_.size(); ++i) {
//     double m = data.compositeMassW[i];
//     Eigen::Vector3d c = data.compositeComW[i];
//     Eigen::Matrix3d Ic = data.compositeInertiaW[i]; // about COM

//     // spatial inertia about *world* origin:
//     // I₀ = Ic + m [   ] shift by c
//     Eigen::Matrix<double, 6, 6> I6;
//     I6.setZero();

//     // top‐left: rotational inertia about COM
//     I6.topLeftCorner<3, 3>() = Ic;
//     // upper‐right & lower‐left: coupling
//     I6.topRightCorner<3, 3>() = m * skew(c);
//     I6.bottomLeftCorner<3, 3>() = -m * skew(c);
//     // bottom‐right: mass matrix
//     I6.bottomRightCorner<3, 3>() = m * Eigen::Matrix3d::Identity();

//     data.compositeSpatialInertia6[i] = I6;
//   }
// }

// inline void crba(const Model &model, Data &data, const Eigen::VectorXd &gc)
// {
//   data.massMatrix = Eigen::MatrixXd::Zero(model.nv_, model.nv_);
//   data.massMatrix.setZero();

//   // update the inertia of each link to be in the world frame
//   inertiaUpdate(model, data, gc);
//   buildCompositeSpatialInertia(model, data);
//   computeMotionSubspaces(model, data);

//   auto skew = [](const Eigen::Vector3d &v) {
//     Eigen::Matrix3d skew;
//     skew.setZero();
//     skew(0, 1) = -v(2);
//     skew(0, 2) = v(1);
//     skew(1, 0) = v(2);
//     skew(1, 2) = -v(0);
//     skew(2, 0) = -v(1);
//     skew(2, 1) = v(0);
//     return skew;
//   };

//   // auto spatialTransform = [skew](Transform T) {
//   //   Eigen::MatrixXd spatialTransform(6, 6);
//   //   spatialTransform.setZero();

//   //   Eigen::Matrix3d R = T.block<3, 3>(0, 0);
//   //   Eigen::Vector3d t = T.block<3, 1>(0, 3);

//   //   spatialTransform.block<3, 3>(0, 0) = R;
//   //   spatialTransform.block<3, 3>(0, 3) = skew(t) * R;
//   //   spatialTransform.block<3, 3>(3, 0) = Eigen::Matrix3d::Zero();
//   //   spatialTransform.block<3, 3>(3, 3) = R;

//   //   return spatialTransform;
//   // };

//   auto motionSubspace = [skew](Eigen::Vector3d axis) {
//     Eigen::VectorXd motionSubspace(6);
//     motionSubspace.setZero();

//     motionSubspace.head<3>() = Eigen::Vector3d::Zero();
//     motionSubspace.tail<3>() = axis;

//     return motionSubspace;
//   };

//   auto spatialInertia = [skew](Eigen::Matrix3d inertia, double mass,
//                                Eigen::Vector3d com) {
//     Eigen::MatrixXd spatialInertia(6, 6);
//     spatialInertia.setZero();
//     spatialInertia.block<3, 3>(0, 0) = mass * Eigen::Matrix3d::Identity();
//     spatialInertia.block<3, 3>(3, 3) = inertia - mass * skew(com) *
//     skew(com); spatialInertia.block<3, 3>(0, 3) = -mass * skew(com);
//     spatialInertia.block<3, 3>(3, 0) = mass * skew(com).transpose();
//     return spatialInertia;
//   };

//   // manually compute one entry of mass matrix
//   int i = model.joints_.size() - 1;
//   int j = model.joints_.size() - 1;
//   auto linkIndex = model.joints_[i]->child->getIndex();
//   auto I6 = data.compositeSpatialInertia6[linkIndex];
//   auto Sj = data.S_W[j];
//   auto Si = data.S_W[i];

//   std::cout << "IA: " << std::endl;
//   std::cout << I6 << std::endl;

//   auto entry = Sj.dot(I6 * Si);
//   std::cout << "entry: " << entry << std::endl;

//   // // iterate from each joint upwards the kinematic chain
//   // for (int i = model.joints_.size() - 1; i >= 0; i--) {
//   //   auto joint = model.joints_[i];
//   //   auto columnI = model.idx_vs_[joint->getIndex()];

//   //   // now iterate updwards the kinematic chain from i-th joint
//   //   auto parentDof = joint->getParentDof();
//   //   auto jointJ = joint;
//   //   while (true) {
//   //     auto columnJ = model.idx_vs_[jointJ->getIndex()];
//   //     // std::cout << "joint name: " << jointJ->getName()
//   //     //           << " parentdof: " << jointJ->getParentDof()
//   //     //           << " column: " << columnIdx << std::endl;

//   //     std::cout << "M_" << columnI << "," << columnJ << " = " <<
//   parentDof
//   //               << std::endl;

//   //     if (parentDof == -1) {
//   //       // std::cout << "Found root: " << jointJ->getName() <<
//   std::endl;
//   //       // std::cout << "check" << std::endl;
//   //       break;
//   //     }

//   //     // iteration
//   //     jointJ = model.joints_[parentDof];
//   //     parentDof = jointJ->getParentDof();
//   //   }
//   // }
// }
}; // namespace algorithms

// #define _MAKE_STR(x) __MAKE_STR(x)
// #define __MAKE_STR(x) #x

// /// do not change the name of the method
// inline Eigen::Vector3d getEndEffectorPosition(const Eigen::VectorXd &gc) {
//   //////////////////////////
//   ///// Your Code Here /////
//   //////////////////////////

//   Model model(std::string(_MAKE_STR(RESOURCE_DIR)) + "/Panda/panda.urdf");
//   Data data(model);

//   algorithms::forwardKinematics(model, data, gc);

//   // get frame with name "panda_finger_joint3"
//   for (size_t i = 0; i < model.joints_.size(); i++) {
//     auto joint = model.joints_[i];
//     if (joint->getName() == "panda_finger_joint3") {
//       auto frame = data.oTj[i];
//       Eigen::Vector3d pos = frame.block<3, 1>(0, 3);
//       return pos;
//     }
//   }

//   // robot.framesForwardKinematics(data, gc);

//   // for (auto t : data.oTf) {
//   //   std::cout << t << std::endl;
//   //   std::cout << "----------------" << std::endl;
//   // }

//   // for (auto joint : robot.getJoints()) {
//   //   std::cout <<  << std::endl;
//   //   // joint->jointPlacement();
//   // }

//   return Eigen::Vector3d::Ones(); /// replace this
// }

// /// do not change the name of the method
// inline Eigen::Vector3d getLinearVelocity(const Eigen::VectorXd &gc,
//                                          const Eigen::VectorXd &gv) {
//   //////////////////////////
//   ///// Your Code Here /////
//   //////////////////////////

//   Model model(
//       "/home/lvjonok/github.com/lvjonok/ME553_2025/resource/Panda/panda.urdf");
//   Data data(model);

//   algorithms::framesForwardKinematics(model, data, gc);

//   // std::cout << "model body names: " << std::endl;
//   // for (auto link : model.links_) {
//   //   std::cout << link->getIndex() << ": " << link->getName() << std::endl;
//   // }

//   // for (size_t bodyId = 0; bodyId < 10; bodyId++) {
//   //   algorithms::frameRotJacobian(model, data, gc, bodyId);
//   //   algorithms::framePosJacobian(model, data, gc, bodyId);

//   //   std::cout << "frame: " << model.joints_[bodyId]->getName() <<
//   std::endl;

//   //   std::cout << "pos jacobian: " << std::endl;
//   //   std::cout << data.posJacobian << std::endl;
//   //   std::cout << "rot jacobian: " << std::endl;
//   //   std::cout << data.rotJacobian << std::endl;
//   //   std::cout << "----------------" << std::endl;
//   // }

//   for (size_t i = 0; i < model.joints_.size(); i++) {
//     auto joint = model.joints_[i];
//     if (joint->getName() == "panda_finger_joint3") {

//       // std::cout << "frame: " << joint->getName() << std::endl;
//       // std::cout << "frame idx: " << i << std::endl;
//       algorithms::framePosJacobian(model, data, gc, i);
//       // std::cout << "pos jacobian: " << std::endl;
//       // std::cout << data.posJacobian << std::endl;

//       Eigen::Vector3d result = data.posJacobian * gv;
//       // std::cout << "result: " << result.transpose() << std::endl;
//       // std::cout << "----------------" << std::endl;

//       return result;
//     }
//   }

//   return Eigen::Vector3d::Ones(); /// replace this
// }

// /// do not change the name of the method
// inline Eigen::Vector3d getAngularVelocity(const Eigen::VectorXd &gc,
//                                           const Eigen::VectorXd &gv) {
//   //////////////////////////
//   ///// Your Code Here /////
//   //////////////////////////

//   Model model(
//       "/home/lvjonok/github.com/lvjonok/ME553_2025/resource/Panda/panda.urdf");
//   Data data(model);

//   algorithms::framesForwardKinematics(model, data, gc);

//   // std::cout << "model body names: " << std::endl;
//   // for (auto link : model.links_) {
//   //   std::cout << link->getIndex() << ": " << link->getName() << std::endl;
//   // }

//   // for (size_t bodyId = 0; bodyId < 10; bodyId++) {
//   //   algorithms::frameRotJacobian(model, data, gc, bodyId);
//   //   algorithms::framePosJacobian(model, data, gc, bodyId);

//   //   std::cout << "frame: " << model.joints_[bodyId]->getName() <<
//   // std::endl;

//   //   std::cout << "pos jacobian: " << std::endl;
//   //   std::cout << data.posJacobian << std::endl;
//   //   std::cout << "rot jacobian: " << std::endl;
//   //   std::cout << data.rotJacobian << std::endl;
//   //   std::cout << "----------------" << std::endl;
//   // }

//   for (size_t i = 0; i < model.joints_.size(); i++) {
//     auto joint = model.joints_[i];
//     if (joint->getName() == "panda_finger_joint3") {

//       algorithms::checkStructure(model, data, i);

//       //   algorithms::frameRotJacobian(model, data, gc, i);

//       std::cout << "frame: " << joint->getName() << std::endl;
//       std::cout << "frame idx: " << i << std::endl;
//       std::cout << "own rot jacobian: " << std::endl;
//       std::cout << data.rotJacobian << std::endl;

//       std::cout << "own pos jacobian: " << std::endl;
//       std::cout << data.posJacobian << std::endl;

//       Eigen::Vector3d result = Eigen::Vector3d::Zero();
//       // std::cout << "result: " << result.transpose() << std::endl;
//       // std::cout << "----------------" << std::endl;

//       return result;
//     }
//   }

//   return Eigen::Vector3d::Ones(); /// replace this
// }

// /// do not change the name of the method
// inline Eigen::MatrixXd getMassMatrix(const Eigen::VectorXd &gc) {

//   /// !!!!!!!!!! NO RAISIM FUNCTIONS HERE !!!!!!!!!!!!!!!!!

//   Model model("/home/lvjonok/github.com/lvjonok/ME553_2025/resource/"
//               "mini_cheetah/urdf/cheetah.urdf");
//   Data data(model);

//   algorithms::framesForwardKinematics(model, data, gc);

//   // // find bodies CoM in the world frame
//   // for (size_t i = 0; i < model.links_.size(); i++) {
//   //   auto link = model.links_[i];
//   //   auto com = link->getTransform();
//   //   auto comW = data.oTb[i] * com;
//   //   std::cout << "body name: " << link->getName() << std::endl;
//   //   std::cout << "body com: " << comW << std::endl;
//   // }

//   // // find bodies Inertia expressed in the world frame
//   // for (size_t i = 0; i < model.links_.size(); i++) {
//   //   auto link = model.links_[i];
//   //   auto inertia = link->getInertia();

//   //   auto worldInertia = data.oTb[i].block<3, 3>(0, 0) * inertia *
//   //                       data.oTb[i].block<3, 3>(0, 0).transpose();

//   //   std::cout << "body name: " << link->getName() << std::endl;
//   //   std::cout << "body inertia: " << worldInertia << std::endl;
//   // }

//   // // iterate from leaves towards the root
//   // // first, find the leaves
//   // std::vector<std::shared_ptr<Joint>> leaves;
//   // std::cout << "leaves: " << std::endl;
//   // for (auto joint : model.leaves_) {
//   //   std::cout << "joint name: " << joint->getName() << std::endl;
//   // }

//   // // now starting from each leaf, go up to the root
//   // for (auto leaf : model.leaves_) {
//   //   std::cout << "Starting from leaf: " << leaf->getName() << std::endl;

//   //   // iterate upwards
//   //   auto joint = leaf;
//   //   auto parentDof = joint->getParentDof();
//   //   while (true) {
//   //     auto columnIdx = model.idx_vs_[joint->getIndex()];
//   //     std::cout << "joint name: " << joint->getName()
//   //               << " parentdof: " << joint->getParentDof()
//   //               << " column: " << columnIdx << std::endl;

//   //     if (parentDof == -1) {
//   //       std::cout << "Found root: " << joint->getName() << std::endl;
//   //       break;
//   //     }

//   //     // iteration
//   //     joint = model.joints_[parentDof];
//   //     parentDof = joint->getParentDof();
//   //   }
//   // }

//   // algorithms::crba(model, data, gc);

//   for (size_t i = 0; i < model.joints_.size(); i++) {
//     auto joint = model.joints_[i];
//     if (joint->getName() == "LF_JOINT3") {

//       algorithms::checkStructure(model, data, i);

//       //   algorithms::frameRotJacobian(model, data, gc, i);

//       std::cout << "frame: " << joint->getName() << std::endl;
//       std::cout << "frame idx: " << i << std::endl;
//       std::cout << "own rot jacobian: " << std::endl;
//       std::cout << data.rotJacobian << std::endl;

//       std::cout << "own pos jacobian: " << std::endl;
//       std::cout << data.posJacobian << std::endl;

//       Eigen::Vector3d result = Eigen::Vector3d::Zero();
//       // std::cout << "result: " << result.transpose() << std::endl;
//       // std::cout << "----------------" << std::endl;

//       return result;
//     }
//   }

//   std::cout << "mass matrix: " << std::endl;

//   // // print all the positions with the name of joints
//   // for (size_t i = 0; i < model.joints_.size(); i++) {
//   //   auto joint = model.joints_[i];
//   //   std::cout << "joint name: " << joint->getName() << std::endl;
//   //   std::cout << "joint position: " << data.oTj[i].block<3, 1>(0, 3)
//   //             << std::endl;
//   // }

//   return Eigen::MatrixXd::Ones(18, 18);
// }

#endif
