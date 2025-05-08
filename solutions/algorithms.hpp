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

inline void forwardPosition(const Model &model, Data &data,
                            const Eigen::VectorXd &gc) {
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
    data.jointAxis_W[i] = data.rot_WB[i] * joint->getAxis();

    if (parentId != -1) {
      data.joint2joint_W[i] = data.jointPos_W[i] - data.jointPos_W[parentId];
    }
  }
}

inline void forwardVelocity(const Model &model, Data &data,
                            const Eigen::VectorXd &gc,
                            const Eigen::VectorXd &gv) {
  data.bodyLinVel_w.resize(model.nbodies_);
  data.bodyAngVel_w.resize(model.nbodies_);
  data.motionSubspace.resize(model.nbodies_);
  data.dMotionSubspace.resize(model.nbodies_);

  // calculate the velocity of the body
  for (size_t i = 0; i < model.nbodies_; i++) {
    Eigen::Vector3d v, w;
    v.setZero();
    w.setZero();

    auto parentId = model.parents_[i];
    if (parentId != -1) {
      w = data.bodyAngVel_w[parentId];
      v = data.bodyLinVel_w[parentId] + w.cross(data.joint2joint_W[i]);
    }

    auto joint = model.actuated_joints_[i];

    Eigen::Vector3d axis = data.jointAxis_W[i];
    Eigen::Vector3d daxis = data.bodyAngVel_w[i].cross(axis);

    switch (joint->getType()) {
    case JointType::FIXED:
      // do nothing
      break;
    case JointType::FLOATING:
      // floating joint, we have to set the velocity
      v = gv.segment(0, 3);
      w = gv.segment(3, 3);
      // TODO: find out what kind of motion subspace we have
      // likely it is simply the identity matrix
      data.motionSubspace[i] = Eigen::MatrixXd::Identity(6, 6);
      data.dMotionSubspace[i] = Eigen::MatrixXd::Zero(6, 6);
      break;
    case JointType::REVOLUTE:
      // revolute joint, we have to set the angular velocity
      w += axis * gv[model.gv_idx_[i]];
      // find the motion subspace in the world frame
      data.motionSubspace[i] = Eigen::VectorXd::Zero(6);
      data.motionSubspace[i].block<3, 1>(3, 0) = axis;

      // find the derivative of the motion subspace
      data.dMotionSubspace[i] = Eigen::VectorXd::Zero(6);
      data.dMotionSubspace[i].block<3, 1>(0, 0) = daxis;
      break;
    case JointType::PRISMATIC:
      // prismatic joint, we have to set the linear velocity
      v += axis * gv[model.gv_idx_[i]];
      // find the motion subspace in the world frame
      data.motionSubspace[i] = Eigen::VectorXd::Zero(6);
      data.motionSubspace[i].block<3, 1>(0, 0) = axis;

      // find the derivative of the motion subspace
      data.dMotionSubspace[i] = Eigen::VectorXd::Zero(6);
      data.dMotionSubspace[i].block<3, 1>(3, 0) = daxis;
      break;
    }

    data.bodyLinVel_w[i] = v;
    data.bodyAngVel_w[i] = w;
    data.motionCross[i] = Eigen::MatrixXd::Zero(6, 6);
    /*
     * motion cross product matrix
     * [wx 0]
     * [vx wx]
     */
    data.motionCross[i].block<3, 3>(0, 0) = skew(w);
    data.motionCross[i].block<3, 3>(0, 3) = Eigen::Matrix3d::Zero();
    data.motionCross[i].block<3, 3>(3, 0) = skew(v);
    data.motionCross[i].block<3, 3>(3, 3) = skew(w);
  }
}

inline void forwardAcceleration(
    const Model &model, Data &data, const Eigen::VectorXd &gc,
    const Eigen::VectorXd &gv, const Eigen::VectorXd &ga,
    const Eigen::Vector3d &a0 = Eigen::Vector3d(0.0, 0.0, -9.81)) {
  data.bodyAngAcc.resize(model.nbodies_);
  data.bodyLinAcc.resize(model.nbodies_);

  // a0 says the zero linear velocity for the first body
  data.bodyLinAcc[0] = -a0;
  data.bodyAngAcc[0] = Eigen::Vector3d::Zero();

  for (size_t i = 1; i < model.nbodies_; i++) {
    auto parentId = model.parents_[i];
    assert(parentId != -1);

    Eigen::Vector3d a = data.bodyLinAcc[parentId];
    Eigen::Vector3d alpha = data.bodyAngAcc[parentId];

    Eigen::Vector3d pW = data.bodyAngVel_w[parentId];
    Eigen::Vector3d r = data.joint2joint_W[i];
    a += alpha.cross(r) + pW.cross(pW.cross(r));

    auto joint = model.actuated_joints_[i];
    auto S = data.motionSubspace[i];
    auto dS = data.dMotionSubspace[i];
    auto axis = data.jointAxis_W[i];

    auto gvi = gv[model.gv_idx_[i]];
    auto gai = ga[model.gv_idx_[i]];

    // // below is way 1
    // // simplified for 3D computations

    // // compute joint contribution
    // // at this point we should have either prismatic or revolute joint
    // if (joint->getType() == JointType::REVOLUTE) {
    //   alpha += axis * gai;
    //   alpha += pW.cross(axis * gvi);
    // } else if (joint->getType() == JointType::PRISMATIC) {
    //   a += axis * gai;
    //   a += pW.cross(axis * gvi);
    // } else {
    //   // do nothing
    //   std::cerr << "Should not be here, joint type is not prismatic or
    //   revolute"
    //             << std::endl;
    // }
    // data.bodyLinAcc[i] = a;
    // data.bodyAngAcc[i] = alpha;

    // below is way 2
    // using spatial relations
    Eigen::Matrix<double, 6, 1> aJ = S * gai;
    aJ += dS * gvi;
    aJ += data.motionCross[parentId] * (S * gvi);

    data.bodyLinAcc[i] = a + aJ.block<3, 1>(0, 0);
    data.bodyAngAcc[i] = alpha + aJ.block<3, 1>(3, 0);
  }
}

inline void compositeInertia(const Model &model, Data &data,
                             const Eigen::VectorXd &gc) {
  // first, update the com and inertia in the world frame
  data.inertiaW.resize(model.nbodies_);
  data.comW.resize(model.nbodies_);

  for (size_t i = 0; i < model.nbodies_; i++) {
    auto link = model.bodies_[i];
    auto inertia = link->getInertia();

    Eigen::Matrix3d R_BC = link->getTransform().block<3, 3>(0, 0);

    Eigen::Matrix3d worldInertia = data.rot_WB[i] * R_BC * inertia *
                                   R_BC.transpose() *
                                   data.rot_WB[i].transpose();

    // TODO: likely there is just a problem with all the loops
    // the iteration is weird and I end up with this stuff
    // just because first body is fixed in the fixed-based robots
    // raisim has zero inertia for them.
    if (i == 0 && model.actuated_joints_[i]->getType() != JointType::FLOATING) {
      worldInertia = Eigen::Matrix3d::Zero();
    }
    data.inertiaW[i] = worldInertia;

    // find the center of mass in the world frame
    Eigen::Vector3d com = link->getTransform().block<3, 1>(0, 3);
    Eigen::Vector3d comW = data.jointPos_W[i] + data.rot_WB[i] * R_BC * com;

    data.comW[i] = comW;
  }

  // now we compute the composite bodies
  data.compositeInertiaW.resize(model.nbodies_);
  data.compositeMassW.resize(model.nbodies_);
  data.compositeComW.resize(model.nbodies_);
  for (int i = model.nbodies_ - 1; i >= 0; i--) {
    auto link = model.bodies_[i];
    // properties of the current body
    data.compositeMassW[i] = link->getMass();
    data.compositeComW[i] = data.comW[i];
    data.compositeInertiaW[i] = data.inertiaW[i];

    // propagate to children
    for (auto childId : model.children_[i]) {
      // get mass, com and inertia of the child subtree
      double massChild = data.compositeMassW[childId];
      Eigen::Vector3d comChild = data.compositeComW[childId];
      Eigen::Matrix3d inertiaChild = data.compositeInertiaW[childId];

      // get mass, com and inertia of the current body
      double mass = data.compositeMassW[i];
      Eigen::Vector3d com = data.compositeComW[i];
      Eigen::Matrix3d inertia = data.compositeInertiaW[i];

      Eigen::Vector3d comNew =
          (mass * com + massChild * comChild) / (mass + massChild);
      Eigen::Vector3d r1 = com - comNew;
      Eigen::Vector3d r2 = comChild - comNew;

      Eigen::Matrix3d I_new = inertia + inertiaChild -
                              mass * skew(r1) * skew(r1) -
                              massChild * skew(r2) * skew(r2);

      // update composite data
      data.compositeMassW[i] = mass + massChild;
      data.compositeComW[i] = comNew;
      data.compositeInertiaW[i] = I_new;
    }
  }

  // for floating body system we shift total inertia from com to joint
  if (model.actuated_joints_[0]->getType() == JointType::FLOATING) {
    double Mtot = data.compositeMassW[0];
    Eigen::Vector3d Ctot = data.compositeComW[0];
    Eigen::Vector3d d0 = Ctot - data.jointPos_W[0]; // world-joint to COM
    data.compositeInertiaW[0] -= Mtot * skew(d0) * skew(d0);
  }

  // for a fixed base, first composite inertia is zero
  if (model.actuated_joints_[0]->getType() != JointType::FLOATING) {
    data.compositeInertiaW[0] = Eigen::Matrix3d::Zero();
    data.compositeComW[0] = Eigen::Vector3d::Zero();
  }

  data.spatialCompositeInertia6.resize(model.nbodies_);
  for (size_t i = 0; i < model.nbodies_; i++) {
    auto m = data.compositeMassW[i];
    Eigen::Matrix3d I = data.compositeInertiaW[i];
    Eigen::Vector3d r_ac = (data.compositeComW[i] - data.jointPos_W[i]);

    // fill the inertia matrix
    Eigen::MatrixXd sI = Eigen::Matrix<double, 6, 6>::Zero();
    sI.block<3, 3>(0, 0) = m * Eigen::Matrix3d::Identity();
    sI.block<3, 3>(0, 3) = -skew(r_ac) * m;
    sI.block<3, 3>(3, 0) = skew(r_ac) * m;
    if (i != 0) {
      sI.block<3, 3>(3, 3) = I - m * skew(r_ac) * skew(r_ac);
    } else {
      sI.block<3, 3>(3, 3) = I;
    }

    data.spatialCompositeInertia6[i] = sI;
  }
}

inline void crba(const Model &model, Data &data, const Eigen::VectorXd &gc) {
  data.massMatrix.resize(model.gv_size_, model.gv_size_);
  data.massMatrix.setZero();

  // for floating base system we iterate to the root (0),
  // for fixed only till (1)
  int root = 0;
  if (model.actuated_joints_[0]->getType() != JointType::FLOATING) {
    root = 1;
  }

  for (int j = model.nbodies_ - 1; j >= root; --j) {
    // get current subspce motion and composite spatial inertia
    auto Sj = data.motionSubspace[j];            // 6xjoint_dof
    auto Isp = data.spatialCompositeInertia6[j]; // 6x6

    // find a block to put the diagonal entry (it may be a matrix for floating
    // base)
    auto j_col = model.gv_idx_[j];
    auto j_len = model.actuated_joints_[j]->gv_length();

    // force from this joint
    Eigen::MatrixXd F = Isp * Sj; // 6xjoint_dof

    // diagonal entry: S_j^T * Isp_j * S_j
    data.massMatrix.block(j_col, j_col, j_len, j_len) =
        Sj.transpose() * F; // joint_dof x joint_dof

    // walk up the tree for off-diagonal entries
    int k = j;
    while (model.parents_[k] != root - 1) {
      int a = model.parents_[k];

      // child to parent vector
      Eigen::Vector3d r = data.joint2joint_W[k];

      // find plucker force motion matrix
      Eigen::MatrixXd aXb = Eigen::MatrixXd::Identity(6, 6);
      aXb.block<3, 3>(3, 0) = skew(r);
      F = aXb * F; // 6xjoint_dof

      // find a block to put the off-diagonal entry
      auto a_row = model.gv_idx_[a];
      auto a_len = model.actuated_joints_[a]->gv_length();

      // current joint subspace motion matrix
      auto Sa = data.motionSubspace[a];           // 6xjoint_dof
      Eigen::MatrixXd entry = Sa.transpose() * F; // joint_dof x joint_dof
      data.massMatrix.block(a_row, j_col, a_len, j_len) = entry;
      // advance up
      k = a;
    }
  }

  // fill the lower triangle
  data.massMatrix.triangularView<Eigen::Lower>() =
      data.massMatrix.transpose().triangularView<Eigen::Lower>();
}

inline void setState(const Model &model, Data &data, const Eigen::VectorXd &gc,
                     const Eigen::VectorXd &gv) {
  forwardPosition(model, data, gc);
  forwardVelocity(model, data, gc, gv);
  compositeInertia(model, data, gc);
}

inline void getBodyPose(const Model &model, Data &data, size_t bodyId,
                        Eigen::Matrix3d &R, Eigen::Vector3d &p) {
  // get the pose of the body in the world frame
  R = data.rot_WB[bodyId];
  p = data.jointPos_W[bodyId];
}

inline void getBodyTwist(const Model &model, Data &data, size_t bodyId,
                         Eigen::Vector3d &linVel, Eigen::Vector3d &angVel) {
  // get the velocity of the body in the world frame
  linVel = data.bodyLinVel_w[bodyId];
  angVel = data.bodyAngVel_w[bodyId];
}

} // namespace algorithms

#endif
