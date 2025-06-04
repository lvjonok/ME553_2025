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
  data.S.resize(model.nbodies_);
  data.dS.resize(model.nbodies_);

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
      data.S[i] = Eigen::MatrixXd::Identity(6, 6);
      data.dS[i] = Eigen::MatrixXd::Zero(6, 6);
      break;
    case JointType::REVOLUTE:
      // revolute joint, we have to set the angular velocity
      w += axis * gv[model.gv_idx_[i]];
      // find the motion subspace in the world frame
      data.S[i] = Eigen::VectorXd::Zero(6);
      data.S[i].block<3, 1>(3, 0) = axis;

      // find the derivative of the motion subspace
      data.dS[i] = Eigen::VectorXd::Zero(6);
      data.dS[i].block<3, 1>(3, 0) = daxis;
      break;
    case JointType::PRISMATIC:
      // prismatic joint, we have to set the linear velocity
      v += axis * gv[model.gv_idx_[i]];
      // find the motion subspace in the world frame
      data.S[i] = Eigen::VectorXd::Zero(6);
      data.S[i].block<3, 1>(0, 0) = axis;

      // find the derivative of the motion subspace
      data.dS[i] = Eigen::VectorXd::Zero(6);
      data.dS[i].block<3, 1>(0, 0) = daxis;
      break;
    }

    data.bodyLinVel_w[i] = v;
    data.bodyAngVel_w[i] = w;
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

  // data.bodyLinAcc[0] = -a0 + ga.segment(0, 3);
  // data.bodyAngAcc[0] = ga.segment(3, 3);

  for (size_t i = 1; i < model.nbodies_; i++) {
    auto parentId = model.parents_[i];
    assert(parentId != -1);

    Eigen::Vector3d a = data.bodyLinAcc[parentId];
    Eigen::Vector3d alpha = data.bodyAngAcc[parentId];

    Eigen::Vector3d pW = data.bodyAngVel_w[parentId];
    Eigen::Vector3d r = data.joint2joint_W[i];
    a += alpha.cross(r) + pW.cross(pW.cross(r));

    auto joint = model.actuated_joints_[i];
    auto S = data.S[i];
    auto dS = data.dS[i];
    auto axis = data.jointAxis_W[i];

    auto gvi = gv[model.gv_idx_[i]];
    auto gai = ga[model.gv_idx_[i]];

    // // below is way 1
    // // simplified for 3D computations

    // compute joint contribution
    // at this point we should have either prismatic or revolute joint
    if (joint->getType() == JointType::REVOLUTE) {
      alpha += axis * gai;
      alpha += pW.cross(axis * gvi);
    } else if (joint->getType() == JointType::PRISMATIC) {
      a += axis * gai;
      // TODO: find where does this 2 come from
      a += 2 * pW.cross(axis * gvi);
    } else {
      // do nothing
      std::cerr << "Should not be here, joint type is not prismatic or revolute"
                << std::endl;
    }
    data.bodyLinAcc[i] = a;
    data.bodyAngAcc[i] = alpha;

    // TODO: spatial version does not match right now
    // // below is way 2
    // // using spatial relations
    // Eigen::Matrix<double, 6, 1> aJ = S * gai;
    // aJ += dS * gvi;
    // aJ += data.motionCross[parentId] * (S * gvi);

    // data.bodyLinAcc[i] = a + aJ.block<3, 1>(0, 0);
    // data.bodyAngAcc[i] = alpha + aJ.block<3, 1>(3, 0);
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
    auto Sj = data.S[j];                         // 6xjoint_dof
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
      auto Sa = data.S[a];                        // 6xjoint_dof
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

inline void nonlinearities(const Model &model, Data &data,
                           const Eigen::VectorXd &gc, const Eigen::VectorXd &gv,
                           const Eigen::VectorXd &ga) {
  // we should have propagated the velocity through the kinematic chain already
  // TODO: maybe don't call at all
  algorithms::forwardVelocity(model, data, gc, gv);
  algorithms::forwardAcceleration(model, data, gc, gv, ga, {0, 0, -9.81});

  // second step
  // make a backward pass
  // find the force wrench acting on each body
  // collect the forces from the children
  // then propagate to the parent

  // we want to store the wrench in the world frame
  // acting on a body
  std::vector<Eigen::Vector3d> force(model.nbodies_);
  std::vector<Eigen::Vector3d> torque(model.nbodies_);

  // the vector to store the nonlinearities vector
  Eigen::VectorXd b = Eigen::VectorXd::Zero(gv.size());
  const int root =
      (model.actuated_joints_[0]->getType() != JointType::FLOATING) ? 1 : 0;
  for (int i = model.nbodies_ - 1; i >= root; --i) {

    // find the current wrench from the relation
    // f = I * a + v x I * v
    auto link = model.bodies_[i];
    auto r = data.comW[i] - data.jointPos_W[i];
    force[i] =
        link->getMass() * Eigen::Matrix3d::Identity() * data.bodyLinAcc[i] -
        link->getMass() * skew(r) * data.bodyAngAcc[i] +
        link->getMass() * skew(data.bodyAngVel_w[i]) *
            skew(data.bodyAngVel_w[i]) * r;

    torque[i] = link->getMass() * skew(r) * data.bodyLinAcc[i] +
                data.inertiaW[i] * data.bodyAngAcc[i] -
                link->getMass() * skew(r) * skew(r) * data.bodyAngAcc[i] +
                skew(data.bodyAngVel_w[i]) *
                    (data.inertiaW[i] - link->getMass() * skew(r) * skew(r)) *
                    data.bodyAngVel_w[i];

    // iterate over the children
    // and add up the forces they have experienced
    for (auto childId : model.children_[i]) {
      // get the force and torque of the child
      auto fChild = force[childId];
      auto tChild = torque[childId];

      // child to parent vector
      Eigen::Vector3d r = data.joint2joint_W[childId];

      force[i] += fChild;
      torque[i] += tChild + skew(r) * fChild;
    }

    // now component of the nonlinearities is S_j^T * f
    int start_idx = model.gv_idx_[i];
    int len = model.actuated_joints_[i]->gv_length();
    auto S = data.S[i];
    Eigen::VectorXd f = Eigen::VectorXd::Zero(6);
    f.segment(0, 3) = force[i];
    f.segment(3, 3) = torque[i];
    b.segment(start_idx, len) = S.transpose() * f; // joint_dof x 1
  }

  data.nonlinearities = b;
  // std::cout << "b: " << b.transpose() << std::endl;

  // after we get the wrench, we can multiply by the motion subspace
  // and get the generalized forces
  // for the case of a0 = -g. it will be fictitious forces or nonlinearities
}

inline void nonlinearities(const Model &model, Data &data,
                           const Eigen::VectorXd &gc,
                           const Eigen::VectorXd &gv) {
  return nonlinearities(model, data, gc, gv, gv * 0.0);
}

inline Eigen::VectorXd motionAct(const Eigen::VectorXd &v1,
                                 const Eigen::VectorXd &v2) {
  // this function computes the motion actuation
  Eigen::VectorXd res = Eigen::VectorXd::Zero(6);

  Eigen::Vector3d v1_lin = v1.segment(0, 3);
  Eigen::Vector3d v1_ang = v1.segment(3, 3);

  Eigen::Vector3d v2_lin = v2.segment(0, 3);
  Eigen::Vector3d v2_ang = v2.segment(3, 3);

  res.segment(0, 3) = v1_lin.cross(v2_ang) + v1_ang.cross(v2_lin);
  res.segment(3, 3) = v1_ang.cross(v2_ang);

  return res;
}

inline Eigen::VectorXd cross(Eigen::VectorXd v1, Eigen::VectorXd v2) {
  // this function computes the cross product of two vectors
  // in the spatial vector form
  Eigen::VectorXd res = Eigen::VectorXd::Zero(6);

  Eigen::Vector3d v1_lin = v1.segment(0, 3);
  Eigen::Vector3d v1_ang = v1.segment(3, 3);

  Eigen::Vector3d v2_lin = v2.segment(0, 3);
  Eigen::Vector3d v2_ang = v2.segment(3, 3);

  res.segment(0, 3) = v1_ang.cross(v2_lin);
  res.segment(3, 3) = v1_ang.cross(v2_ang) + v1_lin.cross(v2_lin);

  return res;
}

namespace aba {
inline void aba_pass1(const Model &model, Data &data, const Eigen::VectorXd &gc,
                      const Eigen::VectorXd &gv, const Eigen::VectorXd &gf) {
  // first, coming from root to the leaves we compute the kinematics and NE
  // terms Ma and ba
  algorithms::forwardVelocity(model, data, gc, gv);
  algorithms::forwardAcceleration(model, data, gc, gv, gv * 0.0, {0, 0, -9.81});
  algorithms::compositeInertia(model, data, gc);

  // we have already computed the forward pass for kinematics
  // now a few more fields specific to aba
  data.aXb.resize(model.nbodies_);
  data.XT.resize(model.nbodies_);
  data.Ma.resize(model.nbodies_);
  data.Pa.resize(model.nbodies_);

  for (size_t i = 0; i < model.nbodies_; i++) {
    auto joint = model.actuated_joints_[i];
    auto parentId = model.parents_[i];

    data.aXb[i] = algorithms::aXb(data.rot_WB[i], data.jointPos_W[i]);
    data.XT[i] =
        algorithms::aXb(Eigen::Matrix3d::Identity(), -data.joint2joint_W[i]);

    // articulated inertia and bias force
    data.Ma[i] = Eigen::MatrixXd::Zero(6, 6);
    {
      auto m = model.bodies_[i]->getMass();
      auto I = data.inertiaW[i];
      auto r = data.comW[i] - data.jointPos_W[i];
      data.Ma[i].block<3, 3>(0, 0) = m * Eigen::Matrix3d::Identity();
      data.Ma[i].block<3, 3>(0, 3) = -skew(r) * m;
      data.Ma[i].block<3, 3>(3, 0) = skew(r) * m;
      data.Ma[i].block<3, 3>(3, 3) = I - m * skew(r) * skew(r);
    }
    data.Pa[i] = Eigen::VectorXd::Zero(6);
    {
      // find the current wrench from the relation
      // f = v x I * v + f_ext
      // TODO: decide whether to use the external forces or not
      // NOTE: we don't include the gravity here, it should be done later
      auto link = model.bodies_[i];
      auto r = data.comW[i] - data.jointPos_W[i];
      data.Pa[i].segment(0, 3) = link->getMass() * skew(data.bodyAngVel_w[i]) *
                                 skew(data.bodyAngVel_w[i]) * r;

      data.Pa[i].segment(3, 3) =
          skew(data.bodyAngVel_w[i]) *
          (data.inertiaW[i] - link->getMass() * skew(r) * skew(r)) *
          data.bodyAngVel_w[i];
    }
  }
}
inline void aba_pass2(const Model &model, Data &data, const Eigen::VectorXd &gc,
                      const Eigen::VectorXd &gv, const Eigen::VectorXd &gf) {
  // from leaves to the root, we should merge the articulated inertia
  data.STMaXT.resize(model.nbodies_);
  data.STMa.resize(model.nbodies_);
  data.STMaSinvSTMaXT.resize(model.nbodies_);
  data.SdotUpXdotTV.resize(model.nbodies_);
  data.Xbp.resize(model.nbodies_);
  data.dXbpT.resize(model.nbodies_);
  data.W.resize(model.nbodies_);
  data.STMaSinv.resize(model.nbodies_);

  for (int i = model.nbodies_ - 1; i >= 1; --i) {
    auto joint = model.actuated_joints_[i];
    auto parentId = model.parents_[i];

    // simply access
    auto &Ma = data.Ma;
    auto &Pa = data.Pa;
    Eigen::MatrixXd S = data.S[i];                // 6xjoint_dof
    Eigen::MatrixXd Xbp = data.XT[i].transpose(); // 6x6
    Eigen::MatrixXd XbpT = Xbp.transpose();
    Eigen::MatrixXd ST = S.transpose(); // joint_dof x 6
    Eigen::MatrixXd dS = data.dS[i];    // 6xjoint_dof

    Eigen::MatrixXd dXbpT = Eigen::MatrixXd::Zero(6, 6);
    {
      // dXbp is the derivative of Xbp with respect to the joint position
      // it is computed as a cross product of the body angular velocity and
      // the joint2joint vector
      dXbpT.block<3, 3>(0, 3) =
          skew(data.bodyAngVel_w[parentId].cross(-data.joint2joint_W[i]));
    }

    data.Xbp[i] = Xbp;
    data.dXbpT[i] = dXbpT;

    // Eigen::MatrixXd dXbpT = data.XT[i].transpose(); // 6x6
    Eigen::MatrixXd gfi = gf.segment(model.gv_idx_[i], joint->gv_length());
    Eigen::MatrixXd gvi = gv.segment(model.gv_idx_[i], joint->gv_length());
    Eigen::VectorXd W = Eigen::VectorXd::Zero(6);
    {
      W.segment(0, 3) = data.bodyLinVel_w[parentId];
      W.segment(3, 3) = data.bodyAngVel_w[parentId];
    }
    data.W[i] = W;

    // simplifications
    data.STMa[i] = ST * Ma[i];
    data.STMaXT[i] = ST * Ma[i] * XbpT;
    // NOTE: we can do it only because we take care of the floating joint
    // separately
    data.STMaSinv[i] = (data.STMa[i] * S).inverse()(0);
    data.STMaSinvSTMaXT[i] = data.STMaSinv[i] * data.STMaXT[i];
    data.SdotUpXdotTV[i] = dS * gvi + dXbpT * W;

    // find the articulated inertia
    Ma[parentId] += Xbp * Ma[i] * (-S * data.STMaSinvSTMaXT[i] + XbpT);
    Pa[parentId] +=
        Xbp *
        (Ma[i] * (S * data.STMaSinv[i] *
                      (gfi - data.STMa[i] * data.SdotUpXdotTV[i] - ST * Pa[i]) +
                  data.SdotUpXdotTV[i]) +
         Pa[i]);
  }
}
inline void aba_pass3(const Model &model, Data &data, const Eigen::VectorXd &gc,
                      const Eigen::VectorXd &gv, const Eigen::VectorXd &gf) {
  // now forward pass where we will update the acceleration
  data.udot = Eigen::VectorXd::Zero(model.gv_size_);
  data.Wdot.resize(model.nbodies_);

  if (model.actuated_joints_[0]->getType() == JointType::FLOATING) {
    // for floating base we have the first 6 entries
    // of the generalized velocity vector
    data.Wdot[0] =
        data.Ma[0].inverse() * (gf.segment(0, 6) - data.Pa[0]); // base is zero
    data.udot.segment(0, 6) = data.Wdot[0];

    // add gravitation
    data.udot(2) -= 9.81; // add gravity to the vertical component
  } else {
    // for fixed base we have the first 6 entries
    // of the generalized velocity vector as zero
    data.Wdot[0] = Eigen::VectorXd::Zero(6);
  }

  for (int i = 1; i < model.nbodies_; ++i) {
    auto joint = model.actuated_joints_[i];
    auto parentId = model.parents_[i];

    // get the generalized velocity vector
    auto v_start = model.gv_idx_[i];
    auto v_len = joint->gv_length();

    // get the generalized force vector
    auto gfi = gf.segment(v_start, v_len);
    auto gvi = gv.segment(v_start, v_len);

    data.udot.segment(v_start, v_len) =
        data.STMaSinv[i] *
        (gfi -
         data.STMa[i] * (data.SdotUpXdotTV[i] +
                         data.Xbp[i].transpose() * data.Wdot[parentId]) -
         data.S[i].transpose() * data.Pa[i]);

    data.Wdot[i] = data.SdotUpXdotTV[i] +
                   data.S[i] * data.udot.segment(v_start, v_len) +
                   data.Xbp[i].transpose() * data.Wdot[parentId];
  }
}
} // namespace aba

inline void articulatedBodyAlgorithm(const Model &model, Data &data,
                                     const Eigen::VectorXd &gc,
                                     const Eigen::VectorXd &gv,
                                     const Eigen::VectorXd &gf) {
  aba::aba_pass1(model, data, gc, gv, gf);
  aba::aba_pass2(model, data, gc, gv, gf);
  aba::aba_pass3(model, data, gc, gv, gf);
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
