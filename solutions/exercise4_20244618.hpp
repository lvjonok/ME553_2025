#pragma once


#include <Eigen/Core>
#include <Eigen/Geometry>
#include <cstddef>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <tinyxml_rai/tinyxml_rai.h>
#include <vector>

using Transform = Eigen::Matrix4d;
using SpatialTransform = Eigen::Matrix<double, 6, 6>;
using Motion = Eigen::Matrix<double, 6, 1>;

inline Eigen::Matrix3d skew(const Eigen::Vector3d &v) {
  Eigen::Matrix3d skew;

  // clang-format off
    skew << 0, -v(2), v(1), 
            v(2), 0, -v(0), 
            -v(1), v(0), 0;
  // clang-format on

  return skew;
}

inline Eigen::Vector3d unskew(const Eigen::Matrix3d &v) {
  Eigen::Vector3d skew;

  skew(0) = v(2, 1);
  skew(1) = v(0, 2);
  skew(2) = v(1, 0);

  return skew;
}

// TODO: check and make in line with Featherstone notation
inline Eigen::Matrix3d RotX(double theta) {
  Eigen::Matrix3d R;
  R << 1, 0, 0, 0, cos(theta), -sin(theta), 0, sin(theta), cos(theta);
  return R;
}

inline Eigen::Matrix3d RotY(double theta) {
  Eigen::Matrix3d R;
  R << cos(theta), 0, sin(theta), 0, 1, 0, -sin(theta), 0, cos(theta);
  return R;
}

inline Eigen::Matrix3d RotZ(double theta) {
  Eigen::Matrix3d R;
  R << cos(theta), -sin(theta), 0, sin(theta), cos(theta), 0, 0, 0, 1;
  return R;
}

/// Returns a 3×3 rotation matrix from roll, pitch, yaw (ZYX convention).
/// roll  = rotation about X
/// pitch = rotation about Y
/// yaw   = rotation about Z
inline Eigen::Matrix3d RotRPY(double roll, double pitch, double yaw) {
  double cr = std::cos(roll);
  double sr = std::sin(roll);
  double cp = std::cos(pitch);
  double sp = std::sin(pitch);
  double cy = std::cos(yaw);
  double sy = std::sin(yaw);

  Eigen::Matrix3d R;
  R(0, 0) = cy * cp;
  R(0, 1) = cy * sp * sr - sy * cr;
  R(0, 2) = cy * sp * cr + sy * sr;

  R(1, 0) = sy * cp;
  R(1, 1) = sy * sp * sr + cy * cr;
  R(1, 2) = sy * sp * cr - cy * sr;

  R(2, 0) = -sp;
  R(2, 1) = cp * sr;
  R(2, 2) = cp * cr;

  return R;
}

inline Eigen::Matrix3d RotAxisAngle(const Eigen::Vector3d &axis, double angle) {
  Eigen::Vector3d k = axis.normalized();
  double c = std::cos(angle);
  double s = std::sin(angle);
  // skew‐symmetric cross‐product matrix of k
  Eigen::Matrix3d K;
  K << 0, -k.z(), k.y(), k.z(), 0, -k.x(), -k.y(), k.x(), 0;
  // Rodrigues' rotation formula
  return c * Eigen::Matrix3d::Identity() + (1 - c) * (k * k.transpose()) +
         s * K;
}

inline SpatialTransform bXa(const Eigen::Matrix3d &R,
                            const Eigen::Vector3d &p) {
  SpatialTransform X;
  X.setZero();
  X.block<3, 3>(0, 0) = R;
  X.block<3, 3>(3, 3) = R;
  X.block<3, 3>(3, 0) = -R * skew(p);
  return X;
}

inline SpatialTransform aXb(const Eigen::Matrix3d &R,
                            const Eigen::Vector3d &p) {
  SpatialTransform X;
  X.setZero();
  X.block<3, 3>(0, 0) = R.transpose();
  X.block<3, 3>(3, 3) = R.transpose();
  X.block<3, 3>(0, 3) = skew(p) * R.transpose();
  return X;
}

class Link {
public:
  Link(std::string name) : name(name) {}

  Link(std::string name, Eigen::Vector3d origin, Eigen::Vector3d rpy,
       double mass, Eigen::VectorXd inertiaFlat)
      : name(name), _origin(origin), _rpy(rpy), _mass(mass) {
    // inertiaFlat contains 6 values: ixx, ixy, ixz, iyy, iyz, izz
    if (inertiaFlat.size() == 6) {
      // clang-format off
      _inertia << inertiaFlat(0), inertiaFlat(1), inertiaFlat(2),
                  inertiaFlat(1), inertiaFlat(3), inertiaFlat(4),
                  inertiaFlat(2), inertiaFlat(4), inertiaFlat(5);
      // clang-format on
    } else {
      std::cerr << "Error: inertiaFlat must have 6 elements." << std::endl;
    }
  }

  std::string getName() const { return name; }

  void setIndex(int idx) { index = idx; }

  int getIndex() const { return index; }

  double getMass() const { return _mass; }
  Eigen::Matrix3d getInertia() const { return _inertia; }

  Transform getTransform() const {
    Transform T = Transform::Identity();
    T.block<3, 3>(0, 0) = RotZ(_rpy[2]) * RotY(_rpy[1]) * RotX(_rpy[0]);
    T.block<3, 1>(0, 3) = _origin;
    return T;
  }

private:
  std::string name;

  int index;

  Eigen::Vector3d _origin = Eigen::Vector3d::Zero();
  Eigen::Vector3d _rpy = Eigen::Vector3d::Zero();

  double _mass = 0.0;
  Eigen::Matrix3d _inertia = Eigen::Matrix3d::Zero();
};

// Joint types
enum class JointType { REVOLUTE, PRISMATIC, FIXED, FLOATING };

class Joint {
public:
  Joint(std::shared_ptr<Link> parent, std::shared_ptr<Link> child,
        Eigen::Vector3d axis, Eigen::Vector3d xyz, Eigen::Vector3d rpy,
        JointType type, std::string jointName)
      : parent(parent), child(child), axis(axis), type_(type),
        joint_name(jointName) {
    if (type_ == JointType::FIXED) {
      // TODO: this was done specifically to match the doubleCartPole model
      this->axis.setZero();
      this->axis(2) = 1;
    }

    this->origin = Transform::Identity();
    this->origin.block<3, 3>(0, 0) = RotZ(rpy[2]) * RotY(rpy[1]) * RotX(rpy[0]);
    this->origin.block<3, 1>(0, 3) = xyz;
  }

  Joint(std::shared_ptr<Link> parent, std::shared_ptr<Link> child,
        Eigen::Vector3d axis, Transform origin, JointType type,
        std::string jointName)
      : parent(parent), child(child), axis(axis), type_(type), origin(origin),
        joint_name(jointName) {}

  // compute the transform from parent to child
  // this is const and defined by the URDF only
  Transform jointPlacement() const { return origin; }

  SpatialTransform spatialJointPlacement() {
    Transform T = jointPlacement();

    auto R = T.block<3, 3>(0, 0);
    auto p = T.block<3, 1>(0, 3);

    // the order of translation and rotation is different to the ordinary
    // homogeneous transformations
    return bXa(R, R.transpose() * p);
  }

  // mergeJoint should be called when we have a structure like
  // link -> joint -> link -> fixed joint -> link -> joint
  // this way we want to merge first and second joint together
  // because the fixed joint does not contribute to the motion
  // then, first joint along can describe the resulting transformation
  Joint *mergeJoint(const Joint &other) {
    auto T1 = jointPlacement();
    auto T2 = other.jointPlacement();

    // we merge <fixed joints> - link - <joint>
    // that means we have to inherit the type of other joint
    // and the axis of the other joint
    return new Joint(parent, other.child, other.axis, T1 * T2, other.type_,
                     joint_name + "_" + other.joint_name);
  }

  std::tuple<SpatialTransform, Motion> jcalc(const Eigen::VectorXd &gc,
                                             int start_idx,
                                             const Eigen::VectorXd &gv,
                                             int gv_start_idx) {
    if (type_ == JointType::FIXED) {
      throw std::runtime_error("Trying to compute "
                               "motion for a fixed "
                               "joint");
    }

    if (type_ == JointType::FLOATING) {
      Eigen::Vector3d xyz = gc.head<3>();

      // by convention the generalized coordinates are
      // x, y, z, qw, qx, qy, qz
      Eigen::Vector4d quat_raw = gc.segment<4>(3);
      Eigen::Quaterniond quat(quat_raw(0), quat_raw(1), quat_raw(2),
                              quat_raw(3));
      Eigen::Matrix3d R = quat.toRotationMatrix();

      SpatialTransform X = bXa(R, R.transpose() * xyz);
      Eigen::MatrixXd S(6, 6);

      return std::make_tuple(X, S * gv.segment<6>(gv_start_idx));
    }

    if (type_ == JointType::PRISMATIC) {
      // For prismatic joint we have one degree of freedom and our
      // transformation is given by the motion along the axis.
      SpatialTransform X =
          bXa(Eigen::Matrix3d::Identity(), axis * gc[start_idx]);

      Eigen::VectorXd S(6);
      S.setZero();
      S.tail<3>() = axis; // v

      return std::make_tuple(X, S * gv[gv_start_idx]);
    }

    if (type_ == JointType::REVOLUTE) {
      // For revolute joint we have one degree of freedom and our
      // transformation is given by the rotation around the axis.
      SpatialTransform X =
          bXa(RotAxisAngle(axis, gc[start_idx]), Eigen::Vector3d::Zero());
      Eigen::VectorXd S(6);
      S.setZero();
      S.head<3>() = axis; // omega

      return std::make_tuple(X, S * gv[gv_start_idx]);
    }

    throw std::runtime_error(
        "Trying to compute motion for an unknown joint: " + joint_name +
        " with index: " + std::to_string(index) + ".");
  }

  // SpatialTransform jointPlacement() {
  //   // Eigen::Matrix3d R = RotZ(rpy[2]) * RotY(rpy[1]) * RotX(rpy[0]);
  //   Eigen::Matrix3d R = RotRPY(rpy[0], rpy[1], rpy[2]);

  //   // Transform T = Transform::Identity();
  //   // T.block<3, 3>(0, 0) = R;
  //   // T.block<3, 1>(0, 3) = xyz;

  //   // Transform T_inv = T.inverse();
  //   // R = T_inv.block<3, 3>(0, 0);
  //   // xyz = T_inv.block<3, 1>(0, 3);

  //   return bXa(R.transpose(), xyz);

  //   // Eigen::Matrix3d R = Eigen::expm1(xyz);

  //   // SpatialTransform XR = SpatialTransform::Identity();
  //   // XR.block<3, 3>(0, 0) = R;
  //   // XR.block<3, 3>(3, 3) = R;
  //   // XR.block<3, 3>(3, 0) = -R * skew(xyz);

  //   // return XR;

  //   // SpatialTransform XP = SpatialTransform::Identity();
  //   // XR.block<3, 3>(3, 0) = skew(-xyz);

  //   // return (XR * XP).inverse();
  // }

  // TODO: add motion subspace matrix here and everything else
  Transform jcalcOld(const Eigen::VectorXd &gc, int start_idx) {
    std::cout << "entering jcalc" << std::endl;
    return motion(gc, start_idx, gc_length());
  }

  // compute the transform due to motion
  Transform motion(double theta) {
    Eigen::Matrix3d R;
    if (type_ == JointType::REVOLUTE) {
      // TODO: what if axis are not z?
      // R = RotZ(theta);
      R = RotAxisAngle(axis, theta);
    } else if (type_ == JointType::PRISMATIC) {
      R = Eigen::Matrix3d::Identity();
    } else if (type_ == JointType::FLOATING) {
      std::cout << "Should have passed a vector of 7 values" << std::endl;
      R = Eigen::Matrix3d::Identity();
    } else {
      R = Eigen::Matrix3d::Identity();
    }

    Eigen::Vector3d t;
    if (type_ == JointType::PRISMATIC) {
      t = axis * theta;
    } else {
      t = Eigen::Vector3d::Zero();
    }

    Transform T = Transform::Identity();
    T.block<3, 3>(0, 0) = R;
    T.block<3, 1>(0, 3) = t;
    return T;
  }

  // Motion transformation is defined as a function of generalized
  // coordinates
  // we pass the start index and length of the generalized coordinates
  Transform motion(const Eigen::VectorXd &gc, size_t start_idx, size_t length) {
    if (start_idx == -1) {
      return Transform::Identity();
    }

    if (length == 1)
      return motion(gc[start_idx]);

    return motion(gc.segment(start_idx, length));
  }

  // compute the transform due to motion (for floating body)
  Transform motion(const Eigen::VectorXd &theta) {
    auto xyz = theta.head<3>();
    auto quat = theta.segment<4>(3);

    // convert manually to rotation matrix
    Eigen::Matrix3d R;
    R(0, 0) = 1 - 2 * (quat(2) * quat(2) + quat(3) * quat(3));
    R(0, 1) = 2 * (quat(1) * quat(2) - quat(0) * quat(3));
    R(0, 2) = 2 * (quat(1) * quat(3) + quat(0) * quat(2));
    R(1, 0) = 2 * (quat(1) * quat(2) + quat(0) * quat(3));
    R(1, 1) = 1 - 2 * (quat(1) * quat(1) + quat(3) * quat(3));
    R(1, 2) = 2 * (quat(2) * quat(3) - quat(0) * quat(1));
    R(2, 0) = 2 * (quat(1) * quat(3) - quat(0) * quat(2));
    R(2, 1) = 2 * (quat(2) * quat(3) + quat(0) * quat(1));
    R(2, 2) = 1 - 2 * (quat(1) * quat(1) + quat(2) * quat(2));
    // std::cout << "R: " << R << std::endl;
    // std::cout << "quat: " << quat.transpose() << std::endl;
    // std::cout << "xyz: " << xyz.transpose() << std::endl;
    Transform T = Transform::Identity();
    T.block<3, 3>(0, 0) = R;
    T.block<3, 1>(0, 3) = xyz;

    return T;
  }

  int gc_length() {
    // return 1;
    if (type_ == JointType::FIXED) {
      return 0;
    } else if (type_ == JointType::FLOATING) {
      return 7;
    } else if (type_ == JointType::REVOLUTE) {
      return 1;
    }

    return 1;
  }

  int gc_start_idx() { return this->gc_start_idx_; }

  int gv_length() {
    // return 1;
    if (type_ == JointType::FIXED) {
      return 0;
    } else if (type_ == JointType::FLOATING) {
      return 6;
    } else if (type_ == JointType::REVOLUTE) {
      return 1;
    }

    return 1;
  }

  int gv_start_idx() { return gv_start_idx_; }

  std::shared_ptr<Link> parent;
  std::shared_ptr<Link> child;

  std::string getName() const { return joint_name; }

  JointType getType() const { return type_; }

  void setIndex(int idx) { index = idx; }

  int getIndex() const { return index; }

  Eigen::Vector3d getAxis() const { return axis; }

  void setParentDof(int parentDof) { parentDof_ = parentDof; }
  int getParentDof() const { return parentDof_; }

private:
  // axis definition, optional, default is [1, 0, 0]
  Eigen::Vector3d axis;

  // // origin definition
  // Eigen::Vector3d xyz; // default is [0, 0, 0]
  // // XYZ Euler angle definition
  // Eigen::Vector3d rpy; // default is [0, 0, 0]

  // origin definition through transform
  Transform origin;

  // in order to apply generalized coordinates and velocities
  // later, we have to keep track of the joint type and starting index
  // in the generalized coordinate and velocity vectors
  // this is due to the fact that the URDF may contain fixed joints
  // and we have to skip them when applying generalized coordinates
  // and velocities
  JointType type_;
  int gc_start_idx_;
  int gv_start_idx_;

  std::string joint_name;
  int index;

  int parentDof_ = -2; // index of the joint that is before this joint
};

class Model {
public:
  Model(std::string urdf_path) {
    raisim::TiXmlDocument urdf;
    urdf.LoadFile(urdf_path);

    parseURDF(urdf);
    // at this points we only filled the links and joints

    parents_.resize(0);

    // now with the parent array
    // first, find the link that is not a child of any joint
    // this is the root link
    std::shared_ptr<Link> rootLink = nullptr;
    {
      for (auto link : links_) {
        bool isChild = false;
        for (auto joint : joints_) {
          if (joint->child == link) {
            isChild = true;
            break;
          }
        }
        if (!isChild) {
          rootLink = link;
          break;
        }
      }
    }

    if (rootLink == nullptr) {
      std::cerr << "Error: no root link found." << std::endl;
      return;
    }

    // std::cout << "root link: " << rootLink->getName() << std::endl;

    // if rootLink is "world", we have a fixed-base robot
    // otherwise, we have a floating base robot and need to manually add
    // the floating joint

    bool isFixedBase = true;

    if (rootLink->getName() != "world") {
      // std::cout << "We have a floating base robot" << std::endl;
      // std::cout << "Adding floating joint and world under the hood"
      //           << std::endl;

      // make a dummy world link
      auto worldLink = std::make_shared<Link>("world");
      // std::cout << "worldLink: " << worldLink->getName() << std::endl;
      auto floatingJoint = std::make_shared<Joint>(
          worldLink, rootLink, Eigen::Vector3d(0, 0, 0),
          Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(0, 0, 0),
          JointType::FLOATING, "ROOT");
      // put worldLink in front of the list
      links_.insert(links_.begin(), worldLink);
      joints_.insert(joints_.begin(), floatingJoint);

      // rootLink = worldLink;
      isFixedBase = false;
    }

    // now we can traverse the children and fill the parent array
    // addLink(rootLink, -1, false, isFixedBase);
    nbodies_ = 0;
    parents_.clear();
    traverse(links_[0], -1);

    // std::cout << "printing bodies and actuated joints in order" << std::endl;
    nbodies_ = bodies_.size();
    for (size_t i = 0; i < nbodies_; i++) {
      auto link = bodies_[i];
      link->setIndex(i);
      // std::cout << "link name: " << link->getName() << std::endl;
    }
    // std::cout << "printing actuated joints in order" << std::endl;
    // for (size_t i = 0; i < actuated_joints_.size(); i++) {
    //   auto joint = actuated_joints_[i];
    //   std::cout << "joint name: " << joint->getName() << std::endl;
    // }
    children_.clear();
    for (size_t i = 0; i < nbodies_; i++) {
      auto bodyi = bodies_[i];
      std::vector<size_t> bodyi_children;

      for (auto joint : actuated_joints_) {
        if (joint->parent != bodyi) {
          continue;
        }

        bodyi_children.push_back(joint->child->getIndex());
      }
      children_.push_back(bodyi_children);
    }

    // now as we have parsed everything, we can set the tree transformation
    // by URDF convention, bodies do not have a transformation (they are just
    // attached to the corresponding joint)
    T_T_.resize(nbodies_);
    X_T_.resize(nbodies_);
    for (size_t i = 0; i < nbodies_; i++) {
      auto joint = actuated_joints_[i];
      T_T_[i] = joint->jointPlacement();
      X_T_[i] = joint->spatialJointPlacement();
    }
  }

  void traverse(std::shared_ptr<Link> link, int currentBodyId,
                Joint *fixed_joints_ = nullptr) {
    // TODO: we should introduce a way to merge bodies that are connected
    // through a fixed joint this is important for dynamics computations, but
    // may be overkill for now

    // this should be reset on every actuated joint, but accumulated for a chain
    // of fixed joints
    // then, this transform will be used to premultiply the transformation for
    // the next actuated joint
    // this way we can merge the fixed joints together
    // Joint *fixed_joints = nullptr;

    for (const auto &joint : joints_) {
      // skip joints that are not related to the current link
      if (joint->parent != link)
        continue;

      // this joint is attached to this link
      if (joint->getType() == JointType::FIXED && currentBodyId == -1) {
        int newBodyId = nbodies_++;
        parents_.push_back(currentBodyId);
        bodies_.push_back(joint->child); // world -> fixed_joint -> child

        gc_idx_.push_back(-1);
        gv_idx_.push_back(-1);

        // TODO: this is a little weird, but even RAILAB does this
        // we consider first fixed joint for fixed-base robots
        actuated_joints_.push_back(joint);

        traverse(joint->child, newBodyId);
        return;
      }

      if (joint->getType() == JointType::FIXED) {
        // this is a fixed joint, it should not be added to the bodies

        Joint *new_fixed_joint = nullptr;

        if (fixed_joints_ == nullptr) {
          new_fixed_joint = joint.get();
        } else {
          // we have to merge the fixed joints together
          new_fixed_joint =
              fixed_joints_->mergeJoint(*joint); // merge the fixed joints
          std::cout << "merging joints, got " << new_fixed_joint->getName()
                    << std::endl;
        }
        traverse(joint->child, currentBodyId, new_fixed_joint);

        // auto last_joint = actuated_joints_.back();
        // auto merged_joint = joint->mergeJoint(*last_joint);
        // // auto merged_joint = last_joint->mergeJoint(*joint);
        // actuated_joints_.pop_back();
        // actuated_joints_.push_back(std::make_shared<Joint>(merged_joint));

      } else {
        // this is a new Body
        int newBodyId = nbodies_++;
        parents_.push_back(currentBodyId);
        bodies_.push_back(joint->child); // link1 -> fixed_joint -> link2

        std::shared_ptr<Joint> newJoint = nullptr;

        if (fixed_joints_ != nullptr) {
          // we have a chain of fixed joints behind,
          // now we have to merge them before the actuated joint
          newJoint = std::make_shared<Joint>(
              *fixed_joints_->mergeJoint(*joint)); // merge the fixed joints
          std::cout << "actuated joint with fixed joints, got "
                    << newJoint->getName() << std::endl;
        } else {
          newJoint = joint;
        }
        // fixed_joints_ = nullptr;

        actuated_joints_.push_back(newJoint);

        gc_idx_.push_back(gc_size_);
        gv_idx_.push_back(gv_size_);
        gc_size_ += joint->gc_length();
        gv_size_ += joint->gv_length();

        traverse(joint->child, newBodyId);
      }
    }
  }

  void parseURDF(const raisim::TiXmlDocument &urdf) {
    auto root = urdf.RootElement();
    if (root == nullptr) {
      std::cerr << "URDF is empty." << std::endl;
      return;
    }

    for (auto link = root->FirstChildElement("link"); link;
         link = link->NextSiblingElement("link")) {
      auto linkObj = parseLink(*link);
      links_.push_back(linkObj);
      linkMap_[linkObj->getName()] = linkObj;
    }

    for (auto joint = root->FirstChildElement("joint"); joint;
         joint = joint->NextSiblingElement("joint")) {
      // std::cout << "joint: " << joint->Attribute("name") << std::endl;

      // parse and create a new joint

      // Required fields
      auto parentName = joint->FirstChildElement("parent")->Attribute("link");
      auto childName = joint->FirstChildElement("child")->Attribute("link");

      // check that the parent and child links exist
      if (linkMap_.find(parentName) == linkMap_.end()) {
        std::cerr << "Parent link " << parentName << " does not exist."
                  << std::endl;
        return;
      }
      if (linkMap_.find(childName) == linkMap_.end()) {
        std::cerr << "Child link " << childName << " does not exist."
                  << std::endl;
        return;
      }

      auto parent = linkMap_[parentName];
      auto child = linkMap_[childName];

      auto jointObj = parseJoint(*joint, parent, child);
      joints_.push_back(jointObj);
    }
  }

  std::shared_ptr<Link> parseLink(const raisim::TiXmlElement &link) {
    auto linkName = link.Attribute("name");
    // return std::make_shared<Link>(linkName);
    auto inertial = link.FirstChildElement("inertial");

    if (inertial == nullptr) {
      std::cout << "Link " << linkName << " does not have an inertial element."
                << std::endl;
      return std::make_shared<Link>(linkName);
    }

    assert(inertial != nullptr);

    auto massEl = inertial->FirstChildElement("mass")->Attribute("value");
    assert(massEl != nullptr);
    auto inertiaEl = inertial->FirstChildElement("inertia");
    assert(inertiaEl != nullptr);
    auto originEl = inertial->FirstChildElement("origin");
    Eigen::Vector3d origin = Eigen::Vector3d::Zero();
    Eigen::Vector3d rpy = Eigen::Vector3d::Zero();

    if (originEl != nullptr) {
      auto originXYZ = originEl->Attribute("xyz");
      auto originRPY = originEl->Attribute("rpy");

      // parse origin
      if (originXYZ != nullptr) {
        std::stringstream ss(originXYZ);
        ss >> origin[0] >> origin[1] >> origin[2];
      }

      // parse rpy
      if (originRPY != nullptr) {
        std::stringstream ss(originRPY);
        ss >> rpy[0] >> rpy[1] >> rpy[2];
      }
    }

    // parse inertia
    auto inertiaXX = inertiaEl->Attribute("ixx");
    auto inertiaXY = inertiaEl->Attribute("ixy");
    auto inertiaXZ = inertiaEl->Attribute("ixz");
    auto inertiaYY = inertiaEl->Attribute("iyy");
    auto inertiaYZ = inertiaEl->Attribute("iyz");
    auto inertiaZZ = inertiaEl->Attribute("izz");

    // parse mass
    double massValue = 0.0;
    if (massEl != nullptr) {
      std::stringstream ss(massEl);
      ss >> massValue;
    }

    // parse inertia
    Eigen::VectorXd inertiaFlat(6);
    std::vector<const char *> inertiaAttrs = {inertiaXX, inertiaXY, inertiaXZ,
                                              inertiaYY, inertiaYZ, inertiaZZ};
    for (size_t i = 0; i < inertiaAttrs.size(); i++) {
      if (inertiaAttrs[i] != nullptr) {
        std::stringstream ss(inertiaAttrs[i]);
        ss >> inertiaFlat[i];
      }
    }

    return std::make_shared<Link>(std::string(linkName), origin, rpy, massValue,
                                  inertiaFlat);
    // return std::make_shared<Link>(linkName);
  }

  std::shared_ptr<Joint> parseJoint(const raisim::TiXmlElement &jointXML,
                                    std::shared_ptr<Link> parent,
                                    std::shared_ptr<Link> child) {
    // Optional fields
    auto axisField = jointXML.FirstChildElement("axis");
    Eigen::Vector3d axis = Eigen::Vector3d(0, 0, 1);
    if (axisField != nullptr) {
      auto axisStr = axisField->Attribute("xyz");
      // parse axis
      if (axisStr != nullptr) {
        std::stringstream ss(axisStr);
        ss >> axis[0] >> axis[1] >> axis[2];
      }
    }

    auto originField = jointXML.FirstChildElement("origin");
    Eigen::Vector3d originXYZ = Eigen::Vector3d::Zero();
    Eigen::Vector3d originRPY = Eigen::Vector3d::Zero();
    if (originField != nullptr) {
      auto originXYZStr = originField->Attribute("xyz");
      // parse origin xyz
      if (originXYZStr != nullptr) {
        std::stringstream ss(originXYZStr);
        ss >> originXYZ[0] >> originXYZ[1] >> originXYZ[2];
      }

      auto originRPYStr = originField->Attribute("rpy");
      // parse origin rpy
      if (originRPYStr != nullptr) {
        std::stringstream ss(originRPYStr);
        ss >> originRPY[0] >> originRPY[1] >> originRPY[2];
      }
    }

    // parse joint type
    auto jointType = jointXML.Attribute("type");
    JointType type = JointType::REVOLUTE;
    if (jointType != nullptr) {
      if (std::string(jointType) == "prismatic") {
        type = JointType::PRISMATIC;
      } else if (std::string(jointType) == "fixed") {
        type = JointType::FIXED;
      }
    }

    // parse joint name
    auto jointName = jointXML.Attribute("name");

    // create the joint
    return std::make_shared<Joint>(parent, child, axis, originXYZ, originRPY,
                                   type, jointName);
  }

  // int nq_; // number of generalized coordinates
  // int nv_; // number of generalized velocities
  // // int nbodies_; // number of bodies
  // int nframes_; // number of frames
  // int njoints_; // number of joints

  // // data due various types of joints and their generalized representation
  // std::vector<int> idx_qs_; // generalized coordinate indices
  // std::vector<int> nqs_;    // generalized coordinate lengths
  // std::vector<int> idx_vs_; // generalized velocity indices
  // std::vector<int> nvs_;    // generalized velocity lengths

  // // kinematic data
  // // joint-placement transformations wrt world frame
  // std::vector<Transform> jointPlacements;
  // std::vector<size_t> parentIdxs;
  // std::vector<size_t> childIdxs;
  // std::vector<std::string> jointNames;
  // std::vector<std::string> linkNames;

  // std::vector<std::shared_ptr<Link>> sortedLinks;
  // std::vector<std::shared_ptr<Joint>> sortedJoints;

  std::vector<std::shared_ptr<Link>> links_;
  std::vector<std::shared_ptr<Joint>> joints_;
  std::map<std::string, std::shared_ptr<Link>> linkMap_;

  // // for given index of generalized velocity, return the previous generalized
  // // velocity source
  // int dof_;
  // std::vector<int> dof_parent;

  // std::string _world_name = "world";
  // // those joints are the leaves of kinematic tree
  // std::vector<std::shared_ptr<Joint>> leaves_;

  // new structures
  int nbodies_;
  // an array of parent indices, size of Nbodies
  std::vector<size_t> parents_;
  std::vector<std::vector<size_t>> children_;
  // an array of joints, size of Nbodies
  std::vector<std::shared_ptr<Joint>> actuated_joints_;
  std::vector<std::shared_ptr<Link>> bodies_;

  // mapping from body index to generalized coordinate and velocity indices,
  // size of NBodies
  int gc_size_ = 0;
  int gv_size_ = 0;
  std::vector<size_t> gc_idx_;
  std::vector<size_t> gv_idx_;
  // tree transform array,
  // the position of the joint relative to the parent link
  std::vector<Transform> T_T_;
  std::vector<SpatialTransform> X_T_;

  // Joint *fixed_joints_ = nullptr;
};

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
    motionSubspace.resize(model.actuated_joints_.size());
    for (size_t i = 0; i < model.actuated_joints_.size(); i++) {
      rot_WB[i] = Eigen::Matrix3d::Identity();
      jointPos_W[i] = Eigen::Vector3d::Zero();
      jointAxis_W[i] = Eigen::Vector3d::Zero();
      bodyLinVel_w[i] = Eigen::Vector3d::Zero();
      bodyAngVel_w[i] = Eigen::Vector3d::Zero();
      joint2joint_W[i] = Eigen::Vector3d::Zero();
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
  std::vector<Eigen::MatrixXd> motionSubspace;

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
  std::vector<Eigen::MatrixXd> dMotionSubspace;
  Eigen::MatrixXd nonlinearities;

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

inline void nonlinearities(const Model &model, Data &data,
                           const Eigen::VectorXd &gc,
                           const Eigen::VectorXd &gv) {
  // we should have propagated the velocity through the kinematic chain already
  // TODO: maybe don't call at all
  algorithms::forwardVelocity(model, data, gc, gv);
  algorithms::forwardAcceleration(model, data, gc, gv, gv * 0.0, {0, 0, -9.81});

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
    auto S = data.motionSubspace[i];
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

#define _MAKE_STR(x) __MAKE_STR(x)
#define __MAKE_STR(x) #x

/// do not change the name of the method
inline Eigen::VectorXd getNonlinearities (const Eigen::VectorXd& gc, const Eigen::VectorXd& gv) {

  /// !!!!!!!!!! NO RAISIM FUNCTIONS HERE !!!!!!!!!!!!!!!!!
  std::string urdf =
      std::string(_MAKE_STR(RESOURCE_DIR)) + "/mini_cheetah/urdf/cheetah.urdf";

  auto model = Model(urdf);
  auto data = Data(model);
  algorithms::setState(model, data, gc, gv);
  algorithms::crba(model, data, gc);
  algorithms::nonlinearities(model, data, gc, gv);

  return data.nonlinearities;
}

