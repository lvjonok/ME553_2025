#ifndef ME553_2025_MODEL_HPP_
#define ME553_2025_MODEL_HPP_

#include <Eigen/Core>
#include <cstddef>
#include <eigen3/Eigen/src/Core/Matrix.h>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <tinyxml_rai/tinyxml_rai.h>
#include <vector>

using Transform = Eigen::Matrix4d;

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
      : parent(parent), child(child), axis(axis), xyz(xyz), rpy(rpy),
        type_(type), joint_name(jointName) {}

  // compute the transform from parent to child
  // this is const and defined by the URDF only
  Transform jointPlacement() {
    // euler angles
    Eigen::Matrix3d R = RotZ(rpy[2]) * RotY(rpy[1]) * RotX(rpy[0]);

    Transform T = Transform::Identity();
    T.block<3, 3>(0, 0) = R;
    T.block<3, 1>(0, 3) = xyz;
    // std::cout << "T: " << T << std::endl;
    // std::cout << "xyz: " << xyz << std::endl;

    return T;
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

  // Motion transformation is defined as a function of generalized coordinates
  // we pass the start index and length of the generalized coordinates
  Transform motion(const Eigen::VectorXd &gc, int start_idx, int length) {
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

  // origin definition
  Eigen::Vector3d xyz; // default is [0, 0, 0]
  // XYZ Euler angle definition
  Eigen::Vector3d rpy; // default is [0, 0, 0]

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

    parents.resize(0);

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

    std::cout << "root link: " << rootLink->getName() << std::endl;

    // if rootLink is "world", we have a fixed-base robot
    // otherwise, we have a floating base robot and need to manually add
    // the floating joint

    bool isFixedBase = true;

    if (rootLink->getName() != "world") {
      std::cout << "We have a floating base robot" << std::endl;
      std::cout << "Adding floating joint and world under the hood"
                << std::endl;

      // make a dummy world link
      auto worldLink = std::make_shared<Link>("world");
      std::cout << "worldLink: " << worldLink->getName() << std::endl;
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
    parents.clear();
    traverse(links_[0], -1);

    std::cout << "printing bodies and actuated joints in order" << std::endl;
    nbodies_ = bodies.size();
    for (size_t i = 0; i < nbodies_; i++) {
      auto link = bodies[i];
      std::cout << "link name: " << link->getName() << std::endl;
    }
    std::cout << "printing actuated joints in order" << std::endl;
    for (size_t i = 0; i < actuated_joints.size(); i++) {
      auto joint = actuated_joints[i];
      std::cout << "joint name: " << joint->getName() << std::endl;
    }

    // // now as we have parsed everything, we can set the tree transformation
    // X_T.resize(nbodies_);
    // for (size_t i = 0; i < nbodies_; i++) {
    //   auto joint = joints_[i];

    //   X_T[i] = joint->jointPlacement();
    // }
  }

  void traverse(std::shared_ptr<Link> link, int currentBodyId) {
    // TODO: we should introduce a way to merge bodies that are connected
    // through a fixed joint this is important for dynamics computations, but
    // may be overkill for now

    for (const auto &joint : joints_) {
      // skip joints that are not related to the current link
      if (joint->parent != link)
        continue;

      // this joint is attached to this link
      if (joint->getType() == JointType::FIXED && currentBodyId == -1) {
        int newBodyId = nbodies_++;
        parents.push_back(currentBodyId);
        bodies.push_back(joint->child); // world -> fixed_joint -> child

        // TODO: this is a little weird, but even RAILAB does this
        // we consider first fixed joint for fixed-base robots
        actuated_joints.push_back(joint);

        traverse(joint->child, newBodyId);
        return;
      }

      if (joint->getType() == JointType::FIXED) {
        // this is a fixed joint, it should not be added to the bodies
        traverse(joint->child, currentBodyId);
      } else {
        // this is a new Body
        int newBodyId = nbodies_++;
        parents.push_back(currentBodyId);
        bodies.push_back(joint->child); // link1 -> fixed_joint -> link2
        actuated_joints.push_back(joint);

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
  std::vector<size_t> parents;
  // an array of joints, size of Nbodies
  std::vector<std::shared_ptr<Joint>> actuated_joints;
  std::vector<std::shared_ptr<Link>> bodies;
  // tree transform array,
  // the position of the joint relative to the parent link
  std::vector<Transform> X_T;
};

#endif