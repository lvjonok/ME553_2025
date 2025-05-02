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

    nbodies_ = links_.size();
    X_T.resize(nbodies_);
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
    // now we can traverse the children and fill the parent array
    // parents.push_back(-1);
    // rootLink->setIndex(0);
    addLink(rootLink, -1);
  }

  void addLink(std::shared_ptr<Link> link, int parentIdx,
               bool isFixed = false) {
    // this function should be called to instantiate the structures and
    // append a new body

    if (!isFixed) {
      link->setIndex(parents.size());
      parents.push_back(parentIdx);
    }

    // find the joints for which this link is a parent
    for (auto joint : joints_) {
      if (joint->parent != link) {
        continue;
      }

      // this joint is attached to this link
      // we have to traverse the structure
      joint->setIndex(link->getIndex());
      if (isFixed) {
        addJoint(joint, parentIdx);
      } else {
        addJoint(joint, link->getIndex());
      }
    }
  }

  void addJoint(std::shared_ptr<Joint> joint, int parentIdx) {
    // this function should be called to instantiate the structures and
    // append a new joint

    // as we optimize fixed joints, we should indicate whether the next link is
    // connected through a fixed joint
    addLink(joint->child, parentIdx, joint->getType() == JointType::FIXED);
  }

  // nq_ = 0;
  // nv_ = 0;
  // nbodies_ = 0;
  // nframes_ = 0;
  // njoints_ = 0;
  // dof_ = -1;

  // // now we have to register the whole structure
  // // first body should be the world

  // // TODO: for mini cheetah the base link is BODY and there is no world
  // auto worldLink = std::make_shared<Link>("world");
  // std::cout << "worldLink: " << worldLink->getName() << std::endl;
  // auto floatingJoint = std::make_shared<Joint>(
  //     worldLink, links_[0], Eigen::Vector3d(0, 0, 0),
  //     Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(0, 0, 0),
  //     JointType::FLOATING, "floating");

  // // put worldLink in front of the list
  // links_.insert(links_.begin(), worldLink);
  // joints_.insert(joints_.begin(), floatingJoint);

  // // assert(world->getName() == "world");
  // _world_name = "world";
  // // auto worldLink = links_[0];
  // addLink(worldLink);

  // // we have parsed everything, now we can set parent dofs for joints for
  // // easier traversal
  // for (auto joint : joints_) {
  //   auto parentLink = joint->parent;
  //   auto parentJoint = linkUpward(parentLink);

  //   if (parentJoint != nullptr) {
  //     joint->setParentDof(parentJoint->getIndex());
  //   } else {
  //     joint->setParentDof(-1);
  //   }

  //   // std::cout << "joint name: " << joint->getName()
  //   //           << " parentDof: " << joint->getParentDof() << std::endl;
  // }

  // void addLink(std::shared_ptr<Link> link, int jointIdx = -1) {
  //   // this function should be called to instantiate the structures and
  //   append a
  //   // new body

  //   // we have to find all the children (joints which have this link as a
  //   // parent) and call add these joints
  //   link->setIndex(nbodies_);
  //   nbodies_++;

  //   auto no_joints = true;
  //   for (auto joint : joints_) {
  //     if (joint->parent == link) {
  //       no_joints = false;
  //       addJoint(joint, link->getIndex());
  //     }
  //   }

  //   if (no_joints) {
  //     // if there are no joints, jointIdx points towards the leaf joint
  //     // we should add the link to the list of joints
  //     leaves_.push_back(joints_[jointIdx]);
  //   }
  // }

  // void addJoint(std::shared_ptr<Joint> joint, size_t parentId) {
  //   // this function should be called to instantiate the structures and
  //   append a
  //   // new joint

  //   joint->setIndex(njoints_);
  //   njoints_++;

  //   // parent id is an index of the parent link we have added
  //   parentIdxs.push_back(parentId);
  //   childIdxs.push_back(joint->child->getIndex());
  //   jointNames.push_back(joint->getName());
  //   linkNames.push_back(joint->parent->getName());

  //   if (joint->getType() == JointType::FIXED) {
  //     // if the joint is fixed, we should not add any generalized coordinates
  //     // and velocities
  //     idx_qs_.push_back(-1);
  //     nqs_.push_back(-1);
  //     idx_vs_.push_back(-1);
  //     nvs_.push_back(-1);
  //   } else {
  //     idx_qs_.push_back(nq_);
  //     nqs_.push_back(joint->gc_length());
  //     idx_vs_.push_back(nv_);
  //     nvs_.push_back(joint->gv_length());
  //     // we should shift accordingly
  //     nq_ += joint->gc_length();
  //     nv_ += joint->gv_length();
  //   }

  //   // now we have to add a child of this link
  //   addLink(joint->child, joint->getIndex());
  // }

  // std::shared_ptr<Joint> linkUpward(std::shared_ptr<Link> link) const {
  //   // this function finds the joint where the link is the child of the joint
  //   // and returns the joint
  //   for (auto joint : joints_) {
  //     if (joint->child == link) {
  //       return joint;
  //     }
  //   }

  //   std::cerr << "Link " << link->getName() << " is not a child of any
  //   joint."
  //             << std::endl;
  //   return nullptr;
  // }

  // bool isLeaf(std::shared_ptr<Link> link) const {
  //   // this function checks if the link is a leaf
  //   for (auto joint : joints_) {
  //     if (joint->parent == link) {
  //       return false;
  //     }
  //   }

  //   return true;
  // }

  // std::vector<std::shared_ptr<Joint>>
  // getAttachedJoints(std::shared_ptr<Link> link) const {
  //   // this function returns all the joints that are attached to the link
  //   std::vector<std::shared_ptr<Joint>> attachedJoints;
  //   for (auto joint : joints_) {
  //     if (joint->parent == link) {
  //       attachedJoints.push_back(joint);
  //     }
  //   }

  //   return attachedJoints;
  // }

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
    assert(originEl != nullptr);

    auto originXYZ = originEl->Attribute("xyz");
    auto originRPY = originEl->Attribute("rpy");

    // parse origin
    Eigen::Vector3d origin = Eigen::Vector3d::Zero();
    if (originXYZ != nullptr) {
      std::stringstream ss(originXYZ);
      ss >> origin[0] >> origin[1] >> origin[2];
    }

    // parse rpy
    Eigen::Vector3d rpy = Eigen::Vector3d::Zero();
    if (originRPY != nullptr) {
      std::stringstream ss(originRPY);
      ss >> rpy[0] >> rpy[1] >> rpy[2];
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
  std::vector<std::shared_ptr<Joint>> joints;
  std::vector<std::shared_ptr<Link>> links;
  // tree transform array,
  // the position of the joint relative to the parent link
  std::vector<Transform> X_T;
};

#endif