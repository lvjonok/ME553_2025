#ifndef ME553_2025_MODEL_HPP_
#define ME553_2025_MODEL_HPP_

#include <Eigen/Core>
#include <cstddef>
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

  std::string getName() const { return name; }

  void setIndex(int idx) { index = idx; }

  int getIndex() const { return index; }

private:
  std::string name;

  int index;
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

    nq_ = 0;
    nv_ = 0;
    nbodies_ = 0;
    nframes_ = 0;
    njoints_ = 0;
    dof_ = -1;

    // now we have to register the whole structure
    // first body should be the world

    // TODO: for mini cheetah the base link is BODY and there is no world
    auto worldLink = std::make_shared<Link>("world");
    std::cout << "worldLink: " << worldLink->getName() << std::endl;
    auto floatingJoint = std::make_shared<Joint>(
        worldLink, links_[0], Eigen::Vector3d(0, 0, 0),
        Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(0, 0, 0), JointType::FLOATING,
        "floating");

    // put worldLink in front of the list
    links_.insert(links_.begin(), worldLink);
    joints_.insert(joints_.begin(), floatingJoint);

    // assert(world->getName() == "world");
    _world_name = "world";

    addLink(worldLink);

    // now each link and joint has index, let's construct vectors
    // with the correct indices
    sortedLinks.resize(links_.size());
    sortedJoints.resize(joints_.size());
    for (auto link : links_) {
      sortedLinks[link->getIndex()] = link;
    }
    for (auto joint : joints_) {
      sortedJoints[joint->getIndex()] = joint;
    }
  }

  void addLink(std::shared_ptr<Link> link, int jointIdx = -1) {
    // this function should be called to instantiate the structures and append a
    // new body

    // we have to find all the children (joints which have this link as a
    // parent) and call add these joints
    link->setIndex(nbodies_);
    nbodies_++;

    for (auto joint : joints_) {
      if (joint->parent == link) {
        addJoint(joint, link->getIndex());
      }
    }
  }

  void addJoint(std::shared_ptr<Joint> joint, size_t parentId) {
    // this function should be called to instantiate the structures and append a
    // new joint

    joint->setIndex(njoints_);
    njoints_++;

    // parent id is an index of the parent link we have added
    parentIdxs.push_back(parentId);
    childIdxs.push_back(joint->child->getIndex());
    jointNames.push_back(joint->getName());
    linkNames.push_back(joint->parent->getName());

    // right now we have to find out, which actuated joint is the parent of this
    auto parentDof = -1;
    std::cout << "Looking for parentDof for joint " << joint->getName()
              << " with index " << joint->getIndex() << std::endl;

    auto iterJoint = joint;
    while (iterJoint->parent->getName() != _world_name) {

      bool foundActuatedParent = false;

      while (!foundActuatedParent) {
        auto parentLink = iterJoint->parent;
        // find for which joint the parent link is the child
        for (auto j : joints_) {
          if (j->child == parentLink) {
            iterJoint = j;
            break;
          }
        }

        // std::cout << "iterJoint: " << iterJoint->getName()
        //           << " parentName: " << iterJoint->parent->getName()
        //           << std::endl;

        // if we are at the world, we should break
        if (iterJoint->parent->getName() == _world_name) {
          break;
        }

        // if we found the actuated joint, we can break
        if (iterJoint->getType() != JointType::FIXED) {
          foundActuatedParent = true;
          parentDof = iterJoint->getIndex();
          break;
        }
      }

      // Set the parentDof for the joint
      joint->setParentDof(parentDof);
      std::cout << "Found parentDof: " << parentDof << std::endl;
      break;
    }

    // // right now we have to find joint for which we are the child
    // auto parentDof = -1;
    // for (auto j : joints_) {
    //   if (j->child == joint->parent) {
    //     parentDof = j->getIndex();
    //     break;
    //   }
    // }
    // joint->setParentDof(parentDof);

    if (joint->getType() == JointType::FIXED) {
      // if the joint is fixed, we should not add any generalized coordinates
      // and velocities
      idx_qs_.push_back(-1);
      nqs_.push_back(-1);
      idx_vs_.push_back(-1);
      nvs_.push_back(-1);
    } else {
      idx_qs_.push_back(nq_);
      nqs_.push_back(joint->gc_length());
      idx_vs_.push_back(nv_);
      nvs_.push_back(joint->gv_length());
      // we should shift accordingly
      nq_ += joint->gc_length();
      nv_ += joint->gv_length();
    }

    // now we have to add a child of this link
    addLink(joint->child, joint->getIndex());
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
    return std::make_shared<Link>(linkName);
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

  int nq_;      // number of generalized coordinates
  int nv_;      // number of generalized velocities
  int nbodies_; // number of bodies
  int nframes_; // number of frames
  int njoints_; // number of joints

  // data due various types of joints and their generalized representation
  std::vector<int> idx_qs_; // generalized coordinate indices
  std::vector<int> nqs_;    // generalized coordinate lengths
  std::vector<int> idx_vs_; // generalized velocity indices
  std::vector<int> nvs_;    // generalized velocity lengths

  // kinematic data
  // joint-placement transformations wrt world frame
  std::vector<Transform> jointPlacements;
  std::vector<size_t> parentIdxs;
  std::vector<size_t> childIdxs;
  std::vector<std::string> jointNames;
  std::vector<std::string> linkNames;

  std::vector<std::shared_ptr<Link>> sortedLinks;
  std::vector<std::shared_ptr<Joint>> sortedJoints;

  std::vector<std::shared_ptr<Link>> links_;
  std::vector<std::shared_ptr<Joint>> joints_;
  std::map<std::string, std::shared_ptr<Link>> linkMap_;

  // for given index of generalized velocity, return the previous generalized
  // velocity source
  int dof_;
  std::vector<int> dof_parent;

  std::string _world_name = "world";
};

#endif