#include "algorithms.hpp"
#include "model.hpp"
#include "raisim/RaisimServer.hpp"
#include <cassert>
#include <cstdlib> // getenv
#include <iostream>
#include <map>
#include <raisim/math.hpp>
#include <string>

// tolerance for approximate comparisons
constexpr double APPROX_THRESHOLD = 1e-10;

// print and abort if x!=y
#define ASSERT_EQ(x, y)                                                        \
  do {                                                                         \
    auto _x = (x);                                                             \
    auto _y = (y);                                                             \
    if (_x != _y) {                                                            \
      std::cerr << "[ASSERT_EQ] " #x " = " << _x << ", " #y " = " << _y        \
                << ", diff = " << (_x - _y) << std::endl;                      \
    }                                                                          \
    assert(_x == _y);                                                          \
  } while (0)

// print and abort if a and b not approx; shows full matrices and difference
#define ASSERT_APPROX(a, b)                                                    \
  do {                                                                         \
    auto _a = (a);                                                             \
    auto _b = (b);                                                             \
    if (!_a.isApprox(_b, APPROX_THRESHOLD)) {                                  \
      std::cerr << "[ASSERT_APPROX] " #a " vs " #b "\n"                        \
                << _a << "\n"                                                  \
                << _b << "\n diff:\n"                                          \
                << (_a - _b) << std::endl;                                     \
    }                                                                          \
    assert(_a.isApprox(_b, APPROX_THRESHOLD));                                 \
  } while (0)

int main(int argc, char *argv[]) {
  // parse -N option
  int N = 1;
  for (int i = 1; i < argc; ++i) {
    if (std::string(argv[i]) == "-N" && i + 1 < argc) {
      N = std::stoi(argv[++i]);
    }
  }

  auto binaryPath = raisim::Path::setFromArgv(argv[0]);

  // create raisim world
  raisim::World world;
  raisim::RaisimServer server(&world);

  const std::string RESOURCE =
      "/home/lvjonok/github.com/lvjonok/ME553_2025/resource/";

  // map of robot name → URDF path; comment out entries to disable
  std::map<std::string, std::string> robots = {
      {"cheetah", RESOURCE + "mini_cheetah/urdf/cheetah.urdf"},
      {"panda", RESOURCE + "Panda/panda.urdf"},
      {"cart_pole", RESOURCE + "cartPole/doubleCartPole.urdf"},
      {"3drobot", RESOURCE + "2DRobotArm/robot_3D.urdf"},
  };

  // preload systems, models, data, and dims
  std::map<std::string, raisim::ArticulatedSystem *> systems;
  std::map<std::string, Model> models;
  std::map<std::string, Data> datas;
  std::map<std::string, std::pair<int, int>> dims;
  for (auto &kv : robots) {
    auto sys = world.addArticulatedSystem(kv.second);
    systems[kv.first] = sys;
    models.emplace(kv.first, Model(kv.second));
    datas.emplace(kv.first, Data(models.at(kv.first)));
    dims[kv.first] = {sys->getGeneralizedCoordinateDim(), sys->getDOF()};
  }

  // generator for random gc/gv per robot (uses preloaded dims)
  auto genState = [&](const std::string &name) {
    auto [nq, nv] = dims[name];
    Eigen::VectorXd gc = 0.3 * Eigen::VectorXd::Random(nq);
    Eigen::VectorXd gv = 0.3 * Eigen::VectorXd::Random(nv);
    if (nq != nv) {
      Eigen::Vector4d quat = Eigen::Vector4d::Random();
      quat.normalize();
      gc.segment<4>(3) = quat;
    }
    return std::make_pair(gc, gv);
  };

  // testRobot now reuses preloaded sys/model/data
  auto testRobot = [&](const std::string &name, const Eigen::VectorXd &gc,
                       const Eigen::VectorXd &gv) {
    auto sys = systems[name];
    auto &model = models.at(name);
    auto &data = datas.at(name);

    sys->setState(gc, gv);
    auto raisimMassMatrix = sys->getMassMatrix().e();
    algorithms::setState(model, data, gc, gv);
    algorithms::crba(model, data, gc);

    // compare parent arrays
    for (size_t i = 1; i < sys->getParentVector().size(); ++i)
      ASSERT_EQ(sys->getParentVector()[i], model.parents_[i]);

    // compare joint counts
    ASSERT_EQ(sys->getNumberOfJoints(), model.actuated_joints_.size());

    {
      // compare children_ array
      ASSERT_EQ(sys->children_.size(), model.children_.size());
      for (size_t i = 0; i < sys->children_.size(); ++i) {
        ASSERT_EQ(sys->children_[i].size(), model.children_[i].size());
        for (size_t j = 0; j < sys->children_[i].size(); ++j) {
          ASSERT_EQ(sys->children_[i][j], model.children_[i][j]);
        }
        for (size_t j = 0; j < model.children_[i].size(); ++j) {
          ASSERT_EQ(sys->children_[i][j], model.children_[i][j]);
        }
      }
    }

    // compare mapping between body and gc, gv indices
    for (size_t i = 0; i < sys->getBodyNames().size(); ++i) {
      auto raisimGcIdx =
          sys->getMappingFromBodyIndexToGeneralizedCoordinateIndex()[i];
      auto raisimGvIdx =
          sys->getMappingFromBodyIndexToGeneralizedVelocityIndex()[i];

      auto modelGcIdx = model.gc_idx_[i];
      auto modelGvIdx = model.gv_idx_[i];

      ASSERT_EQ(raisimGcIdx, modelGcIdx);
      ASSERT_EQ(raisimGvIdx, modelGvIdx);
    }

    // compare body transformations
    for (size_t i = 0; i < sys->getBodyNames().size(); ++i) {
      raisim::Mat<3, 3> raisimR;
      raisim::Vec<3> raisimP, raisimLinVel, raisimAngVel;

      Eigen::Matrix3d modelR;
      Eigen::Vector3d modelP, modelLinVel, modelAngVel;

      // get the pose of the body in the world frame
      sys->getBodyPose(i, raisimR, raisimP);

      // get the pose of the body in the world frame
      algorithms::getBodyPose(model, data, i, modelR, modelP);

      // std::cout << "Body " << i << ": " << sys->getBodyNames()[i] << " "
      //           << model.bodies_[i]->getName() << '\n';

      ASSERT_APPROX(raisimR.e(), modelR);
      ASSERT_APPROX(raisimP.e(), modelP);

      // compare the joint axis in the world frame
      auto raisimAxis = sys->jointAxis_W[i].e();
      auto modelAxis = data.jointAxis_W[i];
      ASSERT_APPROX(raisimAxis, modelAxis);

      // compare joint2joint in the world frame
      {
        auto raisimJoint2Joint = sys->joint2joint_W[i].e();
        ASSERT_APPROX(raisimJoint2Joint, data.joint2joint_W[i]);
      }

      // compare body velocities
      sys->getVelocity(i, raisimLinVel);
      sys->getAngularVelocity(i, raisimAngVel);

      // get the twist of the body in the world frame
      algorithms::getBodyTwist(model, data, i, modelLinVel, modelAngVel);

      ASSERT_APPROX(raisimLinVel.e(), modelLinVel);
      ASSERT_APPROX(raisimAngVel.e(), modelAngVel);

      auto raisimMassMatrix = sys->getMassMatrix();

      {
        // compare that we have parsed the inertia correctly
        auto raisimInertia = sys->getInertia()[i].e();
        auto modelInertia = model.bodies_[i]->getInertia();
      }

      {
        // compare the inertia and com of bodies in the world frame
        auto raisimCom = sys->comPos_W[i].e();
        auto modelCom = data.comW[i];

        ASSERT_APPROX(raisimCom, modelCom);

        auto raisimInertia = sys->inertia_comW[i].e();
        auto modelInertia = data.inertiaW[i];

        ASSERT_APPROX(raisimInertia, modelInertia);
      }

      // for panda there is no sense to compare this
      if (name != "panda") {
        // compare composite mass, com and inertia
        auto raisimM = sys->compositeMass[i];
        auto raisimCom = sys->composite_com_W[i].e();
        auto raisimInertia = sys->compositeInertia_W[i].e();

        auto modelM = data.compositeMassW[i];
        auto modelCom = data.compositeComW[i];
        auto modelInertia = data.compositeInertiaW[i];

        ASSERT_EQ(raisimM, modelM);
        ASSERT_APPROX(raisimCom, modelCom);
        ASSERT_APPROX(raisimInertia, modelInertia);
      }
    }
    if (name != "panda") {
      // compare the mass matrix
      algorithms::crba(model, data, gc);
      ASSERT_APPROX(raisimMassMatrix, data.massMatrix);
    }
  };

  // loop robots and samples
  for (auto &kv : robots) {
    for (int i = 0; i < N; ++i) {
      auto [gc, gv] = genState(kv.first);
      testRobot(kv.first, gc, gv);
    }
  }

  return 0;
}
