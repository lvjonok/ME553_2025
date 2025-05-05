#include "algorithms.hpp"
#include "model.hpp"
#include "raisim/RaisimServer.hpp"
#include <cassert>
#include <iostream>
#include <map>
#include <raisim/math.hpp>
#include <string>

#define _MAKE_STR(x) __MAKE_STR(x)
#define __MAKE_STR(x) #x

int main(int argc, char *argv[]) {
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
  };

  std::map<std::string, std::tuple<Eigen::VectorXd, Eigen::VectorXd>> robotData;
  // robotData["cart_pole"] =
  //     std::make_tuple(Eigen::VectorXd::Zero(3), Eigen::VectorXd::Zero(3));
  {
    // set random gc and gv for cart_pole
    Eigen::VectorXd cart_pole_gc = 0.5 * Eigen::VectorXd::Random(3);
    Eigen::VectorXd cart_pole_gv = 0.5 * Eigen::VectorXd::Random(3);
    robotData["cart_pole"] = std::make_tuple(cart_pole_gc, cart_pole_gv);
  }
  // set random gc and gv for panda
  {
    Eigen::VectorXd panda_gc = 0.3 * Eigen::VectorXd::Random(9);
    Eigen::VectorXd panda_gv = 0.3 * Eigen::VectorXd::Random(9);
    robotData["panda"] = std::make_tuple(panda_gc, panda_gv);
  }

  {
    Eigen::VectorXd cheetah_gc = 0.5 * Eigen::VectorXd::Random(19);
    Eigen::VectorXd cheetah_gv = 0.1 * Eigen::VectorXd::Random(18);

    Eigen::Vector4d random_quat = Eigen::Vector4d::Random();
    random_quat.normalize();
    cheetah_gc.segment<4>(3) = random_quat;
    cheetah_gv.segment<6>(0) = Eigen::VectorXd::Zero(6);
    // simulate velocity only for one leg
    // cheetah_gv.segment<9>(9).setZero();

    std::cout << "Cheetah gc: " << cheetah_gc.transpose() << '\n';
    std::cout << "Cheetah gv: " << cheetah_gv.transpose() << '\n';

    robotData["cheetah"] = std::make_tuple(cheetah_gc, cheetah_gv);
  }

  // lambda to test one robot
  auto testRobot = [&](const std::string &name, const std::string &path) {
    auto sys = world.addArticulatedSystem(path);
    auto model = Model(path);
    auto data = Data(model);
    std::cout << "=== Testing " << name << " ===\n";

    // set system state
    Eigen::VectorXd gc = std::get<0>(robotData[name]);
    Eigen::VectorXd gv = std::get<1>(robotData[name]);
    sys->setState(gc, gv);
    algorithms::setState(model, data, gc, gv);

    // // print Raisim info
    // std::cout << "Raisim parent array:\n";
    // for (auto p : sys->getParentVector())
    //   std::cout << p << ' ';
    // std::cout << "\nModel parent array:\n";
    // for (auto p : model.parents_)
    //   std::cout << p << ' ';
    // std::cout << '\n';

    // compare parent arrays
    for (size_t i = 1; i < sys->getParentVector().size(); ++i)
      assert(sys->getParentVector()[i] == model.parents_[i]);

    // compare joint counts
    assert(sys->getNumberOfJoints() == model.actuated_joints_.size());

    // std::cout << "Model gc indices:\n";
    // for (auto gcIdx :
    //      sys->getMappingFromBodyIndexToGeneralizedCoordinateIndex()) {
    //   std::cout << gcIdx << ' ';
    // }
    // std::cout << "\nModel gv indices:\n";
    // for (auto gvIdx :
    //      sys->getMappingFromBodyIndexToGeneralizedVelocityIndex()) {
    //   std::cout << gvIdx << ' ';
    // }
    // std::cout << '\n';

    // compare mapping between body and gc, gv indices
    for (size_t i = 0; i < sys->getBodyNames().size(); ++i) {
      auto raisimGcIdx =
          sys->getMappingFromBodyIndexToGeneralizedCoordinateIndex()[i];
      auto raisimGvIdx =
          sys->getMappingFromBodyIndexToGeneralizedVelocityIndex()[i];

      auto modelGcIdx = model.gc_idx_[i];
      auto modelGvIdx = model.gv_idx_[i];

      assert(raisimGcIdx == modelGcIdx);
      assert(raisimGvIdx == modelGvIdx);
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

      std::cout << "Body " << i << ": " << sys->getBodyNames()[i] << " "
                << model.bodies_[i]->getName() << '\n';

      // std::cout << "Raisim T:\n"
      //           << raisimR << "\n"
      //           << raisimP << "\n"
      //           << "Model T:\n"
      //           << modelR << "\n"
      //           << modelP << "\n";

      assert(raisimR.e().isApprox(modelR));
      assert(raisimP.e().isApprox(modelP));

      // compare the joint axis in the world frame
      auto raisimAxis = sys->jointAxis_W[i].e();
      auto modelAxis = data.jointAxis_W[i];
      // std::cout << "Raisim axis: " << raisimAxis.transpose() << '\n'
      //           << "Model axis: " << modelAxis.transpose() << '\n';
      // std::cout << "Why axis: " << sys->getJointAxis_P()[i].e().transpose()
      //           << '\n'
      //           << "Model axis: "
      //           << model.actuated_joints_[i]->getAxis().transpose() << '\n';
      assert(raisimAxis.isApprox(modelAxis));

      // compare joint2joint in the world frame
      {
        auto raisimJoint2Joint = sys->joint2joint_W[i].e();
        // std::cout << "Raisim joint2joint: " << raisimJoint2Joint.transpose()
        //           << '\n';
        // std::cout << "Model joint2joint: " <<
        // data.joint2joint_W[i].transpose()
        //           << '\n';
        assert(raisimJoint2Joint.isApprox(data.joint2joint_W[i]));
      }

      // compare body velocities
      sys->getVelocity(i, raisimLinVel);
      sys->getAngularVelocity(i, raisimAngVel);

      // get the twist of the body in the world frame
      algorithms::getBodyTwist(model, data, i, modelLinVel, modelAngVel);

      // std::cout << "Raisim velocity: " << raisimAngVel.e().transpose() << ' '
      //           << raisimLinVel.e().transpose() << '\n'
      //           << "Model velocity: " << modelAngVel.transpose() << ' '
      //           << modelLinVel.transpose() << '\n';
      assert(raisimLinVel.e().isApprox(modelLinVel));
      assert(raisimAngVel.e().isApprox(modelAngVel));

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
        std::cout << "Raisim com: " << raisimCom.transpose() << '\n'
                  << "Model com: " << modelCom.transpose() << '\n';

        assert(raisimCom.isApprox(modelCom));

        auto raisimInertia = sys->inertia_comW[i].e();
        auto modelInertia = data.inertiaW[i];

        std::cout << "Raisim inertia: " << raisimInertia.transpose() << '\n'
                  << "Model inertia: " << modelInertia << '\n';

        assert(raisimInertia.isApprox(modelInertia));
      }
    }
  };

  // loop over enabled robots
  for (auto &kv : robots) {
    testRobot(kv.first, kv.second);
  }

  return 0;
}
