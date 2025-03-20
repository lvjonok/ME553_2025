//
// Created by Jemin Hwangbo on 2022/03/17.
//

#define _MAKE_STR(x) __MAKE_STR(x)
#define __MAKE_STR(x) #x

#include "exercise1_STUDENTID.hpp"
#include "raisim/RaisimServer.hpp"


int main(int argc, char* argv[]) {
  // create raisim world
  raisim::World world; // physics world
  raisim::RaisimServer server(&world); // visualization server
  world.addGround();

  // panda
  auto panda = world.addArticulatedSystem(std::string(_MAKE_STR(RESOURCE_DIR)) + "/Panda/panda.urdf");
  panda->setName("panda");
  server.focusOn(panda);

  // panda configuration
  Eigen::VectorXd jointNominalConfig(panda->getGeneralizedCoordinateDim());
  jointNominalConfig << 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.05, 0.05;
  panda->setGeneralizedCoordinate(jointNominalConfig);
  panda->updateKinematics();

  // debug sphere
  auto debugSphere = server.addVisualSphere("debug_sphere", 0.02);
  debugSphere->setColor(1,0,0,1);
  debugSphere->setPosition(getEndEffectorPosition(jointNominalConfig));

  // solution sphere
  auto answerSphere = server.addVisualSphere("answer_sphere", 0.02);
  answerSphere->setColor(0,1,0,1);
  raisim::Vec<3> pos;
  panda->getFramePosition("panda_finger_joint3", pos);
  answerSphere->setPosition(pos.e());

  // visualization
  server.launchServer();
  for (int i=0; i<2000000; i++)
    std::this_thread::sleep_for(std::chrono::microseconds(1000));

  server.killServer();
}
