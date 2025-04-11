#define _MAKE_STR(x) __MAKE_STR(x)
#define __MAKE_STR(x) #x

#include "raisim/RaisimServer.hpp"
#include "exercise2_STUDENTID.hpp"

int main(int argc, char* argv[]) {
  auto binaryPath = raisim::Path::setFromArgv(argv[0]);

  // create raisim world
  raisim::World world; // physics world
  raisim::RaisimServer server(&world); // visualization server
  world.addGround();
  world.setTimeStep(0.001);

  // panda
  auto panda = world.addArticulatedSystem(std::string(_MAKE_STR(RESOURCE_DIR)) + "/Panda/panda.urdf");
  panda->setName("panda");
  server.focusOn(panda);

  Eigen::VectorXd gc(panda->getGeneralizedCoordinateDim());
  Eigen::VectorXd gv(panda->getDOF());

  gc << 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.05, 0.05;
  gv << 0.12, 0.22, 0.32, 0.12, 0.42, 0.32, 0.13, 0.102, 0.101;
  panda->setState(gc, gv);

  // visualization
  server.launchServer();
  raisim::Vec<3> tipVel, tipAngVel;

  panda->getFrameVelocity("panda_finger_joint3", tipVel);
  panda->getFrameAngularVelocity("panda_finger_joint3", tipAngVel);

  if((tipVel.e() - getLinearVelocity(gc, gv)).norm() < 1e-8) {
    std::cout<<"the linear velocity is correct "<<std::endl;
  } else {
    std::cout<<"the linear velocity is not correct. It should be "<< tipVel.e().transpose() <<std::endl;
  }

  if((tipAngVel.e() - getAngularVelocity(gc, gv)).norm() < 1e-8) {
    std::cout<<"the angular velocity is correct "<<std::endl;
  } else {
    std::cout<<"the angular velocity is not correct. It should be "<< tipAngVel.e().transpose() <<std::endl;
  }

  for (int i=0; i<2000000; i++)
    std::this_thread::sleep_for(std::chrono::microseconds(1000));

  server.killServer();

  return 0;
}
