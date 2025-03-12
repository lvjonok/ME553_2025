//
// Created by jemin on 25. 3. 13.
//


#include "raisim/World.hpp"
#include "raisim/RaisimServer.hpp"

#define _MAKE_STR(x) __MAKE_STR(x)
#define __MAKE_STR(x) #x

int main() {
  raisim::World world;
  world.addGround();

  raisim::RaisimServer server(&world);

  auto panda = world.addArticulatedSystem(std::string(_MAKE_STR(RESOURCE_DIR)) + "/Panda/panda.urdf");

  server.launchServer();

  while(true) {
    raisim::MSLEEP(4);
  }

  server.killServer();

  return 0;
}