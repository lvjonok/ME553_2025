//
// Created by jemin on 25. 3. 13.
//


#include "raisim/World.hpp"
#include "raisim/RaisimServer.hpp"
#include <cstdlib> 

int main() {
  raisim::World world;
  world.addGround();
  raisim::RaisimServer server(&world);

  const char* resourceDir = std::getenv("RESOURCE_DIR");
  if (!resourceDir) {
    std::cerr << "RESOURCE_DIR environment variable not set" << std::endl;
    return 1;
  }
  
  auto panda = world.addArticulatedSystem(std::string(resourceDir) + "/Panda/panda.urdf");

  server.launchServer();

  while(true) {
    raisim::MSLEEP(4);
  }

  server.killServer();

  return 0;
}