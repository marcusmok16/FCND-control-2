#include "Common.h"
#include "QuadControl.h"

#include "Utility/SimpleConfig.h"

#include "Utility/StringUtils.h"
#include "Trajectory.h"
#include "BaseController.h"
#include "Math/Mat3x3F.h"

#include <iostream>
using namespace std;

#ifdef __PX4_NUTTX
#include <systemlib/param/param.h>
#endif

void QuadControl::Init()
{
  BaseController::Init();

  // variables needed for integral control
  integratedAltitudeError = 0;
    
#ifndef __PX4_NUTTX
  // Load params from simulator parameter system
  ParamsHandle config = SimpleConfig::GetInstance();
   
  // Load parameters (default to 0)
  kpPosXY = config->Get(_config+".kpPosXY", 0);
  kpPosZ = config->Get(_config + ".kpPosZ", 0);
  KiPosZ = config->Get(_config + ".KiPosZ", 0);
     
  kpVelXY = config->Get(_config + ".kpVelXY", 0);
  kpVelZ = config->Get(_config + ".kpVelZ", 0);

  kpBank = config->Get(_config + ".kpBank", 0);
  kpYaw = config->Get(_config + ".kpYaw", 0);

  kpPQR = config->Get(_config + ".kpPQR", V3F());

  maxDescentRate = config->Get(_config + ".maxDescentRate", 100);
  maxAscentRate = config->Get(_config + ".maxAscentRate", 100);
  maxSpeedXY = config->Get(_config + ".maxSpeedXY", 100);
  maxAccelXY = config->Get(_config + ".maxHorizAccel", 100);

  maxTiltAngle = config->Get(_config + ".maxTiltAngle", 100);

  minMotorThrust = config->Get(_config + ".minMotorThrust", 0);
  maxMotorThrust = config->Get(_config + ".maxMotorThrust", 100);

  L = config->Get(_config + ".L", 100);
  kappa = config->Get(_config + ".kappa", 100);
  Ixx = config->Get(_config + ".Ixx", 100);
  Iyy = config->Get(_config + ".Iyy", 100);
  Izz = config->Get(_config + ".Izz", 100);
  
  

#else
  // load params from PX4 parameter system
  //TODO
  param_get(param_find("MC_PITCH_P"), &Kp_bank);
  param_get(param_find("MC_YAW_P"), &Kp_yaw);
#endif
}

VehicleCommand QuadControl::GenerateMotorCommands(float collThrustCmd, V3F momentCmd)
{
  // Convert a desired 3-axis moment and collective thrust command to 
  //   individual motor thrust commands
  // INPUTS: 
  //   collThrustCmd: desired collective thrust [N]
  //   momentCmd: desired rotation moment about each axis [N m]
  // OUTPUT:
  //   set class member variable cmd (class variable for graphing) where
  //   cmd.desiredThrustsN[0..3]: motor commands, in [N]

  // HINTS: 
  // - you can access parts of momentCmd via e.g. momentCmd.x
  // You'll need the arm length parameter L, and the drag/thrust ratio kappa

  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////

  //cmd.desiredThrustsN[0] = mass * 9.81f / 4.f; // front left
  //cmd.desiredThrustsN[1] = mass * 9.81f / 4.f; // front right
  //cmd.desiredThrustsN[2] = mass * 9.81f / 4.f; // rear left
  //cmd.desiredThrustsN[3] = mass * 9.81f / 4.f; // rear right
   
  float L2 = L / (2.f* sqrt(2.f));
  cmd.desiredThrustsN[0] = CONSTRAIN((collThrustCmd + momentCmd.x / L2 + momentCmd.y / L2 + momentCmd.z / kappa)/4.f, minMotorThrust, maxMotorThrust); //front left
  cmd.desiredThrustsN[1] = CONSTRAIN((collThrustCmd - momentCmd.x / L2 + momentCmd.y / L2 - momentCmd.z / kappa)/4.f, minMotorThrust, maxMotorThrust); // front right
  cmd.desiredThrustsN[2] = CONSTRAIN((collThrustCmd + momentCmd.x / L2 - momentCmd.y / L2 - momentCmd.z / kappa)/4.f, minMotorThrust, maxMotorThrust); // rear left
  cmd.desiredThrustsN[3] = CONSTRAIN((collThrustCmd - momentCmd.x / L2 - momentCmd.y / L2 + momentCmd.z / kappa)/4.f, minMotorThrust, maxMotorThrust); // rear right
  
  // is Z in the correct sense ?
  /////////////////////////////// END STUDENT CODE ////////////////////////////

  return cmd;
}

V3F QuadControl::BodyRateControl(V3F pqrCmd, V3F pqr)
{
  // Calculate a desired 3-axis moment given a desired and current body rate
  // INPUTS: 
  //   pqrCmd: desired body rates [rad/s]
  //   pqr: current or estimated body rates [rad/s]
  // OUTPUT:
  //   return a V3F containing the desired moments for each of the 3 axes

  // HINTS: 
  //  - you can use V3Fs just like scalars: V3F a(1,1,1), b(2,3,4), c; c=a-b;
  //  - you'll need parameters for moments of inertia Ixx, Iyy, Izz
  //  - you'll also need the gain parameter kpPQR (it's a V3F)

  V3F momentCmd;
    ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////
 
  V3F I = V3F(Ixx, Iyy, Izz);
  V3F pqr_error = pqrCmd - pqr;
  momentCmd = I * kpPQR * pqr_error;    //[M]=[I].[W^2] = [I] x [Kpqr] x [PQRerror]. 
  
 
  return momentCmd;
}

// returns a desired roll and pitch rate 
V3F QuadControl::RollPitchControl(V3F accelCmd, Quaternion<float> attitude, float collThrustCmd)
{
  // Calculate a desired pitch and roll angle rates based on a desired global
  //   lateral acceleration, the current attitude of the quad, and desired
  //   collective thrust command
  // INPUTS: 
  //   accelCmd: desired acceleration in global XY coordinates [m/s2]
  //   attitude: current or estimated attitude of the vehicle
  //   collThrustCmd: desired collective thrust of the quad [N]
  // OUTPUT:
  //   return a V3F containing the desired pitch and roll rates. The Z
  //     element of the V3F should be left at its default value (0)

  // HINTS: 
  //  - we already provide rotation matrix R: to get element R[1,2] (python) use R(1,2) (C++)
  //  - you'll need the roll/pitch gain kpBank
  //  - collThrustCmd is a force in Newtons! You'll likely want to convert it to acceleration first

  V3F pqrCmd;
  Mat3x3F R = attitude.RotationMatrix_IwrtB();

  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////
 float b_x_c, b_y_c, b_x_a, b_y_a, b_x_c_dot, b_y_c_dot;
  if (collThrustCmd > 0.f)
  {
      float c = -collThrustCmd / mass;
      b_x_c = CONSTRAIN(accelCmd.x / c, -maxTiltAngle, maxTiltAngle);
      b_y_c = CONSTRAIN(accelCmd.y / c, -maxTiltAngle, maxTiltAngle);
      b_x_a = R(0, 2);
      b_y_a = R(1, 2);
      b_x_c_dot = (b_x_c - b_x_a);
      b_y_c_dot = (b_y_c - b_y_a);
      
      pqrCmd.x = kpBank * ((R(1, 0)*b_x_c_dot) - (R(0, 0)*b_y_c_dot)) / R(2, 2);
      pqrCmd.y = kpBank * ((R(1, 1)*b_x_c_dot) - (R(0, 1)*b_y_c_dot)) / R(2, 2);
      pqrCmd.z = 0.f;
  }
  
  else 
  {
      pqrCmd.x = 0.f;
      pqrCmd.y = 0.f;
      pqrCmd.z = 0.f;
  }
  
  
  /////////////////////////////// END STUDENT CODE ////////////////////////////
 return pqrCmd;
}

float QuadControl::AltitudeControl(float posZCmd, float velZCmd, float posZ, float velZ, Quaternion<float> attitude, float accelZCmd, float dt)
{
  // Calculate desired quad thrust based on altitude setpoint, actual altitude,
  //   vertical velocity setpoint, actual vertical velocity, and a vertical 
  //   acceleration feed-forward command
  // INPUTS: 
  //   posZCmd, velZCmd: desired vertical position and velocity in NED [m]
  //   posZ, velZ: current vertical position and velocity in NED [m]
  //   accelZCmd: feed-forward vertical acceleration in NED [m/s2]
  //   dt: the time step of the measurements [seconds]
  // OUTPUT:
  //   return a collective thrust command in [N]

  // HINTS: 
  //  - we already provide rotation matrix R: to get element R[1,2] (python) use R(1,2) (C++)
  //  - you'll need the gain parameters kpPosZ and kpVelZ
  //  - maxAscentRate and maxDescentRate are maximum vertical speeds. Note they're both >=0!
  //  - make sure to return a force, not an acceleration
  //  - remember that for an upright quad in NED, thrust should be HIGHER if the desired Z acceleration is LOWER

  Mat3x3F R = attitude.RotationMatrix_IwrtB();
  float thrust = 0;
  float z_error, z_dot_error, u1_bar = 0;
  float g = 9.81;
  float h_dot;
  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////
  /*
  z_error = posZCmd - posZ; //P term
  z_dot_error = velZCmd - velZ; //D term
  h_dot = kpPosZ * z_error + velZCmd;
  h_dot = CONSTRAIN(h_dot, -maxDescentRate, maxAscentRate);

  u1_bar = (kpPosZ * z_error) + (KiPosZ * dt * z_error) + kpVelZ * (h_dot - velZ) + accelZCmd;

  //u1_bar = (kpPosZ * z_error) + (dt * z_error) + (kpVelZ * z_dot_error) + accelZCmd;
  thrust = (u1_bar - g) / R(2, 2) * mass * -1; //Collective thrust
  thrust = CONSTRAIN(thrust, minMotorThrust, maxMotorThrust);
  */
  /////////////////////////////// END STUDENT CODE ////////////////////////////
  
  return thrust;
}

// returns a desired acceleration in global frame
V3F QuadControl::LateralPositionControl(V3F posCmd, V3F velCmd, V3F pos, V3F vel, V3F accelCmdFF)
{
  // Calculate a desired horizontal acceleration based on 
  //  desired lateral position/velocity/acceleration and current pose
  // INPUTS: 
  //   posCmd: desired position, in NED [m]
  //   velCmd: desired velocity, in NED [m/s]
  //   pos: current position, NED [m]
  //   vel: current velocity, NED [m/s]
  //   accelCmdFF: feed-forward acceleration, NED [m/s2]
  // OUTPUT:
  //   return a V3F with desired horizontal accelerations. 
  //     the Z component should be 0
  // HINTS: 
  //  - use the gain parameters kpPosXY and kpVelXY
  //  - make sure you limit the maximum horizontal velocity and acceleration
  //    to maxSpeedXY and maxAccelXY

  // make sure we don't have any incoming z-component
  accelCmdFF.z = 0;
  velCmd.z = 0;
  posCmd.z = pos.z;

  // we initialize the returned desired acceleration to the feed-forward value.
  // Make sure to _add_, not simply replace, the result of your controller
  // to this variable
  V3F accelCmd = accelCmdFF;

  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////
  
  //limit the maximum horizontal velocity and acceleration to maxSpeedXY and maxAccelXY
  //if velCmd.x > maxSpeedXY
  //accelCmdFF.x = kpPosXY * (posCmd.x - pos.x) + kpVelXY * (velCmd.x - vel.x) + accelCmdFF.x; // x_c_d_d
  // Code below written by Cedric Parait
  /*
  V3F kpPos;
  kpPos.x = kpPosXY;
  kpPos.y = kpPosXY;
  kpPos.z = 0.f;
  
  V3F kpVel;
  kpVel.x = kpVelXY;
  kpVel.y = kpVelXY;
  kpVel.z = 0.f;

  V3F cappedVelCmd;
  if (velCmd.mag() > maxSpeedXY) {
      cappedVelCmd = velCmd.norm() * maxSpeedXY;
  }
  else {
      cappedVelCmd = velCmd;
  }
  V3F posErr = posCmd - pos;
  V3F velErr = cappedVelCmd - vel;
  accelCmd = kpPos * posErr + kpVel * velErr + accelCmd;
  
  if (accelCmd.mag() > maxAccelXY) {
      accelCmd = accelCmd.norm() * maxAccelXY;
  }
  */
  /////////////////////////////// END STUDENT CODE ////////////////////////////

  return accelCmd;
}

// returns desired yaw rate
float QuadControl::YawControl(float yawCmd, float yaw)
{
  // Calculate a desired yaw rate to control yaw to yawCmd
  // INPUTS: 
  //   yawCmd: commanded yaw [rad]
  //   yaw: current yaw [rad]
  // OUTPUT:
  //   return a desired yaw rate [rad/s]
  // HINTS: 
  //  - use fmodf(foo,b) to unwrap a radian angle measure float foo to range [0,b]. 
  //  - use the yaw control gain parameter kpYaw

  float yawRateCmd=0;
  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////
  //yawRateCmd = kpYaw * (yawCmd - yaw);

  //float yaw_mod, yawCmd_mod, 
  float yaw_error = 0;
  float pi = 3.142;
  //yaw_mod = fmodf(yaw, 2 * pi);
  //yawCmd_mod = fmodf(yawCmd, 2 * pi);
   //yaw_error = yawCmd_mod - yaw_mod;
  
  /*
  yawCmd = fmodf(yawCmd, 2 * pi);
  yaw_error = yawCmd - yaw;

  if (yaw_error > pi)
  {
      yaw_error -= 2 * pi;
  }
  if (yaw_error < -pi)
  {
      yaw_error += 2 * pi;
  }
  yawRateCmd = kpYaw * yaw_error;
  */
  
  /////////////////////////////// END STUDENT CODE ////////////////////////////

  return yawRateCmd;

}

VehicleCommand QuadControl::RunControl(float dt, float simTime)
{
  curTrajPoint = GetNextTrajectoryPoint(simTime);

  float collThrustCmd = AltitudeControl(curTrajPoint.position.z, curTrajPoint.velocity.z, estPos.z, estVel.z, estAtt, curTrajPoint.accel.z, dt);

  // reserve some thrust margin for angle control
  float thrustMargin = .1f*(maxMotorThrust - minMotorThrust);
  collThrustCmd = CONSTRAIN(collThrustCmd, (minMotorThrust+ thrustMargin)*4.f, (maxMotorThrust-thrustMargin)*4.f);
  
  V3F desAcc = LateralPositionControl(curTrajPoint.position, curTrajPoint.velocity, estPos, estVel, curTrajPoint.accel);
  
  V3F desOmega = RollPitchControl(desAcc, estAtt, collThrustCmd);
  desOmega.z = YawControl(curTrajPoint.attitude.Yaw(), estAtt.Yaw());

  V3F desMoment = BodyRateControl(desOmega, estOmega);

  return GenerateMotorCommands(collThrustCmd, desMoment);
}
