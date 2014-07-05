/* Copyright (c) <2009> <Newton Game Dynamics>
* 
* This software is provided 'as-is', without any express or implied
* warranty. In no event will the authors be held liable for any damages
* arising from the use of this software.
* 
* Permission is granted to anyone to use this software for any purpose,
* including commercial applications, and to alter it and redistribute it
* freely
*/


// NewtonVehicleControllerManager.h: interface for the NewtonVehicleControllerManager class.
//
//////////////////////////////////////////////////////////////////////

#ifndef D_CUSTOM_VEHICLE_CONTROLLER_BODY_STATES_H_
#define D_CUSTOM_VEHICLE_CONTROLLER_BODY_STATES_H_

#include <CustomJointLibraryStdAfx.h>
#include <CustomAlloc.h>
#include <CustomVehicleControllerJoint.h>

class CustomVehicleController;

class CustomVehicleControllerBodyState
{
	public:
	CUSTOM_JOINTS_API CustomVehicleControllerBodyState();
	CUSTOM_JOINTS_API virtual ~CustomVehicleControllerBodyState() {}

	CUSTOM_JOINTS_API void UpdateInertia();
	CUSTOM_JOINTS_API void Init(CustomVehicleController* const controller);
	CUSTOM_JOINTS_API virtual void IntegrateForce (dFloat timestep, const dVector& force, const dVector& torque);
	CUSTOM_JOINTS_API virtual void CalculateAverageAcceleration (dFloat invTimestep, const dVector& veloc, const dVector& omega);

	dMatrix m_matrix;
	dMatrix m_localFrame;
	dMatrix m_inertia;
	dMatrix m_invInertia;

	dVector m_localInertia;
	dVector m_localInvInertia;

	dVector m_veloc;
	dVector m_omega;
	dVector m_externalForce;
	dVector m_externalTorque;
	dVector m_globalCentreOfMass;

	dFloat m_mass;
	dFloat m_invMass;
	int m_myIndex;
	CustomVehicleController* m_controller;
};


class CustomVehicleControllerBodyStateChassis: public CustomVehicleControllerBodyState
{
	public:
	CUSTOM_JOINTS_API void Init (CustomVehicleController* const controller, const dMatrix& localframe);
	CUSTOM_JOINTS_API void UpdateDynamicInputs();
	CUSTOM_JOINTS_API virtual void CalculateAverageAcceleration (dFloat invTimestep, const dVector& veloc, const dVector& omega);

	dVector m_com;
	dVector m_comOffset;
	dVector m_gravity;
};


class CustomVehicleControllerBodyStateEngine: public CustomVehicleControllerBodyState
{
	public:
	CUSTOM_JOINTS_API void Init (CustomVehicleController* const controller);

	CUSTOM_JOINTS_API void Update (dFloat timestep, CustomVehicleController* const controller);
	CUSTOM_JOINTS_API void CalculateAverageAcceleration (dFloat invTimestep, const dVector& veloc, const dVector& omega);
	CUSTOM_JOINTS_API int CalculateActiveJoints (CustomVehicleController* const controller, VehicleJoint** const jointArray);

	EngineGearJoint m_leftTire;
	EngineGearJoint m_rightTire;
	EngineIdleJoint	m_idleFriction;
	dFloat m_radianPerSecund;
};


class CustomVehicleControllerBodyStateTire: public CustomVehicleControllerBodyState
{
	public:
	class TireCreationInfo
	{
		public:
		dVector m_location;
		dFloat m_mass;
		dFloat m_radio;
		dFloat m_width;
		dFloat m_dampingRatio;
		dFloat m_springStrength;
		dFloat m_suspesionlenght;
		void* m_userData;
	};

	CUSTOM_JOINTS_API void Init (CustomVehicleController* const controller, const TireCreationInfo& tireInfo);


	CUSTOM_JOINTS_API dMatrix CalculateSteeringMatrix () const;
	CUSTOM_JOINTS_API dMatrix CalculateSuspensionMatrix () const;

	CUSTOM_JOINTS_API void UpdateDynamicInputs(dFloat timestep);
	CUSTOM_JOINTS_API void Collide (CustomControllerConvexCastPreFilter& filter, dFloat timestepInv);
/*
	CUSTOM_JOINTS_API dFloat GetAdhesionCoefficient() const;
	CUSTOM_JOINTS_API void SetAdhesionCoefficient(dFloat Coefficient);

	CUSTOM_JOINTS_API dMatrix CalculateMatrix () const;
	


	

	CUSTOM_JOINTS_API void UpdateTransform ();
	CUSTOM_JOINTS_API virtual void IntegrateForce (dFloat timestep, const dVector& force, const dVector& torque);
	CUSTOM_JOINTS_API virtual void CalculateAverageAcceleration (dFloat invTimestep, const dVector& veloc, const dVector& omega);
*/
	dVector m_tireLoad;
	dVector m_lateralForce;
	dVector m_longitidinalForce;
	dFloat m_radio;
	dFloat m_width;
	dFloat m_posit;
	dFloat m_speed;
	dFloat m_breakTorque;
	dFloat m_engineTorque;
	dFloat m_rotatonSpeed;
	dFloat m_rotationAngle;
	dFloat m_steeringAngle;
	dFloat m_dampingRatio;
	dFloat m_springStrength;
	dFloat m_suspensionlenght;
	dFloat m_adhesionCoefficient; 
	dFloat m_idleRollingResistance;
	//dFloat m_engineTorqueResistance;

	void* m_userData;
	NewtonCollision* m_shape;
	TireJoint m_chassisJoint;
	ContactJoint m_contactJoint;
};



#endif 

