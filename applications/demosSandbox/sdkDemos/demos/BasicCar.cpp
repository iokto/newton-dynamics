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

#include <toolbox_stdafx.h>
#include "SkyBox.h"
#include "NewtonDemos.h"
#include "PhysicsUtils.h"
#include "TargaToOpenGl.h"
#include "dSoundManager.h"
#include "HeightFieldPrimitive.h"
#include "DemoMesh.h"
#include "DemoCamera.h"
#include "DemoEntityManager.h"
#include "DebugDisplay.h"
#include "UserPlaneCollision.h"
#include "CustomVehicleControllerManager.h"


// some figures form the 2000 SRT Viper data sheet: http://www.vipercentral.com/specifications/
//the 2000 Vipers� 8.4-liter engine generates
// max speed: 164 miles per hours					= 73.0f meter per seconds		
// horse power: 450 hp @ 5,200 rpm                    
// peak torque: 490 lb-ft @ 3700 rpm				
// read line:   6000 rpm : 6200 rpm fuel cut-off
// Curb weight: 3460 lb								= 1560 kilogram

// gear Box:
// 1st Gear 2.66:1 
// 2nd Gear 1.78:1 
// 3rd Gear 1.30:1 
// 4th Gear 1.00:1 
// 5th Gear 0.74:1 
// 6th Gear 0.50:1 
// Reverse 2.90:1 

// vehicle definition for a Muscle car
//#define VIPER_IDLE_TORQUE					 40.0f
#define VIPER_IDLE_TORQUE					300.0f
#define VIPER_IDLE_TORQUE_RPM				500.0f

#define VIPER_ENGINE_MOMENT_OF_INERTIA		10.0f

#define VIPER_PEAK_TORQUE					490.0f
#define VIPER_PEAK_TORQUE_RPM				3700.0f

#define VIPER_PEAK_HP						450.0f
#define VIPER_PEAK_HP_RPM					5200.0f


#define VIPER_REDLINE_TORQUE				30.0f
#define VIPER_REDLINE_TORQUE_RPM			6000.0f

#define VIPER_MASS							1500.0f
#define VIPER_TIRE_STEER_ANGLE				35.0f

// note: tire unsprung mass (note the engine impose a limit of 1.0 / 50.0 mass ratio between connected bodies 
//#define VIPER_TIRE_MASS					(VIPER_MASS / 50.0f)  
#define VIPER_TIRE_MASS						40.0f  

//#define VIPER_TIRE_TOP_SPEED				164 mile / hours
#define VIPER_TIRE_TOP_SPEED_KMH			264.0f			 

#define VIPER_TIRE_LATERAL_STIFFNESS		20.0f
#define VIPER_TIRE_LONGITUDINAL_STIFFNESS	1000.0f
#define VIPER_TIRE_ALIGNING_MOMENT_TRAIL	0.5f
#define VIPER_TIRE_SUSPENSION_SPRING		15000.0f
#define VIPER_TIRE_SUSPENSION_DAMPER		600.0f
#define VIPER_TIRE_SUSPENSION_LENGTH		0.20f
#define VIPER_TIRE_BRAKE_TORQUE				2000.0f

#define VIPER_TIRE_GEAR_1					2.66f
#define VIPER_TIRE_GEAR_2					1.78f
#define VIPER_TIRE_GEAR_3 					1.30f
#define VIPER_TIRE_GEAR_4 					1.00f
#define VIPER_TIRE_GEAR_5 					0.74f
#define VIPER_TIRE_GEAR_6 					0.50f
#define VIPER_TIRE_GEAR_REVERSE				2.90f
#define VIPER_COM_Y_OFFSET					-0.30f

#define VEHICLE_THIRD_PERSON_VIEW_HIGHT		2.0f
#define VEHICLE_THIRD_PERSON_VIEW_DIST		10.0f
#define VEHICLE_THIRD_PERSON_VIEW_FILTER	0.125f

class BasicVehicleEntity: public DemoEntity
{
	public: 
	class TireAligmentTransform: public UserData
	{
		public: 
		TireAligmentTransform (const dMatrix& matrix)
			:UserData()
			,m_matrix( matrix)
		{
		}

		virtual void OnRender (dFloat timestep) const{};
		virtual void OnInterpolateMatrix (DemoEntityManager& world, dFloat param) const{};

		dMatrix m_matrix;
	};


	BasicVehicleEntity (DemoEntityManager* const scene, CustomVehicleControllerManager* const manager, const dMatrix& location, const char* const filename)
		:DemoEntity (GetIdentityMatrix(), NULL)
		,m_controller(NULL)
		,m_helpKey (true)
		,m_gearUpKey (false)
		,m_gearDownKey (false)
		,m_engineKeySwitch(false)
		,m_automaticTransmission(true)
		,m_engineKeySwitchCounter(0)
		,m_engineOldKeyState(false)
		,m_engineRPMOn(false)
	{
		// add this entity to the scene for rendering
		scene->Append(this);

		// load the visual mesh 
		LoadNGD_mesh (filename, scene->GetNewton());

		// place entity in the world
		ResetMatrix (*scene, location);

		NewtonWorld* const world = scene->GetNewton();

		// create the car rigid body
		//NewtonBody* const body = CreateChassisBody (carEntity, world, materialID, params);
		NewtonCollision* const chassisCollision = CreateChassisCollision (world);

		// create the vehicle controller
		dMatrix chassisMatrix;
		chassisMatrix.m_front = dVector (1.0f, 0.0f, 0.0f, 0.0f);			// this is the vehicle direction of travel
		chassisMatrix.m_up	  = dVector (0.0f, 1.0f, 0.0f, 0.0f);			// this is the downward vehicle direction
		chassisMatrix.m_right = chassisMatrix.m_front * chassisMatrix.m_up;	// this is in the side vehicle direction (the plane of the wheels)
		chassisMatrix.m_posit = dVector (0.0f, 0.0f, 0.0f, 1.0f);

		// create a default vehicle 
		m_controller = manager->CreateVehicle (chassisCollision, chassisMatrix, VIPER_MASS, dVector (0.0f, DEMO_GRAVITY, 0.0f, 0.0f));

		// get body from player
		NewtonBody* const body = m_controller->GetBody();

		// set the user data
		NewtonBodySetUserData(body, this);

		// set the transform callback
		NewtonBodySetTransformCallback (body, DemoEntity::TransformCallback);

		// set the standard force and torque call back
		NewtonBodySetForceAndTorqueCallback(body, PhysicsApplyGravityForce);

		// set the player matrix 
		NewtonBodySetMatrix(body, &location[0][0]);

		// destroy the collision helper shape 
		NewtonDestroyCollision(chassisCollision);

		for (int i = 0; i < int ((sizeof (m_gearMap) / sizeof (m_gearMap[0]))); i ++) {
			m_gearMap[i] = i;
		}
		m_gearMap[0] = 1;
		m_gearMap[1] = 0;
	}

	~BasicVehicleEntity ()
	{
		((CustomVehicleControllerManager*)m_controller->GetManager())->DestroyController(m_controller);
	}

	NewtonCollision* CreateChassisCollision (NewtonWorld* const world) const
	{
		DemoEntity* const chassis = dHierarchy<DemoEntity>::Find("car_body");
		dAssert (chassis);

		DemoMesh* const mesh = chassis->GetMesh();
		//dAssert (chassis->GetMeshMatrix().TestIdentity());
		const dMatrix& meshMatrix = chassis->GetMeshMatrix();
		float* const vertexList = mesh->m_vertex;
		return NewtonCreateConvexHull(world, mesh->m_vertexCount, vertexList, 3 * sizeof (float), 0.001f, 0, &meshMatrix[0][0]);
	}

	// interpolate all skeleton transform 
	virtual void InterpolateMatrix (DemoEntityManager& world, dFloat param)
	{
		DemoEntity::InterpolateMatrix (world, param);
		if (m_controller){
			for (CustomVehicleControllerBodyStateTire* tire = m_controller->GetFirstTire(); tire; tire = m_controller->GetNextTire(tire)) {
				DemoEntity* const tirePart = (DemoEntity*) tire->GetUserData();
				tirePart->InterpolateMatrix (world, param);
			}
		}
	}

	void CalculateTireDimensions (const char* const tireName, dFloat& width, dFloat& radius) const
	{
		NewtonBody* const body = m_controller->GetBody();
		NewtonWorld* const world = NewtonBodyGetWorld(body);

		// find the the tire visual mesh 
		DemoEntity* const entity = (DemoEntity*) NewtonBodyGetUserData(body);
		DemoEntity* const tirePart = entity->Find (tireName);
		dAssert (tirePart);

		// make a convex hull collision shape to assist in calculation of the tire shape size
		DemoMesh* const tireMesh = tirePart->GetMesh();
		//dAssert (tirePart->GetMeshMatrix().TestIdentity());
		const dMatrix& meshMatrix = tirePart->GetMeshMatrix();
		dVector* const temp = new dVector [tireMesh->m_vertexCount];
		meshMatrix.TransformTriplex (&temp[0].m_x, sizeof (dVector), tireMesh->m_vertex, 3 * sizeof (dFloat), tireMesh->m_vertexCount);
		NewtonCollision* const collision = NewtonCreateConvexHull(world, tireMesh->m_vertexCount, &temp[0].m_x, sizeof (dVector), 0, 0, NULL);
		delete[] temp;

		// get the location of this tire relative to the car chassis
		dMatrix tireLocalMatrix (entity->GetNextMatrix());
		dMatrix tireMatrix (tirePart->CalculateGlobalMatrix() * tireLocalMatrix);

		// find the support points that can be used to define the with and high of the tire collision mesh
		dVector extremePoint;
		dVector upDir (tireMatrix.UnrotateVector(dVector (0.0f, 1.0f, 0.0f, 0.0f)));
		NewtonCollisionSupportVertex (collision, &upDir[0], &extremePoint[0]);
		radius = dAbs (upDir % extremePoint);

		//dVector widthDir (tireMatrix.UnrotateVector(tireLocalMatrix.m_right));
		dVector widthDir (tireMatrix.UnrotateVector(tireLocalMatrix.m_front));
		NewtonCollisionSupportVertex (collision, &widthDir[0], &extremePoint[0]);
		width = widthDir % extremePoint;

		widthDir = widthDir.Scale (-1.0f);
		NewtonCollisionSupportVertex (collision, &widthDir[0], &extremePoint[0]);
		width += widthDir % extremePoint;

		// destroy the auxiliary collision
		NewtonDestroyCollision (collision);	
	}


	CustomVehicleControllerBodyStateTire* AddTire (const char* const tireName, const dVector& offset, dFloat width, dFloat radius, dFloat mass, dFloat suspensionLength, dFloat suspensionSpring, dFloat suspensionDamper, dFloat lateralStiffness, dFloat longitudinalStiffness, dFloat aligningMOmentTrail) 
	{
		NewtonBody* const body = m_controller->GetBody();
		DemoEntity* const entity = (DemoEntity*) NewtonBodyGetUserData(body);
		DemoEntity* const tirePart = entity->Find (tireName);

		// add this tire, get local position and rise it by the suspension length 
		dMatrix tireMatrix (tirePart->CalculateGlobalMatrix(entity));

		// add the offset location
		tireMatrix.m_posit += offset;

		//lower the tire base position by some distance
		tireMatrix.m_posit.m_y -= suspensionLength * 0.25f;

		// add and alignment matrix,to match visual mesh to physics collision shape
		dMatrix aligmentMatrix (GetIdentityMatrix());
		if (tireMatrix[0][2] < 0.0f) {
			aligmentMatrix = dYawMatrix(3.141592f);
		}
		TireAligmentTransform* const m_ligmentMatrix = new TireAligmentTransform(aligmentMatrix);
		tirePart->SetUserData(m_ligmentMatrix);

		// add the tire to the vehicle
		CustomVehicleControllerBodyStateTire::TireCreationInfo tireInfo;
		tireInfo.m_location = tireMatrix.m_posit;
		tireInfo.m_mass = mass;
		tireInfo.m_radio = radius;
		tireInfo.m_width = width;
		tireInfo.m_dampingRatio = suspensionDamper;
		tireInfo.m_springStrength = suspensionSpring;
		tireInfo.m_suspesionlenght = suspensionLength;
		tireInfo.m_lateralStiffness = lateralStiffness;
		tireInfo.m_longitudialStiffness = longitudinalStiffness;
		tireInfo.m_aligningMomentTrail =  aligningMOmentTrail;
		tireInfo.m_userData = tirePart;

		return m_controller->AddTire (tireInfo);
	}


	void UpdateTireTransforms()
	{
		NewtonBody* const body = m_controller->GetBody();
		DemoEntityManager* const scene = (DemoEntityManager*) NewtonWorldGetUserData(NewtonBodyGetWorld(body));
		
#if 0
		// this is the general way for getting the tire matrices
		dMatrix rootMatrixInv (GetNextMatrix().Inverse());
		for (CustomVehicleControllerBodyStateTire* tire = m_controller->GetFirstTire(); tire; tire = m_controller->GetNextTire(tire)) {
			DemoEntity* const tirePart = (DemoEntity*) tire->GetUserData();
			TireAligmentTransform* const aligmentMatrix = (TireAligmentTransform*)tirePart->GetUserData();
			dMatrix matrix (aligmentMatrix->m_matrix * tire->CalculateGlobalMatrix() * rootMatrixInv);
			dQuaternion rot (matrix);
			tirePart->SetMatrix(*scene, rot, matrix.m_posit);
		}
#else
		// this saves some calculation since it get the tire local to the chassis
		for (CustomVehicleControllerBodyStateTire* tire = m_controller->GetFirstTire(); tire; tire = m_controller->GetNextTire(tire)) {
			DemoEntity* const tirePart = (DemoEntity*) tire->GetUserData();
			TireAligmentTransform* const aligmentMatrix = (TireAligmentTransform*)tirePart->GetUserData();
			dMatrix matrix (aligmentMatrix->m_matrix * tire->CalculateLocalMatrix());
			dQuaternion rot (matrix);
			tirePart->SetMatrix(*scene, rot, matrix.m_posit);
		}
#endif
	}


	// this function is an example of how to make a high performance super car
	void BuildMuscleCar ()
	{
		// step one: find the location of each tire, in the visual mesh and add them one by one to the vehicle controller 
		dFloat width;
		dFloat radius;

		// Muscle cars have the front engine, we need to shift the center of mass to the front to represent that
		m_controller->SetCenterOfGravity (dVector (0.35f, VIPER_COM_Y_OFFSET, 0.0f, 0.0f)); 

		// a car may have different size front an rear tire, therefore we do this separate for front and rear tires
		CalculateTireDimensions ("fl_tire", width, radius);
		dVector offset (0.0f, 0.0f, 0.0f, 0.0f);
		CustomVehicleControllerBodyStateTire* const leftFrontTire = AddTire ("fl_tire", offset, width, radius, VIPER_TIRE_MASS, VIPER_TIRE_SUSPENSION_LENGTH, VIPER_TIRE_SUSPENSION_SPRING, VIPER_TIRE_SUSPENSION_DAMPER, VIPER_TIRE_LATERAL_STIFFNESS, VIPER_TIRE_LONGITUDINAL_STIFFNESS, VIPER_TIRE_ALIGNING_MOMENT_TRAIL);
		CustomVehicleControllerBodyStateTire* const rightFrontTire = AddTire ("fr_tire", offset, width, radius, VIPER_TIRE_MASS, VIPER_TIRE_SUSPENSION_LENGTH, VIPER_TIRE_SUSPENSION_SPRING, VIPER_TIRE_SUSPENSION_DAMPER, VIPER_TIRE_LATERAL_STIFFNESS, VIPER_TIRE_LONGITUDINAL_STIFFNESS, VIPER_TIRE_ALIGNING_MOMENT_TRAIL);

		// add real tires
		CalculateTireDimensions ("rl_tire", width, radius);
		dVector offset1 (0.0f, 0.05f, 0.0f, 0.0f);
		CustomVehicleControllerBodyStateTire* const leftRearTire = AddTire ("rl_tire", offset1, width, radius, VIPER_TIRE_MASS, VIPER_TIRE_SUSPENSION_LENGTH, VIPER_TIRE_SUSPENSION_SPRING, VIPER_TIRE_SUSPENSION_DAMPER, VIPER_TIRE_LATERAL_STIFFNESS, VIPER_TIRE_LONGITUDINAL_STIFFNESS, VIPER_TIRE_ALIGNING_MOMENT_TRAIL);
		CustomVehicleControllerBodyStateTire* const rightRearTire = AddTire ("rr_tire", offset1, width, radius, VIPER_TIRE_MASS, VIPER_TIRE_SUSPENSION_LENGTH, VIPER_TIRE_SUSPENSION_SPRING, VIPER_TIRE_SUSPENSION_DAMPER, VIPER_TIRE_LATERAL_STIFFNESS, VIPER_TIRE_LONGITUDINAL_STIFFNESS, VIPER_TIRE_ALIGNING_MOMENT_TRAIL);

		// add an steering Wheel
		CustomVehicleControllerComponentSteering* const steering = new CustomVehicleControllerComponentSteering (m_controller, VIPER_TIRE_STEER_ANGLE * 3.141592f / 180.0f);
		steering->AddSteeringTire(leftFrontTire, -1.0f);
		steering->AddSteeringTire(rightFrontTire, -1.0f);
		m_controller->SetSteering(steering);

		// add vehicle brakes
		CustomVehicleControllerComponentBrake* const brakes = new CustomVehicleControllerComponentBrake (m_controller, VIPER_TIRE_BRAKE_TORQUE * 0.5f);
		brakes->AddBrakeTire (leftFrontTire);
		brakes->AddBrakeTire (rightFrontTire);
		brakes->AddBrakeTire (leftRearTire);
		brakes->AddBrakeTire (rightRearTire);
		m_controller->SetBrakes(brakes);

		// add vehicle hand brakes
		CustomVehicleControllerComponentBrake* const handBrakes = new CustomVehicleControllerComponentBrake (m_controller, VIPER_TIRE_BRAKE_TORQUE);
		handBrakes->AddBrakeTire (leftRearTire);
		handBrakes->AddBrakeTire (rightRearTire);
		m_controller->SetHandBrakes(handBrakes);

		// add an engine
		// first make the gear Box
		dFloat fowardSpeedGearsBoxRatios[] = {VIPER_TIRE_GEAR_1, VIPER_TIRE_GEAR_2, VIPER_TIRE_GEAR_3, VIPER_TIRE_GEAR_4, VIPER_TIRE_GEAR_5, VIPER_TIRE_GEAR_6};
		CustomVehicleControllerComponentEngine::dGearBox* const gearBox = new CustomVehicleControllerComponentEngine::dGearBox(m_controller, VIPER_TIRE_GEAR_REVERSE, sizeof (fowardSpeedGearsBoxRatios) / sizeof (fowardSpeedGearsBoxRatios[0]), fowardSpeedGearsBoxRatios); 
		CustomVehicleControllerComponentEngine* const engine = new CustomVehicleControllerComponentEngine (m_controller, gearBox, leftRearTire, rightRearTire);

		// the the default transmission type
		engine->SetTransmissionMode(m_automaticTransmission.GetPushButtonState());

		m_controller->SetEngine(engine);


		dFloat viperIdleRPM = VIPER_IDLE_TORQUE_RPM;
		dFloat viperIdleTorquePoundPerFoot = VIPER_IDLE_TORQUE;

		dFloat viperPeakTorqueRPM = VIPER_PEAK_TORQUE_RPM;
		dFloat viperPeakTorquePoundPerFoot = VIPER_PEAK_TORQUE;

		dFloat viperPeakHorsePowerRPM = VIPER_PEAK_HP_RPM;
		dFloat viperPeakHorsePower = VIPER_PEAK_HP;
		
		dFloat viperRedLineRPM = VIPER_REDLINE_TORQUE_RPM;
		dFloat viperRedLineTorquePoundPerFoot = VIPER_REDLINE_TORQUE;

		dFloat vehicleTopSpeedKPH = VIPER_TIRE_TOP_SPEED_KMH;
        dFloat vehicleMomentOfInteria  = VIPER_ENGINE_MOMENT_OF_INERTIA;
		engine->InitEngineTorqueCurve (vehicleTopSpeedKPH, vehicleMomentOfInteria, viperIdleTorquePoundPerFoot, viperIdleRPM, viperPeakTorquePoundPerFoot, viperPeakTorqueRPM, viperPeakHorsePower, viperPeakHorsePowerRPM, viperRedLineTorquePoundPerFoot, viperRedLineRPM);

		// do not forget to call finalize after all components are added or after any change is made to the vehicle
		m_controller->Finalize();
	}


	void BuildRaceCar ()
	{
/*
		// step one: find the location of each tire, in the visual mesh and add them one by one to the vehicle controller 
		dFloat width;
		dFloat radius;

		// Muscle cars have the front engine, we need to shift the center of mass to the front to represent that
		m_controller->SetCenterOfGravity (dVector (0.0f, 0.0f, 0.0f, 0.0f)); 

		// a car may have different size front an rear tire, therefore we do this separate for front and rear tires
		CalculateTireDimensions ("fl_tire", width, radius);
		dVector offset (0.0f, 0.0f, 0.0f, 0.0f);
		CustomVehicleControllerBodyStateTire* const leftFrontTireHandle = AddTire ("fl_tire", offset, width, radius, VIPER_TIRE_MASS, VIPER_TIRE_SUSPENSION_LENGTH, VIPER_TIRE_SUSPENSION_SPRING, VIPER_TIRE_SUSPENSION_DAMPER);
		CustomVehicleControllerBodyStateTire* const rightFrontTireHandle = AddTire ("fr_tire", offset, width, radius, VIPER_TIRE_MASS, VIPER_TIRE_SUSPENSION_LENGTH, VIPER_TIRE_SUSPENSION_SPRING, VIPER_TIRE_SUSPENSION_DAMPER);

		// add real tires
		CalculateTireDimensions ("rl_tire", width, radius);
		CustomVehicleControllerBodyStateTire* const leftRearTireHandle = AddTire ("rl_tire", offset, width, radius, VIPER_TIRE_MASS, VIPER_TIRE_SUSPENSION_LENGTH, VIPER_TIRE_SUSPENSION_SPRING, VIPER_TIRE_SUSPENSION_DAMPER);
		CustomVehicleControllerBodyStateTire* const rightRearTireHandle = AddTire ("rr_tire", offset, width, radius, VIPER_TIRE_MASS, VIPER_TIRE_SUSPENSION_LENGTH, VIPER_TIRE_SUSPENSION_SPRING, VIPER_TIRE_SUSPENSION_DAMPER);

		// add an engine
		// first make the gear Box
		dFloat gearsSpeed[] = {VIPER_TIRE_GEAR_1, VIPER_TIRE_GEAR_2, VIPER_TIRE_GEAR_3, VIPER_TIRE_GEAR_4, VIPER_TIRE_GEAR_5, VIPER_TIRE_GEAR_6};
		//CustomVehicleController::EngineComponent::GearBox* const gearBox = new CustomVehicleController::EngineComponent::GearBox(m_controller, VIPER_TIRE_GEAR_REVERSE, sizeof (gearsSpeed) / sizeof (gearsSpeed[0]), gearsSpeed); 
		CustomVehicleController::EngineComponent::GearBox* const gearBox = new CustomVehicleController::EngineComponent::GearBox(m_controller, VIPER_TIRE_GEAR_REVERSE, 2, gearsSpeed); 

		CustomVehicleController::EngineComponent* const engine = new CustomVehicleController::EngineComponent (m_controller, gearBox, leftRearTireHandle, rightRearTireHandle);

		dFloat viperIdleRPM = VIPER_IDLE_TORQUE_RPM;
		dFloat viperIdleTorquePoundPerFoot = VIPER_IDLE_TORQUE;

		dFloat viperPeakTorqueRPM = VIPER_PEAK_TORQUE_RPM;
		dFloat viperPeakTorquePoundPerFoot = VIPER_PEAK_TORQUE;

		dFloat viperPeakHorsePowerRPM = VIPER_PEAK_HP_RPM;
		dFloat viperPeakHorsePower = VIPER_PEAK_HP;

		dFloat viperRedLineRPM = VIPER_REDLINE_TORQUE_RPM;
		dFloat viperRedLineTorquePoundPerFoot = VIPER_REDLINE_TORQUE;

		dFloat vehicleTopSpeedKPH = VIPER_TIRE_TOP_SPEED_KMH;
        dFloat vehicleMomentOfInteria = VIPER_ENGINE_MOMENT_OF_INERTIA;
		engine->InitEngineTorqueCurve (vehicleTopSpeedKPH, vehicleMomentOfInteria, viperIdleTorquePoundPerFoot, viperIdleRPM, viperPeakTorquePoundPerFoot, viperPeakTorqueRPM, viperPeakHorsePower, viperPeakHorsePowerRPM, viperRedLineTorquePoundPerFoot, viperRedLineRPM);
		m_controller->SetEngine(engine);


		// add an steering Wheel
		CustomVehicleController::SteeringComponent* const steering = new CustomVehicleController::SteeringComponent (m_controller, VIPER_TIRE_STEER_ANGLE * 3.141592f / 180.0f);
		steering->AddSteeringTire(leftFrontTireHandle, -1.0f);
		steering->AddSteeringTire(rightFrontTireHandle, -1.0f);
		m_controller->SetSteering(steering);


		// add vehicle brakes
		CustomVehicleController::BrakeComponent* const brakes = new CustomVehicleController::BrakeComponent (m_controller, VIPER_TIRE_BRAKE_TORQUE);
		brakes->AddBrakeTire (leftFrontTireHandle);
		brakes->AddBrakeTire (rightFrontTireHandle);
		brakes->AddBrakeTire (leftRearTireHandle);
		brakes->AddBrakeTire (rightRearTireHandle);
		m_controller->SetBrakes(brakes);


		// add vehicle hand brakes
		CustomVehicleController::BrakeComponent* const handBrakes = new CustomVehicleController::BrakeComponent (m_controller, VIPER_TIRE_BRAKE_TORQUE);
		handBrakes->AddBrakeTire (leftRearTireHandle);
		handBrakes->AddBrakeTire (rightRearTireHandle);
		m_controller->SetHandBrakes(handBrakes);
*/
	}


	void ApplyPlayerControl ()
	{
		NewtonBody* const body = m_controller->GetBody();
		NewtonWorld* const world = NewtonBodyGetWorld(body);
		DemoEntityManager* const scene = (DemoEntityManager*) NewtonWorldGetUserData(world);
		NewtonDemos* const mainWindow = scene->GetRootWindow();

		CustomVehicleControllerComponentEngine* const engine = m_controller->GetEngine();
		CustomVehicleControllerComponentSteering* const steering = m_controller->GetSteering();
		CustomVehicleControllerComponentBrake* const brakes = m_controller->GetBrakes();
		CustomVehicleControllerComponentBrake* const handBrakes = m_controller->GetHandBrakes();

		// get the throttler input
		dFloat joyPosX;
		dFloat joyPosY;
		int joyButtons;

		int gear = engine ? engine->GetGear() : CustomVehicleControllerComponentEngine::dGearBox::m_newtralGear;
		dFloat steeringVal = 0.0f;
		dFloat engineGasPedal = 0.0f;
		dFloat brakePedal = 0.0f;
		dFloat handBrakePedal = 0.0f;
	
		bool hasJopytick = mainWindow->GetJoytickPosition (joyPosX, joyPosY, joyButtons);
		if (hasJopytick) {
			dAssert (0);
/*
			// apply a cubic attenuation to the joystick inputs
			joyPosX = joyPosX * joyPosX * joyPosX;
			joyPosY = joyPosY * joyPosY * joyPosY;

			steeringVal = joyPosX;
			brakePedal = (joyPosY < 0.0f) ? -joyPosY: 0.0f;
			engineGasPedal = (joyPosY >= 0.0f) ? joyPosY: 0.0f;

			gear += int (m_gearUpKey.UpdateTriggerJoystick(mainWindow, joyButtons & 2)) - int (m_gearDownKey.UpdateTriggerJoystick(mainWindow, joyButtons & 4));
			handBrakePedal = (joyButtons & 1) ? 1.0f : 0.0f;
*/			
		} else {

	
			// get keyboard controls
			if (mainWindow->GetKeyState ('W')) {
				engineGasPedal = 1.0f;
			}

			// see if HandBreale is on
			if (mainWindow->GetKeyState ('S')) {
				brakePedal = 1.0f;
			}

			// see if HandBreale is on
			if (mainWindow->GetKeyState (' ')) {
				handBrakePedal = 1.0f;
			}

			// get the steering input
			steeringVal = (dFloat (mainWindow->GetKeyState ('D')) - dFloat (mainWindow->GetKeyState ('A')));

/*
			// check for gear change (note key code for '>' = '.' and key code for '<' == ',')
			gear += int (m_gearUpKey.UpdateTriggerButton(mainWindow, '.')) - int (m_gearDownKey.UpdateTriggerButton(mainWindow, ','));
			// do driving heuristic for automatic transmission
			if (engine->GetTransmissionMode()) {
				dFloat speed = engine->GetSpeed();
				// check if vehicle is parked
				if ((dAbs (speed) < 1.0f) && !engineGasPedal && !brakePedal && !handBrakePedal) {
					handBrakePedal = 0.5f;
				}
			}
*/
		}
				
		// set the help key
		m_helpKey.UpdatePushButton (mainWindow, 'H');

		// count the engine key switch
		m_engineKeySwitchCounter += m_engineKeySwitch.UpdateTriggerButton(mainWindow, 'Y');
		

		// check transmission type
		int toggleTransmission = m_automaticTransmission.UpdateTriggerButton (mainWindow, 0x0d) ? 1 : 0;

#if 0
	#if 0
		static FILE* file = fopen ("log.bin", "wb");
		if (file) {
			fwrite (&m_engineKeySwitchCounter, sizeof (int), 1, file);
			fwrite (&toggleTransmission, sizeof (int), 1, file);
			fwrite (&gear, sizeof (int), 1, file);
			fwrite (&steeringVal, sizeof (dFloat), 1, file);
			fwrite (&engineGasPedal, sizeof (dFloat), 1, file);
			fwrite (&handBrakePedal, sizeof (dFloat), 1, file);
			fwrite (&brakePedal, sizeof (dFloat), 1, file);
			fflush(file);
		}
	#else 
		static FILE* file = fopen ("log.bin", "rb");
		if (file) {		
			fread (&m_engineKeySwitchCounter, sizeof (int), 1, file);
			fread (&toggleTransmission, sizeof (int), 1, file);
			fread (&gear, sizeof (int), 1, file);
			fread (&steeringVal, sizeof (dFloat), 1, file);
			fread (&engineGasPedal, sizeof (dFloat), 1, file);
			fread (&handBrakePedal, sizeof (dFloat), 1, file);
			fread (&brakePedal, sizeof (dFloat), 1, file);
		}
	#endif
#endif

		if (engine) {
			m_engineOldKeyState = engine->GetKey();
			engine->SetKey ((m_engineKeySwitchCounter & 1) ? true : false);

			if (toggleTransmission) {
				engine->SetTransmissionMode (!engine->GetTransmissionMode());
			}
			if (!engine->GetTransmissionMode()) {
				engine->SetGear(gear);
			}
			engine->SetParam(engineGasPedal);
		}
		if (steering) {
			steering->SetParam(steeringVal);
		}
		if (brakes) {
			brakes->SetParam(brakePedal);
		}
		if (handBrakes) {
			handBrakes->SetParam(handBrakePedal);
		}
		
	}

	void ApplyNPCControl ()
	{
		dAssert (0);
	}


	void Debug () const 
	{
		NewtonBody* const body = m_controller->GetBody();
		const CustomVehicleControllerBodyStateChassis& chassis = m_controller->GetChassisState ();

		dFloat scale = -4.0f / (chassis.GetMass() * DEMO_GRAVITY);
		dVector p0 (chassis.GetCenterOfMass());

		glDisable (GL_LIGHTING);
		glDisable(GL_TEXTURE_2D);

		glLineWidth(3.0f);
		glBegin(GL_LINES);


		// draw vehicle weight at the center of mass
		dFloat lenght = scale * chassis.GetMass() * DEMO_GRAVITY;
		glColor3f(0.0f, 0.0f, 1.0f);
		glVertex3f (p0.m_x, p0.m_y, p0.m_z);
		glVertex3f (p0.m_x, p0.m_y - lenght, p0.m_z);

		// draw vehicle front dir
		glColor3f(1.0f, 1.0f, 1.0f);
		dVector r0 (p0 + chassis.GetMatrix()[1].Scale (0.5f));
		dVector r1 (r0 + chassis.GetMatrix()[0].Scale (1.0f));
		glVertex3f (r0.m_x, r0.m_y, r0.m_z);
		glVertex3f (r1.m_x, r1.m_y, r1.m_z);


		// draw the velocity vector, a little higher so that is not hidden by the vehicle mesh 
		dVector veloc;
		NewtonBodyGetVelocity(body, &veloc[0]);
		dVector q0 (p0 + chassis.GetMatrix()[1].Scale (1.0f));
		dVector q1 (q0 + veloc.Scale (0.25f));
		glColor3f(1.0f, 1.0f, 0.0f);
		glVertex3f (q0.m_x, q0.m_y, q0.m_z);
		glVertex3f (q1.m_x, q1.m_y, q1.m_z);
		
		for (CustomVehicleControllerBodyStateTire* node = m_controller->GetFirstTire(); node; node = m_controller->GetNextTire(node)) {
			const CustomVehicleControllerBodyStateTire& tire = *node;
			dVector p0 (tire.GetCenterOfMass());

const dMatrix& tireMatrix = tire.GetLocalMatrix ();
p0 += chassis.GetMatrix()[2].Scale ((tireMatrix.m_posit.m_z > 0.0f ? 1.0f : -1.0f) * 0.5f);

			// draw the tire load 
			dVector p1 (p0 + tire.GetTireLoad().Scale (scale));
			glColor3f (0.0f, 0.0f, 1.0f);
			glVertex3f (p0.m_x, p0.m_y, p0.m_z);
			glVertex3f (p1.m_x, p1.m_y, p1.m_z);

			// show tire lateral force
			dVector p2 (p0 - tire.GetLateralForce().Scale (scale));
			glColor3f(1.0f, 0.0f, 0.0f);
			glVertex3f (p0.m_x, p0.m_y, p0.m_z);
			glVertex3f (p2.m_x, p2.m_y, p2.m_z);

			// show tire longitudinal force
			dVector p3 (p0 - tire.GetLongitudinalForce().Scale (scale));
			glVertex3f (p0.m_x, p0.m_y, p0.m_z);
			glVertex3f (p3.m_x, p3.m_y, p3.m_z);
		}

		glEnd();

		glLineWidth(1.0f);

	}


	CustomVehicleController* m_controller;
	DemoEntityManager::ButtonKey m_helpKey;
	DemoEntityManager::ButtonKey m_gearUpKey;
	DemoEntityManager::ButtonKey m_gearDownKey;
	DemoEntityManager::ButtonKey m_engineKeySwitch;
	DemoEntityManager::ButtonKey m_automaticTransmission;
	int m_gearMap[10];
	int m_engineKeySwitchCounter;
	bool m_engineOldKeyState;
	bool m_engineRPMOn;
};


class BasicVehicleControllerManager: public CustomVehicleControllerManager
{
	public:

	BasicVehicleControllerManager (NewtonWorld* const world)
		:CustomVehicleControllerManager (world)
		,m_externalView(true)
		,m_player (NULL) 
		,m_soundsCount(0)
	{
		// hook a callback for 2d help display
		DemoEntityManager* const scene = (DemoEntityManager*) NewtonWorldGetUserData(world);
		scene->Set2DDisplayRenderFunction (RenderVehicleHud, this);

		// load 2d display assets
		m_gears = LoadTexture ("gears_font.tga");
		m_tachometer = LoadTexture ("rpm_dial.tga");
		m_odometer = LoadTexture ("kmh_dial.tga");
		m_redNeedle = LoadTexture ("needle_red.tga");
		m_greenNeedle = LoadTexture ("needle_green.tga");

		// create the vehicle sound 
		const char* engineSounds[] = {"starter.wav", "tire_skid.wav", "revLow.wav", "revMiddle.wav", "revHigh.wav"};
		dSoundManager* const soundManager = scene->GetSoundManager();
		for (int i = 0; i < int (sizeof (engineSounds) / sizeof (engineSounds[0])); i ++) {
			void* const sound = soundManager->CreateSound(engineSounds[i]);
			void* const channel = soundManager->CreatePlayChannel(sound);
			m_engineSounds[i] = channel;
			m_soundsCount ++;
		}
	}

	~BasicVehicleControllerManager ()
	{
		ReleaseTexture (m_gears);
		ReleaseTexture (m_tachometer);
		ReleaseTexture (m_odometer);
		ReleaseTexture (m_redNeedle);
		ReleaseTexture (m_greenNeedle);
	}
	

	static void RenderVehicleHud (DemoEntityManager* const scene, void* const context, int lineNumber)
	{
		BasicVehicleControllerManager* const me = (BasicVehicleControllerManager*) context;
		me->RenderVehicleHud (scene, lineNumber);
	}

	void DrawGage(GLuint gage, GLuint needle, dFloat param, dFloat origin_x, dFloat origin_y, dFloat size) const
	{
		size *= 0.5f;
		dMatrix origin (GetIdentityMatrix());
		origin.m_posit = dVector(origin_x, origin_y, 0.0f, 1.0f);

		// render dial
		glPushMatrix();
		glMultMatrix (&origin[0][0]);
		glBindTexture(GL_TEXTURE_2D, gage);
		glBegin(GL_QUADS);
		glTexCoord2f(0.0f, 1.0f); glVertex3f(-size,  size, 0.0f);
		glTexCoord2f(0.0f, 0.0f); glVertex3f(-size, -size, 0.0f);
		glTexCoord2f(1.0f, 0.0f); glVertex3f( size, -size, 0.0f);
		glTexCoord2f(1.0f, 1.0f); glVertex3f( size,  size, 0.0f);
		glEnd();

		// render needle
		const dFloat minAngle = 180.0f * 3.141592f / 180.0f;
		const dFloat maxAngle = -90.0f * 3.141592f / 180.0f;
		dFloat angle = minAngle + (maxAngle - minAngle) * param;
		dMatrix needleMatrix (dRollMatrix (angle));

		float x = size * 0.7f;
		float y = size * 0.7f;

		glPushMatrix();
		glMultMatrix (&needleMatrix[0][0]);
		glBindTexture(GL_TEXTURE_2D, needle);
		glBegin(GL_QUADS);
		glTexCoord2f(0.0f, 1.0f); glVertex3f(-x,  y, 0.0f);
		glTexCoord2f(0.0f, 0.0f); glVertex3f(-x, -y, 0.0f);
		glTexCoord2f(1.0f, 0.0f); glVertex3f( x, -y, 0.0f);
		glTexCoord2f(1.0f, 1.0f); glVertex3f( x,  y, 0.0f);
		glEnd();

		glPopMatrix();
		glPopMatrix();
	}

	void DrawGear(dFloat param, dFloat origin_x, dFloat origin_y, int gear, dFloat size) const
	{
		dMatrix origin (GetIdentityMatrix());
		origin.m_posit = dVector(origin_x + size * 0.3f, origin_y - size * 0.25f, 0.0f, 1.0f);

		glPushMatrix();
		glMultMatrix (&origin[0][0]);

		dFloat uwith = 0.1f;
		float u0 = uwith * gear;
		float u1 = u0 + uwith;

		dFloat x1 = 10.0f;
		dFloat y1 = 10.0f;
		glColor4f(1, 1, 0, 1);
		glBindTexture(GL_TEXTURE_2D, m_gears);
		glBegin(GL_QUADS);
		glTexCoord2f(u0, 1.0f); glVertex3f(-x1,  y1, 0.0f);
		glTexCoord2f(u0, 0.0f); glVertex3f(-x1, -y1, 0.0f);
		glTexCoord2f(u1, 0.0f); glVertex3f( x1, -y1, 0.0f);
		glTexCoord2f(u1, 1.0f); glVertex3f( x1,  y1, 0.0f);
		glEnd();

		glPopMatrix();
	}

	void DrawHelp(DemoEntityManager* const scene, int lineNumber) const
	{
		if (m_player->m_helpKey.GetPushButtonState()) {
			dVector color(1.0f, 1.0f, 0.0f, 0.0f);
			lineNumber = scene->Print (color, 10, lineNumber + 20, "Vehicle driving keyboard control:   Joystick control");
			lineNumber = scene->Print (color, 10, lineNumber + 20, "keySwich			: 'Y'           start engine");
			lineNumber = scene->Print (color, 10, lineNumber + 20, "accelerator         : 'W'           stick forward");
			lineNumber = scene->Print (color, 10, lineNumber + 20, "brakes              : 'S'           stick back");
			lineNumber = scene->Print (color, 10, lineNumber + 20, "turn right          : 'D'           stick right");
			lineNumber = scene->Print (color, 10, lineNumber + 20, "turn right          : 'S'           stick left");
			lineNumber = scene->Print (color, 10, lineNumber + 20, "gear up             : '>'           button 2");
			lineNumber = scene->Print (color, 10, lineNumber + 20, "gear down           : '<'           button 4");
			lineNumber = scene->Print (color, 10, lineNumber + 20, "manual transmission : enter         button 4");
			lineNumber = scene->Print (color, 10, lineNumber + 20, "hand brakes         : space         button 1");
			lineNumber = scene->Print (color, 10, lineNumber + 20, "hide help           : H");
		}
	}


	void RenderVehicleHud (DemoEntityManager* const scene, int lineNumber) const
	{
		// set to transparent color
		glEnable (GL_BLEND);
		glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        
		CustomVehicleControllerComponentEngine* const engine = m_player->m_controller->GetEngine();
		if (engine) {
			glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
			dFloat width = scene->GetWidth();
			//dFloat height = scene->getHeight();
			dFloat gageSize = 200.0f;

			//dFloat y = -height / 2.0f + gageSize / 2.0f + 60.0f;
			dFloat y = gageSize / 2.0f + 60.0f;

			// draw the tachometer
			dFloat x = gageSize / 2 + 80.0f;
			dFloat rpm = engine->GetRPM () / engine->GetRedLineRPM();
			DrawGage(m_tachometer, m_redNeedle, rpm, x, y, gageSize);

			// draw the odometer
			x = width - gageSize / 2 - 80.0f;
			dFloat speed = dAbs(engine->GetSpeed()) * 3.6f / 340.0f;
			DrawGage(m_odometer, m_greenNeedle, speed, x, y, gageSize);
			DrawGear(speed, x, y, m_player->m_gearMap[engine->GetGear()], gageSize);
		}

		//print controllers help 
		DrawHelp(scene, lineNumber);

		// restore color and blend mode
		glDisable (GL_BLEND);
	}


	void SetAsPlayer (BasicVehicleEntity* const player)
	{
		m_player = player;
	}


	virtual void PreUpdate (dFloat timestep)
	{
		// apply the vehicle controls, and all simulation time effect
		NewtonWorld* const world = GetWorld(); 
		DemoEntityManager* const scene = (DemoEntityManager*) NewtonWorldGetUserData(world);
		dSoundManager* const soundManager = scene->GetSoundManager();
		for (dListNode* ptr = GetFirst(); ptr; ptr = ptr->GetNext()) {
			CustomVehicleController* const controller = &ptr->GetInfo();
		
			NewtonBody* const body = controller->GetBody();
			BasicVehicleEntity* const vehicleEntity = (BasicVehicleEntity*) NewtonBodyGetUserData(body);

			if (vehicleEntity == m_player) {
				// do player control
				vehicleEntity->ApplyPlayerControl ();
			} else {
				// do no player control
				vehicleEntity->ApplyNPCControl ();
			}

			CustomVehicleControllerComponentEngine* const engine = vehicleEntity->m_controller->GetEngine();
			if (engine->GetKey()) {
				// see if the player started the engine
				void* const startEngine = m_engineSounds[0];

				if (!vehicleEntity->m_engineOldKeyState) {
					// play engine start sound;
					soundManager->PlayChannel(startEngine);

					for (int i = 2; i < m_soundsCount; i ++) {
						void* const rpmEngine = m_engineSounds[i];
						soundManager->PlayChannel(rpmEngine);
						soundManager->SetChannelPitch (rpmEngine, 1.0f);
						soundManager->SetChannelVolume(rpmEngine, 0.0f);
						soundManager->SetChannelLoopMode(rpmEngine, true);
					}
				} 

				// player need to check if the start engine sound is still on
				void* const startEngineSoundAsset = soundManager->GetAsset(startEngine);
				dFloat length = soundManager->GetSoundlength (startEngineSoundAsset);
				dFloat posit = soundManager->GetChannelSecPosition(startEngine);
				
				if (posit >= length * 0.5f) {
					// attenuate star engine volume 
					dFloat volume = soundManager->GetChannelVolume(startEngine);
					if (volume < 1.0e-3f) {
						volume = 0.0f;
					}
					soundManager->SetChannelVolume(startEngine, volume * 0.98f);
					vehicleEntity->m_engineRPMOn = true;
				}

				// update engine rpm sound for all vehicles
				if (vehicleEntity->m_engineRPMOn) {
					int count = m_soundsCount - 3;
					//dFloat rpm = engine->GetRPM () / engine->GetRedLineRPM();
					dFloat rpm = count * (dClamp (engine->GetRPM () / engine->GetRedLineRPM(), 0.0f, 0.999f));

					int index = dFloor(rpm);
					rpm = dClamp (rpm - dFloat(index), 0.0f, 1.0f);
					dAssert (index >= 0);
					dAssert (index < count);
					soundManager->SetChannelVolume(m_engineSounds[index + 2], 1.0f - rpm);
					soundManager->SetChannelPitch (m_engineSounds[index + 2], rpm + 1.0f);

					soundManager->SetChannelVolume(m_engineSounds[index + 3], rpm);
					soundManager->SetChannelPitch (m_engineSounds[index + 3], rpm * 0.5f + 0.5f);
				}
			} else {
				vehicleEntity->m_engineRPMOn = false;
				if (vehicleEntity->m_engineOldKeyState) {
					for (int i = 0; i < m_soundsCount; i ++) {
						soundManager->SetChannelVolume(m_engineSounds[i], 0.0f);
						soundManager->StopChannel(m_engineSounds[i]);
					}
				}
			}

		}

		// do the base class post update
		CustomVehicleControllerManager::PreUpdate(timestep);
	}


	virtual void PostUpdate (dFloat timestep)
	{
		// do the base class post update
		CustomVehicleControllerManager::PostUpdate(timestep);
		
		// update the visual transformation matrices for all vehicle tires
		for (dListNode* node = GetFirst(); node; node = node->GetNext()) {
			CustomVehicleController* const controller = &node->GetInfo();
			BasicVehicleEntity* const vehicleEntity = (BasicVehicleEntity*)NewtonBodyGetUserData (controller->GetBody());
			vehicleEntity->UpdateTireTransforms();
		}

		UpdateCamera (timestep);
	}


	void UpdateCamera (dFloat timestep)
	{
		DemoEntityManager* const scene = (DemoEntityManager*) NewtonWorldGetUserData(GetWorld());
		DemoCamera* const camera = scene->GetCamera();
		dMatrix camMatrix (camera->GetNextMatrix ());
		dMatrix playerMatrix (m_player->GetNextMatrix());

		dVector frontDir (camMatrix[0]);
		dVector camOrigin; 
		if (m_externalView) {
			camOrigin = playerMatrix.m_posit + dVector(0.0f, VEHICLE_THIRD_PERSON_VIEW_HIGHT, 0.0f, 0.0f);
			camOrigin -= frontDir.Scale(VEHICLE_THIRD_PERSON_VIEW_DIST);
		} else {
			dAssert (0);
			//            camMatrix = camMatrix * playerMatrix;
			//            camOrigin = playerMatrix.TransformVector(dVector(-0.8f, ARTICULATED_VEHICLE_CAMERA_EYEPOINT, 0.0f, 0.0f));
		}

		camera->SetNextMatrix (*scene, camMatrix, camOrigin);
	}

	// use this to display debug information about vehicle 
	void Debug () const
	{
		for (dListNode* ptr = GetFirst(); ptr; ptr = ptr->GetNext()) {
			CustomVehicleController* const controller = &ptr->GetInfo();
			BasicVehicleEntity* const vehicleEntity = (BasicVehicleEntity*)NewtonBodyGetUserData (controller->GetBody());
			vehicleEntity->Debug();
		}
	}

	bool m_externalView;
	BasicVehicleEntity* m_player;

	GLuint m_gears;
	GLuint m_redNeedle;
	GLuint m_greenNeedle;
	GLuint m_odometer;
	GLuint m_tachometer;
	GLuint m_needle;
	int m_soundsCount;
	void* m_engineSounds[16];
};



// *************************************************************************************************
// 
//  create a simple racing game with a simple controlled by a Newton AI engine and newton physics
//
// *************************************************************************************************
void BasicCar (DemoEntityManager* const scene)
{

	// load the sky box
	scene->CreateSkyBox();
	
	CreateLevelMesh (scene, "flatPlane.ngd", 1);
	//CreateLevelMesh (scene, "raceTrack2.ngd", 0);
	//CreateLevelMesh (scene, "raceTrack2.ngd", 1);
	//CreateLevelMesh (scene, "raceTrack1.ngd", 0);
	//CreateLevelMesh (scene, "raceTrack1.ngd", 1);
	//CreateHeightFieldTerrain (scene, 10, 8.0f, 1.5f, 0.2f, 200.0f, -50.0f);
	//CreatePlaneCollision (scene, dVector (0.0f, 1.0f, 0.0f, 0.0f));


	NewtonWorld* const world = scene->GetNewton();

	// create a vehicle controller manager
	BasicVehicleControllerManager* const manager = new BasicVehicleControllerManager (world);

	// create a simple vehicle
//	dMatrix location (dYawMatrix(90.0f * 3.1416f/ 180.0f));
//	location.m_posit.m_x = 126.0f;
//	location.m_posit.m_y = 50.0f;
//	location.m_posit.m_z = 50.0f;

dMatrix location (GetIdentityMatrix());
location.m_posit.m_y = 50.0f;
location.m_posit.m_x = -150.0f;

	location.m_posit = FindFloor (scene->GetNewton(), location.m_posit, 100.0f);
	location.m_posit.m_y += 0.5f;

	// make a vehicle entity shell
	//BasicVehicleEntity* const vehicle = new BasicVehicleEntity (scene, manager, location, "f1.ngd");
	BasicVehicleEntity* const vehicle = new BasicVehicleEntity (scene, manager, location, "viper.ngd");

	// build a muscle car from this vehicle controller
	vehicle->BuildMuscleCar();

	// set this vehicle as the player
	manager->SetAsPlayer(vehicle);

	// set the camera matrix, we only care the initial direction since it will be following the player vehicle
	dMatrix camMatrix (manager->m_player->GetNextMatrix());
	scene->SetCameraMouseLock (true);
	scene->SetCameraMatrix(camMatrix, camMatrix.m_posit);


//	int defaultMaterialID = NewtonMaterialGetDefaultGroupID (scene->GetNewton());
//	dVector location (origin);
//	location.m_x += 20.0f;
//	location.m_z += 20.0f;
//	dVector size (0.5f, 0.5f, 0.75f, 0.0f);

//	int count = 5;
//	AddPrimitiveArray(scene, 10.0f, location, size, count, count, 5.0f, _SPHERE_PRIMITIVE, defaultMaterialID, shapeOffsetMatrix);
//	AddPrimitiveArray(scene, 10.0f, location, size, count, count, 5.0f, _BOX_PRIMITIVE, defaultMaterialID, shapeOffsetMatrix);
//	AddPrimitiveArray(scene, 10.0f, location, size, count, count, 5.0f, _CAPSULE_PRIMITIVE, defaultMaterialID, shapeOffsetMatrix);
//	AddPrimitiveArray(scene, 10.0f, location, size, count, count, 5.0f, _CYLINDER_PRIMITIVE, defaultMaterialID, shapeOffsetMatrix);
//	AddPrimitiveArray(scene, 10.0f, location, size, count, count, 5.0f, _TAPERED_CAPSULE_PRIMITIVE, defaultMaterialID, shapeOffsetMatrix);
//	AddPrimitiveArray(scene, 10.0f, location, size, count, count, 5.0f, _TAPERED_CYLINDER_PRIMITIVE, defaultMaterialID, shapeOffsetMatrix);
//	AddPrimitiveArray(scene, 10.0f, location, size, count, count, 5.0f, _CHAMFER_CYLINDER_PRIMITIVE, defaultMaterialID, shapeOffsetMatrix);
//	AddPrimitiveArray(scene, 10.0f, location, size, count, count, 5.0f, _CONE_PRIMITIVE, defaultMaterialID, shapeOffsetMatrix);
//	AddPrimitiveArray(scene, 10.0f, location, size, count, count, 5.0f, _REGULAR_CONVEX_HULL_PRIMITIVE, defaultMaterialID, shapeOffsetMatrix);
//	AddPrimitiveArray(scene, 10.0f, location, size, count, count, 5.0f, _RANDOM_CONVEX_HULL_PRIMITIVE, defaultMaterialID, shapeOffsetMatrix);

//	NewtonSerializeToFile (scene->GetNewton(), "C:/Users/Julio/Desktop/newton-dynamics/applications/media/xxxxx.bin");

}

