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
#include "TargaToOpenGl.h"
#include "DemoMesh.h"
#include "DemoEntityManager.h"
#include "DemoCamera.h"
#include "PhysicsUtils.h"
#include "HeightFieldPrimitive.h"
#include "DebugDisplay.h"


struct BasciCarParameters
{
	enum DifferentialType
	{
		m_4WD,
		m_RWD,
		m_FWD,
	};

	dFloat MASS;
	dFloat TIRE_MASS;
	dFloat STEER_ANGLE;
	dFloat BRAKE_TORQUE;
	dFloat COM_Y_OFFSET;
	dFloat TIRE_TOP_SPEED_KMH;

	dFloat IDLE_TORQUE;
	dFloat IDLE_TORQUE_RPM;

	dFloat PEAK_TORQUE;
	dFloat PEAK_TORQUE_RPM;

	dFloat PEAK_HP;
	dFloat PEAK_HP_RPM;

	dFloat REDLINE_TORQUE;
	dFloat REDLINE_TORQUE_RPM;

	dFloat GEAR_1;
	dFloat GEAR_2;
	dFloat GEAR_3;
	dFloat REVERSE_GEAR;

	dFloat SUSPENSION_LENGTH;
	dFloat SUSPENSION_SPRING;
	dFloat SUSPENSION_DAMPER;
	dFloat LATERAL_STIFFNESS;
	dFloat LONGITUDINAL_STIFFNESS;
	dFloat ALIGNING_MOMENT_TRAIL;

	DifferentialType m_differentialType;
	dMatrix m_tireaLigment;
};

static BasciCarParameters basicCarParameters = 
{
	2500.0f,	// MASS
	 100.0f,	// TIRE_MASS
	  25.0f,	// STEER_ANGLE
	10000.0f,	// BRAKE_TORQUE
	  -0.6f,	// COM_Y_OFFSET
	120.0f,		// TIRE_TOP_SPEED_KMH
	 400.0f,	// IDLE_TORQUE
	500.0f,		// IDLE_TORQUE_RPM
	 500.0f,	// PEAK_TORQUE
	3000.0f,	// PEAK_TORQUE_RPM
	 300.0f,	// PEAK_HP
	4000.0f,	// PEAK_HP_RPM
	 50.0f,		// REDLINE_TORQUE
	4500.0f,	// REDLINE_TORQUE_RPM
		2.5f,	// GEAR_1
		2.0f,	// GEAR_2
		1.5f,	// GEAR_3
		2.9f,	// REVERSE_GEAR
	   0.7f,	// SUSPENSION_LENGTH
	 700.0f,	// SUSPENSION_SPRING
	  40.0f,	// SUSPENSION_DAMPER
	  20.0f,	// LATERAL_STIFFNESS
   10000.0f,	// LONGITUDINAL_STIFFNESS
	   1.5f,	// ALIGNING_MOMENT_TRAIL
	   BasciCarParameters::m_4WD,

	   dGetIdentityMatrix(),
};


#define VEHICLE_THIRD_PERSON_VIEW_HIGHT		2.0f
#define VEHICLE_THIRD_PERSON_VIEW_DIST		10.0f
#define VEHICLE_THIRD_PERSON_VIEW_FILTER		0.125f


static float VehicleHullShape0[][3] =  
{
	{-2.3f, 0.0f, -0.9f}, {-2.3f, 0.0f, 0.9f}, {2.3f, 0.0f, -0.9f}, {2.3f, 0.0f, 0.9f},
	{-2.1f, 0.7f, -0.9f}, {-2.1f, 0.7f, 0.9f}, {2.1f, 0.7f, -0.9f}, {2.1f, 0.7f, 0.9f},
};

static float VehicleHullShape1[][3] =  
{
	{-1.5f, 0.0f, -0.9f}, {-1.5f, 0.0f, 0.9f}, {1.2f, 0.0f, -0.9f}, {1.2f, 0.0f, 0.9f},
	{-1.1f, 0.7f, -0.9f}, {-1.1f, 0.7f, 0.9f}, {0.8f, 0.7f, -0.9f}, {0.8f, 0.7f, 0.9f},
};


class BasicCarEntity: public DemoEntity
{
	public: 
	BasicCarEntity (DemoEntityManager* const scene, CustomVehicleControllerManager* const manager, const dMatrix& location, const BasciCarParameters& parameters)
		:DemoEntity (dGetIdentityMatrix(), NULL)
		,m_tireaLigmentMatrix (dYawMatrix(3.141592f * 90.0f / 180.0f))
		,m_controller(NULL)
		,m_helpKey (true)
		,m_gearUpKey (false)
		,m_gearDownKey (false)
		,m_reverseGear (false)
		,m_engineKeySwitch(false)
		,m_automaticTransmission(true)
		,m_engineKeySwitchCounter(0)
		,m_engineOldKeyState(false)
		,m_engineRPMOn(false)
	{
		// add this entity to the scene for rendering
		scene->Append(this);

		// place entity in the world
		ResetMatrix (*scene, location);

		NewtonWorld* const world = scene->GetNewton();

		// create the vehicle collision shape
		NewtonCollision* const chassisCollision = CreateChassisCollision (world);

		// caret the visual mesh form the collision shape
		DemoMesh* const visualMesh = new DemoMesh ("vehicle chassis", chassisCollision, "metal_30.tga", "metal_30.tga", "metal_30.tga");
		SetMesh (visualMesh, dGetIdentityMatrix());
		visualMesh->Release();
		
		// create the coordinate system 
		dMatrix chassisMatrix;
		chassisMatrix.m_front = dVector (1.0f, 0.0f, 0.0f, 0.0f);			// this is the vehicle direction of travel
		chassisMatrix.m_up	  = dVector (0.0f, 1.0f, 0.0f, 0.0f);			// this is the downward vehicle direction
		chassisMatrix.m_right = chassisMatrix.m_front * chassisMatrix.m_up;	// this is in the side vehicle direction (the plane of the wheels)
		chassisMatrix.m_posit = dVector (0.0f, 0.0f, 0.0f, 1.0f);

		// create a default vehicle controller
		m_controller = manager->CreateVehicle (chassisCollision, chassisMatrix, parameters.MASS, dVector (0.0f, DEMO_GRAVITY, 0.0f, 0.0f));

		// get body the vehicle rigid body and set the Newton rigid body physics properties
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

		// map the gear to a look up table: gear 0 is reverse, gea 1 is neutral, gear 1 is first, gear 2 is second and so on
		for (int i = 0; i < int ((sizeof (m_gearMap) / sizeof (m_gearMap[0]))); i ++) {
			m_gearMap[i] = i;
		}
		m_gearMap[0] = 1;
		m_gearMap[1] = 0;
	}

	~BasicCarEntity ()
	{
		((CustomVehicleControllerManager*)m_controller->GetManager())->DestroyController(m_controller);
	}

	NewtonCollision* CreateChassisCollision (NewtonWorld* const world) const
	{
		dMatrix offset (dGetIdentityMatrix());
		offset.m_posit.m_y = 0.7f;

		NewtonCollision* const convex0 = NewtonCreateConvexHull (world, 8, &VehicleHullShape0[0][0], 3 * sizeof (dFloat), 0.001f, 0, NULL);
		NewtonCollision* const convex1 = NewtonCreateConvexHull (world, 8, &VehicleHullShape1[0][0], 3 * sizeof (dFloat), 0.001f, 0, &offset[0][0]);
		
		NewtonCollision* const collision = NewtonCreateCompoundCollision(world, 0);

		NewtonCompoundCollisionBeginAddRemove(collision);
		NewtonCompoundCollisionAddSubCollision (collision, convex0);
		NewtonCompoundCollisionAddSubCollision (collision, convex1);
		NewtonCompoundCollisionEndAddRemove (collision);	

		NewtonDestroyCollision (convex0);
		NewtonDestroyCollision (convex1);

		return collision;
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


	CustomVehicleControllerBodyStateTire* AddTire (const dVector& offset, dFloat width, dFloat radius, dFloat mass, dFloat suspensionLength, dFloat suspensionSpring, dFloat suspensionDamper, dFloat lateralStiffness, dFloat longitudinalStiffness, dFloat aligningMomentTrail, const dMatrix& tireAligmentMatrix) 
	{
		NewtonBody* const body = m_controller->GetBody();

		// make the tire matrix from the offset and the body matrix
		dMatrix tireMatrix (GetNextMatrix());
		tireMatrix.m_posit = offset;

		// add the visual representation of the is tire to as a child of the vehicle model 
		NewtonCollision*  const tireMeshGenerator = NewtonCreateChamferCylinder (NewtonBodyGetWorld(body), 0.5f, 1.0f, 0, NULL);
		NewtonCollisionSetScale (tireMeshGenerator, width, radius, radius);

		DemoEntity* const tireEntity = new DemoEntity (tireMatrix, this);
		DemoMesh* const visualMesh = new DemoMesh ("tireMesh", tireMeshGenerator, "smilli.tga", "smilli.tga", "smilli.tga");
		tireEntity->SetMesh (visualMesh, dYawMatrix(3.141592f * 90.0f / 180.0f));
		visualMesh->Release();
		NewtonDestroyCollision (tireMeshGenerator);

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
		tireInfo.m_aligningMomentTrail = aligningMomentTrail;
		tireInfo.m_userData = tireEntity;

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
			dMatrix matrix (m_tireaLigmentMatrix * tire->CalculateGlobalMatrix() * rootMatrixInv);
			dQuaternion rot (matrix);
			tirePart->SetMatrix(*scene, rot, matrix.m_posit);
		}
#else
		// this saves some calculation since it get the tire local to the chassis
		for (CustomVehicleControllerBodyStateTire* tire = m_controller->GetFirstTire(); tire; tire = m_controller->GetNextTire(tire)) {
			DemoEntity* const tirePart = (DemoEntity*) tire->GetUserData();
			dMatrix matrix (m_tireaLigmentMatrix * tire->CalculateLocalMatrix());
			dQuaternion rot (matrix);
			tirePart->SetMatrix(*scene, rot, matrix.m_posit);
		}
#endif

	}

	// create a simple vehicle 
	void BuidlBasicCar (const BasciCarParameters& parameters)
	{
		// step one: find the location of each tire, in the visual mesh and add them one by one to the vehicle controller 
		dFloat width = 0.35f;
		dFloat radius = 0.5f;

		// Muscle cars have the front engine, we need to shift the center of mass to the front to represent that
		m_controller->SetCenterOfGravity (dVector (0.0f, parameters.COM_Y_OFFSET, 0.0f, 0.0f)); 

		// add left tires
		CustomVehicleControllerBodyStateTire* leftTire[2]; 
		dVector offset (1.5f, 0.0f, -1.0f, 1.0f);
		leftTire[0] = AddTire (offset, width, radius, parameters.TIRE_MASS, parameters.SUSPENSION_LENGTH, parameters.SUSPENSION_SPRING, parameters.SUSPENSION_DAMPER, parameters.LATERAL_STIFFNESS, parameters.LONGITUDINAL_STIFFNESS, parameters.ALIGNING_MOMENT_TRAIL, parameters.m_tireaLigment);
		offset = dVector (-1.7f, 0.0f, -1.0f, 1.0f);
		leftTire[1] = AddTire (offset, width, radius, parameters.TIRE_MASS, parameters.SUSPENSION_LENGTH, parameters.SUSPENSION_SPRING, parameters.SUSPENSION_DAMPER, parameters.LATERAL_STIFFNESS, parameters.LONGITUDINAL_STIFFNESS, parameters.ALIGNING_MOMENT_TRAIL, parameters.m_tireaLigment);

		// add right tires
		CustomVehicleControllerBodyStateTire* rightTire[2];
		offset = dVector (1.5f, 0.0f, 1.0f, 1.0f);		
		rightTire[0] = AddTire (offset, width, radius, parameters.TIRE_MASS, parameters.SUSPENSION_LENGTH, parameters.SUSPENSION_SPRING, parameters.SUSPENSION_DAMPER, parameters.LATERAL_STIFFNESS, parameters.LONGITUDINAL_STIFFNESS, parameters.ALIGNING_MOMENT_TRAIL, parameters.m_tireaLigment);
		offset = dVector (-1.7f, 0.0f, 1.0f, 1.0f);
		rightTire[1] = AddTire (offset, width, radius, parameters.TIRE_MASS, parameters.SUSPENSION_LENGTH, parameters.SUSPENSION_SPRING, parameters.SUSPENSION_DAMPER, parameters.LATERAL_STIFFNESS, parameters.LONGITUDINAL_STIFFNESS, parameters.ALIGNING_MOMENT_TRAIL, parameters.m_tireaLigment);

		// add a steering Wheel
		CustomVehicleControllerComponentSteering* const steering = new CustomVehicleControllerComponentSteering (m_controller, parameters.STEER_ANGLE * 3.141592f / 180.0f);
		steering->AddSteeringTire (leftTire[0]);
		steering->AddSteeringTire (rightTire[0]);
		m_controller->SetSteering(steering);

		// add all wheels brakes
		CustomVehicleControllerComponentBrake* const brakes = new CustomVehicleControllerComponentBrake (m_controller, parameters.BRAKE_TORQUE);
		for (int i = 0; i < 2; i ++) {
			brakes->AddBrakeTire (leftTire[i]);
			brakes->AddBrakeTire (rightTire[i]);
		}
		m_controller->SetBrakes(brakes);

		// add hand brakes
		CustomVehicleControllerComponentBrake* const handBrakes = new CustomVehicleControllerComponentBrake (m_controller, parameters.BRAKE_TORQUE);
		handBrakes->AddBrakeTire (leftTire[1]);
		handBrakes->AddBrakeTire (rightTire[1]);
		m_controller->SetHandBrakes (handBrakes);

		// add an engine
		CustomVehicleControllerComponentEngine::dMultiAxelDifferential* differencial = NULL;
		switch (parameters.m_differentialType)
		{
			case BasciCarParameters::m_4WD:
			{
				CustomVehicleControllerComponentEngine::dSingleAxelDifferential* axles[2];
				for (int i = 0; i < 2; i ++) {
					axles[i] = new CustomVehicleControllerComponentEngine::dSingleAxelDifferential (m_controller, leftTire[i], rightTire[i]);
				}
				differencial = new CustomVehicleControllerComponentEngine::dMultiAxelDifferential (m_controller, 2, axles);
				break;
			}

			case BasciCarParameters::m_FWD:
			{
				CustomVehicleControllerComponentEngine::dSingleAxelDifferential* axles[1];
				axles[0] = new CustomVehicleControllerComponentEngine::dSingleAxelDifferential (m_controller, leftTire[0], rightTire[0]);
				differencial = new CustomVehicleControllerComponentEngine::dMultiAxelDifferential (m_controller, 1, axles);
				break;
			}

			case BasciCarParameters::m_RWD:
			default:
			{
				CustomVehicleControllerComponentEngine::dSingleAxelDifferential* axles[1];
				axles[0] = new CustomVehicleControllerComponentEngine::dSingleAxelDifferential (m_controller, leftTire[1], rightTire[1]);
				differencial = new CustomVehicleControllerComponentEngine::dMultiAxelDifferential (m_controller, 1, axles);
				break;
			}
		}
		

		dFloat fowardSpeedGearsBoxRatios[] = {parameters.GEAR_1, parameters.GEAR_1, parameters.GEAR_3};
		CustomVehicleControllerComponentEngine::dGearBox* const gearBox = new CustomVehicleControllerComponentEngine::dGearBox(m_controller, parameters.REVERSE_GEAR, sizeof (fowardSpeedGearsBoxRatios) / sizeof (fowardSpeedGearsBoxRatios[0]), fowardSpeedGearsBoxRatios); 
		CustomVehicleControllerComponentEngine* const engine = new CustomVehicleControllerComponentEngine (m_controller, gearBox, differencial);

		// the the default transmission type
		engine->SetTransmissionMode(m_automaticTransmission.GetPushButtonState());
		m_controller->SetEngine(engine);


		dFloat viperIdleRPM = parameters.IDLE_TORQUE_RPM;
		dFloat viperIdleTorquePoundPerFoot = parameters.IDLE_TORQUE;

		dFloat viperPeakTorqueRPM = parameters.PEAK_TORQUE_RPM;
		dFloat viperPeakTorquePoundPerFoot = parameters.PEAK_TORQUE;

		dFloat viperPeakHorsePowerRPM = parameters.PEAK_HP_RPM;
		dFloat viperPeakHorsePower = parameters.PEAK_HP;

		dFloat viperRedLineRPM = parameters.REDLINE_TORQUE_RPM;
		dFloat viperRedLineTorquePoundPerFoot = parameters.REDLINE_TORQUE;

		dFloat vehicleTopSpeedKPH = parameters.TIRE_TOP_SPEED_KMH;
		engine->InitEngineTorqueCurve (vehicleTopSpeedKPH, viperIdleTorquePoundPerFoot, viperIdleRPM, viperPeakTorquePoundPerFoot, viperPeakTorqueRPM, viperPeakHorsePower, viperPeakHorsePowerRPM, viperRedLineTorquePoundPerFoot, viperRedLineRPM);

		// do not forget to call finalize after all components are added or after any change is made to the vehicle
		m_controller->Finalize();


/*
		// test tire update interface
		for (CustomVehicleControllerBodyStateTire* node = m_controller->GetFirstTire(); node; node = m_controller->GetNextTire(node)) {
			CustomVehicleControllerBodyStateTire& tire = *node;
			CustomVehicleControllerBodyStateTire::TireCreationInfo tireInfo;
			tire.GetInfo(tireInfo);
			tireInfo.m_location.m_y += 0.1f;
			tire.UpdateInfo(tireInfo);
			tire.GetInfo(tireInfo);
			tire.GetInfo(tireInfo);
		}

		m_controller->Finalize();

		for (CustomVehicleControllerBodyStateTire* node = m_controller->GetFirstTire(); node; node = m_controller->GetNextTire(node)) {
			CustomVehicleControllerBodyStateTire& tire = *node;
			CustomVehicleControllerBodyStateTire::TireCreationInfo tireInfo;
			tire.GetInfo(tireInfo);
			tire.GetInfo(tireInfo);
		}
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

		int gear = engine ? engine->GetGear() : CustomVehicleControllerComponentEngine::dGearBox::m_newtralGear;
		dFloat steeringVal = 0.0f;
		dFloat engineGasPedal = 0.0f;
		dFloat brakePedal = 0.0f;
		dFloat handBrakePedal = 0.0f;

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
		steeringVal = (dFloat (mainWindow->GetKeyState ('A')) - dFloat (mainWindow->GetKeyState ('D')));

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
		// set the help key
		m_helpKey.UpdatePushButton (mainWindow, 'H');

		int reverseGear = m_reverseGear.UpdateTriggerButton (mainWindow, 'R') ? 1 : 0;

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
			fwrite (&reverseGear, sizeof (int), 1, file);
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
			fread (&reverseGear, sizeof (int), 1, file);
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
			if (reverseGear) {
				if (engine->GetGear() == CustomVehicleControllerComponentEngine::dGearBox::m_reverseGear) {
					engine->SetGear(CustomVehicleControllerComponentEngine::dGearBox::m_newtralGear);
				} else {
					engine->SetGear(CustomVehicleControllerComponentEngine::dGearBox::m_reverseGear);
				}
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
//		dAssert (0);
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
p0 += chassis.GetMatrix()[2].Scale ((tireMatrix.m_posit.m_z > 0.0f ? 1.0f : -1.0f) * 0.25f);

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
			glColor3f (0.0f, 1.0f, 0.0f);
			glVertex3f (p0.m_x, p0.m_y, p0.m_z);
			glVertex3f (p3.m_x, p3.m_y, p3.m_z);
		}

		glEnd();

		glLineWidth(1.0f);

	}

	dMatrix m_tireaLigmentMatrix;
	CustomVehicleController* m_controller;
	DemoEntityManager::ButtonKey m_helpKey;
	DemoEntityManager::ButtonKey m_gearUpKey;
	DemoEntityManager::ButtonKey m_gearDownKey;
	DemoEntityManager::ButtonKey m_reverseGear;
	DemoEntityManager::ButtonKey m_engineKeySwitch;
	DemoEntityManager::ButtonKey m_automaticTransmission;
	int m_gearMap[10];
	int m_engineKeySwitchCounter;
	bool m_engineOldKeyState;
	bool m_engineRPMOn;
};

class BasicCarControllerManager: public CustomVehicleControllerManager
{
	public:
	BasicCarControllerManager (NewtonWorld* const world)
		:CustomVehicleControllerManager (world)
		,m_externalView(true)
		,m_player (NULL) 
	{
		// hook a callback for 2d help display
		DemoEntityManager* const scene = (DemoEntityManager*) NewtonWorldGetUserData(world);
		scene->Set2DDisplayRenderFunction (RenderVehicleHud, this);
	}

	~BasicCarControllerManager ()
	{
	}

	static void RenderVehicleHud (DemoEntityManager* const scene, void* const context, int lineNumber)
	{
		BasicCarControllerManager* const me = (BasicCarControllerManager*) context;
		me->RenderVehicleHud (scene, lineNumber);
	}

	void DrawHelp(DemoEntityManager* const scene, int lineNumber) const
	{
		if (m_player->m_helpKey.GetPushButtonState()) {
			dVector color(1.0f, 1.0f, 0.0f, 0.0f);
			lineNumber = scene->Print (color, 10, lineNumber + 20, "Vehicle driving keyboard control");
			lineNumber = scene->Print (color, 10, lineNumber + 20, "key switch          : 'Y'");
			lineNumber = scene->Print (color, 10, lineNumber + 20, "accelerator         : 'W'");
			lineNumber = scene->Print (color, 10, lineNumber + 20, "brakes              : 'S'");
			lineNumber = scene->Print (color, 10, lineNumber + 20, "had brakes          : 'space'");
			lineNumber = scene->Print (color, 10, lineNumber + 20, "turn right          : 'D'");
			lineNumber = scene->Print (color, 10, lineNumber + 20, "turn right          : 'S'");
			lineNumber = scene->Print (color, 10, lineNumber + 20, "toggle Reverse Gear : 'R'");
			lineNumber = scene->Print (color, 10, lineNumber + 20, "hide help           : 'H'");
		}
	}


	void RenderVehicleHud (DemoEntityManager* const scene, int lineNumber) const
	{
		if (m_player) {
			// set to transparent color
			glEnable (GL_BLEND);
			glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

			//print controllers help 
			DrawHelp(scene, lineNumber);

			// restore color and blend mode
			glDisable (GL_BLEND);
		}
	}


	void SetAsPlayer (BasicCarEntity* const player)
	{
		m_player = player;
	}


	virtual void PreUpdate (dFloat timestep)
	{
		// apply the vehicle controls, and all simulation time effect
		for (dListNode* ptr = GetFirst(); ptr; ptr = ptr->GetNext()) {
			CustomVehicleController* const controller = &ptr->GetInfo();

			NewtonBody* const body = controller->GetBody();
			BasicCarEntity* const vehicleEntity = (BasicCarEntity*) NewtonBodyGetUserData(body);

			if (vehicleEntity == m_player) {
				// do player control
				vehicleEntity->ApplyPlayerControl ();
			} else {
				// do no player control
				vehicleEntity->ApplyNPCControl ();
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
			BasicCarEntity* const vehicleEntity = (BasicCarEntity*)NewtonBodyGetUserData (controller->GetBody());
			vehicleEntity->UpdateTireTransforms();
		}

		UpdateCamera (timestep);
	}


	void UpdateCamera (dFloat timestep)
	{
		if (m_player) {
			DemoEntityManager* const scene = (DemoEntityManager*) NewtonWorldGetUserData(GetWorld());
			DemoCamera* const camera = scene->GetCamera();
			dMatrix camMatrix (camera->GetNextMatrix ());
			dMatrix playerMatrix (m_player->GetNextMatrix());

			dVector frontDir (camMatrix[0]);
			dVector camOrigin; 
			if (m_externalView) {
				camOrigin = playerMatrix.m_posit + dVector(0.0f, VEHICLE_THIRD_PERSON_VIEW_HIGHT, 0.0f, 0.0f);
				camOrigin -= frontDir.Scale (VEHICLE_THIRD_PERSON_VIEW_DIST);
			} else {
				dAssert (0);
				//            camMatrix = camMatrix * playerMatrix;
				//            camOrigin = playerMatrix.TransformVector(dVector(-0.8f, ARTICULATED_VEHICLE_CAMERA_EYEPOINT, 0.0f, 0.0f));
			}

			camera->SetNextMatrix (*scene, camMatrix, camOrigin);
		}
	}

	// use this to display debug information about vehicle 
	void Debug () const
	{
		for (dListNode* ptr = GetFirst(); ptr; ptr = ptr->GetNext()) {
			CustomVehicleController* const controller = &ptr->GetInfo();
			BasicCarEntity* const vehicleEntity = (BasicCarEntity*)NewtonBodyGetUserData (controller->GetBody());
			vehicleEntity->Debug();
		}
	}

	bool m_externalView;
	BasicCarEntity* m_player;
};


void BasicCar (DemoEntityManager* const scene)
{
	// load the sky box
	scene->CreateSkyBox();

	CreateLevelMesh (scene, "flatPlane.ngd", 1);
//	CreateHeightFieldTerrain (scene, 10, 8.0f, 5.0f, 0.2f, 200.0f, -50.0f);


	dMatrix location (dGetIdentityMatrix());
	location.m_posit = dVector (0.0f, 10.0f, 0.0f, 1.0f);

	location.m_posit = FindFloor (scene->GetNewton(), location.m_posit, 100.0f);
	location.m_posit.m_y += 2.0f;

	NewtonWorld* const world = scene->GetNewton();

	// create a vehicle controller manager
	BasicCarControllerManager* const manager = new BasicCarControllerManager (world);
	
	// load 
	basicCarParameters.m_differentialType = BasciCarParameters::m_RWD;
	BasicCarEntity* const heavyVehicle = new BasicCarEntity (scene, manager, location, basicCarParameters);
	heavyVehicle->BuidlBasicCar (basicCarParameters);

	// set this vehicle as the player
	manager->SetAsPlayer(heavyVehicle);

	dMatrix camMatrix (manager->m_player->GetNextMatrix());
	scene->SetCameraMouseLock (true);
	scene->SetCameraMatrix(camMatrix, camMatrix.m_posit);

/*
	int count = 1;
	dMatrix shapeOffsetMatrix (dGetIdentityMatrix());
	int defaultMaterialID = NewtonMaterialGetDefaultGroupID (scene->GetNewton());

	dVector size (3.0f, 0.125f, 3.0f, 0.0f);
	//AddPrimitiveArray(scene, 100.0f, location.m_posit, size, count, count, 5.0f, _BOX_PRIMITIVE, defaultMaterialID, shapeOffsetMatrix);

	size = dVector(1.0f, 0.5f, 1.0f, 0.0f);
	//	AddPrimitiveArray(scene, 10.0f, location.m_posit, size, count, count, 5.0f, _SPHERE_PRIMITIVE, defaultMaterialID, shapeOffsetMatrix);
	//	AddPrimitiveArray(scene, 10.0f, location.m_posit, size, count, count, 5.0f, _BOX_PRIMITIVE, defaultMaterialID, shapeOffsetMatrix);
	//	AddPrimitiveArray(scene, 10.0f, location.m_posit, size, count, count, 5.0f, _CAPSULE_PRIMITIVE, defaultMaterialID, shapeOffsetMatrix);
		AddPrimitiveArray(scene, 10.0f, location.m_posit, size, count, count, 5.0f, _CYLINDER_PRIMITIVE, defaultMaterialID, shapeOffsetMatrix);
	//	AddPrimitiveArray(scene, 10.0f, location.m_posit, size, count, count, 5.0f, _TAPERED_CAPSULE_PRIMITIVE, defaultMaterialID, shapeOffsetMatrix);
	//	AddPrimitiveArray(scene, 10.0f, location.m_posit, size, count, count, 5.0f, _TAPERED_CYLINDER_PRIMITIVE, defaultMaterialID, shapeOffsetMatrix);
	//	AddPrimitiveArray(scene, 10.0f, location.m_posit, size, count, count, 5.0f, _CHAMFER_CYLINDER_PRIMITIVE, defaultMaterialID, shapeOffsetMatrix);
	//	AddPrimitiveArray(scene, 10.0f, location.m_posit, size, count, count, 5.0f, _CONE_PRIMITIVE, defaultMaterialID, shapeOffsetMatrix);
	//	AddPrimitiveArray(scene, 10.0f, location.m_posit, size, count, count, 5.0f, _REGULAR_CONVEX_HULL_PRIMITIVE, defaultMaterialID, shapeOffsetMatrix);
	//	AddPrimitiveArray(scene, 10.0f, location.m_posit, size, count, count, 5.0f, _RANDOM_CONVEX_HULL_PRIMITIVE, defaultMaterialID, shapeOffsetMatrix);
	//	NewtonSerializeToFile (scene->GetNewton(), "C:/Users/Julio/Desktop/newton-dynamics/applications/media/xxxxx.bin");
*/		
}

