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
#include "DemoEntityManager.h"
#include "DemoCamera.h"
#include "DemoMesh.h"
#include "PhysicsUtils.h"


class TornadoTwiester: public DemoEntity
{
	public:
	class PerThreadData
	{
		public:
		NewtonBody* m_body;
		TornadoTwiester* m_me;
	};

	TornadoTwiester (DemoEntityManager* const scene, const dVector& origin, const dVector& size, const dVector& omega)
		:DemoEntity (GetIdentityMatrix(), NULL)
		,m_omega (omega)
		,m_veloc (dFloat (0.0f), dFloat (0.0f), dFloat (0.0f), dFloat (0.0f))
		,m_prevPosit (dFloat (0.0f), dFloat (0.0f), dFloat (0.0f), dFloat (0.0f))
	{
		dMatrix matrix (GetIdentityMatrix());
		matrix.m_posit = origin; 
		matrix.m_posit.m_w = 1.0f;
		m_prevPosit = matrix.m_posit;

		dMatrix aligment (dRollMatrix(3.141592f * 0.5f));
		aligment.m_posit.m_w = 1.0f;
		aligment.m_posit.m_y += size.m_x / 2.0f;

		
		NewtonWorld* const world = scene->GetNewton();
		// create a large body to represent the volume of a tornado
		//NewtonCollision* const collision = CreateConvexCollision (world, &aligment[0][0], size, _CAPSULE_PRIMITIVE, 0);
		NewtonCollision* const collision = CreateConvexCollision (world, &aligment[0][0], size, _BOX_PRIMITIVE, 0);
		
		// the tornado is represented by a kinematic body
		m_tornadoBody = NewtonCreateKinematicBody(world, collision, &matrix[0][0]);
		NewtonBodySetUserData(m_tornadoBody, this);
		NewtonBodySetTransformCallback(m_tornadoBody, DemoEntity::TransformCallback);
		

		DemoMesh* const geometry = new DemoMesh("capsule", collision, "smilli.tga", "smilli.tga", "smilli.tga");
		// I need to make this geometry transparent
//		SetMesh(geometry);
		// do not forget to release the assets	
		geometry->Release(); 
		NewtonDestroyCollision (collision);

		scene->Append (this);
	}

	~TornadoTwiester ()
	{
	}

/*
	virtual void SimulationPreListener(DemoEntityManager* const scene, DemoEntityManager::dListNode* const mynode, dFloat timestep)
	{
		dAssert (0);

		NewtonWorld* const world = scene->GetNewton();

		// calculate the tornado linear velocity
		dMatrix matrix;
		NewtonBodyGetMatrix(m_tornadoBody, &matrix[0][0]);
		m_veloc = (matrix.m_posit - m_prevPosit).Scale(1.0f/ m_timestep);

		// save the position and current positions
		m_origin = matrix.m_posit;
		m_prevPosit = matrix.m_posit;

		// make sure this body is set to never sleep
		NewtonBodySetAutoSleep(m_tornadoBody, 0);

		// submit all bodies in parallel
		int count = 0;
		m_timestep = timestep;
		for (NewtonJoint* joint = NewtonBodyGetFirstContactJoint(m_tornadoBody); joint; joint = NewtonBodyGetNextContactJoint (m_tornadoBody, joint)) {
			dFloat Ixx;
			dFloat Iyy;
			dFloat Izz;
			dFloat mass;

			NewtonBody* const body0 = NewtonJointGetBody0(joint);
			NewtonBody* const body1 = NewtonJointGetBody1(joint);
			NewtonBody* const body = (NewtonBodyGetUserData(body0) == this) ? body1 : body0;	
			NewtonBodyGetMassMatrix(body, &mass, &Ixx, &Iyy, &Izz);
			if (mass > 0.0f) {
				m_perThreadData[count].m_me = this;
				m_perThreadData[count].m_body = (NewtonBodyGetUserData(body0) == this) ? body1 : body0;
				NewtonDispachThreadJob(world, ApplyTornadoForcesKernel, &m_perThreadData[count]);
				count ++;
				if (count >= int ((sizeof (m_perThreadData) / sizeof (m_perThreadData[0])))) {
					NewtonSyncThreadJobs(world);
					count = 0;
				}
			}
		}
		if (count) {
			NewtonSyncThreadJobs(world);
		}
	}
*/

	static void ApplyTornadoForcesKernel (void* const context, int threadIndex)
	{
		PerThreadData* const data = (PerThreadData*) context;
		data->m_me->ApplyTornadoForces (data->m_body);
		
	}
	
	void ApplyTornadoForces (NewtonBody* const body) const
	{
		// a tornado is volume of space with a low pleasure center that point in the direction of the center of the tornado in the vertical direction.
		// the motion of a body generated a Coriolis's forces generated by the a time variant distance from the body to the center of the tornado 
		// plus centripetal acceleration generated but the relative angle of the radio between a body and the tornado center.

		//this  can be derive as follow
		// let the tornado has a location p and every body inside the tornado be a a position r relative to the tornado center.
		// the tornado center will have a velocity v, an each body inside the tornado volume will has a instance angular velocity relave to the tornado origin

		//this can be written as
		//pb = p + r
		//where pb is the position of a particle in tornado volume
		//the time derivation of this is given by 
		//vp = v + w x r + r'

		//where vp is the velocity of a body in the tornado volume
		//v is the velocity of the tornado center
		//w is the angular velocity of radio a in the tornado filed
		//r' is the velocity of the body in the tornado filed alone the ratios r, this velocity is just the velocy that 
		//is cerated by the low pleasure center.
		//the acceleration is given by the second derivative of vp
		//ap = a + w x (w x r + r') + alpha x r + w x r'

		//the tornado is not lineally accelerating  nor it is angularly acceleration, therefore  a = 0 and alpha x r'= 0
		// we also assume for simplicity that that that the tornado si moving and rotation with constant velocities

		//with these assumption we get
		 //ap = w x w x r + 2.0f * w x r';

		//this is the acceleration generated by the tornado which is a low pressure center at the origin of the volume
		// w x w x r  which is the centripetal acceleration  generated by a rotating frame and point tornado the origin on the tornado, and
		// 2.0f * w x r' which ios eh Coriolis acceleration generation by a time varying radios on a rotation frame, and point in a directions perpendicular to the radio
		// so let us do this


		dFloat Ixx;
		dFloat Iyy;
		dFloat Izz;
		dFloat mass;

		dVector veloc;
		dMatrix matrix;

		NewtonBodyGetVelocity(body, &veloc[0]);
		NewtonBodyGetMatrix(body, &matrix[0][0]);

		dVector r (matrix.m_posit - m_origin);
		// calculate centripetal acceleration  
		dVector centripetal (m_omega * (m_omega * r));

		//calculate Coriolis
		dVector coriolis (0.0, 0.0, 0.0, 0.0);
		dFloat mag2 = r % r;
		if (mag2 > 0.0f) {
			dVector dir (r.Scale (1.0f / dSqrt(mag2)));

			dVector rr (veloc - m_veloc);
			rr = dir.Scale(2.0f * (rr % dir));
			coriolis = m_omega * rr;
		}

		// wake up the body for this frame while it is inside the tornado volume
		NewtonBodySetSleepState (body, 0);
		NewtonBodyGetMassMatrix(body, &mass, &Ixx, &Iyy, &Izz);
		dVector tornadoForce ((coriolis + centripetal).Scale(mass)) ;
		NewtonBodyAddForce(body, &tornadoForce[0]);
	}

	dVector m_omega;
	dVector m_veloc;
	dVector m_origin;
	dVector m_prevPosit;
	dFloat m_timestep;
	NewtonBody* m_tornadoBody;	
	PerThreadData m_perThreadData[1024];
};




void KinematicBodies (DemoEntityManager* const scene)
{
	// load the skybox
	scene->CreateSkyBox();

	// load the scene from a ngd file format
	CreateLevelMesh (scene, "flatPlane.ngd", 1);

	// build a large kinematic capsule that represent the volume of a tornado
	new TornadoTwiester (scene, dVector(0.0f, 0.0f, 0.0f, 0.0f), dVector (100, 100, 100.0f, 0.0), dVector (0, 5, 0, 0));
	

	// add some bodies 
	int defaultMaterialID = NewtonMaterialGetDefaultGroupID (scene->GetNewton());
	dVector location (0,0,0,0);
	location.m_x += 0.0f;
	location.m_z += 0.0f;
	dVector size (0.5f, 0.5f, 0.75f, 0.0f);

	int count = 1;
	dMatrix shapeOffsetMatrix (GetIdentityMatrix());
//	AddPrimitiveArray(scene, 10.0f, location, size, count, count, 5.0f, _SPHERE_PRIMITIVE, defaultMaterialID, shapeOffsetMatrix);
	AddPrimitiveArray(scene, 10.0f, location, size, count, count, 5.0f, _BOX_PRIMITIVE, defaultMaterialID, shapeOffsetMatrix);
//	AddPrimitiveArray(scene, 10.0f, location, size, count, count, 5.0f, _CAPSULE_PRIMITIVE, defaultMaterialID, shapeOffsetMatrix);
//	AddPrimitiveArray(scene, 10.0f, location, size, count, count, 5.0f, _CYLINDER_PRIMITIVE, defaultMaterialID, shapeOffsetMatrix);
//	AddPrimitiveArray(scene, 10.0f, location, size, count, count, 5.0f, _CONE_PRIMITIVE, defaultMaterialID, shapeOffsetMatrix);
//	AddPrimitiveArray(scene, 10.0f, location, size, count, count, 5.0f, _TAPERED_CAPSULE_PRIMITIVE, defaultMaterialID, shapeOffsetMatrix);
//	AddPrimitiveArray(scene, 10.0f, location, size, count, count, 5.0f, _TAPERED_CYLINDER_PRIMITIVE, defaultMaterialID, shapeOffsetMatrix);
//	AddPrimitiveArray(scene, 10.0f, location, size, count, count, 5.0f, _CHAMFER_CYLINDER_PRIMITIVE, defaultMaterialID, shapeOffsetMatrix);
//	AddPrimitiveArray(scene, 10.0f, location, size, count, count, 5.0f, _REGULAR_CONVEX_HULL_PRIMITIVE, defaultMaterialID, shapeOffsetMatrix);
//	AddPrimitiveArray(scene, 10.0f, location, size, count, count, 5.0f, _RANDOM_CONVEX_HULL_PRIMITIVE, defaultMaterialID, shapeOffsetMatrix);
//	AddPrimitiveArray(scene, 10.0f, location, size, count, count, 5.0f, _COMPOUND_CONVEX_CRUZ_PRIMITIVE, defaultMaterialID, shapeOffsetMatrix);


	// place camera into position
	dQuaternion rot;
	dVector origin (-100.0f, 40.0f, 0.0f, 0.0f);
	scene->SetCameraMatrix(rot, origin);



//	ExportScene (scene->GetNewton(), "../../../media/test1.ngd");

}

