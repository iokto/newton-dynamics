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
#include "TargaToOpenGl.h"
#include "DemoMesh.h"
#include "DemoEntityManager.h"
#include "DemoCamera.h"
#include "PhysicsUtils.h"



#define PARALLET_RAYS_COUNT 1000

class dRayCastRecord
{
	CUSTOM_CONTROLLER_GLUE(dRayCastRecord);

	public:

	void PreUpdate(dFloat timestep, int threadIndex)
	{
	}


	static dFloat RayCast (const NewtonBody* const body, const NewtonCollision* const collisionHit, const dFloat* const contact, const dFloat* const normal, dLong collisionID, void* const userData, dFloat intersetParam)
	{
		dFloat* const paramPtr = (dFloat*)userData;
		if (intersetParam < paramPtr[0]) {
			paramPtr[0] = intersetParam;
		}
		return paramPtr[0];
	}

	void PostUpdate(dFloat timestep, int threadIndex)
	{
		DemoEntity* const targetEnt = (DemoEntity*) NewtonBodyGetUserData (m_target);
		const dMatrix& matrix = targetEnt->GetRenderMatrix(); 

		float parameter = 1.1f;
		NewtonWorld* const world = NewtonBodyGetWorld (m_target);
		NewtonWorldRayCast(world, &m_p0[0], &matrix.m_posit[0], RayCast, &parameter, NULL, threadIndex);

		dAssert (parameter <= 1.0f);
		m_p1 = m_p0 + (matrix.m_posit - m_p0).Scale (parameter);
	}

	dVector m_p0;
	dVector m_p1;
	NewtonBody* m_target;
};



class LineOfSightRayCastEntity: public DemoEntity
{
	public:
	LineOfSightRayCastEntity (DemoEntityManager* const scene, CustomControllerManager<dRayCastRecord>* casterManager)
		:DemoEntity (GetIdentityMatrix(), NULL)
		,m_casterManager(casterManager)
	{
		scene->Append(this);
	}

	~LineOfSightRayCastEntity ()
	{
	}

	void Render(dFloat timeStep) const
	{
		glDisable (GL_LIGHTING);
		glDisable(GL_TEXTURE_2D);

	
		glColor3f(0.0f, 0.5f, 0.5f);
		glBegin(GL_LINES);
		//for (void* node = m_caterManager->GetFirstController(); node; node = m_caterManager->GetNextController(node)) {
		for (CustomControllerManager<dRayCastRecord>::CustomController* ray = m_casterManager->GetFirstController(); ray; ray = m_casterManager->GetNextController(ray)) {
			//dRayCastRecord* const ray = (dRayCastRecord*) m_caterManager->GetControllerFromNode(node);
			glVertex3f(ray->m_p0.m_x, ray->m_p0.m_y, ray->m_p0.m_z);
			glVertex3f(ray->m_p1.m_x, ray->m_p1.m_y, ray->m_p1.m_z);
		}
		glEnd();


		glColor3f(1.0f, 0.0f, 0.0f);
		glPointSize(6.0f);
		glBegin(GL_POINTS);
		for (CustomControllerManager<dRayCastRecord>::CustomController* ray = m_casterManager->GetFirstController(); ray; ray = m_casterManager->GetNextController(ray)) {
//			dRayCastRecord* const ray = (dRayCastRecord*) m_caterManager->GetControllerFromNode(node);
			glVertex3f(ray->m_p1.m_x, ray->m_p1.m_y, ray->m_p1.m_z);
		}
		glEnd();
		glPointSize(1.0f);

		glColor3f(1.0f, 1.0f, 1.0f);
	}

	CustomControllerManager<dRayCastRecord>* m_casterManager; 
};


class dRayCasterManager: public CustomControllerManager<dRayCastRecord> 
{
	public:
	dRayCasterManager(DemoEntityManager* const scene, NewtonBody* const skipLevelMesh)
		:CustomControllerManager<dRayCastRecord>(scene->GetNewton(), "dRayCasterManager")
	{
		// make 16 casting center
		dVector location[16];
		for (int i = 0; i < 4; i ++) {
			float x = (i - 2) * 4.0f;
			for (int j = 0; j < 4; j ++) {
				float z = (j - 2) * 4.0f;
				location [i * 4 + j] = dVector (x, 30.0f, z, 0.0f);
			}
		}

		// each cast center will cast a group of bodies
		NewtonBody* body = NULL;
		NewtonWorld* const world = scene->GetNewton();
		for (int index = 0; index < PARALLET_RAYS_COUNT; index ++) {
			int locationIndex = dRand() & 15;

			while (!body || (body == skipLevelMesh)) {
				if (!body) {
					body = NewtonWorldGetFirstBody(world);
				}else if (body == skipLevelMesh) {
					body = NewtonWorldGetNextBody(world, body);
				}
			}

			CreateCaster (location[locationIndex], body);
			body = NewtonWorldGetNextBody(world, body);
		}

		// add a visual entity so that we can render the ray cast
		new LineOfSightRayCastEntity (scene, this);
	}

	virtual ~dRayCasterManager()
	{
	}

	dRayCastRecord* CreateCaster (const dVector& origin, NewtonBody* const targetBody)
	{
		dRayCastRecord* const caster = (dRayCastRecord*) CreateController();
		caster->m_p0 = origin;
		caster->m_p1 = dVector (0.0f, 0.0f, 0.0f, 0.0f);
		caster->m_target = targetBody;
		return caster;
	}

	virtual void Debug () const {};

};


void MultiRayCast (DemoEntityManager* const scene)
{
	// load the skybox
	scene->CreateSkyBox();

	// load the scene from a ngd file format
	NewtonBody* const levelBody = CreateLevelMesh (scene, "flatPlane.ngd", true);
//	CreateLevelMesh (scene, "playground.ngd", true);
//	CreateLevelMesh (scene, "sponza.ngd", true);


	NewtonWorld* const world = scene->GetNewton();
	int defaultMaterialID = NewtonMaterialGetDefaultGroupID (world);
	dVector location0 (0.0f, 0.0f, 0.0f, 0.0f);
	dVector location1 (0.2f, 0.0f, 0.0f, 0.0f);
	dVector location2 (0.0f, 0.0f, 0.2f, 0.0f);
	dVector location3 (0.2f, 0.0f, 0.2f, 0.0f);
	dVector size (1.0f, 0.5f, 0.5f, 0.0f);

	dMatrix shapeOffsetMatrix (GetIdentityMatrix());

	// populate the world with few objects
	int count = 8;
	float separation = 4.0f;
	AddPrimitiveArray(scene, 10.0f, location0, size, count, count, separation, _SPHERE_PRIMITIVE, defaultMaterialID, shapeOffsetMatrix);
	AddPrimitiveArray(scene, 10.0f, location1, size, count, count, separation, _BOX_PRIMITIVE, defaultMaterialID, shapeOffsetMatrix);
	AddPrimitiveArray(scene, 10.0f, location2, size, count, count, separation, _CAPSULE_PRIMITIVE, defaultMaterialID, shapeOffsetMatrix);
	AddPrimitiveArray(scene, 10.0f, location3, size, count, count, separation, _CYLINDER_PRIMITIVE, defaultMaterialID, shapeOffsetMatrix);
	AddPrimitiveArray(scene, 10.0f, location0, size, count, count, separation, _CHAMFER_CYLINDER_PRIMITIVE, defaultMaterialID, shapeOffsetMatrix);
	AddPrimitiveArray(scene, 10.0f, location1, size, count, count, separation, _BOX_PRIMITIVE, defaultMaterialID, shapeOffsetMatrix);
	AddPrimitiveArray(scene, 10.0f, location2, size, count, count, separation, _CAPSULE_PRIMITIVE, defaultMaterialID, shapeOffsetMatrix);
	AddPrimitiveArray(scene, 10.0f, location3, size, count, count, separation, _REGULAR_CONVEX_HULL_PRIMITIVE, defaultMaterialID, shapeOffsetMatrix);
	AddPrimitiveArray(scene, 10.0f, location0, size, count, count, separation, _RANDOM_CONVEX_HULL_PRIMITIVE, defaultMaterialID, shapeOffsetMatrix);


	// add special Entity to ray cast these object 
	new dRayCasterManager (scene, levelBody);

	// place camera into position
	dQuaternion rot;
//	dVector origin (-40.0f, 10.0f, 0.0f, 0.0f);
	dVector origin (-30.0f, 5.0f, 0.0f, 0.0f);
	scene->SetCameraMatrix(rot, origin);

	//ExportScene (scene->GetNewton(), "../../../media/test1.ngd");
}


