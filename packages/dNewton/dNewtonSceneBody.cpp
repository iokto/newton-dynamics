/* Copyright (c) <2003-2011> <Julio Jerez, Newton Game Dynamics>
* 
* This software is provided 'as-is', without any express or implied
* warranty. In no event will the authors be held liable for any damages
* arising from the use of this software.
* 
* Permission is granted to anyone to use this software for any purpose,
* including commercial applications, and to alter it and redistribute it
* freely, subject to the following restrictions:
* 
* 1. The origin of this software must not be misrepresented; you must not
* claim that you wrote the original software. If you use this software
* in a product, an acknowledgment in the product documentation would be
* appreciated but is not required.
* 
* 2. Altered source versions must be plainly marked as such, and must not be
* misrepresented as being the original software.
* 
* 3. This notice may not be removed or altered from any source distribution.
*/

#include "dStdAfxNewton.h"
#include "dNewton.h"
#include "dNewtonCollision.h"
#include "dNewtonSceneBody.h"


/*
dNewtonSceneBody::dNewtonSceneBody(dNewton* const dWorld)
	:dNewtonBody()
{
	NewtonWorld* const world = dWorld->GetNewton ();

	dNewtonCollisionScene collision (dWorld);

	dMatrix matrix (GetIdentityMatrix());
	SetBody (NewtonCreateDynamicBody (world, collision.GetShape(), &matrix[0][0]));
}

dNewtonSceneBody::~dNewtonSceneBody()
{
}

void dNewtonSceneBody::BeginAddRemoveCollision()
{
	dNewtonCollisionScene* const scene = (dNewtonCollisionScene*) GetCollision();
	scene->BeginAddRemoveCollision();
}

void dNewtonSceneBody::EndAddRemoveCollision()
{
	
	dNewtonCollisionScene* const scene = (dNewtonCollisionScene*) GetCollision();
	scene->EndAddRemoveCollision();

	// need to update the aabb in the broad phase, for this we call set matrix
	dMatrix matrix;
	NewtonBody* const body = GetNewtonBody();
	NewtonBodyGetMatrix(body, &matrix[0][0]);
	NewtonBodySetMatrix(body, &matrix[0][0]);
}

void* dNewtonSceneBody::AddCollision(const dNewtonCollision* const collision)
{
	dNewtonCollisionScene* const scene = (dNewtonCollisionScene*) GetCollision();
	return scene->AddCollision(collision);
}

void dNewtonSceneBody::RemoveCollision (void* const handle)
{
	dNewtonCollisionScene* const scene = (dNewtonCollisionScene*) GetCollision();
	scene->RemoveCollision(handle);
}
*/

