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
#include "dNewtonMesh.h"
#include "dNewtonCollision.h"


dNewtonMesh::dNewtonMesh(dNewton* const world)
	:m_mesh (NewtonMeshCreate (world->GetNewton()))
{
}

dNewtonMesh::dNewtonMesh(const dNewtonMesh& clone)
	:m_mesh (NewtonMeshCreateFromMesh (clone.m_mesh))
{
}

dNewtonMesh::dNewtonMesh(const dNewtonCollision& collision)
	:m_mesh (NewtonMeshCreateFromCollision(collision.GetShape()))
{
}

dNewtonMesh::dNewtonMesh(dNewton* const world, int pointCount, const dFloat* const vertexCloud, int strideInBytes, dFloat tolerance)
	:m_mesh (NewtonMeshCreateConvexHull (world->GetNewton(), pointCount, vertexCloud, strideInBytes, tolerance))
{
}

void dNewtonMesh::CreateVoronoiConvexDecomposition (const dNewtonMesh& contexMesh)
{
	NewtonMeshDestroy (m_mesh);
	dAssert (0);
//	m_mesh = NewtonMeshCreateVoronoiConvexDecomposition (const NewtonWorld* const newtonWorld, int pointCount, const dFloat* const vertexCloud, int strideInBytes, int materialID, const dFloat* const textureMatrix);
}

void dNewtonMesh::CreateApproximateConvexDecomposition (const dNewtonMesh& mesh, dFloat maxConcavity, dFloat backFaceDistanceFactor, int maxCount, int maxVertexPerHull)
{
	NewtonMeshDestroy (m_mesh);
	m_mesh = NewtonMeshApproximateConvexDecomposition (mesh.m_mesh, maxConcavity, backFaceDistanceFactor, maxCount, maxVertexPerHull, NULL);
}



dNewtonMesh::~dNewtonMesh()
{
	NewtonMeshDestroy (m_mesh);
}



