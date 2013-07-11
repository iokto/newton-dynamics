/* Copyright (c) <2003-2013> <Julio Jerez, Newton Game Dynamics>
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

#ifndef _D_NEWTON_MESH_H_
#define _D_NEWTON_MESH_H_

#include "dStdAfxNewton.h"
#include "dNewtonAlloc.h"

class dNewton;
class dNewtonCollision;

class dNewtonMesh: public dNewtonAlloc
{
	public:
	class dUV
	{
		public:
		dFloat m_u;
		dFloat m_v;
	};

	class dPoint
	{
		public:
		dFloat m_x;
		dFloat m_y;
		dFloat m_z;
	};

	CNEWTON_API dNewtonMesh(dNewton* const world);
	CNEWTON_API dNewtonMesh(const dNewtonMesh& clone);
	CNEWTON_API dNewtonMesh(const dNewtonCollision& collision);
	CNEWTON_API dNewtonMesh(dNewton* const world, int pointCount, const dFloat* const vertexCloud, int strideInBytes, dFloat tolerance);
	CNEWTON_API virtual ~dNewtonMesh();

	CNEWTON_API void CreateVoronoiConvexDecomposition (const dNewtonMesh& convexMesh);
	CNEWTON_API void CreateApproximateConvexDecomposition (const dNewtonMesh& mesh, dFloat maxConcavity, dFloat backFaceDistanceFactor, int maxCount, int maxVertexPerHull);

	CNEWTON_API int GetPointCount() const;
	CNEWTON_API void GetVertexStreams(dPoint* const posit, dPoint* const normal, dUV* const uv0, dUV* const uv1) const;

	CNEWTON_API int GetTotalIndexCount() const;
	CNEWTON_API int GetTotalFaceCount() const;

	CNEWTON_API void* BeginMaterialHandle () const; 
	CNEWTON_API void EndMaterialHandle (void* const materialHandle) const; 

	CNEWTON_API void ApplySphericalMapping (int matId); 
	CNEWTON_API void ApplyCylindricalMapping (int cylinderMatId, int capMatId); 
	CNEWTON_API void ApplyBoxMapping (int topMatId, int sideMatId, int frontMatId); 
	

	CNEWTON_API int GetMaterialIndex (void* const materialHandle) const; 
	CNEWTON_API int GetNextMaterialIndex (void* const materialHandle, int materialIndex) const; 

	int MaterialGetMaterial (void* const materialHandle, int materialIndex) const; 
	int MaterialGetIndexCount (void* const materialHandle, int materialIndex) const; 
	void MaterialGetIndexStream (void* const materialHandle, int materialIndex, int* const indexes) const; 

	CNEWTON_API void Polygonize ();
	CNEWTON_API void Triangulate ();

	protected:
	NewtonMesh* m_mesh;
};

#endif
