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

#ifndef __DGCOLLISION_HEIGHT_FIELD__
#define __DGCOLLISION_HEIGHT_FIELD__

#include "dgCollision.h"
#include "dgCollisionMesh.h"

class dgCollisionHeightField;
typedef dgFloat32 (*dgCollisionHeightFieldRayCastCallback) (const dgBody* const body, const dgCollisionHeightField* const heightFieldCollision, dgFloat32 interception, dgInt32 row, dgInt32 col, dgVector* const normal, int faceId, void* const usedData);


class dgCollisionHeightField: public dgCollisionMesh
{
	public:
	dgCollisionHeightField (dgWorld* const world, dgInt32 width, dgInt32 height, dgInt32 contructionMode, 
							const dgFloat32* const elevationMap, const dgInt8* const atributeMap, dgFloat32 horizontalScale);

	dgCollisionHeightField (dgWorld* const world, dgDeserialize deserialization, void* const userData);

	virtual ~dgCollisionHeightField(void);

	void SetCollisionRayCastCallback (dgCollisionHeightFieldRayCastCallback rayCastCallback);
	dgCollisionHeightFieldRayCastCallback GetDebugRayCastCallback() const { return m_userRayCastCallback;} 


	private:
	
	void AllocateVertex(dgWorld* const world, dgInt32 thread) const;
	void CalculateMinExtend2d (const dgVector& p0, const dgVector& p1, dgVector& boxP0, dgVector& boxP1) const;
	void CalculateMinExtend3d (const dgVector& p0, const dgVector& p1, dgVector& boxP0, dgVector& boxP1) const;
	dgFloat32 RayCastCell (const dgFastRayTest& ray, dgInt32 xIndex0, dgInt32 zIndex0, dgVector& normalOut) const;

	virtual void Serialize(dgSerialize callback, void* const userData) const;
	virtual dgFloat32 RayCast (const dgVector& localP0, const dgVector& localP1, dgContactPoint& contactOut, const dgBody* const body, void* const userData) const;
	virtual void GetCollidingFaces (dgPolygonMeshDesc* const data) const;

	virtual void GetCollisionInfo(dgCollisionInfo* const info) const;
	virtual void DebugCollision (const dgMatrix& matrixPtr, OnDebugCollisionMeshCallback callback, void* const userData) const;

	void GetVertexListIndexList (const dgVector& p0, const dgVector& p1, dgMeshVertexListIndexList &data) const;

	void GetLocalAABB (const dgVector& p0, const dgVector& p1, dgVector& boxP0, dgVector& boxP1) const;


	dgInt32 m_width;
	dgInt32 m_height;
	dgInt32 m_diagonalMode;
	dgInt8* m_atributeMap;
	dgFloat32* m_elevationMap;
	dgFloat32 m_horizontalScale;
	dgFloat32 m_horizontalScaleInv;
	dgCollisionHeightFieldRayCastCallback m_userRayCastCallback;

	dgVector m_minBox;
	dgVector m_maxBox;

	struct dgPerIntanceData
	{
		dgWorld* m_world;
		dgInt32 m_refCount;
		dgInt32 m_vertexCount[DG_MAX_THREADS_HIVE_COUNT];
		dgVector *m_vertex[DG_MAX_THREADS_HIVE_COUNT];
	};

	dgPerIntanceData* m_instanceData;
	friend class dgCollisionCompound;
};


#endif
