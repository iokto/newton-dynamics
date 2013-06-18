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

/****************************************************************************
*
*  File Name  : Bitmap.C
*  Visual C++ 4.0 base by Julio Jerez
*
****************************************************************************/
#ifndef __DG_AABB_POLYGON_SOUP_H_
#define __DG_AABB_POLYGON_SOUP_H_

#include "dgStdafx.h"
#include "dgPolygonSoupDatabase.h"


class dgPolygonSoupDatabaseBuilder;

//#define DG_NEW_AABB_TREE

#ifdef DG_NEW_AABB_TREE

class dgAABBPolygonSoup: public dgPolygonSoupDatabase
{
	public:
	DG_MSC_VECTOR_AVX_ALIGMENT
	class dgNode
	{
		public:
		dgNode ()
			:m_indexBox0(0)
			,m_indexBox1(0)
			,m_left(NULL)
			,m_right(NULL)
		{
		}

		dgInt32 m_indexBox0;
		dgInt32 m_indexBox1;
		dgNode* m_left;
		dgNode* m_right;
	} DG_GCC_VECTOR_AVX_ALIGMENT;

	class dgNodeBuilder;

//	dgInt32 GetIndexCount() const
//	{
//		return m_indexCount;
//	}

//	dgInt32* GetIndexPool() const
//	{
//		return m_indices;
//	}

//	virtual void GetAABB (dgVector& p0, dgVector& p1) const;
	virtual void Serialize (dgSerialize callback, void* const userData) const {dgAssert (0);}
	virtual void Deserialize (dgDeserialize callback, void* const userData) {dgAssert (0);}

	protected:
	dgAABBPolygonSoup ();
	virtual ~dgAABBPolygonSoup ();

	void Create (dgMemoryAllocator* const allocator, const dgPolygonSoupDatabaseBuilder& builder, bool optimizedBuild);
	virtual void ForAllSectors (const dgVector& minBox, const dgVector& maxBox, const dgVector& boxDistanceTravel, dgFloat32 m_maxT, dgAABBIntersectCallback callback, void* const context) const {dgAssert (0);}

	void* GetRootNode() const {dgAssert (0);return NULL;}
	void* GetBackNode(const void* const root) const {dgAssert (0);return NULL;}
	void* GetFrontNode(const void* const root) const {dgAssert (0);return NULL;}
	void GetNodeAABB(const void* const root, dgVector& p0, dgVector& p1) const {dgAssert (0);}
	virtual dgVector ForAllSectorsSupportVectex (const dgVector& dir) const {dgAssert (0); return dgVector(0.0f);}	
/*
	void CalculateAdjacendy ();
	
	dgFloat32 CalculateFaceMaxSize (dgTriplex* const vertex, dgInt32 indexCount, const dgInt32* const indexArray) const;

	virtual void ForAllSectorsRayHit (const dgFastRayTest& ray, dgRayIntersectCallback callback, void* const context) const;
*/	

	private:
//	static dgIntersectStatus CalculateManifoldFaceEdgeNormals (void* const context, const dgFloat32* const polygon, dgInt32 strideInBytes, const dgInt32* const indexArray, dgInt32 indexCount);
//	static dgIntersectStatus CalculateDisjointedFaceEdgeNormals (void* const context, const dgFloat32* const polygon, dgInt32 strideInBytes, const dgInt32* const indexArray, dgInt32 indexCount, dgFloat32 hitDistance);
//	static dgIntersectStatus CalculateAllFaceEdgeNormals (void* const context, const dgFloat32* const polygon, dgInt32 strideInBytes, const dgInt32* const indexArray, dgInt32 indexCount, dgFloat32 hitDistance);
//	dgInt32 m_nodesCount;
//	dgInt32 m_indexCount;
//	dgInt32 *m_indices;
//	void* m_aabb;
	void ImproveNodeFitness (dgNodeBuilder* const node) const;

	dgInt32 m_nodesCount;
	dgInt32 m_indexCount;
	dgNode* m_aabb;
	dgInt32* m_indices;
};

#else
class dgAABBPolygonSoup: public dgPolygonSoupDatabase
{
	public:
	dgInt32 GetIndexCount() const
	{
		return m_indexCount;
	}

	dgInt32* GetIndexPool() const
	{
		return m_indices;
	}

	virtual void GetAABB (dgVector& p0, dgVector& p1) const;
	virtual void Serialize (dgSerialize callback, void* const userData) const;
	virtual void Deserialize (dgDeserialize callback, void* const userData);

	protected:
	dgAABBPolygonSoup ();
	~dgAABBPolygonSoup ();

	void* GetRootNode() const;
	void* GetBackNode(const void* const root) const;
	void* GetFrontNode(const void* const root) const;
	void GetNodeAABB(const void* const root, dgVector& p0, dgVector& p1) const;

	void CalculateAdjacendy ();
	void Create (dgMemoryAllocator* const allocator, const dgPolygonSoupDatabaseBuilder& builder, bool optimizedBuild);
	dgFloat32 CalculateFaceMaxSize (dgTriplex* const vertex, dgInt32 indexCount, const dgInt32* const indexArray) const;

	virtual void ForAllSectors (const dgVector& minBox, const dgVector& maxBox, const dgVector& boxDistanceTravel, dgFloat32 m_maxT, dgAABBIntersectCallback callback, void* const context) const;
	virtual void ForAllSectorsRayHit (const dgFastRayTest& ray, dgRayIntersectCallback callback, void* const context) const;
	virtual dgVector ForAllSectorsSupportVectex (const dgVector& dir) const;	

	private:
	static dgIntersectStatus CalculateManifoldFaceEdgeNormals (void* const context, const dgFloat32* const polygon, dgInt32 strideInBytes, const dgInt32* const indexArray, dgInt32 indexCount);
	static dgIntersectStatus CalculateDisjointedFaceEdgeNormals (void* const context, const dgFloat32* const polygon, dgInt32 strideInBytes, const dgInt32* const indexArray, dgInt32 indexCount, dgFloat32 hitDistance);
	static dgIntersectStatus CalculateAllFaceEdgeNormals (void* const context, const dgFloat32* const polygon, dgInt32 strideInBytes, const dgInt32* const indexArray, dgInt32 indexCount, dgFloat32 hitDistance);

	dgInt32 m_nodesCount;
	dgInt32 m_indexCount;
	dgInt32 *m_indices;
	void* m_aabb;
};

#endif

#endif


