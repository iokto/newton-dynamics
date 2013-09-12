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

#ifndef __dgPolyhedra__
#define __dgPolyhedra__

#include "dgStdafx.h"
#include "dgList.h"
#include "dgTree.h"
#include "dgHeap.h"
#include "dgDebug.h"




class dgEdge;
class dgPlane;
class dgSphere;
class dgMatrix;
class dgPolyhedra;
class dgVertexCollapseVertexMetric;

typedef dgInt64 dgEdgeKey;


DG_MSC_VECTOR_ALIGMENT
class dgEdge
{
	public:
	dgEdge ();
	dgEdge (dgInt32 vertex, dgInt32 face, dgUnsigned64 userdata = 0);

	dgInt32 m_incidentVertex;
	dgInt32 m_incidentFace;
	dgUnsigned64 m_userData;
	dgEdge* m_next;
	dgEdge* m_prev;
	dgEdge* m_twin;
	dgInt32 m_mark;
} DG_GCC_VECTOR_ALIGMENT;


class dgPolyhedra: public dgTree <dgEdge, dgEdgeKey>
{
	public:

	struct dgPairKey
	{
		dgPairKey ();
		dgPairKey (dgInt64 val);
		dgPairKey (dgInt32 v0, dgInt32 v1);
		dgInt64 GetVal () const; 
		dgInt32 GetLowKey () const;
		dgInt32 GetHighKey () const;

		private:
		dgUnsigned64 m_key;
	};

	dgPolyhedra (dgMemoryAllocator* const allocator);
	dgPolyhedra (const dgPolyhedra &polyhedra);
	virtual ~dgPolyhedra();

	void BeginFace();
	dgEdge* AddFace (dgInt32 v0, dgInt32 v1, dgInt32 v2);
	dgEdge* AddFace (dgInt32 count, const dgInt32* const index);
	dgEdge* AddFace (dgInt32 count, const dgInt32* const index, const dgInt64* const userdata);
	void EndFace ();
	void DeleteFace(dgEdge* const edge);

	dgInt32 GetFaceCount() const;
	dgInt32 GetEdgeCount() const;
	dgInt32 GetLastVertexIndex() const;

	dgInt32 IncLRU() const;
	void SetLRU(dgInt32 lru) const;

	dgEdge* FindEdge (dgInt32 v0, dgInt32 v1) const;
	dgTreeNode* FindEdgeNode (dgInt32 v0, dgInt32 v1) const;

	dgEdge* AddHalfEdge (dgInt32 v0, dgInt32 v1);
	void DeleteEdge (dgEdge* const edge);
	void DeleteEdge (dgInt32 v0, dgInt32 v1);

	dgEdge* ConnectVertex (dgEdge* const e0, dgEdge* const e1);
	
	bool FlipEdge (dgEdge* const edge);
	dgEdge* SpliteEdge (dgInt32 newIndex, dgEdge* const edge);
	dgBigVector FaceNormal (const dgEdge* const face, const dgFloat64* const vertex, dgInt32 strideInBytes) const;

	void BeginConectedSurface() const;
	bool GetConectedSurface (dgPolyhedra &polyhedra) const;
	void EndConectedSurface() const;

	dgSphere CalculateSphere (const dgFloat64* const vertex, dgInt32 strideInBytes, const dgMatrix* const basis = NULL) const;
	void ChangeEdgeIncidentVertex (dgEdge* const edge, dgInt32 newIndex);	
	void DeleteDegenerateFaces (const dgFloat64* const pool, dgInt32 dstStrideInBytes, dgFloat64 minArea);

	void Optimize (const dgFloat64* const pool, dgInt32 strideInBytes, dgFloat64 tol, dgInt32 maxFaceCount = 1<<28);
	void Triangulate (const dgFloat64* const vertex, dgInt32 strideInBytes, dgPolyhedra* const leftOversOut);
	void ConvexPartition (const dgFloat64* const vertex, dgInt32 strideInBytes, dgPolyhedra* const leftOversOut);
	dgEdge* CollapseEdge(dgEdge* const edge);
	void PolygonizeFace (dgEdge* face, const dgFloat64* const pool, dgInt32 stride);

	private:
	void RefineTriangulation (const dgFloat64* const vertex, dgInt32 stride);
	void RefineTriangulation (const dgFloat64* const vertex, dgInt32 stride, dgBigVector* const normal, dgInt32 perimeterCount, dgEdge** const perimeter);
	void OptimizeTriangulation (const dgFloat64* const vertex, dgInt32 strideInBytes);
	void MarkAdjacentCoplanarFaces (dgPolyhedra& polyhedraOut, dgEdge* const face, const dgFloat64* const pool, dgInt32 strideInBytes);
	dgEdge* FindEarTip (dgEdge* const face, const dgFloat64* const pool, dgInt32 stride, dgDownHeap<dgEdge*, dgFloat64>& heap, const dgBigVector &normal) const;
	dgEdge* TriangulateFace (dgEdge* const face, const dgFloat64* const pool, dgInt32 stride, dgDownHeap<dgEdge*, dgFloat64>& heap, dgBigVector* const faceNormalOut);
	

	void RemoveHalfEdge (dgEdge* const edge);
	dgEdge* OptimizeCollapseEdge (dgEdge* const edge);
	bool IsOkToCollapse (const dgBigVector* const pool, dgEdge* const edge) const;
	dgFloat64 EdgePenalty (const dgBigVector* const pool, dgEdge* const edge, dgFloat64 dist) const;
	dgBigPlane EdgePlane (dgInt32 i0, dgInt32 i1, dgInt32 i2, const dgBigVector* const pool) const;
	void CalculateAllMetrics (dgVertexCollapseVertexMetric* const table, const dgBigVector* const pool) const;
	void CalculateVertexMetrics (dgVertexCollapseVertexMetric* const table, const dgBigVector* const pool, dgEdge* const edge) const;
	

	mutable dgInt32 m_baseMark;
	mutable dgInt32 m_edgeMark;
	mutable dgInt32 m_faceSecuence;

	friend class dgPolyhedraDescriptor;
	
};



inline dgEdge::dgEdge ()
{
}

inline dgEdge::dgEdge (dgInt32 vertex, dgInt32 face, dgUnsigned64 userdata)
	:m_incidentVertex(vertex)
	,m_incidentFace(face)
	,m_userData(userdata)
	,m_next(NULL)
	,m_prev(NULL)
	,m_twin(NULL)
	,m_mark(0)
{
}


inline dgPolyhedra::dgPairKey::dgPairKey ()
{
}

inline dgPolyhedra::dgPairKey::dgPairKey (dgInt64 val)
	:m_key(dgUnsigned64 (val))
{
}

inline dgPolyhedra::dgPairKey::dgPairKey (dgInt32 v0, dgInt32 v1)
	:m_key (dgUnsigned64 ((dgInt64 (v0) << 32) | v1))
{
}

inline dgInt64 dgPolyhedra::dgPairKey::GetVal () const 
{
	return dgInt64 (m_key);
}

inline dgInt32 dgPolyhedra::dgPairKey::GetLowKey () const 
{
	return dgInt32 (m_key>>32);
}

inline dgInt32 dgPolyhedra::dgPairKey::GetHighKey () const 
{
	return dgInt32 (m_key & 0xffffffff);
}

inline void dgPolyhedra::BeginFace ()
{
}

inline dgEdge* dgPolyhedra::AddFace (dgInt32 count, const dgInt32* const index) 
{
	return AddFace (count, index, NULL);
}

inline dgEdge* dgPolyhedra::AddFace (dgInt32 v0, dgInt32 v1, dgInt32 v2)
{
	dgInt32 vertex[3];

	vertex [0] = v0;
	vertex [1] = v1;
	vertex [2] = v2;
	return AddFace (3, vertex, NULL);
}

inline dgInt32 dgPolyhedra::GetEdgeCount() const
{
#ifdef _DEBUG
	dgInt32 edgeCount = 0;
	Iterator iter(*this);
	for (iter.Begin(); iter; iter ++) {
		edgeCount ++;
	}
	dgAssert (edgeCount == GetCount());;
#endif
	return GetCount();
}

inline dgInt32 dgPolyhedra::GetLastVertexIndex() const
{
	dgInt32 maxVertexIndex = -1;
	Iterator iter(*this);
	for (iter.Begin(); iter; iter ++) {
		const dgEdge* const edge = &(*iter);
		if (edge->m_incidentVertex > maxVertexIndex) {
			maxVertexIndex = edge->m_incidentVertex;
		}
	}
	return maxVertexIndex + 1;
}


inline dgInt32 dgPolyhedra::IncLRU() const
{	
	m_edgeMark ++;
	dgAssert (m_edgeMark < 0x7fffffff);
	return m_edgeMark;
}

inline void dgPolyhedra::SetLRU(dgInt32 lru) const
{
	if (lru > m_edgeMark) {
		m_edgeMark = lru;
	}
}

inline void dgPolyhedra::BeginConectedSurface() const
{
	m_baseMark = IncLRU();
}

inline void dgPolyhedra::EndConectedSurface() const
{
}


inline dgPolyhedra::dgTreeNode* dgPolyhedra::FindEdgeNode (dgInt32 i0, dgInt32 i1) const
{
	dgPairKey key (i0, i1);
	return Find (key.GetVal());
}


inline dgEdge *dgPolyhedra::FindEdge (dgInt32 i0, dgInt32 i1) const
{
	//	dgTreeNode *node;
	//	dgPairKey key (i0, i1);
	//	node = Find (key.GetVal());
	//	return node ? &node->GetInfo() : NULL;
	dgTreeNode* const node = FindEdgeNode (i0, i1);
	return node ? &node->GetInfo() : NULL;
}

inline void dgPolyhedra::DeleteEdge (dgInt32 v0, dgInt32 v1)
{
	dgPairKey pairKey (v0, v1);
	dgTreeNode* const node = Find(pairKey.GetVal());
	dgEdge* const edge = node ? &node->GetInfo() : NULL;
	if (!edge) {
		return;
	}
	DeleteEdge (edge);
}


#endif

