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

#ifndef __AFX_BROADPHASE_H_
#define __AFX_BROADPHASE_H_

#include "dgPhysicsStdafx.h"

class dgBody;
class dgWorld;
class dgContact;
class dgCollision;
class dgCollisionInstance;
class dgBroadphaseSyncDescriptor;

typedef void (dgApi *OnBodiesInAABB) (dgBody* body, void* const userData);
typedef void (dgApi *OnLeavingWorldAction) (dgBody* body, dgInt32 threadIndex);
typedef dgUnsigned32 (dgApi *OnRayPrecastAction) (const dgBody* const body, const dgCollisionInstance* const collision, void* const userData);
typedef dgFloat32 (dgApi *OnRayCastAction) (const dgBody* const body, const dgCollisionInstance* const collision, const dgVector& normal, void* collisionID, void* const userData, dgFloat32 intersetParam);


DG_MSC_VECTOR_ALIGMENT
struct dgLineBox
{
	dgVector m_l0;
	dgVector m_l1;
	dgVector m_boxL0;
	dgVector m_boxL1;
} DG_GCC_VECTOR_ALIGMENT;


class dgConvexCastReturnInfo
{
	public:
	dgFloat32 m_point[4];					// collision point in global space
	dgFloat32 m_normal[4];					// surface normal at collision point in global space
	dgFloat32 m_normalOnHitPoint[4];		// surface normal at the surface of the hit body, 
											// is the same as the normal calculate by a raycast passing by the hit point in the direction of the cast
	dgInt64  m_contaID;	                // collision ID at contact point
	const dgBody* m_hitBody;				// body hit at contact point
	dgFloat32 m_penetration;                // contact penetration at collision point
};


class dgBroadPhase
{
	public:
	DG_CLASS_ALLOCATOR(allocator);

	enum dgType
	{
		m_dynamic,
		m_static,
		m_hybrid,
	};

	class dgNode;
	

	dgBroadPhase(dgWorld* const world);
	~dgBroadPhase();

	void GetWorldSize (dgVector& p0, dgVector& p1) const;
//	void SetWorldSize (const dgVector& min, const dgVector& max);
	void RayCast (const dgVector& p0, const dgVector& p1, OnRayCastAction filter, OnRayPrecastAction prefilter, void* const userData) const;
	dgInt32 ConvexCast (dgCollisionInstance* const shape, const dgMatrix& p0, const dgVector& p1, dgFloat32& timetoImpact, OnRayPrecastAction prefilter, void* const userData, dgConvexCastReturnInfo* const info, dgInt32 maxContacts, dgInt32 threadIndex) const;
	void ForEachBodyInAABB (const dgVector& q0, const dgVector& q1, OnBodiesInAABB callback, void* const userData) const;

	dgInt32 GetBroadPhaseType () const;
	void SelectBroadPhaseType (dgInt32 algorthmType);

	protected:
	class dgFitnessList: public dgList <dgNode*>
	{
		public:
		dgFitnessList (dgMemoryAllocator* const allocator);
		dgFloat64 TotalCost () const;
	};

	void Add (dgBody* const body);
	void Remove (dgBody* const body);
	void InvalidateCache ();
	void UpdateContacts (dgFloat32 timestep);
	void UpdateBodyBroadphase(dgBody* const body, dgInt32 threadIndex);
	void UpdateBodyBroadphaseSimd(dgBody* const body, dgInt32 threadIndex);

	void ImproveFitness();
	void ImproveNodeFitness (dgNode* const node);
	void ImproveNodeFitnessSimd (dgNode* const node);
	dgNode* InsertNode (dgNode* const node);
	dgFloat32 CalculateSurfaceArea (const dgNode* const node0, const dgNode* const node1, dgVector& minBox, dgVector& maxBox) const;
	dgFloat32 CalculateSurfaceAreaSimd (const dgNode* const node0, const dgNode* const node1, dgSimd& minBox, dgSimd& maxBox) const;

	void AddPair (dgBody* const body0, dgBody* const body1, dgInt32 threadID);

	static void ForceAndToqueKernel (void* const descriptor, void* const worldContext, dgInt32 threadID);
	static void CollidingPairsKernel (void* const descriptor, void* const worldContext, dgInt32 threadID);
	static void UpdateContactsKernel (void* const descriptor, void* const worldContext, dgInt32 threadID);
	static void UpdateSoftBodyForcesKernel (void* const descriptor, void* const worldContext, dgInt32 threadID);
	
	void UpdateContactsBroadPhaseEnd ();
	void ApplyForceAndtorque (dgBroadphaseSyncDescriptor* const desctiptor, dgInt32 threadID);
	void ApplyDeformableForceAndtorque (dgBroadphaseSyncDescriptor* const desctiptor, dgInt32 threadID);
	void CalculatePairContacts (dgBroadphaseSyncDescriptor* const descriptor, dgInt32 threadID);
	void UpdateSoftBodyForcesKernel (dgBroadphaseSyncDescriptor* const descriptor, dgInt32 threadID);
	
	void SubmitPairsStatic (dgNode* const body, dgNode* const node, dgInt32 threadID);
	void FindCollidingPairsDynamics (dgBroadphaseSyncDescriptor* const desctiptor, dgInt32 threadID);
	void FindCollidingPairsStatic (dgBroadphaseSyncDescriptor* const desctiptor, dgInt32 threadID);
	void FindCollidingPairsHybrid (dgBroadphaseSyncDescriptor* const desctiptor, dgInt32 threadID);

	void KinematicBodyActivation (dgContact* const contatJoint) const;

	bool TestOverlaping (const dgBody* const body0, const dgBody* const body1) const;

	dgWorld* m_world;
	dgNode* m_rootNode;
	dgUnsigned32 m_lru;
	dgFitnessList m_fitness;
	dgType m_broadPhaseType;
	dgThread::dgCriticalSection m_contacJointLock;
	dgThread::dgCriticalSection m_criticalSectionLock;

	friend class dgBody;
	friend class dgWorld;
	friend class dgWorldDynamicUpdate;
};
#endif
