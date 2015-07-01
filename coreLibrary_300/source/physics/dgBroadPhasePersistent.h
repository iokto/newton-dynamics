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

#ifndef __AFX_BROADPHASE_PERSINTENT_H_
#define __AFX_BROADPHASE_PERSINTENT_H_

#include "dgPhysicsStdafx.h"
#include "dgBroadPhase.h"


class dgBroadPhasePersistent: public dgBroadPhase
{
	public:
	DG_CLASS_ALLOCATOR(allocator);

	dgBroadPhasePersistent(dgWorld* const world);
	virtual ~dgBroadPhasePersistent();

	protected:
	virtual dgInt32 GetType() const;
	virtual void Add(dgBody* const body);
	virtual void Remove(dgBody* const body);
	virtual void InvalidateCache();
	virtual dgBroadPhaseNodeAggegate* CreateAggegate();
	virtual void DestroyAggregate(dgBroadPhaseNodeAggegate* const aggregate);

	virtual void ResetEntropy();
	virtual void UpdateFitness();
	virtual void ForEachBodyInAABB(const dgVector& minBox, const dgVector& maxBox, OnBodiesInAABB callback, void* const userData) const;
	virtual void RayCast(const dgVector& p0, const dgVector& p1, OnRayCastAction filter, OnRayPrecastAction prefilter, void* const userData) const;
	virtual void ConvexRayCast(dgCollisionInstance* const shape, const dgMatrix& matrix, const dgVector& target, OnRayCastAction filter, OnRayPrecastAction prefilter, void* const userData, dgInt32 threadId) const;
	virtual dgInt32 ConvexCast(dgCollisionInstance* const shape, const dgMatrix& matrix, const dgVector& target, dgFloat32& timeToImpact, OnRayPrecastAction prefilter, void* const userData, dgConvexCastReturnInfo* const info, dgInt32 maxContacts, dgInt32 threadIndex) const;

	virtual void CheckStaticDynamic(dgBody* const body, dgFloat32 mass);
	virtual void ScanForContactJoints(dgBroadphaseSyncDescriptor& syncPoints);
	virtual void FindCollidingPairs (dgBroadphaseSyncDescriptor* const descriptor, dgBodyMasterList::dgListNode* node, dgInt32 threadID);

	

	dgFloat64 m_staticEntropy;
	dgFloat64 m_dynamicsEntropy;
	dgFitnessList m_staticFitness;
	dgFitnessList m_dynamicsFitness;
	bool m_staticNeedsUpdate;
};
#endif
