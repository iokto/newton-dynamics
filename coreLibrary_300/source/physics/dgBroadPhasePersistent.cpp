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

#include "dgPhysicsStdafx.h"
#include "dgBody.h"
#include "dgWorld.h"
#include "dgCollisionInstance.h"
#include "dgBroadPhasePersistent.h"

dgBroadPhasePersistent::dgBroadPhasePersistent(dgWorld* const world)
	:dgBroadPhase(world)
	,m_staticEntropy(dgFloat32 (0.0f))
	,m_dynamicsEntropy(dgFloat32 (0.0f))
	,m_staticFitness(world->GetAllocator())
	,m_dynamicsFitness(world->GetAllocator())
	,m_staticNeedsUpdate(true)
{
	dgAssert (0);
//	m_rootNode = new (world->GetAllocator()) dgBroadPhaseNode(NULL);
}

dgBroadPhasePersistent::~dgBroadPhasePersistent()
{
	dgAssert (0);
	delete m_rootNode;
}

dgInt32 dgBroadPhasePersistent::GetType() const
{
	return dgWorld::m_persistentBroadphase;
}

void dgBroadPhasePersistent::CheckStaticDynamic(dgBody* const body, dgFloat32 mass)
{
	dgBroadPhaseNode* const node = body->GetBroadPhase();
	if (node) {
		dgVector temp (body->GetInvMass());
		if (((mass != dgFloat32 (0.0f)) && (temp.m_w == dgFloat32 (0.0f))) || ((mass == dgFloat32 (0.0f)) && (temp.m_w != dgFloat32 (0.0f)))) {
			Remove(body);
			body->SetInvMass (dgVector (dgFloat32 (1.0f)));			
			Add(body);
			body->SetInvMass (temp);
		}
	}
}

void dgBroadPhasePersistent::Add(dgBody* const body)
{
dgAssert (0);
/*
	dgAssert (!body->GetCollision()->IsType (dgCollision::dgCollisionNull_RTTI));
	if (body->GetInvMass().m_w == dgFloat32(0.0f)) {
		m_staticNeedsUpdate = true;
		if (m_rootNode->m_right) {
			dgBroadPhaseNode* const node = InsertNode(m_rootNode->m_right, new (m_world->GetAllocator()) dgBroadPhaseNode(body));
			node->m_fitnessNode = m_staticFitness.Append(node);
		} else {
			m_rootNode->m_right = new (m_world->GetAllocator()) dgBroadPhaseNode(body);
			m_rootNode->m_right->m_parent = m_rootNode;
		}
	} else {
		if (m_rootNode->m_left) {
			dgBroadPhaseNode* const node = InsertNode(m_rootNode->m_left, new (m_world->GetAllocator()) dgBroadPhaseNode(body));
			node->m_fitnessNode = m_dynamicsFitness.Append(node);
		} else {
			m_rootNode->m_left = new (m_world->GetAllocator()) dgBroadPhaseNode(body);
			m_rootNode->m_left->m_parent = m_rootNode;
		}
	}
*/
}


void dgBroadPhasePersistent::Remove(dgBody* const body)
{
dgAssert (0);
/*
	dgBroadPhaseNode* const node = body->GetBroadPhase();
	if (node) {
		dgAssert(node->m_parent);
		dgAssert(!node->m_fitnessNode);

		dgBroadPhaseNode* const grandParent = node->m_parent->m_parent;
		if (grandParent) {
			if (grandParent->m_left == node->m_parent) {
				if (node->m_parent->m_right == node) {
					grandParent->m_left = node->m_parent->m_left;
					node->m_parent->m_left->m_parent = grandParent;
					node->m_parent->m_left = NULL;
					node->m_parent->m_parent = NULL;
				} else {
					grandParent->m_left = node->m_parent->m_right;
					node->m_parent->m_right->m_parent = grandParent;
					node->m_parent->m_right = NULL;
					node->m_parent->m_parent = NULL;
				}
			} else {
				if (node->m_parent->m_right == node) {
					grandParent->m_right = node->m_parent->m_left;
					node->m_parent->m_left->m_parent = grandParent;
					node->m_parent->m_left = NULL;
					node->m_parent->m_parent = NULL;
				} else {
					grandParent->m_right = node->m_parent->m_right;
					node->m_parent->m_right->m_parent = grandParent;
					node->m_parent->m_right = NULL;
					node->m_parent->m_parent = NULL;
				}
			}

			dgAssert(node->m_parent->m_fitnessNode);
			if (body->GetInvMass().m_w == dgFloat32(0.0f)) {
				m_staticNeedsUpdate = true;
				m_staticFitness.Remove(node->m_parent->m_fitnessNode);
			} else {
				m_dynamicsFitness.Remove(node->m_parent->m_fitnessNode);
			}
			delete node->m_parent;
		} else {
			if (node->m_parent->m_right == node) {
				m_rootNode->m_right = NULL;
			} else {
				m_rootNode->m_left = NULL;
			}
			node->m_parent = NULL;
			delete node;
		}
	}
*/
}


dgBroadPhaseNodeAggegate* dgBroadPhasePersistent::CreateAggegate()
{
	dgAssert (0);
	return NULL;
}

void dgBroadPhasePersistent::DestroyAggregate(dgBroadPhaseNodeAggegate* const aggregate)
{
	dgAssert (0);
}


void dgBroadPhasePersistent::ResetEntropy()
{
	m_staticNeedsUpdate = true;
	m_staticEntropy = dgFloat32(0.0f);
	m_dynamicsEntropy = dgFloat32(0.0f);
}


void dgBroadPhasePersistent::InvalidateCache()
{
dgAssert (0);
/*
	ResetEntropy();
	m_staticNeedsUpdate = false;
	ImproveFitness(m_staticFitness, m_staticEntropy, &m_rootNode->m_right);
	ImproveFitness(m_dynamicsFitness, m_dynamicsEntropy, &m_rootNode->m_left);
*/
}

void dgBroadPhasePersistent::UpdateFitness()
{
dgAssert (0);
/*
	if (m_staticNeedsUpdate) {
		m_staticNeedsUpdate = false;
		ImproveFitness(m_staticFitness, m_staticEntropy, &m_rootNode->m_right);
	}
	ImproveFitness(m_dynamicsFitness, m_dynamicsEntropy, &m_rootNode->m_left);
*/
}

void dgBroadPhasePersistent::ForEachBodyInAABB(const dgVector& minBox, const dgVector& maxBox, OnBodiesInAABB callback, void* const userData) const
{
dgAssert (0);
/*
	const dgBroadPhaseNode* stackPool[DG_BROADPHASE_MAX_STACK_DEPTH];

	dgInt32 stack = 0;
	if (m_rootNode->m_left) {
		stackPool[stack] = m_rootNode->m_left;
		stack++;
	}

	if (m_rootNode->m_right) {
		stackPool[stack] = m_rootNode->m_right;
		stack++;
	}
	dgBroadPhase::ForEachBodyInAABB(stackPool, stack, minBox, maxBox, callback, userData);
*/
}

void dgBroadPhasePersistent::RayCast(const dgVector& l0, const dgVector& l1, OnRayCastAction filter, OnRayPrecastAction prefilter, void* const userData) const
{
dgAssert (0);
/*
	if (filter && (m_rootNode->m_left || m_rootNode->m_right)) {
		dgVector segment(l1 - l0);
		dgFloat32 dist2 = segment % segment;
		if (dist2 > dgFloat32(1.0e-8f)) {

			dgFloat32 distance[DG_BROADPHASE_MAX_STACK_DEPTH];
			const dgBroadPhaseNode* stackPool[DG_BROADPHASE_MAX_STACK_DEPTH];

			dgFastRayTest ray(l0, l1);

			dgInt32 stack = 0;
			if (m_rootNode->m_left) {
				stackPool[stack] = m_rootNode->m_left;
				distance[stack] = ray.BoxIntersect(m_rootNode->m_left->m_minBox, m_rootNode->m_left->m_maxBox);
				stack++;
			}
			if (m_rootNode->m_right) {
				stackPool[stack] = m_rootNode->m_right;
				distance[stack] = ray.BoxIntersect(m_rootNode->m_right->m_minBox, m_rootNode->m_right->m_maxBox);
				stack++;
			}
			if (stack == 2) {
				if (distance[0] < distance[1]) {
					dgSwap(distance[0], distance[1]);
					dgSwap(stackPool[0], stackPool[1]);
				}
			}

			dgBroadPhase::RayCast(stackPool, distance, stack, l0, l1, ray, filter, prefilter, userData);
		}
	}
*/
}

void dgBroadPhasePersistent::ConvexRayCast(dgCollisionInstance* const shape, const dgMatrix& matrix, const dgVector& target, OnRayCastAction filter, OnRayPrecastAction prefilter, void* const userData, dgInt32 threadId) const
{
dgAssert (0);
/*
	if (filter && m_rootNode && shape->IsType(dgCollision::dgCollisionConvexShape_RTTI)) {
		dgVector boxP0;
		dgVector boxP1;
		shape->CalcAABB(shape->GetLocalMatrix() * matrix, boxP0, boxP1);

		//dgInt32 stack = 1;
		dgFloat32 distance[DG_COMPOUND_STACK_DEPTH];
		const dgBroadPhaseNode* stackPool[DG_BROADPHASE_MAX_STACK_DEPTH];

		dgVector velocA((target - matrix.m_posit) & dgVector::m_triplexMask);
		dgFastRayTest ray(dgVector(dgFloat32(0.0f)), velocA);

		//dgVector minBox(m_rootNode->m_minBox - boxP1);
		//dgVector maxBox(m_rootNode->m_maxBox - boxP0);
		//stackPool[0] = m_rootNode;
		//distance[0] = ray.BoxIntersect(minBox, maxBox);

		dgInt32 stack = 0;
		if (m_rootNode->m_left) {
			dgVector minBox(m_rootNode->m_left->m_minBox - boxP1);
			dgVector maxBox(m_rootNode->m_left->m_maxBox - boxP0);
			stackPool[stack] = m_rootNode->m_left;
			distance[stack] = ray.BoxIntersect(minBox, maxBox);
			stack++;
		}
		if (m_rootNode->m_right) {
			dgVector minBox(m_rootNode->m_right->m_minBox - boxP1);
			dgVector maxBox(m_rootNode->m_right->m_maxBox - boxP0);

			stackPool[stack] = m_rootNode->m_right;
			distance[stack] = ray.BoxIntersect(minBox, maxBox);
			stack++;
		}
		if (stack == 2) {
			if (distance[0] < distance[1]) {
				dgSwap(distance[0], distance[1]);
				dgSwap(stackPool[0], stackPool[1]);
			}
		}

		dgBroadPhase::ConvexRayCast(stackPool, distance, stack, velocA, ray, shape, matrix, target, filter, prefilter, userData, threadId);
	}
*/
}

dgInt32 dgBroadPhasePersistent::ConvexCast(dgCollisionInstance* const shape, const dgMatrix& matrix, const dgVector& target, dgFloat32& timeToImpact, OnRayPrecastAction prefilter, void* const userData, dgConvexCastReturnInfo* const info, dgInt32 maxContacts, dgInt32 threadIndex) const
{
dgAssert (0);
return 0;
/*
	dgInt32 totalCount = 0;
	if (m_rootNode) {
		dgVector boxP0;
		dgVector boxP1;
		dgAssert(matrix.TestOrthogonal());
		shape->CalcAABB(matrix, boxP0, boxP1);

		dgFloat32 distance[DG_BROADPHASE_MAX_STACK_DEPTH];
		const dgBroadPhaseNode* stackPool[DG_BROADPHASE_MAX_STACK_DEPTH];

		dgVector velocA((target - matrix.m_posit) & dgVector::m_triplexMask);
		dgVector velocB(dgFloat32(0.0f));
		dgFastRayTest ray(dgVector(dgFloat32(0.0f)), velocA);

		//dgVector minBox(m_rootNode->m_minBox - boxP1);
		//dgVector maxBox(m_rootNode->m_maxBox - boxP0);
		//stackPool[0] = m_rootNode;
		//distance[0] = ray.BoxIntersect(minBox, maxBox);

		dgInt32 stack = 0;
		if (m_rootNode->m_left) {
			dgVector minBox(m_rootNode->m_left->m_minBox - boxP1);
			dgVector maxBox(m_rootNode->m_left->m_maxBox - boxP0);
			stackPool[stack] = m_rootNode->m_left;
			distance[stack] = ray.BoxIntersect(minBox, maxBox);
			stack++;
		}
		if (m_rootNode->m_right) {
			dgVector minBox(m_rootNode->m_right->m_minBox - boxP1);
			dgVector maxBox(m_rootNode->m_right->m_maxBox - boxP0);

			stackPool[stack] = m_rootNode->m_right;
			distance[stack] = ray.BoxIntersect(minBox, maxBox);
			stack++;
		}
		if (stack == 2) {
			if (distance[0] < distance[1]) {
				dgSwap(distance[0], distance[1]);
				dgSwap(stackPool[0], stackPool[1]);
			}
		}

		totalCount = dgBroadPhase::ConvexCast(stackPool, distance, 2, velocA, velocB, ray, shape, matrix, target, timeToImpact, prefilter, userData, info, maxContacts, threadIndex);
	}

	return totalCount;
*/
}

//void dgBroadPhasePersistent::FindCollidingPairs(dgBroadphaseSyncDescriptor* const descriptor, dgBodyMasterList::dgListNode* node, dgInt32 threadID)
void dgBroadPhasePersistent::FindCollidingPairs (dgBroadphaseSyncDescriptor* const descriptor, dgList<dgBroadPhaseNode*>::dgListNode* const nodePtr, dgInt32 threadID)
{
dgAssert (0);
/*
	const dgFloat32 timestep = descriptor->m_timestep;

	const dgInt32 threadCount = descriptor->m_world->GetThreadCount();
	while (node) {
		dgBody* const body = node->GetInfo().GetBody();
		dgBroadPhaseNode* const broadPhaseNode = body->GetBroadPhase();
		if (broadPhaseNode) {
			for (dgBroadPhaseNode* ptr = broadPhaseNode; ptr->m_parent; ptr = ptr->m_parent) {
				dgBroadPhaseNode* const sibling = ptr->m_parent->m_right;
				if (sibling && (sibling != ptr)) {
					SubmitPairs(broadPhaseNode, sibling, timestep, threadID);
				}
			}
		}
		for (dgInt32 i = 0; i < threadCount; i++) {
			node = (node && (node->GetPrev()->GetInfo().GetBody()->GetInvMass().m_w != dgFloat32(0.0f))) ? node->GetPrev() : NULL;
		}
	}
*/
}


void dgBroadPhasePersistent::ScanForContactJoints(dgBroadphaseSyncDescriptor& syncPoints)
{
	dgInt32 threadsCount = m_world->GetThreadCount();
	const dgBodyMasterList* const masterList = m_world;
	dgBodyMasterList::dgListNode* node = (masterList->GetLast()->GetInfo().GetBody()->GetInvMass().m_w != dgFloat32(0.0f)) ? masterList->GetLast() : NULL;
	for (dgInt32 i = 0; i < threadsCount; i++) {
		m_world->QueueJob(CollidingPairsKernel, &syncPoints, node);
		node = (node && (node->GetPrev()->GetInfo().GetBody()->GetInvMass().m_w != NULL)) ? node->GetPrev() : NULL;
	}
	m_world->SynchronizationBarrier();
}