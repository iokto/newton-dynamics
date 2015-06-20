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
#include "dgConstraint.h"
#include "dgDynamicBody.h"
#include "dgAcyclicContainer.h"
#include "dgWorldDynamicUpdate.h"
#include "dgBilateralConstraint.h"


#define DG_ACYCLIC_STACK_SIZE		512


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////// 
dgInt32 dgAcyclicContainer::m_uniqueID = 10;

dgAcyclicContainer::dgAcyclicGraph::dgAcyclicGraph ()
	:m_body(NULL)
	,m_parent(NULL)
	,m_children(NULL)
{
	dgAssert (0);
}

dgAcyclicContainer::dgAcyclicGraph::dgAcyclicGraph (dgMemoryAllocator* const allocator, dgDynamicBody* const body, dgDynamicBody* const parent)
	:m_body(body)
	,m_parent(parent)
	,m_children (allocator)
{
}

dgAcyclicContainer::dgAcyclicGraph::~dgAcyclicGraph()
{
	for (dgList<dgAcyclicGraph*>::dgListNode* ptr = m_children.GetFirst(); ptr; ptr = ptr->GetNext()) {
		delete ptr->GetInfo();
	}
}

dgAcyclicContainer::dgAcyclicContainer (dgDynamicBody* const rootBody)
	:m_skeleton(rootBody->GetWorld()->GetAllocator(), rootBody, NULL)
	,m_upDownOrder(NULL)
	,m_downUpOrder(NULL)
	,m_id(m_uniqueID)
	,m_jointCount(0)
{
	m_uniqueID ++;
}

dgAcyclicContainer::~dgAcyclicContainer ()
{
	dgMemoryAllocator* const allocator = m_skeleton.m_body->GetWorld()->GetAllocator();
	if (m_upDownOrder) {
		allocator->Free(m_upDownOrder);
		allocator->Free(m_downUpOrder);
	}
}

dgInt32 dgAcyclicContainer::NodeCount () const
{
	dgInt32 count = 0;
	dgInt32 stack = 1;
	dgAcyclicGraph* stackPool[DG_ACYCLIC_STACK_SIZE];

	stackPool[0] = (dgAcyclicGraph*)&m_skeleton;
	while (stack) {
		stack--;
		dgAcyclicGraph* const node = stackPool[stack];
		count ++;
		for (dgList<dgAcyclicGraph*>::dgListNode* ptr = node->m_children.GetFirst(); ptr; ptr = ptr->GetNext()) {
			stackPool[stack] = ptr->GetInfo();
			stack++;
			dgAssert(stack < dgInt32(sizeof (stackPool) / sizeof (stackPool[0])));
		}
	}
	return count;
}

dgAcyclicContainer::dgAcyclicGraph* dgAcyclicContainer::FindNode (dgDynamicBody* const body) const
{
	dgInt32 stack = 1;
	dgAcyclicGraph* stackPool[DG_ACYCLIC_STACK_SIZE];

	stackPool[0] = (dgAcyclicGraph*) &m_skeleton;
	while (stack) {
		stack --;
		dgAcyclicGraph* const node = stackPool[stack];
		if (node->m_body == body) {
			return node;
		}

		for (dgList<dgAcyclicGraph*>::dgListNode* ptr = node->m_children.GetFirst(); ptr; ptr = ptr->GetNext()) {
			stackPool[stack] = ptr->GetInfo();
			stack ++;
			dgAssert (stack < dgInt32 (sizeof (stackPool) / sizeof (stackPool[0])));
		}
	}
	return NULL;
}

void dgAcyclicContainer::AddChild (dgBody* const parent, dgBody* const child)
{
	dgAssert (child);
	dgBody* const parent1 = parent ? parent : m_skeleton.m_body;
	dgAssert (child->GetType() == dgBody::m_dynamicBody);
	dgAssert (parent1->GetType() == dgBody::m_dynamicBody);
	AddChild ((dgDynamicBody*) parent1, (dgDynamicBody*) child);
}

void dgAcyclicContainer::AddChild (dgDynamicBody* const parent, dgDynamicBody* const child)
{
	dgAcyclicGraph* const parentNode = FindNode (parent);
	dgAssert (parentNode);
	dgMemoryAllocator* const allocator = child->GetWorld()->GetAllocator();
	dgAcyclicGraph* const node = new (allocator) dgAcyclicGraph (allocator, child, parent);
	parentNode->m_children.Append (node);
}


void dgAcyclicContainer::SortGraph (dgAcyclicGraph* const root, dgAcyclicGraph* const parent, const dgInt32 count, dgInt32& index)
{
	dgWorld* const world = m_skeleton.m_body->GetWorld();
	for (dgList<dgAcyclicGraph*>::dgListNode* node = root->m_children.GetFirst(); node; node = node->GetNext()) {
		SortGraph (node->GetInfo(), root, count, index);
	}

	dgBilateralConstraint* const joint = parent ? world->FindBilateralJoint (root->m_body, parent->m_body) : NULL;

	if (joint) {
		joint->m_priority = (m_id << DG_ACYCLIC_BIT_SHIFT_KEY) + index;
	}

	dgAssert ((count - index - 1) >= 0);
	m_upDownOrder[count - index - 1].m_Body = root->m_body;
	m_upDownOrder[count - index - 1].m_joint = joint;

	m_downUpOrder[index].m_Body = root->m_body;
	m_downUpOrder[index].m_joint = joint;
	index ++;
	dgAssert (index <= count);
}

void dgAcyclicContainer::Finalize ()
{
	dgInt32 count = NodeCount ();
	dgAssert (count >= 2);
	m_jointCount = count / 2; 

	dgMemoryAllocator* const allocator = m_skeleton.m_body->GetWorld()->GetAllocator();
	m_upDownOrder = (dgAcyclicJointBodyPair*) allocator->Malloc(count * sizeof (dgAcyclicJointBodyPair));
	m_downUpOrder = (dgAcyclicJointBodyPair*) allocator->Malloc(count * sizeof (dgAcyclicJointBodyPair));

	dgInt32 index = 0;
	SortGraph (&m_skeleton, NULL, count, index);
}

void dgAcyclicContainer::CalculateJointForce (dgJointInfo* const jointInfoArray, const dgBodyInfo* const bodyArray, dgJacobian* const internalForces, dgJacobianMatrixElement* const matrixRow) const
{

}