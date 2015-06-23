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



class dgAcyclicContainer::dgAcyclicMatrix
{
	public:
	dgMatrix m_d[2][2];
};


class dgAcyclicContainer::dgAcyclicNodeInfo
{
	public:
	dgAcyclicMatrix m_matrix;
	//dgAcyclicMatrix m_invMatrix;
	dgAcyclicMatrix m_jacobian;
	dgJacobian m_vector;
};


class dgAcyclicContainer::dgAcyclicGraph
{
	public:
	DG_CLASS_ALLOCATOR(allocator)

	dgAcyclicGraph(dgMemoryAllocator* const allocator, dgAcyclicGraph* const parent)
		:m_parent(parent)
		,m_block(NULL)
		,m_children(allocator)
		,m_index(0)
	{
	}

	virtual ~dgAcyclicGraph()
	{
		for (dgList<dgAcyclicGraph*>::dgListNode* ptr = m_children.GetFirst(); ptr; ptr = ptr->GetNext()) {
			delete ptr->GetInfo();
		}
	}

	void AddChild(dgAcyclicGraph* const child)
	{
		m_children.Append(child);
	}

	virtual dgDynamicBody* GetBody() const 
	{
		return NULL;
	}

	virtual dgBilateralConstraint* GetJoint() const 
	{
		return NULL;
	}

	virtual void SetPriority(dgUnsigned32 priority) const 
	{
	}


	virtual void Init(dgAcyclicNodeInfo* const buffer, dgJointInfo* const jointInfoArray, dgJacobianMatrixElement* const matrixRow) 
	{
		m_block = buffer;
	}

	void AccumulateMatrix()
	{
		dgAssert (0);
	}

	void CalculateInverse()
	{
		dgAssert (0);
	}

	void AccumulateJacobian()
	{
		dgAssert (0);
	}

	dgAcyclicGraph* m_parent;
	dgAcyclicNodeInfo* m_block;
	dgList<dgAcyclicGraph*> m_children;
	dgInt32 m_index;
};

class dgAcyclicContainer::dgAcyclicGraphMassNode: public dgAcyclicGraph
{
	public:
	dgAcyclicGraphMassNode(dgMemoryAllocator* const allocator, dgDynamicBody* const body, dgAcyclicGraph* const parent)
		:dgAcyclicGraph (allocator, parent)
		,m_body(body)
	{
		if (parent) {
			parent->AddChild (this);
		}
	}

	virtual dgDynamicBody* GetBody() const 
	{
		return m_body;
	}

	virtual void Init(dgAcyclicNodeInfo* const buffer, dgJointInfo* const jointInfoArray, dgJacobianMatrixElement* const matrixRow)
	{
		dgAcyclicGraph::Init(buffer, jointInfoArray, matrixRow);
		dgVector mass = m_body->GetMass();
		m_block->m_matrix.m_d[0][0][0] = mass & dgVector::m_xMask;
		m_block->m_matrix.m_d[0][0][1] = mass & dgVector::m_xMask;
		m_block->m_matrix.m_d[0][0][2] = mass & dgVector::m_xMask;
		m_block->m_matrix.m_d[0][0][3] = dgVector::m_wOne;
		m_block->m_matrix.m_d[0][1] = dgGetZeroMatrix();
		m_block->m_matrix.m_d[1][0] = dgGetZeroMatrix();
		m_block->m_matrix.m_d[2][2] = m_body->CalculateInertiaMatrix();

		if (m_parent) {
			dgBilateralConstraint* const joint = m_parent->GetJoint();
			dgJointInfo* const jointInfo = &jointInfoArray[joint->m_index];
			dgAssert (jointInfo->m_joint == joint);
			dgInt32 rowsCount = jointInfo->m_pairCount;
			dgInt32 rowsCount1 = rowsCount > 4 ? 4 : rowsCount;

			m_block->m_matrix.m_d[0][0] = dgGetIdentityMatrix();
			m_block->m_matrix.m_d[0][1] = dgGetIdentityMatrix();
			m_block->m_matrix.m_d[1][0] = dgGetIdentityMatrix();
			m_block->m_matrix.m_d[1][1] = dgGetIdentityMatrix();
			for (dgInt32 i = 0; i < rowsCount1; i ++) {
				dgInt32 index = jointInfo->m_pairStart + i;
				m_block->m_matrix.m_d[0][0][i] = matrixRow[index].m_JMinv.m_jacobianM0.m_linear;
				m_block->m_matrix.m_d[1][0][i] = matrixRow[index].m_JMinv.m_jacobianM0.m_angular;
			}
			//m_block->m_matrix.m_d[0][0] = m_block->m_matrix.m_d[0][0].Transpose4X4();
			//m_block->m_matrix.m_d[1][0] = m_block->m_matrix.m_d[1][0].Transpose4X4();

			for (dgInt32 i = rowsCount1; i < rowsCount; i++) {
				dgInt32 index = jointInfo->m_pairStart + i;
				m_block->m_matrix.m_d[0][1][i] = matrixRow[index].m_JMinv.m_jacobianM0.m_linear;
				m_block->m_matrix.m_d[1][1][i] = matrixRow[index].m_JMinv.m_jacobianM0.m_angular;
			}
			//m_block->m_matrix.m_d[0][1] = m_block->m_matrix.m_d[0][0].Transpose4X4();
			//m_block->m_matrix.m_d[1][1] = m_block->m_matrix.m_d[1][0].Transpose4X4();
		}
	}

	dgDynamicBody* m_body;
};

class dgAcyclicContainer::dgAcyclicGraphJointNode: public dgAcyclicGraph
{
	public:
	dgAcyclicGraphJointNode(dgMemoryAllocator* const allocator, dgBilateralConstraint* const Joint, dgAcyclicGraph* const parent)
		:dgAcyclicGraph(allocator, parent)
		,m_joint(Joint)
	{
		parent->AddChild (this);
	}

	virtual dgBilateralConstraint* GetJoint() const
	{
		return m_joint;
	}

	virtual void SetPriority(dgUnsigned32 priority) const
	{
		m_joint->m_priority = priority;
	}

	virtual void Init(dgAcyclicNodeInfo* const buffer, dgJointInfo* const jointInfoArray, dgJacobianMatrixElement* const matrixRow) 
	{
		dgAcyclicGraph::Init(buffer, jointInfoArray, matrixRow);

		dgJointInfo* const jointInfo = &jointInfoArray[m_joint->m_index];
		dgAssert (jointInfo->m_joint == m_joint);

		m_block->m_matrix.m_d[0][0] = dgGetZeroMatrix();
		m_block->m_matrix.m_d[0][1] = dgGetZeroMatrix();
		m_block->m_matrix.m_d[1][0] = dgGetZeroMatrix();
		m_block->m_matrix.m_d[1][1] = dgGetZeroMatrix();

		dgAssert (m_parent);
		dgInt32 rowsCount = jointInfo->m_pairCount;
		dgInt32 rowsCount1 = rowsCount > 4 ? 4 : rowsCount;

		m_block->m_matrix.m_d[0][0] = dgGetIdentityMatrix();
		m_block->m_matrix.m_d[0][1] = dgGetIdentityMatrix();
		m_block->m_matrix.m_d[1][0] = dgGetIdentityMatrix();
		m_block->m_matrix.m_d[1][1] = dgGetIdentityMatrix();
		for (dgInt32 i = 0; i < rowsCount1; i++) {
			dgInt32 index = jointInfo->m_pairStart + i;
			m_block->m_matrix.m_d[0][0][i] = matrixRow[index].m_JMinv.m_jacobianM0.m_linear;
			m_block->m_matrix.m_d[1][0][i] = matrixRow[index].m_JMinv.m_jacobianM0.m_angular;
		}
		for (dgInt32 i = rowsCount1; i < rowsCount; i++) {
			dgInt32 index = jointInfo->m_pairStart + i;
			m_block->m_matrix.m_d[0][1][i] = matrixRow[index].m_JMinv.m_jacobianM0.m_linear;
			m_block->m_matrix.m_d[1][1][i] = matrixRow[index].m_JMinv.m_jacobianM0.m_angular;
		}
	}

	dgBilateralConstraint* m_joint;
};



dgAcyclicContainer::dgAcyclicContainer (dgDynamicBody* const rootBody)
	:m_skeleton(new (rootBody->GetWorld()->GetAllocator()) dgAcyclicGraphMassNode (rootBody->GetWorld()->GetAllocator(), rootBody, NULL))
	,m_topDownOrder(NULL)
	,m_downTopOrder(NULL)
	,m_id(m_uniqueID)
	,m_nodeCount(1)
{
	m_uniqueID ++;
}

dgAcyclicContainer::~dgAcyclicContainer ()
{
	dgMemoryAllocator* const allocator = m_skeleton->m_body->GetWorld()->GetAllocator();
	if (m_topDownOrder) {
		allocator->Free(m_topDownOrder);
		allocator->Free(m_downTopOrder);
	}
	delete m_skeleton;
}


dgAcyclicContainer::dgAcyclicGraph* dgAcyclicContainer::FindNode (dgDynamicBody* const body) const
{
	dgInt32 stack = 1;
	dgAcyclicGraph* stackPool[DG_ACYCLIC_STACK_SIZE];

	stackPool[0] = m_skeleton;
	while (stack) {
		stack --;
		dgAcyclicGraph* const node = stackPool[stack];
		if (node->GetBody() == body) {
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
	dgBody* const parent1 = parent ? parent : m_skeleton->m_body;
	dgAssert (child->GetType() == dgBody::m_dynamicBody);
	dgAssert (parent1->GetType() == dgBody::m_dynamicBody);
	AddChild ((dgDynamicBody*) parent1, (dgDynamicBody*) child);
}

void dgAcyclicContainer::AddChild (dgDynamicBody* const parent, dgDynamicBody* const child)
{
	dgAcyclicGraph* const parentNode = FindNode (parent);
	dgAssert (parentNode);

	dgWorld* const world = m_skeleton->m_body->GetWorld();
	dgMemoryAllocator* const allocator = world->GetAllocator();
	dgBilateralConstraint* const joint = world->FindBilateralJoint (parent, child);
	dgAssert (joint);
	dgAcyclicGraph* const massParent = new (allocator) dgAcyclicGraphJointNode (allocator, joint, parentNode);
	new (allocator) dgAcyclicGraphMassNode (allocator, child, massParent);
	m_nodeCount += 2;
}


dgInt32 dgAcyclicContainer::GetBufferSize () const
{
	dgInt32 blocksize = sizeof (dgAcyclicNodeInfo);
	return m_nodeCount * blocksize;
}

void dgAcyclicContainer::SortGraph (dgAcyclicGraph* const root, dgAcyclicGraph* const parent, dgInt32& index)
{
	for (dgList<dgAcyclicGraph*>::dgListNode* node = root->m_children.GetFirst(); node; node = node->GetNext()) {
		SortGraph (node->GetInfo(), root, index);
	}

	root->SetPriority((m_id << DG_ACYCLIC_BIT_SHIFT_KEY) + index);
	dgAssert ((m_nodeCount - index - 1) >= 0);
	m_downTopOrder[index] = root;
	m_topDownOrder[m_nodeCount - index - 1] = root;
	root->m_index = index;
	index ++;
	dgAssert (index <= m_nodeCount);
}

void dgAcyclicContainer::Finalize ()
{
	dgAssert (m_nodeCount >= 1);
	dgMemoryAllocator* const allocator = m_skeleton->m_body->GetWorld()->GetAllocator();
	m_topDownOrder = (dgAcyclicGraph**) allocator->Malloc(m_nodeCount * sizeof (dgAcyclicGraph*));
	m_downTopOrder = (dgAcyclicGraph**) allocator->Malloc(m_nodeCount * sizeof (dgAcyclicGraph*));

	dgInt32 index = 0;
	SortGraph (m_skeleton, NULL, index);

/*
	dgInt32 stack = 1;
	dgAcyclicGraph* stackPool[DG_ACYCLIC_STACK_SIZE];
	index = 0;
	stackPool[0] = m_skeleton;
	while (stack) {
		stack--;
		dgAcyclicGraph* const node = stackPool[stack];
		node->m_index = index;
		index ++;
		for (dgList<dgAcyclicGraph*>::dgListNode* ptr = node->m_children.GetLast(); ptr; ptr = ptr->GetPrev()) {
			stackPool[stack] = ptr->GetInfo();
			stack++;
			dgAssert(stack < dgInt32(sizeof (stackPool) / sizeof (stackPool[0])));
		}
	}
*/
}


dgFloat32 dgAcyclicContainer::CalculateJointForce (dgJointInfo* const jointInfoArray, const dgBodyInfo* const bodyArray, dgJacobian* const internalForces, dgJacobianMatrixElement* const matrixRow) const
{
dgAssert (0);
/*
	for (dgInt32 i = 0; i < m_jointCount; i ++) {
		dgAcyclicGraph* const node = m_topDownOrder[i];

		dgInt32 index = node->m_index;
		dgAcyclicElement& elemBody = array[index * 2];
		dgBody* const body = node->m_body;
		dgVector mass (body->GetMass());
		dgMatrix inertia (body->CalculateInertiaMatrix());
		elemBody.m_d.m[0].m_linear = mass & dgVector::m_xMask;
		elemBody.m_d.m[0].m_angular = dgVector::m_zero;
		elemBody.m_d.m[1].m_linear = mass & dgVector::m_yMask;
		elemBody.m_d.m[1].m_angular = dgVector::m_zero;
		elemBody.m_d.m[2].m_linear = mass & dgVector::m_zMask;
		elemBody.m_d.m[2].m_angular = dgVector::m_zero;
		elemBody.m_d.m[3].m_linear = dgVector::m_zero;
		elemBody.m_d.m[3].m_angular = inertia[0];
		elemBody.m_d.m[4].m_linear = dgVector::m_zero;
		elemBody.m_d.m[4].m_angular = inertia[1];
		elemBody.m_d.m[5].m_linear = dgVector::m_zero;
		elemBody.m_d.m[5].m_angular = inertia[2];


		dgAcyclicElement& elemJoint = array[index * 2 + 1];
		//dgAcyclicGraph* const node1 = m_downTopOrder___[i];

		elemBody.m_d.m[0].m_linear = dgVector::m_zero;
		elemBody.m_d.m[0].m_angular = dgVector::m_zero;
		elemBody.m_d.m[1].m_linear = dgVector::m_zero;
		elemBody.m_d.m[1].m_angular = dgVector::m_zero;
		elemBody.m_d.m[2].m_linear = dgVector::m_zero;
		elemBody.m_d.m[2].m_angular = dgVector::m_zero;
		elemBody.m_d.m[3].m_linear = dgVector::m_zero;
		elemBody.m_d.m[3].m_angular = dgVector::m_zero;
		elemBody.m_d.m[4].m_linear = dgVector::m_zero;
		elemBody.m_d.m[4].m_angular = dgVector::m_zero;
		elemBody.m_d.m[5].m_linear = dgVector::m_zero;
		elemBody.m_d.m[5].m_angular = dgVector::m_zero;
	}
*/
	return dgFloat32 (0.0f);
}


void dgAcyclicContainer::InitMassMatrix (void* const buffer, dgJointInfo* const jointInfoArray, dgJacobianMatrixElement* const matrixRow)
{
	dgAcyclicNodeInfo* const block = (dgAcyclicNodeInfo*) buffer;

	for (dgInt32 i = 0; i < m_nodeCount; i++) {
		dgAcyclicGraph* const node = m_topDownOrder[i];
		node->Init(block + i, jointInfoArray, matrixRow);
	}

	for (dgInt32 i = 0; i < m_nodeCount; i++) {
		dgAcyclicGraph* const node = m_topDownOrder[i];
		for (dgList<dgAcyclicGraph*>::dgListNode* child = node->m_children.GetFirst(); child; child = child->GetNext()) {
			node->AccumulateMatrix();
		}
		node->CalculateInverse();
		if (node->m_parent) {
			node->AccumulateJacobian();
		}
	}

}