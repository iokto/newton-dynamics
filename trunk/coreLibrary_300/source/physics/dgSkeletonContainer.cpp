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
#include "dgSkeletonContainer.h"
#include "dgWorldDynamicUpdate.h"
#include "dgBilateralConstraint.h"


#define DG_SKELETON_STACK_SIZE		512


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////// 
dgInt32 dgSkeletonContainer::m_uniqueID = 10;



class dgSkeletonContainer::dgSkeletonMatrix
{
	public:
	dgMatrix m_d[2][2];
};


class dgSkeletonContainer::dgSkeletonNodeInfo
{
	public:
	dgSkeletonMatrix m_matrix;
	//dgAcyclicMatrix m_invMatrix;
	dgSkeletonMatrix m_jacobian;
	dgJacobian m_vector;
};


class dgSkeletonContainer::dgSkeletonGraph
{
	public:
	DG_CLASS_ALLOCATOR(allocator)

	dgSkeletonGraph(dgMemoryAllocator* const allocator, dgSkeletonGraph* const parent)
		:m_parent(parent)
		,m_block(NULL)
		,m_children(allocator)
		,m_index(0)
	{
	}

	virtual ~dgSkeletonGraph()
	{
		for (dgList<dgSkeletonGraph*>::dgListNode* ptr = m_children.GetFirst(); ptr; ptr = ptr->GetNext()) {
			delete ptr->GetInfo();
		}
	}

	void AddChild(dgSkeletonGraph* const child)
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


	virtual void Init(dgSkeletonNodeInfo* const buffer, dgJointInfo* const jointInfoArray, dgJacobianMatrixElement* const matrixRow) 
	{
		m_block = buffer;
	}

	void AccumulateMatrix(dgSkeletonGraph* const child)
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

	dgSkeletonGraph* m_parent;
	dgSkeletonNodeInfo* m_block;
	dgList<dgSkeletonGraph*> m_children;
	dgInt32 m_index;
};

class dgSkeletonContainer::dgSkeletonGraphMassNode: public dgSkeletonGraph
{
	public:
	dgSkeletonGraphMassNode(dgMemoryAllocator* const allocator, dgDynamicBody* const body, dgSkeletonGraph* const parent)
		:dgSkeletonGraph (allocator, parent)
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

	virtual void Init(dgSkeletonNodeInfo* const buffer, dgJointInfo* const jointInfoArray, dgJacobianMatrixElement* const matrixRow)
	{
		dgSkeletonGraph::Init(buffer, jointInfoArray, matrixRow);
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

class dgSkeletonContainer::dgSkeletonGraphJointNode: public dgSkeletonGraph
{
	public:
	dgSkeletonGraphJointNode(dgMemoryAllocator* const allocator, dgBilateralConstraint* const Joint, dgSkeletonGraph* const parent)
		:dgSkeletonGraph(allocator, parent)
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

	virtual void Init(dgSkeletonNodeInfo* const buffer, dgJointInfo* const jointInfoArray, dgJacobianMatrixElement* const matrixRow) 
	{
		dgSkeletonGraph::Init(buffer, jointInfoArray, matrixRow);

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



dgSkeletonContainer::dgSkeletonContainer (dgDynamicBody* const rootBody)
	:m_skeleton(new (rootBody->GetWorld()->GetAllocator()) dgSkeletonGraphMassNode (rootBody->GetWorld()->GetAllocator(), rootBody, NULL))
	,m_topDownOrder(NULL)
	,m_downTopOrder(NULL)
	,m_id(m_uniqueID)
	,m_nodeCount(1)
{
	m_uniqueID ++;
}

dgSkeletonContainer::~dgSkeletonContainer ()
{
	dgMemoryAllocator* const allocator = m_skeleton->m_body->GetWorld()->GetAllocator();
	if (m_topDownOrder) {
		allocator->Free(m_topDownOrder);
		allocator->Free(m_downTopOrder);
	}
	delete m_skeleton;
}


dgSkeletonContainer::dgSkeletonGraph* dgSkeletonContainer::FindNode (dgDynamicBody* const body) const
{
	dgInt32 stack = 1;
	dgSkeletonGraph* stackPool[DG_SKELETON_STACK_SIZE];

	stackPool[0] = m_skeleton;
	while (stack) {
		stack --;
		dgSkeletonGraph* const node = stackPool[stack];
		if (node->GetBody() == body) {
			return node;
		}

		for (dgList<dgSkeletonGraph*>::dgListNode* ptr = node->m_children.GetFirst(); ptr; ptr = ptr->GetNext()) {
			stackPool[stack] = ptr->GetInfo();
			stack ++;
			dgAssert (stack < dgInt32 (sizeof (stackPool) / sizeof (stackPool[0])));
		}
	}
	return NULL;
}

void dgSkeletonContainer::AddChild (dgBody* const parent, dgBody* const child)
{
	dgAssert (child);
	dgBody* const parent1 = parent ? parent : m_skeleton->m_body;
	dgAssert (child->GetType() == dgBody::m_dynamicBody);
	dgAssert (parent1->GetType() == dgBody::m_dynamicBody);
	AddChild ((dgDynamicBody*) parent1, (dgDynamicBody*) child);
}

void dgSkeletonContainer::AddChild (dgDynamicBody* const parent, dgDynamicBody* const child)
{
	dgSkeletonGraph* const parentNode = FindNode (parent);
	dgAssert (parentNode);

	dgWorld* const world = m_skeleton->m_body->GetWorld();
	dgMemoryAllocator* const allocator = world->GetAllocator();
	dgBilateralConstraint* const joint = world->FindBilateralJoint (parent, child);
	dgAssert (joint);
	dgSkeletonGraph* const massParent = new (allocator) dgSkeletonGraphJointNode (allocator, joint, parentNode);
	new (allocator) dgSkeletonGraphMassNode (allocator, child, massParent);
	m_nodeCount += 2;
}


dgInt32 dgSkeletonContainer::GetBufferSize () const
{
	dgInt32 blocksize = sizeof (dgSkeletonNodeInfo);
	return m_nodeCount * blocksize;
}

void dgSkeletonContainer::SortGraph (dgSkeletonGraph* const root, dgSkeletonGraph* const parent, dgInt32& index)
{
	for (dgList<dgSkeletonGraph*>::dgListNode* node = root->m_children.GetFirst(); node; node = node->GetNext()) {
		SortGraph (node->GetInfo(), root, index);
	}

	root->SetPriority((m_id << DG_SKELETON_BIT_SHIFT_KEY) + index);
	dgAssert ((m_nodeCount - index - 1) >= 0);
	m_downTopOrder[index] = root;
	m_topDownOrder[m_nodeCount - index - 1] = root;
	root->m_index = index;
	index ++;
	dgAssert (index <= m_nodeCount);
}

void dgSkeletonContainer::Finalize ()
{
	dgAssert (m_nodeCount >= 1);
	dgMemoryAllocator* const allocator = m_skeleton->m_body->GetWorld()->GetAllocator();
	m_topDownOrder = (dgSkeletonGraph**) allocator->Malloc(m_nodeCount * sizeof (dgSkeletonGraph*));
	m_downTopOrder = (dgSkeletonGraph**) allocator->Malloc(m_nodeCount * sizeof (dgSkeletonGraph*));

	dgInt32 index = 0;
	SortGraph (m_skeleton, NULL, index);
}


dgFloat32 dgSkeletonContainer::CalculateJointForce (dgJointInfo* const jointInfoArray, const dgBodyInfo* const bodyArray, dgJacobian* const internalForces, dgJacobianMatrixElement* const matrixRow) const
{
dgAssert (0);
/*
	for (dgInt32 i = 0; i < m_jointCount; i ++) {
		dgSkeletonGraph* const node = m_topDownOrder[i];

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
		//dgSkeletonGraph* const node1 = m_downTopOrder___[i];

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


void dgSkeletonContainer::InitMassMatrix (void* const buffer, dgJointInfo* const jointInfoArray, dgJacobianMatrixElement* const matrixRow)
{
	dgSkeletonNodeInfo* const block = (dgSkeletonNodeInfo*) buffer;

	for (dgInt32 i = 0; i < m_nodeCount; i++) {
		dgSkeletonGraph* const node = m_topDownOrder[i];
		node->Init(block + i, jointInfoArray, matrixRow);
	}

	for (dgInt32 i = 0; i < m_nodeCount; i++) {
		dgSkeletonGraph* const node = m_topDownOrder[i];
		for (dgList<dgSkeletonGraph*>::dgListNode* child = node->m_children.GetFirst(); child; child = child->GetNext()) {
			node->AccumulateMatrix(child->GetInfo());
		}
		node->CalculateInverse();
		if (node->m_parent) {
			node->AccumulateJacobian();
		}
	}

}