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


class dgSkeletonContainer::dgSkeletonGraph
{
	public:
	DG_CLASS_ALLOCATOR(allocator)

	class dgData
	{
		public:
		dgJacobian m_jacobian[6];
		dgJacobian m_vector;
	};

	class dgBodyData: public dgData
	{
		public:
		dgMatrix m_inertia;
		dgMatrix m_invInertia;
		dgVector m_mass;
		dgVector m_invMass;
	};

	class dgJointData: public dgData
	{
		public:
		dgFloat32 m_JJ[6][6];
		dgFloat32 m_invJJ[6][6];
	};

	dgSkeletonGraph(dgMemoryAllocator* const allocator, dgSkeletonGraph* const parent)
		:m_parent(parent)
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

	virtual bool IsMass() const
	{
		return false;
	}

	virtual bool IsJoint() const
	{
		return false;
	}

	virtual dgDynamicBody* GetBody() const
	{
		return NULL;
	}

	virtual dgBilateralConstraint* GetJoint() const
	{
		return NULL;
	}

	virtual dgBodyData* GetBodyData() const 
	{
		dgAssert (0);
		return NULL;
	}

	virtual dgJointData* GetJointData() const
	{
		dgAssert(0);
		return NULL;
	}


	void AddChild(dgSkeletonGraph* const child)
	{
		m_children.Append(child);
	}


	virtual void SetPriority(dgUnsigned32 priority) const 
	{
	}

	virtual void Init (char** const buffer, dgJointInfo* const jointInfoArray, dgJacobianMatrixElement* const matrixRow) 
	{
		dgAssert (0);
	}

	virtual void CalculateDiagonal(dgSkeletonGraph* const child, dgJointInfo* const jointInfoArray)
	{
		dgAssert (0);
	}

	virtual void CalculateInverse(dgJointInfo* const jointInfoArray)
	{
		dgAssert (0);
	}

	virtual void CalculateOffDiagonalBlock(dgJointInfo* const jointInfoArray)
	{
		dgAssert (0);
	}

	dgSkeletonGraph* m_parent;
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

	virtual bool IsMass() const
	{
		return true;
	}

	virtual dgBodyData* GetBodyData() const
	{
		return m_data;
	}

	virtual dgDynamicBody* GetBody() const 
	{
		return m_body;
	}

	virtual void Init(char** const buffer, dgJointInfo* const jointInfoArray, dgJacobianMatrixElement* const matrixRow)
	{
		m_data = (dgBodyData*) *buffer;
		*buffer += sizeof (dgBodyData);
		dgAssert( (dgUnsigned64(m_data) & 0x0f) == 0);

		dgVector mass (m_body->GetMass().m_w);
		m_data->m_mass = m_body->GetMass();
		m_data->m_mass.m_w = dgFloat32 (1.0f);
		m_data->m_inertia = m_body->CalculateInertiaMatrix();

		if (m_parent) {
			dgBilateralConstraint* const joint = m_parent->GetJoint();
			dgAssert (joint->GetBody0() == m_body);
			dgAssert (joint->GetBody1() == m_parent->m_parent->GetBody());
			dgJointInfo* const jointInfo = &jointInfoArray[joint->m_index];
			dgAssert (jointInfo->m_joint == joint);
			dgInt32 rowsCount = jointInfo->m_pairCount;

			for (dgInt32 i = 0; i < rowsCount; i ++) {
				dgInt32 index = jointInfo->m_pairStart + i;
				m_data->m_jacobian[i] = matrixRow[index].m_JMinv.m_jacobianM0;
			}
		}
	}

	virtual void CalculateDiagonal(dgSkeletonGraph* const child, dgJointInfo* const jointInfoArray)
	{
		dgAssert(0);
/*
		dgJacobian mJtArray[6];

		dgAssert(child->IsMass());
		dgJointInfo* const jointInfo = &jointInfoArray[m_joint->m_index];
		const dgInt32 rows = jointInfo->m_pairCount;
		dgBodyData* const bodyData = child->GetBodyData();

		const dgVector& mass = bodyData->m_mass;
		const dgMatrix& inertiaMatrix = bodyData->m_inertia;
		dgJacobian* const jacobians = bodyData->m_jacobian;
		for (dgInt32 i = 0; i < rows; i++) {
			mJtArray[i].m_linear = mass.CompProduct4(jacobians[i].m_linear);
			mJtArray[i].m_angular = inertiaMatrix.RotateVector(jacobians[i].m_angular);
		}

		for (dgInt32 i = 0; i < rows; i++) {
			m_data->m_JJ[i][i] -= (jacobians[i].m_linear.DotProduct4(mJtArray[i].m_linear) + jacobians[i].m_angular.DotProduct4(mJtArray[i].m_angular)).GetScalar();
			for (dgInt32 j = 0; j < i; j++) {
				dgFloat32 a = (jacobians[i].m_linear.DotProduct4(mJtArray[j].m_linear) + jacobians[i].m_angular.DotProduct4(mJtArray[j].m_angular)).GetScalar();
				dgAssert((a - (jacobians[j].m_linear.DotProduct4(mJtArray[i].m_linear) + jacobians[j].m_angular.DotProduct4(mJtArray[i].m_angular)).GetScalar()) < dgFloat32(1.0e-4f));
				m_data->m_JJ[i][j] -= a;
				m_data->m_JJ[j][i] -= a;
			}
		}
*/
	}

	virtual void CalculateInverse(dgJointInfo* const jointInfoArray)
	{
		m_data->m_invMass = dgVector (dgFloat32 (1.0f) / m_data->m_mass[0], dgFloat32 (1.0f) / m_data->m_mass[1], dgFloat32 (1.0f) / m_data->m_mass[2], dgFloat32 (1.0f));
		m_data->m_invInertia = m_data->m_inertia.Symetric3by3Inverse();
	}

	virtual void CalculateOffDiagonalBlock(dgJointInfo* const jointInfoArray)
	{
		dgAssert (m_parent && m_parent->GetJoint());
		dgBilateralConstraint* const joint = m_parent->GetJoint();	
		dgJointInfo* const jointInfo = &jointInfoArray[joint->m_index];
		const dgInt32 rows = jointInfo->m_pairCount;

		const dgVector& invMass = m_data->m_invMass;
		const dgMatrix& invInertia = m_data->m_invInertia;
		for (dgInt32 i = 0; i < rows; i ++) {
			m_data->m_jacobian[i].m_linear = m_data->m_jacobian[i].m_linear.CompProduct4(invMass);
			m_data->m_jacobian[i].m_angular = invInertia.RotateVector(m_data->m_jacobian[i].m_angular);
		}
	}

	dgDynamicBody* m_body;
	dgBodyData* m_data;
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

	virtual void Init (char** const buffer, dgJointInfo* const jointInfoArray, dgJacobianMatrixElement* const matrixRow) 
	{
		m_data = (dgJointData*)*buffer;
		*buffer += sizeof (dgJointData);
		dgAssert( (dgUnsigned64(m_data) & 0x0f) == 0);

		dgJointInfo* const jointInfo = &jointInfoArray[m_joint->m_index];

		dgAssert (m_parent);
		dgAssert (jointInfo->m_joint == m_joint);
		dgAssert (jointInfo->m_joint->GetBody1() == m_parent->GetBody());
		dgInt32 rowsCount = jointInfo->m_pairCount;
		for (dgInt32 i = 0; i < rowsCount; i++) {
			dgInt32 index = jointInfo->m_pairStart + i;
			m_data->m_jacobian[i] = matrixRow[index].m_JMinv.m_jacobianM1;
			memset (&m_data->m_JJ[i][0], 0, rowsCount * sizeof (dgFloat32));
			memset (&m_data->m_invJJ[i][0], 0, rowsCount * sizeof (dgFloat32));
			m_data->m_invJJ[i][i] = dgFloat32 (1.0f);
		}
	}

	void CalculateDiagonal(dgSkeletonGraph* const child, dgJointInfo* const jointInfoArray)
	{
		dgJacobian mJtArray[6];

		dgAssert (child->IsMass());
		dgJointInfo* const jointInfo = &jointInfoArray[m_joint->m_index];
		const dgInt32 rows = jointInfo->m_pairCount;
		dgBodyData* const bodyData = child->GetBodyData();

		const dgVector& mass = bodyData->m_mass;
		const dgMatrix& inertiaMatrix = bodyData->m_inertia;
		dgJacobian* const jacobians =  bodyData->m_jacobian;
		for (dgInt32 i = 0; i < rows; i ++) {
			mJtArray[i].m_linear = mass.CompProduct4(jacobians[i].m_linear);
			mJtArray[i].m_angular = inertiaMatrix.RotateVector(jacobians[i].m_angular);
		}

		for (dgInt32 i = 0; i < rows; i ++) {
			m_data->m_JJ[i][i] -= (jacobians[i].m_linear.DotProduct4(mJtArray[i].m_linear) + jacobians[i].m_angular.DotProduct4(mJtArray[i].m_angular)).GetScalar();
			for (dgInt32 j = 0; j < i; j ++) {
				dgFloat32 a = (jacobians[i].m_linear.DotProduct4(mJtArray[j].m_linear) + jacobians[i].m_angular.DotProduct4(mJtArray[j].m_angular)).GetScalar();
				dgAssert ((a - (jacobians[j].m_linear.DotProduct4(mJtArray[i].m_linear) + jacobians[j].m_angular.DotProduct4(mJtArray[i].m_angular)).GetScalar()) < dgFloat32 (1.0e-4f));
				m_data->m_JJ[i][j] -= a;
				m_data->m_JJ[j][i] -= a;
			}
		}
	}

	virtual void CalculateInverse(dgJointInfo* const jointInfoArray)
	{
		dgJointInfo* const jointInfo = &jointInfoArray[m_joint->m_index];
		const dgInt32 rows = jointInfo->m_pairCount;

		dgFloat32 copy[6][6];
		memcpy (copy, m_data->m_JJ, sizeof(copy));

		for (dgInt32 i = 0; i < rows; i ++) {
			dgAssert (dgAbsf(copy[i][i]) > dgFloat32 (1.0e-5f));
			dgFloat32  den = dgFloat32(1.0f) / copy[i][i];
			for(dgInt32 j = 0; j < rows; j ++) {
				copy[i][j] *= den;
				m_data->m_invJJ[i][j] *= den;
			}
			copy[i][i] = dgFloat32 (1.0f);

			for (dgInt32 j = 0; j < rows; j ++) {
				if (j != i) {
					dgFloat32 pivot = copy[j][i];
					for (dgInt32 k = 0; k < rows; k ++) {
						copy[j][k] -= pivot * copy[i][k];
						m_data->m_invJJ[j][k] -= pivot * m_data->m_invJJ[i][k];
					}
					copy[j][i] = dgFloat32 (0.0f);
				}
			}
		}
	}

	virtual void CalculateOffDiagonalBlock(dgJointInfo* const jointInfoArray)
	{
		dgJointInfo* const jointInfo = &jointInfoArray[m_joint->m_index];
		const dgInt32 rows = jointInfo->m_pairCount;

		dgJacobian copy[6];
		dgJacobian* const jacobians = m_data->m_jacobian;
		for (dgInt32 i = 0; i < rows; i ++) {
			copy[i] = jacobians[i];
		}

		for (dgInt32 i = 0; i < rows; i ++) {
			for (dgInt32 j = 0; j < 3; j ++) {
				dgFloat32 linear = dgFloat32 (0.0f);
				dgFloat32 angular = dgFloat32 (0.0f);
				for (dgInt32 k = 0; k < rows; k ++) {
					linear += m_data->m_invJJ[i][k] * copy[k].m_linear[j];
					angular += m_data->m_invJJ[i][k] * copy[k].m_angular[j];
				}
				jacobians[i].m_linear[j] = linear;
				jacobians[i].m_angular[j] = angular;
			}
		}
	}


	dgBilateralConstraint* m_joint;
	dgJointData* m_data;
};



dgSkeletonContainer::dgSkeletonContainer (dgDynamicBody* const rootBody)
	:m_skeleton(new (rootBody->GetWorld()->GetAllocator()) dgSkeletonGraphMassNode (rootBody->GetWorld()->GetAllocator(), rootBody, NULL))
	,m_bottomTopOrder(NULL)
	,m_topBottomOrder(NULL)
	,m_id(m_uniqueID)
	,m_nodeCount(1)
{
	m_uniqueID ++;
}

dgSkeletonContainer::~dgSkeletonContainer ()
{
	dgMemoryAllocator* const allocator = m_skeleton->m_body->GetWorld()->GetAllocator();
	if (m_bottomTopOrder) {
		allocator->Free(m_bottomTopOrder);
		allocator->Free(m_topBottomOrder);
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

void dgSkeletonContainer::AddChild (dgBody* const child, dgBody* const parent)
{
	dgAssert (child);
	dgBody* const parent1 = parent ? parent : m_skeleton->m_body;
	dgAssert (child->GetType() == dgBody::m_dynamicBody);
	dgAssert (parent1->GetType() == dgBody::m_dynamicBody);
	AddChild ((dgDynamicBody*) child, (dgDynamicBody*) parent1);
}

void dgSkeletonContainer::AddChild (dgDynamicBody* const child, dgDynamicBody* const parent)
{
	dgSkeletonGraph* const parentNode = FindNode (parent);
	dgAssert (parentNode);
	dgWorld* const world = m_skeleton->m_body->GetWorld();
	dgMemoryAllocator* const allocator = world->GetAllocator();
	dgBilateralConstraint* const joint = world->FindBilateralJoint (child, parent);
	dgAssert (joint);

	dgAssert (joint->GetBody0() == child);
	dgAssert (joint->GetBody1() == parent);
	dgSkeletonGraph* const massParent = new (allocator) dgSkeletonGraphJointNode (allocator, joint, parentNode);
	new (allocator) dgSkeletonGraphMassNode (allocator, child, massParent);
	m_nodeCount += 2;
}


dgInt32 dgSkeletonContainer::GetBufferSize () const
{
	dgInt32 blocksize = sizeof (dgSkeletonGraphMassNode::dgBodyData) + sizeof (dgSkeletonGraphJointNode::dgJointData);
	return blocksize * (m_nodeCount  + 1) / 2;
}

void dgSkeletonContainer::SortGraph (dgSkeletonGraph* const root, dgSkeletonGraph* const parent, dgInt32& index)
{
	for (dgList<dgSkeletonGraph*>::dgListNode* node = root->m_children.GetFirst(); node; node = node->GetNext()) {
		SortGraph (node->GetInfo(), root, index);
	}

	root->SetPriority((m_id << DG_SKELETON_BIT_SHIFT_KEY) + index);
	dgAssert ((m_nodeCount - index - 1) >= 0);
	m_topBottomOrder[m_nodeCount - index - 1] = root;
	m_bottomTopOrder[index] = root;
	root->m_index = index;
	index ++;
	dgAssert (index <= m_nodeCount);
}

void dgSkeletonContainer::Finalize ()
{
	dgAssert (m_nodeCount >= 1);
	dgMemoryAllocator* const allocator = m_skeleton->m_body->GetWorld()->GetAllocator();
	m_bottomTopOrder = (dgSkeletonGraph**) allocator->Malloc(m_nodeCount * sizeof (dgSkeletonGraph*));
	m_topBottomOrder = (dgSkeletonGraph**) allocator->Malloc(m_nodeCount * sizeof (dgSkeletonGraph*));

	dgInt32 index = 0;
	SortGraph (m_skeleton, NULL, index);
}


dgFloat32 dgSkeletonContainer::CalculateJointForce (dgJointInfo* const jointInfoArray, const dgBodyInfo* const bodyArray, dgJacobian* const internalForces, dgJacobianMatrixElement* const matrixRow) const
{
dgAssert (0);
	return dgFloat32 (0.0f);
}


void dgSkeletonContainer::InitMassMatrix (char* const buffer, dgJointInfo* const jointInfoArray, dgJacobianMatrixElement* const matrixRow)
{
	char* ptr = buffer;
	dgAssert( (dgUnsigned64(buffer) & 0x0f) == 0);
	for (dgInt32 i = 0; i < m_nodeCount; i++) {
		dgSkeletonGraph* const node = m_bottomTopOrder[i];
		node->Init(&ptr, jointInfoArray, matrixRow);
	}

	for (dgInt32 i = 0; i < m_nodeCount; i++) {
		dgSkeletonGraph* const node = m_bottomTopOrder[i];
		for (dgList<dgSkeletonGraph*>::dgListNode* child = node->m_children.GetFirst(); child; child = child->GetNext()) {
			node->CalculateDiagonal(child->GetInfo(), jointInfoArray);
		}
		node->CalculateInverse(jointInfoArray);
		if (node->m_parent) {
			node->CalculateOffDiagonalBlock(jointInfoArray);
		}
	}

}