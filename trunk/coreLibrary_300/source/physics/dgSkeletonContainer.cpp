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

	class dgSpacialVector: public dgJacobian
	{
		public:
		
	};

	class dgSpacialMatrix
	{
		public:
		void SetZero (dgInt32 rows)
		{
			for (dgInt32 i = 0; i < rows; i ++) {
				m_rows[i].m_linear = dgVector::m_zero;
				m_rows[i].m_angular = dgVector::m_zero;
			}
		}

		void SetIdentity (dgInt32 rows)
		{
			if (rows <= 3) {
				for (dgInt32 i = 0; i < rows; i ++) {
					m_rows[i].m_linear = dgVector::m_zero;
					m_rows[i].m_angular = dgVector::m_zero;
					m_rows[i].m_linear[i] = dgFloat32 (1.0f);
				}
			} else {
				for (dgInt32 i = 0; i < 3; i ++) {
					m_rows[i].m_linear = dgVector::m_zero;
					m_rows[i].m_angular = dgVector::m_zero;
					m_rows[i].m_linear[i] = dgFloat32 (1.0f);
				}

				for (dgInt32 i = 3; i < rows; i ++) {
					m_rows[i].m_linear = dgVector::m_zero;
					m_rows[i].m_angular = dgVector::m_zero;
					m_rows[i].m_angular[i - 3] = dgFloat32 (1.0f);
				}
			}
		}

		void Inverse (const dgSpacialMatrix& src, dgInt32 rows)
		{
			dgSpacialMatrix copy;
			for (dgInt32 i = 0; i < rows; i ++) {
				copy.m_rows[i] = src.m_rows[i];
			}
			SetIdentity(rows);
			if (rows < 3) {
				for (dgInt32 i = 0; i < rows; i++) {
					dgFloat32 val = copy.m_rows[i].m_linear[i];
					dgAssert(dgAbsf(val) > dgFloat32(1.0e-4f));
					dgVector den(dgFloat32(1.0f) / val);
					copy.m_rows[i].m_linear = copy.m_rows[i].m_linear.CompProduct4(den);
					m_rows[i].m_linear = m_rows[i].m_linear.CompProduct4(den);

					for (dgInt32 j = 0; j < rows; j++) {
						if (j != i) {
							dgVector pivot(copy.m_rows[j].m_linear[i]);
							copy.m_rows[j].m_linear -= copy.m_rows[i].m_linear.CompProduct4(pivot);
							m_rows[j].m_linear -= m_rows[i].m_linear.CompProduct4(pivot);
						}
					}
				}

			} else {
				for (dgInt32 i = 0; i < 3; i++) {
					dgFloat32 val = copy.m_rows[i].m_linear[i];
					dgAssert (dgAbsf (val) > dgFloat32 (1.0e-4f));
					dgVector den (dgFloat32(1.0f) / val);
					copy.m_rows[i].m_linear = copy.m_rows[i].m_linear.CompProduct4(den);
					copy.m_rows[i].m_angular = copy.m_rows[i].m_angular.CompProduct4(den);
					m_rows[i].m_linear = m_rows[i].m_linear.CompProduct4(den);
					m_rows[i].m_angular = m_rows[i].m_angular.CompProduct4(den);

					for (dgInt32 j = 0; j < rows; j++) {
						if (j != i) {
							dgVector pivot(copy.m_rows[j].m_linear[i]);
							copy.m_rows[j].m_linear -= copy.m_rows[i].m_linear.CompProduct4(pivot);
							copy.m_rows[j].m_angular -= copy.m_rows[i].m_angular.CompProduct4(pivot);
							m_rows[j].m_linear -= m_rows[i].m_linear.CompProduct4(pivot);
							m_rows[j].m_angular -= m_rows[i].m_angular.CompProduct4(pivot);
						}
					}
				}

				for (dgInt32 i = 3; i < rows; i++) {
					dgFloat32 val = copy.m_rows[i].m_angular[i - 3];
					dgAssert(dgAbsf(val) > dgFloat32(1.0e-4f));
					dgVector den(dgFloat32(1.0f) / val);
					copy.m_rows[i].m_linear = copy.m_rows[i].m_linear.CompProduct4(den);
					copy.m_rows[i].m_angular = copy.m_rows[i].m_angular.CompProduct4(den);
					m_rows[i].m_linear = m_rows[i].m_linear.CompProduct4(den);
					m_rows[i].m_angular = m_rows[i].m_angular.CompProduct4(den);

					for (dgInt32 j = 0; j < rows; j++) {
						if (j != i) {
							dgVector pivot(copy.m_rows[j].m_angular[i - 3]);
							copy.m_rows[j].m_linear -= copy.m_rows[i].m_linear.CompProduct4(pivot);
							copy.m_rows[j].m_angular -= copy.m_rows[i].m_angular.CompProduct4(pivot);
							m_rows[j].m_linear -= m_rows[i].m_linear.CompProduct4(pivot);
							m_rows[j].m_angular -= m_rows[i].m_angular.CompProduct4(pivot);
						}
					}
				}
			}
		}

		void MultiplyMatrix6x6TimeJacobianTransposed (dgSpacialVector& jacobian) const
		{
			dgSpacialVector tmp (jacobian);
			for (dgInt32 i = 0; i < 3; i ++) {
				jacobian.m_linear[i] = (m_rows[i].m_linear.DotProduct4(tmp.m_linear) + m_rows[i].m_angular.DotProduct4(tmp.m_angular)).GetScalar();
				jacobian.m_angular[i] = (m_rows[i + 3].m_linear.DotProduct4(tmp.m_linear) + m_rows[i + 3].m_angular.DotProduct4(tmp.m_angular)).GetScalar();
			}
		}

		void SubstractCovariance (const dgSpacialVector& jt, const dgSpacialVector& j)
		{
			for (dgInt32 i = 0; i < 3; i ++) {
				dgVector l (jt.m_linear[i]);
				m_rows[i].m_linear -= j.m_linear.CompProduct4(l);
				m_rows[i].m_angular -= j.m_angular.CompProduct4(l);

				dgVector a(jt.m_angular[i]);
				m_rows[i + 3].m_linear -= j.m_linear.CompProduct4(a);
				m_rows[i + 3].m_angular -= j.m_angular.CompProduct4(a);
			}

			#ifdef _DEBUG
			for (dgInt32 i = 0; i < 6; i ++) {
				for (dgInt32 j = 0; j < 6; j ++) {
					dgFloat32 a = (j <= 3) ? m_rows[i].m_linear[j] : m_rows[i].m_angular[j - 3];
					dgFloat32 b = (i <= 3) ? m_rows[j].m_linear[i] : m_rows[j].m_angular[i - 3];
					dgAssert (dgAbsf(a - b) < dgFloat32 (1.0e-4f));
				}
			}
			#endif
		}

		dgSpacialVector m_rows[6];
	};

	class dgData
	{
		public:
		dgSpacialMatrix m_diagonal;
		dgSpacialMatrix m_invDiagonal;
		dgSpacialMatrix m_offDiagonal;
	};


	dgSkeletonGraph(dgMemoryAllocator* const allocator, dgDynamicBody* const root)
		:m_parent(NULL)
		,m_body(root)
		,m_joint(NULL)
		,m_data(NULL)
		,m_children(allocator)
		,m_index(0)
		,m_diaginalDof(0)
		,m_jacobialDof(0)
	{
	}

	dgSkeletonGraph(dgMemoryAllocator* const allocator, dgDynamicBody* const child, dgSkeletonGraph* const parent)
		:m_parent(parent)
		,m_body(child)
		,m_joint(NULL)
		,m_data(NULL)
		,m_children(allocator)
		,m_index(0)
		, m_diaginalDof(0)
		, m_jacobialDof(0)
	{
		m_parent->m_children.Append(this);
	}

	dgSkeletonGraph(dgMemoryAllocator* const allocator, dgBilateralConstraint* const Joint, dgSkeletonGraph* const parent)
		:m_parent(parent)
		,m_body(NULL)
		,m_joint(Joint)
		,m_data(NULL)
		,m_children(allocator)
		,m_index(0)
		, m_diaginalDof(0)
		, m_jacobialDof(0)
	{
		m_parent->m_children.Append(this);
	}


	~dgSkeletonGraph()
	{
		for (dgList<dgSkeletonGraph*>::dgListNode* ptr = m_children.GetFirst(); ptr; ptr = ptr->GetNext()) {
			delete ptr->GetInfo();
		}
	}

	bool IsMass() const
	{
		return false;
	}

	bool IsJoint() const
	{
		return false;
	}

	virtual void SetPriority(dgUnsigned32 priority) const
	{
		if (m_joint) {
			m_joint->m_priority = priority;
		}
	}

	void Init (dgData* const buffer, dgJointInfo* const jointInfoArray, dgJacobianMatrixElement* const matrixRow) 
	{
		m_data = buffer;
		if (m_body) {
			m_diaginalDof = 6;
			dgFloat32 mass = m_body->GetMass().m_w;
			dgMatrix inertia (m_body->CalculateInertiaMatrix());
			for (dgInt32 i = 0; i < 3; i++) {
				m_data->m_diagonal.m_rows[i].m_linear = dgVector::m_zero;
				m_data->m_diagonal.m_rows[i].m_angular = dgVector::m_zero;
				m_data->m_diagonal.m_rows[i].m_linear[i] = mass;

				m_data->m_diagonal.m_rows[i + 3].m_linear = dgVector::m_zero;
				m_data->m_diagonal.m_rows[i + 3].m_angular = inertia[i];
			}

			if (m_parent) {
				dgBilateralConstraint* const joint = m_parent->m_joint;
				dgAssert (joint);
				dgAssert (joint->GetBody0() == m_body);
				dgAssert (joint->GetBody1() == m_parent->m_parent->m_body);
				dgJointInfo* const jointInfo = &jointInfoArray[joint->m_index];
				dgAssert (jointInfo->m_joint == joint);
				m_jacobialDof = jointInfo->m_pairCount;

				for (dgInt32 i = 0; i < m_jacobialDof; i ++) {
					dgInt32 index = jointInfo->m_pairStart + i;
					m_data->m_offDiagonal.m_rows[i].m_linear = matrixRow[index].m_JMinv.m_jacobianM0.m_linear;
					m_data->m_offDiagonal.m_rows[i].m_angular = matrixRow[index].m_JMinv.m_jacobianM0.m_angular;
				}
			}
		} else {
			dgJointInfo* const jointInfo = &jointInfoArray[m_joint->m_index];

			dgAssert (m_parent);
			dgAssert (jointInfo->m_joint == m_joint);
			dgAssert (jointInfo->m_joint->GetBody1() == m_parent->m_body);
			m_diaginalDof = jointInfo->m_pairCount;
			m_jacobialDof = jointInfo->m_pairCount;

			m_data->m_diagonal.SetZero(m_diaginalDof);
			for (dgInt32 i = 0; i < m_diaginalDof; i++) {
				dgInt32 index = jointInfo->m_pairStart + i;
				m_data->m_offDiagonal.m_rows[i].m_linear = matrixRow[index].m_JMinv.m_jacobianM1.m_linear;
				m_data->m_offDiagonal.m_rows[i].m_angular = matrixRow[index].m_JMinv.m_jacobianM1.m_angular;
			}
		}
	}

	void CalculateOffDiagonalBlock ()
	{
		if (m_body) {
			for (dgInt32 i = 0; i < m_jacobialDof; i ++) {
				m_data->m_invDiagonal.MultiplyMatrix6x6TimeJacobianTransposed (m_data->m_offDiagonal.m_rows[i]);
			}
		} else {
			dgSpacialMatrix copy;
			for (dgInt32 i = 0; i < m_jacobialDof; i ++) {
				copy.m_rows[i] = m_data->m_offDiagonal.m_rows[i];
			}
			m_data->m_offDiagonal.SetZero(m_jacobialDof);
			for (dgInt32 i = 0; i < m_jacobialDof; i ++) {
				const dgFloat32* const line = &m_data->m_invDiagonal.m_rows[i].m_linear.m_x; 
				const dgSpacialVector& jacovian = copy.m_rows[i];
				for (dgInt32 j = 0; j < m_jacobialDof; j ++) {
					dgVector val (line[i]);
					m_data->m_offDiagonal.m_rows[j].m_linear += jacovian.m_linear.CompProduct4(val);
					m_data->m_offDiagonal.m_rows[j].m_angular += jacovian.m_angular.CompProduct4(val);
				}
			}
		}
	}

	virtual void CalculateDiagonal(dgSkeletonGraph* const child, dgJointInfo* const jointInfoArray)
	{
		dgSpacialMatrix tmp;
		const dgSpacialMatrix& childDiagonal = child->m_data->m_diagonal;
		const dgSpacialMatrix& jacobianTransposed = child->m_data->m_offDiagonal;
		if (m_body) {
			dgSpacialMatrix copy;
			copy.SetZero(child->m_jacobialDof);

			dgAssert (child->m_jacobialDof  <= 3);
			for (dgInt32 i = 0; i < child->m_jacobialDof; i++) {
				const dgFloat32* const row = &childDiagonal.m_rows[i].m_linear.m_x;
				dgAssert (row[3] == dgFloat32 (0.0f));
				dgAssert (row[7] == dgFloat32 (0.0f));
				const dgSpacialVector& jacovian = jacobianTransposed.m_rows[i];
				for (dgInt32 j = 0; j < child->m_jacobialDof; j++) {
					dgVector val0(row[i]);
					copy.m_rows[j].m_linear += jacovian.m_linear.CompProduct4(val0);
					copy.m_rows[j].m_angular += jacovian.m_angular.CompProduct4(val0);
				}
			}
			
			dgAssert (m_diaginalDof == 6);
			dgSpacialMatrix& diagonal = m_data->m_diagonal;
			for (dgInt32 i = 0; i < child->m_jacobialDof; i++) {
				const dgSpacialVector& jt = jacobianTransposed.m_rows[i];
				for (dgInt32 j = 0; j < child->m_jacobialDof; j ++) {
					diagonal.SubstractCovariance (jt, copy.m_rows[j]);
				}
			}

		} else {
			dgAssert (child->m_body);
			for (dgInt32 i = 0; i < m_jacobialDof; i++) {
				tmp.m_rows[i] = jacobianTransposed.m_rows[i];
				childDiagonal.MultiplyMatrix6x6TimeJacobianTransposed(tmp.m_rows[i]);
			}

			dgSpacialMatrix& diagonal = m_data->m_diagonal;
			if (m_jacobialDof <= 3) {
				for (dgInt32 i = 0; i < m_jacobialDof; i++) {
					dgFloat32 a = (jacobianTransposed.m_rows[i].m_linear.DotProduct4(tmp.m_rows[i].m_linear) + jacobianTransposed.m_rows[i].m_angular.DotProduct4(tmp.m_rows[i].m_angular)).GetScalar();
					diagonal.m_rows[i].m_linear[i] -= a;
					for (dgInt32 j = 0; j < i; j ++) {
						dgFloat32 a = (jacobianTransposed.m_rows[i].m_linear.DotProduct4(tmp.m_rows[j].m_linear) + jacobianTransposed.m_rows[i].m_angular.DotProduct4(tmp.m_rows[j].m_angular)).GetScalar();
						dgAssert (dgAbsf ((jacobianTransposed.m_rows[j].m_linear.DotProduct4(tmp.m_rows[i].m_linear) + jacobianTransposed.m_rows[j].m_angular.DotProduct4(tmp.m_rows[i].m_angular)).GetScalar() - a) < dgFloat32 (1.0e-5f));
						diagonal.m_rows[i].m_linear[j] -= a;
						diagonal.m_rows[j].m_linear[i] -= a;
					}
				}
			} else {
				dgAssert(0);
			}
		}
	}

	virtual void CalculateDiagonalInverse ()
	{
		m_data->m_invDiagonal.Inverse (m_data->m_diagonal, m_diaginalDof);
	}

	dgSkeletonGraph* m_parent;
	dgDynamicBody* m_body;
	dgBilateralConstraint* m_joint;
	dgData* m_data;
	dgList<dgSkeletonGraph*> m_children;
	dgInt32 m_index;
	dgInt16 m_diaginalDof;
	dgInt16 m_jacobialDof;
};


dgSkeletonContainer::dgSkeletonContainer (dgDynamicBody* const rootBody)
	:m_skeleton(new (rootBody->GetWorld()->GetAllocator()) dgSkeletonGraph (rootBody->GetWorld()->GetAllocator(), rootBody))
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
		if (node->m_body == body) {
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
	dgSkeletonGraph* const massParent = new (allocator) dgSkeletonGraph (allocator, joint, parentNode);
	new (allocator) dgSkeletonGraph (allocator, child, massParent);
	m_nodeCount += 2;
}


dgInt32 dgSkeletonContainer::GetBufferSize () const
{
	dgInt32 blocksize = sizeof(dgSkeletonGraph::dgData);
	return blocksize * m_nodeCount;
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
	dgSkeletonGraph::dgData* const ptr = (dgSkeletonGraph::dgData*) buffer;
	dgAssert( (dgUnsigned64(buffer) & 0x0f) == 0);
	for (dgInt32 i = 0; i < m_nodeCount; i++) {
		dgSkeletonGraph* const node = m_bottomTopOrder[i];
		node->Init(&ptr[i], jointInfoArray, matrixRow);
	}

	for (dgInt32 i = 0; i < m_nodeCount; i++) {
		dgSkeletonGraph* const node = m_bottomTopOrder[i];
		for (dgList<dgSkeletonGraph*>::dgListNode* child = node->m_children.GetFirst(); child; child = child->GetNext()) {
			node->CalculateDiagonal(child->GetInfo(), jointInfoArray);
		}
		node->CalculateDiagonalInverse();
		if (node->m_parent) {
			node->CalculateOffDiagonalBlock();
		}
	}
}