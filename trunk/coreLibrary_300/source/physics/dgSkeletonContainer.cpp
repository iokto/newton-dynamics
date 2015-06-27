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
			dgAssert(0);
		}
	}

	virtual void CalculateDiagonal(dgSkeletonGraph* const child, dgJointInfo* const jointInfoArray)
	{
		if (m_body) {
			dgAssert(0);
		} else {
			dgAssert (child->m_body);
			const dgSpacialMatrix& childDiagonal = child->m_data->m_diagonal;
			const dgSpacialMatrix& jacobianTransposed = child->m_data->m_offDiagonal;

			dgSpacialMatrix tmp;
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


/*
	void CalculateDiagonal(dgSkeletonGraph* const child, dgJointInfo* const jointInfoArray)
	{
		dgAssert (0);
	}

	void CalculateInverse(dgJointInfo* const jointInfoArray)
	{
		dgAssert (0);
	}

	void CalculateOffDiagonalBlock(dgJointInfo* const jointInfoArray)
	{
		dgAssert (0);
	}
*/
	dgSkeletonGraph* m_parent;
	dgDynamicBody* m_body;
	dgBilateralConstraint* m_joint;
	dgData* m_data;
	dgList<dgSkeletonGraph*> m_children;
	dgInt32 m_index;
	dgInt16 m_diaginalDof;
	dgInt16 m_jacobialDof;
};

/*
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

*/

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