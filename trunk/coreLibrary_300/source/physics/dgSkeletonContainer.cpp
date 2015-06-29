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
/*
	class dgSpacialVector
	{
		public:
		dgVector m_v[2];
		
		DG_INLINE dgFloat32& operator[] (dgInt32 i)
		{
			dgAssert(i < 8);
			dgAssert(i >= 0);
			dgFloat32* const ptr = &m_v[0].m_x;
			return ptr[i];
		}

		DG_INLINE const dgFloat32& operator[] (dgInt32 i) const
		{
			dgAssert(i < 8);
			dgAssert(i >= 0);
			const dgFloat32* const ptr = &m_v[0].m_x;
			return ptr[i];
		}

		DG_INLINE dgFloat32 DotProduct (const dgSpacialVector& v) const
		{
			return (m_v[0].DotProduct4(v.m_v[0]) + m_v[1].DotProduct4(v.m_v[1])).GetScalar();
		}

		DG_INLINE void Scale (const dgVector& s, dgSpacialVector& dst) const
		{
			dst.m_v[0] = m_v[0].CompProduct4(s);
			dst.m_v[1] = m_v[1].CompProduct4(s);
		}

		DG_INLINE void ScaleAdd(const dgVector& s, const dgSpacialVector& b, dgSpacialVector& dst) const
		{
			dst.m_v[0] = b.m_v[0] + m_v[0].CompProduct4(s);
			dst.m_v[1] = b.m_v[1] + m_v[1].CompProduct4(s);
		}
	};

	class dgSpacialMatrix
	{
		public:
		DG_INLINE dgSpacialVector& operator[] (dgInt32 i)
		{
			dgAssert(i < 6);
			dgAssert(i >= 0);
			return m_rows[i];
		}

		DG_INLINE const dgSpacialVector& operator[] (dgInt32 i) const
		{
			dgAssert(i < 6);
			dgAssert(i >= 0);
			return m_rows[i];
		}

		DG_INLINE void SetZero ()
		{
			memset (m_rows, 0, sizeof(m_rows));
		}

		DG_INLINE void SetIdentity (dgInt32 rows)
		{
			for (dgInt32 i = 0; i < rows; i ++) {
				m_rows[i].m_v[0] = dgVector::m_zero;
				m_rows[i].m_v[1] = dgVector::m_zero;
				m_rows[i][i] = dgFloat32 (1.0f);
			}
		}

		DG_INLINE void Inverse (const dgSpacialMatrix& src, dgInt32 rows)
		{
			dgSpacialMatrix copy;
			for (dgInt32 i = 0; i < rows; i ++) {
				copy[i] = src[i];
			}
			SetIdentity(rows);

			for (dgInt32 i = 0; i < rows; i++) {
				dgFloat32 val = copy.m_rows[i][i];
				dgAssert(dgAbsf(val) > dgFloat32(1.0e-12f));
				dgVector den(dgFloat32(1.0f) / val);

				m_rows[i].Scale(den, m_rows[i]);
				copy[i].Scale (den, copy[i]);

				for (dgInt32 j = 0; j < rows; j++) {
					if (j != i) {
						dgVector pivot(-copy[j][i]);
						m_rows[i].ScaleAdd (pivot, m_rows[j], m_rows[j]);
						copy[i].ScaleAdd (pivot, copy[j], copy[j]);
					}
				}
			}
		}

		DG_INLINE void MultiplyMatrix6x6TimeJacobianTransposed (dgSpacialVector& jacobian) const
		{
			dgSpacialVector tmp (jacobian);
			dgAssert (tmp[6] == dgFloat32 (0.0f));
			dgAssert (tmp[7] == dgFloat32 (0.0f));
			for (dgInt32 i = 0; i < 6; i ++) {
				jacobian[i] = m_rows[i].DotProduct(tmp);
			}
		}

		DG_INLINE void Multiply (const dgSpacialMatrix& B, dgSpacialMatrix& dst) const
		{
			const dgSpacialMatrix& A = *this;
			for (dgInt32 i = 0; i < 6; i ++) {
				for (dgInt32 j = 0; j < 6; j ++) {
					dgFloat32 acc = dgFloat32 (0.0f);
					for (dgInt32 k = 0; k < 6; k ++) {
						dgFloat32 a = A[i][k];
						dgFloat32 b = B[k][j];
						acc += a * b;
					}
					dst[i][j] = acc;
				}
			}
		}

		void Trace (dgInt32 n) const
		{
			dgTrace(("\n"));
			for (dgInt32 i = 0; i < n; i ++) {
				for (dgInt32 j = 0; j < 6; j ++) {
					dgTrace(("%f, ", m_rows[i][j]));
				}
				dgTrace(("\n"));
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
*/

	dgSkeletonGraph(dgMemoryAllocator* const allocator, dgDynamicBody* const root)
		:m_parent(NULL)
		,m_body(root)
		,m_joint(NULL)
		//,m_data(NULL)
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
		//,m_data(NULL)
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
		//,m_data(NULL)
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

	virtual void SetPriority(dgUnsigned32 priority) const
	{
		if (m_joint) {
			m_joint->m_priority = priority;
		}
	}

/*
	void Init (dgData* const buffer, dgJointInfo* const jointInfoArray, dgJacobianMatrixElement* const matrixRow) 
	{
		m_data = buffer;

		m_data->m_diagonal.SetZero();
		m_data->m_offDiagonal.SetZero();
		if (m_body) {
			m_diaginalDof = 6;
			dgFloat32 mass = m_body->GetMass().m_w;
			dgMatrix inertia (m_body->CalculateInertiaMatrix());
			
			for (dgInt32 i = 0; i < 3; i++) {
				m_data->m_diagonal[i][i] = mass;
				for (dgInt32 j = 0; j < 3; j++) {
					m_data->m_diagonal[i + 3][j + 3] = inertia[i][j];
				}
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
					for (dgInt32 j = 0; j < 3; j ++) {
						m_data->m_offDiagonal[i][j + 0] = matrixRow[index].m_JMinv.m_jacobianM0.m_linear[j];
						m_data->m_offDiagonal[i][j + 3] = matrixRow[index].m_JMinv.m_jacobianM0.m_angular[j];
					}
				}
			}
		} else {
			dgJointInfo* const jointInfo = &jointInfoArray[m_joint->m_index];

			dgAssert (m_parent);
			dgAssert (jointInfo->m_joint == m_joint);
			dgAssert (jointInfo->m_joint->GetBody1() == m_parent->m_body);
			m_diaginalDof = jointInfo->m_pairCount;
			m_jacobialDof = jointInfo->m_pairCount;
		
			for (dgInt32 i = 0; i < m_diaginalDof; i++) {
				dgInt32 index = jointInfo->m_pairStart + i;
				for (dgInt32 j = 0; j < 3; j ++) {
					m_data->m_offDiagonal[i][j + 0] = matrixRow[index].m_JMinv.m_jacobianM1.m_linear[j];
					m_data->m_offDiagonal[i][j + 3] = matrixRow[index].m_JMinv.m_jacobianM1.m_angular[j];
				}
			}
		}
	}

	void CalculateOffDiagonalBlock ()
	{
		if (m_body) {
			for (dgInt32 i = 0; i < m_jacobialDof; i ++) {
				m_data->m_invDiagonal.MultiplyMatrix6x6TimeJacobianTransposed (m_data->m_offDiagonal[i]);
			}
		} else {
			dgSpacialMatrix copy;
			for (dgInt32 i = 0; i < m_jacobialDof; i ++) {
				copy[i] = m_data->m_offDiagonal[i];
			}
			m_data->m_offDiagonal.SetZero();
			for (dgInt32 i = 0; i < m_jacobialDof; i ++) {
				const dgSpacialVector& jacobian = copy[i];
				const dgSpacialVector& invDiagonalRow = m_data->m_invDiagonal[i];
				for (dgInt32 j = 0; j < m_jacobialDof; j ++) {
					dgVector val (invDiagonalRow[j]);
					jacobian.ScaleAdd(val, m_data->m_offDiagonal[j], m_data->m_offDiagonal[j]);
				}
			}
		}
	}

	virtual void CalculateDiagonal(dgSkeletonGraph* const child, dgJointInfo* const jointInfoArray)
	{
		dgSpacialMatrix tmp;
		const dgSpacialMatrix& childDiagonal = child->m_data->m_diagonal;
		if (m_body) {
			dgSpacialMatrix copy;
			copy.SetZero();

			const dgSpacialMatrix& jacobianMatrix = child->m_data->m_offDiagonal;
			for (dgInt32 i = 0; i < child->m_jacobialDof; i++) {
				const dgSpacialVector& jacobian = jacobianMatrix[i];
				for (dgInt32 j = 0; j < child->m_jacobialDof; j++) {
					dgAssert (dgAbsf (childDiagonal[i][j] - childDiagonal[j][i]) < dgFloat32 (1.0e-5f));
					dgVector val(childDiagonal[i][j]);
					jacobian.ScaleAdd (val, copy[j], copy[j]);
				}
			}

//childDiagonal.Trace(child->m_jacobialDof);
//copy.Trace (child->m_jacobialDof);
//jacobianMatrix.Trace (child->m_jacobialDof);

			dgAssert (m_diaginalDof == 6);
			dgSpacialMatrix& diagonal = m_data->m_diagonal;
			for (dgInt32 i = 0; i < child->m_jacobialDof; i++) {
				const dgSpacialVector& Jacobian = copy[i];
				const dgSpacialVector& JacobianTranspose = jacobianMatrix[i];
				for (dgInt32 j = 0; j < 6; j ++) {
					dgFloat32 val (-Jacobian[j]);
					JacobianTranspose.ScaleAdd (val, diagonal[j], diagonal[j]);
				}
			}
//diagonal.Trace(6);

			#ifdef _DEBUG
 			for (dgInt32 i = 0; i < 6; i++) {
				for (dgInt32 k = 0; k < i; k++) {
					dgFloat32 a = diagonal[i][k];
					dgFloat32 b = diagonal[k][i];
					dgAssert(dgAbsf(a - b) < dgFloat32(1.0e-5f));
				}
			}
			#endif
		} else {
			dgAssert (child->m_body);
			const dgSpacialMatrix& jacobianTransposed = child->m_data->m_offDiagonal;
			for (dgInt32 i = 0; i < m_jacobialDof; i++) {
				tmp[i] = jacobianTransposed[i];
				childDiagonal.MultiplyMatrix6x6TimeJacobianTransposed(tmp[i]);
			}

			dgSpacialMatrix& diagonal = m_data->m_diagonal;
			for (dgInt32 i = 0; i < m_jacobialDof; i++) {
				diagonal[i][i] -= jacobianTransposed[i].DotProduct(tmp[i]);
				for (dgInt32 j = 0; j < i; j ++) {
					dgFloat32 a = jacobianTransposed[i].DotProduct(tmp[j]);
					diagonal[i][j] -= a;
					diagonal[j][i] -= a;
				}
			}
		}
	}

	virtual void CalculateDiagonalInverse ()
	{
		m_data->m_invDiagonal.Inverse (m_data->m_diagonal, m_diaginalDof);
	}
*/
	dgSkeletonGraph* m_parent;
	dgDynamicBody* m_body;
	dgBilateralConstraint* m_joint;
	//dgData* m_data;
	dgList<dgSkeletonGraph*> m_children;
	dgInt32 m_index;
	dgInt16 m_diaginalDof;
	dgInt16 m_jacobialDof;
};


dgSkeletonContainer::dgSkeletonContainer (dgDynamicBody* const rootBody)
	:m_solverData(NULL)
	,m_skeleton(new (rootBody->GetWorld()->GetAllocator()) dgSkeletonGraph (rootBody->GetWorld()->GetAllocator(), rootBody))
	,m_jointArray(NULL)
	,m_bottomTopOrder(NULL)
//	,m_topBottomOrder(NULL)
	,m_id(m_uniqueID)
	,m_nodeCount(1)
{
	m_uniqueID ++;
}

dgSkeletonContainer::~dgSkeletonContainer ()
{
	dgMemoryAllocator* const allocator = m_skeleton->m_body->GetWorld()->GetAllocator();
	if (m_jointArray) {
		allocator->Free(m_jointArray);
	}
	delete m_skeleton;
}

void dgSkeletonContainer::ResetUniqueId(dgInt32 id)
{
	m_uniqueID = 10;
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
	dgAssert (0);
	return 0;
//	dgInt32 blocksize = sizeof(dgSkeletonGraph::dgData);
//	return blocksize * m_nodeCount;
}

void dgSkeletonContainer::SortGraph (dgSkeletonGraph* const root, dgSkeletonGraph* const parent, dgInt32& index)
{
	dgAssert (0);
/*
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
*/
}

void dgSkeletonContainer::Finalize ()
{
	dgAssert (m_nodeCount >= 1);

	dgInt32 jointCount = (m_nodeCount - 1) / 2;
	dgMemoryAllocator* const allocator = m_skeleton->m_body->GetWorld()->GetAllocator();
	m_jointArray = (dgSkeletonGraph**) allocator->Malloc ((2 * m_nodeCount + jointCount) * sizeof (dgSkeletonGraph*));
	m_bottomTopOrder = &m_jointArray[jointCount];
//	m_topBottomOrder = &m_bottomTopOrder[m_nodeCount];
	
	dgAssert (0);

	dgInt32 index = 0;
	SortGraph (m_skeleton, NULL, index);

	dgAssert (index == m_nodeCount);
	index = 0;
	for (dgInt32 i = 0; i < m_nodeCount; i ++) {
		dgSkeletonGraph* const node = m_bottomTopOrder[i];
		if (node->m_joint) {
			m_jointArray[index] = node;
			index ++;
		}
	}
}



bool dgSkeletonContainer::Sanity() const
{
/*
	for (dgInt32 i = 0; i < m_nodeCount; i++) {
		dgSkeletonGraph* const node = m_bottomTopOrder[i];
		dgInt32 n = node->m_diaginalDof;

		dgSkeletonGraph::dgSpacialMatrix identity;
		node->m_data->m_diagonal.Multiply(node->m_data->m_invDiagonal, identity);

		for (dgInt32 j = 0; j < n; j++) {
			for (dgInt32 k = 0; k < n; k++) {
				if (dgAbsf (node->m_data->m_diagonal[j][k] - node->m_data->m_diagonal[k][j]) > dgFloat32 (1.0e-5f)) {
					return false;
				}
				dgTrace (("%f, ", node->m_data->m_diagonal[j][k]));
			}

			dgTrace (("    "));
			for (dgInt32 k = 0; k < n; k++) {
				dgTrace (("%f, ", node->m_data->m_invDiagonal[j][k]));
			}

			dgTrace (("    "));
			for (dgInt32 k = 0; k < n; k++) {
				dgTrace(("%f, ", identity[j][k]));
			}
			dgTrace (("\n"));

		}
		dgTrace (("\n"));
	}
*/
	return true;
}


void dgSkeletonContainer::InitMassMatrix (char* const buffer, dgJointInfo* const jointInfoArray, dgJacobianMatrixElement* const matrixRow)
{
/*
	dgSkeletonGraph::dgData* const ptr = (dgSkeletonGraph::dgData*) buffer;
	dgAssert( (dgUnsigned64(buffer) & 0x0f) == 0);
	for (dgInt32 i = 0; i < m_nodeCount; i++) {
		dgSkeletonGraph* const node = m_bottomTopOrder[i];
		node->Init(&ptr[i], jointInfoArray, matrixRow);
	}

	dgAssert (Sanity ());
	for (dgInt32 i = 0; i < m_nodeCount; i++) {
		dgSkeletonGraph* const node = m_bottomTopOrder[i];
		for (dgList<dgSkeletonGraph*>::dgListNode* child = node->m_children.GetFirst(); child; child = child->GetNext()) {
			node->CalculateDiagonal(child->GetInfo(), jointInfoArray);
			dgAssert (Sanity ());
		}
		node->CalculateDiagonalInverse();
		dgAssert (Sanity ());
		if (node->m_parent) {
			node->CalculateOffDiagonalBlock();
			dgAssert (Sanity ());
		}
	}
dgAssert (0);
*/
}

dgFloat32 dgSkeletonContainer::CalculateJointForce(dgJointInfo* const jointInfoArray, const dgBodyInfo* const bodyArray, dgJacobian* const internalForces, dgJacobianMatrixElement* const matrixRow) const
{
	dgAssert (0);
/*
	dgInt32 jointCount = (m_nodeCount - 1) / 2;
//	dgJointInfo* const constraintArrayPtr = (dgJointInfo*)&world->m_jointsMemory[0];
//	dgJointInfo* const constraintArray = &constraintArrayPtr[island->m_jointStart];
	
	for (dgInt32 i = 0; i < jointCount; i++) {
		dgJointInfo* const jointInfo = &jointInfoArray[i];
		const dgInt32 m0 = jointInfo->m_m0;
		const dgInt32 m1 = jointInfo->m_m1;
		internalForces[m0].m_linear = dgVector::m_zero;
		internalForces[m1].m_angular = dgVector::m_zero;
	}

	for (dgInt32 i = 0; i < jointCount; i++) {
		dgJacobian y0;
		dgJacobian y1;
		y0.m_linear = dgVector::m_zero;
		y0.m_angular = dgVector::m_zero;
		y1.m_linear = dgVector::m_zero;
		y1.m_angular = dgVector::m_zero;

		dgJointInfo* const jointInfo = &jointInfoArray[i];
		const dgInt32 first = jointInfo->m_pairStart;
		const dgInt32 count = jointInfo->m_pairActiveCount;
		for (dgInt32 j = 0; j < count; j++) {
			dgJacobianMatrixElement* const row = &matrixRow[j + first];
			dgVector val(row->m_force);
			dgAssert(dgCheckFloat(row->m_force));
			y0.m_linear += row->m_Jt.m_jacobianM0.m_linear.CompProduct4(val);
			y0.m_angular += row->m_Jt.m_jacobianM0.m_angular.CompProduct4(val);
			y1.m_linear += row->m_Jt.m_jacobianM1.m_linear.CompProduct4(val);
			y1.m_angular += row->m_Jt.m_jacobianM1.m_angular.CompProduct4(val);
		}

		const dgInt32 m0 = jointInfo->m_m0;
		const dgInt32 m1 = jointInfo->m_m1;
		internalForces[m0].m_linear += y0.m_linear;
		internalForces[m0].m_angular += y0.m_angular;
		internalForces[m1].m_linear += y1.m_linear;
		internalForces[m1].m_angular += y1.m_angular;
	}

	dgFloat32 accNorm = dgFloat32(0.0f);
	dgFloat32 akNum = dgFloat32(0.0f);
	for (dgInt32 i = 0; i < jointCount; i++) {
		dgJointInfo* const jointInfo = &jointInfoArray[i];
		const dgInt32 first = jointInfo[i].m_pairStart;
		const dgInt32 count = jointInfo[i].m_pairCount;
		const dgInt32 m0 = jointInfo[i].m_m0;
		const dgInt32 m1 = jointInfo[i].m_m1;
		const dgJacobian& y0 = internalForces[m0];
		const dgJacobian& y1 = internalForces[m1];
		for (dgInt32 j = 0; j < count; j++) {
			dgJacobianMatrixElement* const row = &matrixRow[j + first];
			dgVector acc (row->m_JMinv.m_jacobianM0.m_linear.CompProduct4(y0.m_linear) + row->m_JMinv.m_jacobianM0.m_angular.CompProduct4(y0.m_angular) +
						  row->m_JMinv.m_jacobianM1.m_linear.CompProduct4(y1.m_linear) + row->m_JMinv.m_jacobianM1.m_angular.CompProduct4(y1.m_angular));
			acc = dgVector(row->m_coordenateAccel) - acc.AddHorizontal();
			//acc.StoreScalar(&row->m_accel);
			row->m_accel = acc.GetScalar();
			row->m_deltaAccel = acc.GetScalar() * row->m_invDJMinvJt;
			accNorm += acc.Abs().GetScalar();
			akNum += acc.GetScalar() * acc.GetScalar() * row->m_invDJMinvJt;
		}
	}

	dgFloat32 maxAccel = accNorm;
	const dgInt32 maxPasses = 16;
	for (dgInt32 passes = 0; (passes < maxPasses) && (accNorm > dgFloat32 (1.0e-1f)); passes ++) {
		for (dgInt32 i = 0; i < jointCount; i++) {
			dgJointInfo* const jointInfo = &jointInfoArray[i];
			const dgInt32 m0 = jointInfo->m_m0;
			const dgInt32 m1 = jointInfo->m_m1;
			internalForces[m0].m_linear = dgVector::m_zero;
			internalForces[m1].m_angular = dgVector::m_zero;
		}

		for (dgInt32 i = 0; i < jointCount; i++) {
			dgJacobian y0;
			dgJacobian y1;
			y0.m_linear = dgVector::m_zero;
			y0.m_angular = dgVector::m_zero;
			y1.m_linear = dgVector::m_zero;
			y1.m_angular = dgVector::m_zero;

			dgJointInfo* const jointInfo = &jointInfoArray[i];
			const dgInt32 first = jointInfo->m_pairStart;
			const dgInt32 count = jointInfo->m_pairActiveCount;
			for (dgInt32 j = 0; j < count; j++) {
				dgJacobianMatrixElement* const row = &matrixRow[j + first];
				dgVector val(row->m_deltaAccel);
				dgAssert(dgCheckFloat(row->m_deltaAccel));
				y0.m_linear += row->m_Jt.m_jacobianM0.m_linear.CompProduct4(val);
				y0.m_angular += row->m_Jt.m_jacobianM0.m_angular.CompProduct4(val);
				y1.m_linear += row->m_Jt.m_jacobianM1.m_linear.CompProduct4(val);
				y1.m_angular += row->m_Jt.m_jacobianM1.m_angular.CompProduct4(val);
			}

			const dgInt32 m0 = jointInfo->m_m0;
			const dgInt32 m1 = jointInfo->m_m1;
			internalForces[m0].m_linear += y0.m_linear;
			internalForces[m0].m_angular += y0.m_angular;
			internalForces[m1].m_linear += y1.m_linear;
			internalForces[m1].m_angular += y1.m_angular;
		}

		dgFloat32 akDen = dgFloat32 (0.0f);
		for (dgInt32 i = 0; i < jointCount; i++) {
			dgJointInfo* const jointInfo = &jointInfoArray[i];
			const dgInt32 first = jointInfo[i].m_pairStart;
			const dgInt32 count = jointInfo[i].m_pairCount;
			const dgInt32 m0 = jointInfo[i].m_m0;
			const dgInt32 m1 = jointInfo[i].m_m1;
			const dgJacobian& y0 = internalForces[m0];
			const dgJacobian& y1 = internalForces[m1];
			for (dgInt32 j = 0; j < count; j++) {
				dgJacobianMatrixElement* const row = &matrixRow[j + first];
				dgVector acc(row->m_JMinv.m_jacobianM0.m_linear.CompProduct4(y0.m_linear) + row->m_JMinv.m_jacobianM0.m_angular.CompProduct4(y0.m_angular) +
							 row->m_JMinv.m_jacobianM1.m_linear.CompProduct4(y1.m_linear) + row->m_JMinv.m_jacobianM1.m_angular.CompProduct4(y1.m_angular));

				row->m_deltaForce = acc.AddHorizontal().GetScalar();
				akDen += row->m_deltaForce  * row->m_deltaAccel;
			}
		}
		dgFloat32 alpha = akNum / akDen;

		dgVector accelMag (dgVector::m_zero);
		for (dgInt32 i = 0; i < jointCount; i++) {
			dgJointInfo* const jointInfo = &jointInfoArray[i];
			const dgInt32 first = jointInfo[i].m_pairStart;
			const dgInt32 count = jointInfo[i].m_pairCount;
			for (dgInt32 j = 0; j < count; j++) {
				dgJacobianMatrixElement* const row = &matrixRow[j + first];
				row->m_force += alpha * row->m_deltaAccel;
				row->m_accel -= alpha * row->m_deltaForce;
				accelMag += dgVector(row->m_accel).Abs();  
			}
		}
		accNorm = accelMag.GetScalar();
		if (accNorm > dgFloat32(1.0e-1f)) {
			dgAssert(0);
		}
	}
	return maxAccel;
*/
	return 0;
}

