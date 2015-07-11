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

#ifndef _DG_SKELETON_CONTAINER_H__
#define _DG_SKELETON_CONTAINER_H__

#define DG_SKELETON_BIT_SHIFT_KEY	16
#define DG_SKELETON_BASEW_UNIQUE_ID	10

#include "dgConstraint.h"

class dgDynamicBody;
class dgSkeletonContainer;
typedef void (dgApi *dgOnSkeletonContainerDestroyCallback) (dgSkeletonContainer* const me);

class dgSkeletonContainer
{
	public:
	class dgSkeletonGraph;

	DG_CLASS_ALLOCATOR(allocator)
	dgSkeletonContainer(dgWorld* const world, dgDynamicBody* const rootBody);
	~dgSkeletonContainer();

	dgWorld* GetWorld() const; 
	void SetDestructorCallback (dgOnSkeletonContainerDestroyCallback destructor);

	dgInt32 GetId () const {return m_id;}
	dgSkeletonGraph* AddChild (dgBody* const child, dgBody* const parent);
	void AddJointList (dgInt32 count, dgBilateralConstraint** const array);
	
	void Finalize ();
	dgInt32 GetJointCount () const {return (m_nodeCount - 1) / 2;}

	void InitMassMatrix (dgJointInfo* const jointInfoArray, dgJacobianMatrixElement* const matrixRow);
	dgFloat32 CalculateJointForce (dgJointInfo* const jointInfo, const dgBodyInfo* const bodyArray, dgJacobian* const internalForces, dgJacobianMatrixElement* const matrixRow) const;

	private:
	dgSkeletonGraph* FindNode (dgDynamicBody* const node) const;
	dgSkeletonGraph* AddChild (dgDynamicBody* const child, dgDynamicBody* const parent);
	
	void SolveFoward () const;
	void SolveBackward () const;
	dgFloat32 SolveUnilaterals (dgJointInfo* const jointInfoArray, const dgBodyInfo* const bodyArray, dgJacobian* const internalForces, dgJacobianMatrixElement* const matrixRow) const;

	void SortGraph (dgSkeletonGraph* const root, dgSkeletonGraph* const parent, dgInt32& index);
	static void ResetUniqueId(dgInt32 id);

	dgWorld* m_world;
	dgSkeletonGraph* m_skeleton;
	dgSkeletonGraph** m_nodesOrder;
	dgOnSkeletonContainerDestroyCallback m_destructor;
	dgInt32 m_id;
	dgInt32 m_nodeCount;
	static dgInt32 m_uniqueID;

	friend class dgWorld;
};

#endif

