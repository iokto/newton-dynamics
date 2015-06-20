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

#ifndef _DG_ACYCLIC_ARTCULATED_CONTAINER_CONSTRAINT_H__
#define _DG_ACYCLIC_ARTCULATED_CONTAINER_CONSTRAINT_H__

#include "dgConstraint.h"

class dgDynamicBody;

class dgAcyclicContainer
{
	public:
	class dgAcyclicGraph
	{
		public: 
		DG_CLASS_ALLOCATOR(allocator)
		dgAcyclicGraph ();
		dgAcyclicGraph (dgMemoryAllocator* const allocator, dgDynamicBody* const body, dgDynamicBody* const parent);
		~dgAcyclicGraph();

		dgDynamicBody* m_body;
		dgDynamicBody* m_parent;
		dgList<dgAcyclicGraph*> m_children;
	};

	class dgAcyclicJointBodyPair
	{
		public:
		dgDynamicBody* m_Body;
		dgBilateralConstraint* m_joint;
	};

	DG_CLASS_ALLOCATOR(allocator)
	dgAcyclicContainer(dgDynamicBody* const rootBody);
	~dgAcyclicContainer();

	dgInt32 GetId () const {return m_id;}
	void AddChild (dgBody* const parent, dgBody* const child);
	void AddChild (dgDynamicBody* const parent, dgDynamicBody* const child);

	void Finalize ();
	
	protected:
	
	dgInt32 NodeCount () const;
	dgAcyclicGraph* FindNode (dgDynamicBody* const node) const;
	void SortGraph (dgAcyclicGraph* const root, dgAcyclicGraph* const parent, const dgInt32 count, dgInt32& index);
	virtual void SetDestructorCallback(OnConstraintDestroy destructor) {dgAssert (0);}

	virtual dgUnsigned32 JacobianDerivative(dgContraintDescritor& params) {return 0;}
	virtual void JointAccelerations(dgJointAccelerationDecriptor* const params) {}
	virtual void JointVelocityCorrection(dgJointAccelerationDecriptor* const params) {dgAssert (0);}

	dgAcyclicGraph m_skeleton;
	dgAcyclicJointBodyPair* m_upDownOrder;
	dgAcyclicJointBodyPair* m_downUpOrder;
	dgInt32 m_id;
	static dgInt32 m_uniqueID;
};

#endif

