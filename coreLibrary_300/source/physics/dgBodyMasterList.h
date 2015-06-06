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

#ifndef __DGBODYMASTER_LIST__
#define __DGBODYMASTER_LIST__


class dgBody;
class dgContact;
class dgConstraint;
class dgBilateralConstraint;

class dgBodyMasterListCell
{
	public:
	dgBodyMasterListCell ()
	{
		dgAssert (0);
	}
	dgBodyMasterListCell (dgConstraint* const joint, dgBody* const body)
		:m_joint(joint)
		,m_bodyNode(body)
		,m_lru(0)
	{
	}
	dgConstraint* m_joint;
	dgBody* m_bodyNode;
	dgInt32 m_lru;
};

class dgBodyMasterListRow: public dgList<dgBodyMasterListCell>
{
	public:
	dgBodyMasterListRow ();
	~dgBodyMasterListRow ();

	dgBody* GetBody() const	{ return m_body;}

	private:
	void SortList();
	void SetBody(dgBody* const body) {m_body = body;}

	dgContact* FindContactJoint (const dgBody* const body) const;
	dgBilateralConstraint* FindBilateralJoint (const dgBody* const body) const;

	dgListNode* AddContactJoint (dgConstraint* const joint, dgBody* const body);
	dgListNode* AddBilateralJoint (dgConstraint* const joint, dgBody* const body);

	void RemoveContactJoint (dgListNode* const link);
	void RemoveBilateralJoint (dgListNode* const link);

	void RemoveAllJoints ();
	
	void SetAcceleratedSearch();

	dgBody* m_body;
	dgListNode* m_acceleratedSearch[3];
	dgInt32 m_contactCount;
	dgInt32 m_lru;
	static dgInt32 m_contactCountReversal[];
	friend class dgBodyMasterList;
	friend class dgWorldDynamicUpdate;
};

class dgBodyMasterList: public dgList<dgBodyMasterListRow>
{
	public:
	dgBodyMasterList (dgMemoryAllocator* const allocator);
	~dgBodyMasterList ();

	void AddBody (dgBody* const body);
	void RemoveBody (dgBody* const body);
	void RemoveConstraint (dgConstraint* const constraint);
	void AttachConstraint (dgConstraint* const constraint, dgBody* const body0, dgBody* const body1);

	dgContact* FindContactJoint (const dgBody* body0, const dgBody* body1) const;
	dgBilateralConstraint* FindBilateralJoint (const dgBody* body0, const dgBody* body1) const;

	dgBodyMasterListRow::dgListNode* FindConstraintLink (const dgBody* const body0, const dgBody* const body1) const;
	//dgBodyMasterListRow::dgListNode* FindConstraintLinkNext (const dgBodyMasterListRow::dgListNode* const me, const dgBody* const body) const;
	dgUnsigned32 MakeSortMask(const dgBody* const body) const;
	void SortMasterList();

	void ResetColor ();

	public:
	dgTree<int, dgBody*> m_disableBodies;
	dgInt32 m_lru;
	dgInt32 m_deformableCount;
	dgInt32 m_constraintCount;
};

#endif
