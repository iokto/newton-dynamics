/* Copyright (c) <2003-2013> <Julio Jerez, Newton Game Dynamics>
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

#ifndef _D_NEWTON_JOINT_H_
#define _D_NEWTON_JOINT_H_

#include "dStdAfxNewton.h"
#include "dNewtonBody.h"
#include "dNewtonAlloc.h"
#include "dNewtonDynamicBody.h"


class dNewtonJoint: public dNewtonAlloc
{
	public:
	enum dJointType
	{
		m_ballAndSocket,
		m_hinge,
		m_slider,
		m_universal,

		m_hingeActuator,
		m_sliderActuator,
		m_universalActuator,
		m_unknown,
	};
	
	CNEWTON_API virtual ~dNewtonJoint();

	CNEWTON_API dNewtonDynamicBody* GetBody0 () const;
	CNEWTON_API dNewtonDynamicBody* GetBody1 () const;


	protected:
	CNEWTON_API dNewtonJoint(dJointType type)
		:m_type(type) 
		,m_joint(NULL)
	{
	}

	CNEWTON_API void SetJoint(CustomJoint* const joint);
	CNEWTON_API virtual void OnSubmitConstraint (dFloat timestep, int threadIndex)
	{
	}

	private:
	CNEWTON_API static void OnJointDestroyCallback (const NewtonUserJoint* const me);
	CNEWTON_API static void OnSubmitConstraintCallback (const NewtonUserJoint* const me, dFloat timestep, int threadIndex);

	protected:
	dJointType m_type;
	CustomJoint* m_joint;
};

class dNewtonBallAndSocketJoint: public dNewtonJoint 
{
	public:
	CNEWTON_API dNewtonBallAndSocketJoint(const dFloat* const pinAndPivotFrame, dNewtonDynamicBody* const body0, dNewtonDynamicBody* const body1 = NULL)
		:dNewtonJoint(m_ballAndSocket)
	{
		SetJoint (new CustomBallAndSocket (dMatrix(pinAndPivotFrame), body0->GetNewtonBody(), body1 ? body1->GetNewtonBody() : NULL));
	}
};

class dNewtonHingeJoint: public dNewtonJoint 
{
	public:
	CNEWTON_API dNewtonHingeJoint(const dFloat* const pinAndPivotFrame, dNewtonDynamicBody* const body0, dNewtonDynamicBody* const body1 = NULL)
		:dNewtonJoint(m_hinge)
	{
		SetJoint (new CustomHinge (dMatrix(pinAndPivotFrame), body0->GetNewtonBody(), body1 ? body1->GetNewtonBody() : NULL));
	}

	CNEWTON_API dFloat GetFriction () const
	{
		return ((CustomHinge*)m_joint)->GetFriction();
	}

	CNEWTON_API void SetFriction (dFloat friction)
	{
		((CustomHinge*)m_joint)->SetFriction(friction);
	}
};

class dNewtonSliderJoint: public dNewtonJoint 
{
	public:
	CNEWTON_API dNewtonSliderJoint(const dFloat* const pinAndPivotFrame, dNewtonDynamicBody* const body0, dNewtonDynamicBody* const body1 = NULL)
		:dNewtonJoint(m_slider)
	{
		SetJoint (new CustomSlider (dMatrix(pinAndPivotFrame), body0->GetNewtonBody(), body1 ? body1->GetNewtonBody() : NULL));
	}
};


class dNewtonUniversalJoint: public dNewtonJoint 
{
	public:
	CNEWTON_API dNewtonUniversalJoint(const dFloat* const pinAndPivotFrame, dNewtonDynamicBody* const body0, dNewtonDynamicBody* const body1 = NULL)
		:dNewtonJoint(m_universal)
	{
		SetJoint (new CustomUniversal (dMatrix(pinAndPivotFrame), body0->GetNewtonBody(), body1 ? body1->GetNewtonBody() : NULL));
	}

	NEWTON_API void EnableLimit_0(bool state)
	{
		((CustomUniversal*) m_joint)->EnableLimit_0(state);
	}

	NEWTON_API void EnableLimit_1(bool state)
	{
		((CustomUniversal*) m_joint)->EnableLimit_1(state);
	}


	NEWTON_API void SetLimis_0(dFloat minAngle, dFloat maxAngle)
	{
		((CustomUniversal*) m_joint)->SetLimis_0 (minAngle, maxAngle);
	}
	NEWTON_API void SetLimis_1(dFloat minAngle, dFloat maxAngle)
	{
		((CustomUniversal*) m_joint)->SetLimis_1 (minAngle, maxAngle);
	}
};



#endif