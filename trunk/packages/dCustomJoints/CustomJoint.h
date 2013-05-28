/* Copyright (c) <2009> <Newton Game Dynamics>
* 
* This software is provided 'as-is', without any express or implied
* warranty. In no event will the authors be held liable for any damages
* arising from the use of this software.
* 
* Permission is granted to anyone to use this software for any purpose,
* including commercial applications, and to alter it and redistribute it
* freely
*/


// CustomJoint.h: interface for the NewtonCustomJoint class.
//
//////////////////////////////////////////////////////////////////////

#ifndef D_NEWTON_CUSTOM_JOINT_H_
#define D_NEWTON_CUSTOM_JOINT_H_

#include "CustomJointLibraryStdAfx.h"
#include "dMathDefines.h"
#include "dVector.h"
#include "dMatrix.h"
#include "dQuaternion.h"
#include "dTree.h"
#include "dList.h"
#include "dRtti.h"
#include "dCRC.h"


struct NewtonUserJoint;
typedef void (*JointUserDestructorCallback) (const NewtonUserJoint* const me);	
typedef void (*JointUserSubmitConstraintCallback) (const NewtonUserJoint* const me, dFloat timestep, int threadIndex);


// this is the base class to implement custom joints, it is not a joint it just provide functionality
// for the user to implement it own joints
class NEWTON_API CustomJoint  
{
	public:
	struct AngularIntegration
	{
		AngularIntegration()
		{
			m_angle = 0.0f;
		}

		dFloat CalculateJointAngle (dFloat newAngleCos, dFloat newAngleSin)
		{
			// the joint angle can be determine by getting the angle between any two non parallel vectors
			dFloat sinJointAngle = dSin (m_angle);
			dFloat cosJointAngle = dCos (m_angle);

			dFloat sin_da = newAngleSin * cosJointAngle - newAngleCos * sinJointAngle; 
			dFloat cos_da = newAngleCos * cosJointAngle + newAngleSin * sinJointAngle; 

			m_angle += dAtan2 (sin_da, cos_da);
			return  m_angle;
		}
		dFloat m_angle;
	};

	CustomJoint();
	CustomJoint(int maxDOF, NewtonBody* const body0, NewtonBody* const body1);
	virtual ~CustomJoint();


	void *operator new (size_t size);
	void operator delete (void *ptr);


	void SetBodiesCollisionState (int state);
	int GetBodiesCollisionState () const;

	NewtonBody* GetBody0 () const;
	NewtonBody* GetBody1 () const;
	NewtonJoint* GetJoint () const;


	// the application needs to implement this function for serialization
	virtual void GetInfo (NewtonJointRecord* const info) const;


	// these member function are only used by the C interface or for hooking callback to customize a particular 
	// joint without deriving a new one
	// note: this is not a extension of a virtual function, DO NOT CALL the base class SubmitConstraints!! 
	void SetUserData (void* userData) {m_userData = userData;}
	void* GetUserData () const {return m_userData;}
	void SetUserDestructorCallback (JointUserDestructorCallback callback) {m_userDestructor = callback;}
	void SetUserSubmintConstraintCallback (JointUserSubmitConstraintCallback callback) {m_userConstrationCallback = callback;}

	private:
	// this are the callback needed to have transparent c++ method interfaces 
	static void Destructor (const NewtonJoint* me);	
	static void SubmitConstraints (const NewtonJoint* const me, dFloat timestep, int threadIndex);
	static void GetInfo (const NewtonJoint* const me, NewtonJointRecord* info);
	

	protected:
	void Init (int maxDOF, NewtonBody* const body0, NewtonBody* const body1);
	// the application needs to implement this function for each derived joint. See examples for more detail
	virtual void SubmitConstraints (dFloat timestep, int threadIndex);
	void CalculateGlobalMatrix (const dMatrix& localMatrix0, const dMatrix& localMatrix1, dMatrix& matrix0, dMatrix& matrix1) const;
	void CalculateLocalMatrix (const dMatrix& pinsAndPivotFrame, dMatrix& localMatrix0, dMatrix& localMatrix1) const;


	void* m_userData;
	NewtonBody* m_body0;
	NewtonBody* m_body1;
	NewtonJoint* m_joint;
	NewtonWorld* m_world;
	JointUserDestructorCallback m_userDestructor;
	JointUserSubmitConstraintCallback m_userConstrationCallback;
	int m_maxDof;
	int m_autoDestroy;

	dRttiRootClassSupportDeclare(CustomJoint);
};



#endif 

