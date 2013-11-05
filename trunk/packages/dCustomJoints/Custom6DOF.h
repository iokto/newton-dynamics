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



// Custom6DOF.h: interface for the Custom6DOF class.
//
//////////////////////////////////////////////////////////////////////

#ifndef _CUSTOM_6DOF_H_
#define _CUSTOM_6DOF_H_

#include "CustomJoint.h"

class Custom6DOF: public CustomJoint  
{
	public:
	CUSTOM_JOINTS_API Custom6DOF (const dMatrix& pinsAndPivotChildFrame, const dMatrix& pinsAndPivotParentFrame, NewtonBody* const child, NewtonBody* const parent = NULL);
	CUSTOM_JOINTS_API virtual ~Custom6DOF();

	CUSTOM_JOINTS_API void SetLinearLimits (const dVector& minLinearLimits, const dVector& maxLinearLimits);
	CUSTOM_JOINTS_API void SetAngularLimits (const dVector& minAngularLimits, const dVector& maxAngularLimits);
	CUSTOM_JOINTS_API void GetLinearLimits (dVector& minLinearLimits, dVector& maxLinearLimits);
	CUSTOM_JOINTS_API void GetAngularLimits (dVector& minAngularLimits, dVector& maxAngularLimits);


	CUSTOM_JOINTS_API void SetReverserUniversal (int order);

	protected:
	CUSTOM_JOINTS_API virtual void SubmitConstraints (dFloat timestep, int threadIndex);
	CUSTOM_JOINTS_API virtual void GetInfo (NewtonJointRecord* const info) const;

	CUSTOM_JOINTS_API void SubmitConstraints (const dMatrix& matrix0, const dMatrix& matrix1, dFloat timestep);

	protected:
	CUSTOM_JOINTS_API dMatrix CalculateBasisAndJointAngle (const dMatrix& matrix0, const dMatrix& matrix1) const;
	CUSTOM_JOINTS_API dMatrix CalculateHinge_Angles (const dMatrix& matrix0, const dMatrix& matrix1, int x, int y, int z) const;
	CUSTOM_JOINTS_API dMatrix CalculateUniversal_Angles (const dMatrix& matrix0, const dMatrix& matrix1, int x, int y, int z) const;

	dMatrix m_localMatrix0;
	dMatrix m_localMatrix1;
	dVector m_minLinearLimits;
	dVector m_maxLinearLimits;
	dVector m_minAngularLimits;
	dVector m_maxAngularLimits;
	dVector m_maxMaxLinearErrorRamp;
	dVector m_maxMaxAngularErrorRamp;
	bool m_reverseUniversal;
};

#endif // !defined(AFX_Custom6DOF_H__B631F556_B7D7_F85ECF3E9ADE_H)

