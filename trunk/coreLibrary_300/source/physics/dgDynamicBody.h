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

#ifndef _DG_DYNAMIC_BODY_H_
#define _DG_DYNAMIC_BODY_H_

#include "dgPhysicsStdafx.h"
#include "dgBody.h"

#define DG_MIN_SPEED_ATT	dgFloat32(0.0f)
#define DG_MAX_SPEED_ATT	dgFloat32(0.02f)
#define DG_FREEZE_MAG		dgFloat32(0.1f)
#define DG_FREEZE_MAG2		dgFloat32(DG_FREEZE_MAG * DG_FREEZE_MAG)

#define DG_ErrTolerance		(1.0e-2f)
#define DG_ErrTolerance2	(DG_ErrTolerance * DG_ErrTolerance)



//DG_MSC_VECTOR_ALIGMENT
DG_MSC_VECTOR_AVX_ALIGMENT
class dgDynamicBody : public dgBody 
{
	public:
	dgDynamicBody();
	dgDynamicBody (dgWorld* const world, const dgTree<const dgCollision*, dgInt32>* const collisionNode, dgDeserialize serializeCallback, void* const userData);
	virtual ~dgDynamicBody ();

	
	virtual const dgVector& GetForce() const;
	virtual const dgVector& GetTorque() const;
	virtual const dgVector& GetNetForce() const;
	virtual const dgVector& GetNetTorque() const;
	
	virtual void AddForce (const dgVector& force);
	virtual void AddTorque (const dgVector& torque);
	virtual void SetForce (const dgVector& force);
	virtual void SetTorque (const dgVector& torque);

	virtual dgFloat32 GetLinearDamping () const;
	virtual dgVector GetAngularDamping () const;
	virtual void SetLinearDamping (dgFloat32 linearDamp);
	virtual void SetAngularDamping (const dgVector& angularDamp);

	virtual const dgMatrix& GetInertiaMatrix () const;
	virtual dgMatrix CalculateInertiaMatrix () const;
	virtual dgMatrix CalculateInvInertiaMatrix () const;

	virtual void SetMassMatrix (dgFloat32 mass, dgFloat32 Ix, dgFloat32 Iy, dgFloat32 Iz);
	virtual void SetMassProperties (dgFloat32 mass, const dgCollisionInstance* const collision);

	virtual void AttachCollision (dgCollisionInstance* const collision);
	virtual dgVector PredictLinearVelocity(dgFloat32 timestep) const;
	virtual dgVector PredictAngularVelocity(dgFloat32 timestep) const;

	virtual void InvalidateCache();
	virtual dgVector GetMass() const;
	virtual dgVector GetInvMass() const;
	virtual void SetMatrixIgnoreSleep(const dgMatrix& matrix);

	virtual bool IsInEquilibrium () const;
	virtual void SetCollidable (bool state) {}

	virtual void ApplyExtenalForces (dgFloat32 timestep, dgInt32 threadIndex);
	virtual OnApplyExtForceAndTorque GetExtForceAndTorqueCallback () const;
	virtual void SetExtForceAndTorqueCallback (OnApplyExtForceAndTorque callback);

	
	virtual void AddImpulse (const dgVector& pointVeloc, const dgVector& pointPosit);
	virtual void ApplyImpulsePair (const dgVector& linearImpulse, const dgVector& angularImpulse);
	virtual void ApplyImpulsesAtPoint (dgInt32 count, dgInt32 strideInBytes, const dgFloat32* const impulseArray, const dgFloat32* const pointArray);

	virtual void Serialize (const dgTree<dgInt32, const dgCollision*>* const collisionNode, dgSerialize serializeCallback, void* const userData);
/*
	virtual dgConstraint* GetFirstJoint() const;
	virtual dgConstraint* GetNextJoint(dgConstraint* const joint) const;
	virtual dgConstraint* GetFirstContact() const;
	virtual dgConstraint* GetNextContact(dgConstraint* const joint) const;
	virtual void AddBuoyancyForce (dgFloat32 fluidDensity, dgFloat32 fluidLinearViscousity, dgFloat32 fluidAngularViscousity,  const dgVector& gravityVector, GetBuoyancyPlane buoyancyPlane, void* const context);
	virtual dgVector CalculateInverseDynamicForce (const dgVector& desiredVeloc, dgFloat32 timestep) const;
*/

	private:
	
	dgVector GetApparentMass() const;
	void SetAparentMassMatrix (const dgVector& massMatrix);
	void AddDampingAcceleration();
	void AddDampingAccelerationSimd();
	void CalcInvInertiaMatrix ();
	void CalcInvInertiaMatrixSimd ();


	dgMatrix m_invWorldInertiaMatrix;
	dgVector m_mass;
	dgVector m_invMass;
	dgVector m_accel;
	dgVector m_alpha;
	dgVector m_netForce;
	dgVector m_netTorque;
	dgVector m_prevExternalForce;
	dgVector m_prevExternalTorque;
	dgVector m_dampCoef;
	dgVector m_aparentMass;

	dgInt32 m_sleepingCounter;
	dgUnsigned32 m_dynamicsLru;	
	dgUnsigned32 m_isInDestructionArrayLRU;

	OnApplyExtForceAndTorque m_applyExtForces;

	friend class dgWorld;
	friend class dgWorldDynamicUpdate;
	friend class dgBroadPhase;
} DG_GCC_VECTOR_AVX_ALIGMENT;

// *****************************************************************************
// 
//	 Implementation	
// 
// *****************************************************************************
//DG_INLINE void dgDynamicBody::SetFreezeTreshhold (dgFloat32 freezeAccel2, dgFloat32 freezeAlpha2, dgFloat32 freezeSpeed2, dgFloat32 freezeOmega2)
//{
//	m_freezeAccel2 = GetMax (freezeAccel2, dgFloat32(DG_FREEZE_MAG2));
//	m_freezeAlpha2 = GetMax (freezeAlpha2, dgFloat32(DG_FREEZE_MAG2));
//	m_freezeSpeed2 = GetMax (freezeSpeed2, dgFloat32(DG_FREEZE_MAG2));
//	m_freezeOmega2 = GetMax (freezeOmega2, dgFloat32(DG_FREEZE_MAG2));
//}

//DG_INLINE void dgDynamicBody::GetFreezeTreshhold (dgFloat32& freezeAccel2, dgFloat32& freezeAlpha2, dgFloat32& freezeSpeed2, dgFloat32& freezeOmega2) const
//{
//	freezeAccel2 = m_freezeAccel2;
//	freezeAlpha2 = m_freezeAlpha2;
//	freezeSpeed2 = m_freezeSpeed2;
//	freezeOmega2 = m_freezeOmega2;
//}

DG_INLINE const dgVector& dgDynamicBody::GetForce() const
{
	return m_accel; 
}

DG_INLINE const dgVector& dgDynamicBody::GetTorque() const
{
	return m_alpha;
}

DG_INLINE const dgVector& dgDynamicBody::GetNetForce() const
{
	return m_netForce; 
}

DG_INLINE const dgVector& dgDynamicBody::GetNetTorque() const
{
	return m_netTorque;
}




DG_INLINE const dgMatrix& dgDynamicBody::GetInertiaMatrix () const 
{
	return m_invWorldInertiaMatrix;
}


DG_INLINE dgFloat32 dgDynamicBody::GetLinearDamping () const
{
	return (m_dampCoef.m_w - DG_MIN_SPEED_ATT) / (DG_MAX_SPEED_ATT - DG_MIN_SPEED_ATT);
}

DG_INLINE dgVector dgDynamicBody::GetAngularDamping () const
{
	return dgVector ((m_dampCoef.m_x - DG_MIN_SPEED_ATT) / (DG_MAX_SPEED_ATT - DG_MIN_SPEED_ATT),
					 (m_dampCoef.m_y - DG_MIN_SPEED_ATT) / (DG_MAX_SPEED_ATT - DG_MIN_SPEED_ATT),
					 (m_dampCoef.m_z - DG_MIN_SPEED_ATT) / (DG_MAX_SPEED_ATT - DG_MIN_SPEED_ATT), dgFloat32 (0.0f));
}

DG_INLINE void dgDynamicBody::SetLinearDamping (dgFloat32 linearDamp)
{
	linearDamp = dgClamp (linearDamp, dgFloat32(0.0f), dgFloat32(1.0f));
	m_dampCoef.m_w = DG_MIN_SPEED_ATT + (DG_MAX_SPEED_ATT - DG_MIN_SPEED_ATT) * linearDamp;
}

DG_INLINE void dgDynamicBody::SetAngularDamping (const dgVector& angularDamp)
{
	dgFloat32 tmp = dgClamp (angularDamp.m_x, dgFloat32(0.0f), dgFloat32(1.0f));
	m_dampCoef.m_x = DG_MIN_SPEED_ATT + (DG_MAX_SPEED_ATT - DG_MIN_SPEED_ATT) * tmp;

	tmp = dgClamp (angularDamp.m_y, dgFloat32(0.0f), dgFloat32(1.0f));
	m_dampCoef.m_y = DG_MIN_SPEED_ATT + (DG_MAX_SPEED_ATT - DG_MIN_SPEED_ATT) * tmp;

	tmp = dgClamp (angularDamp.m_z, dgFloat32(0.0f), dgFloat32(1.0f));
	m_dampCoef.m_z = DG_MIN_SPEED_ATT + (DG_MAX_SPEED_ATT - DG_MIN_SPEED_ATT) * tmp;
}

DG_INLINE dgVector dgDynamicBody::GetApparentMass() const
{
	return m_aparentMass;
}



DG_INLINE dgVector dgDynamicBody::GetInvMass() const
{
	return m_invMass;
}

DG_INLINE dgVector dgDynamicBody::GetMass() const
{
	return m_mass;
}


DG_INLINE void dgDynamicBody::AddForce (const dgVector& force)
{
	SetForce (m_accel + force);
}

DG_INLINE void dgDynamicBody::AddTorque (const dgVector& torque)
{
	SetTorque (torque + m_alpha);
}


DG_INLINE void dgDynamicBody::SetForce (const dgVector& force)
{
	m_accel = force;
	dgVector error (m_accel - m_prevExternalForce);
	dgFloat32 errMag2 = (error % error) * m_invMass[3] * m_invMass[3];
	if (errMag2 > DG_ErrTolerance2) {
		m_sleepingCounter = 0;
		m_equilibrium = false;
	}
}

DG_INLINE void dgDynamicBody::SetTorque (const dgVector& torque)
{
	m_alpha = torque;
	dgVector error (m_alpha - m_prevExternalTorque);
	dgFloat32 errMag2 = (error % error) * m_invMass[3] * m_invMass[3];
	if (errMag2 > DG_ErrTolerance2) {
		m_sleepingCounter = 0;
		m_equilibrium = false;
	}
}

DG_INLINE void dgDynamicBody::AddDampingAcceleration()
{
	m_veloc -= m_veloc.Scale (m_dampCoef.m_w);
	dgVector omega (m_matrix.UnrotateVector (m_omega));
	omega -= omega.CompProduct (m_dampCoef);
	m_omega = m_matrix.RotateVector (omega);
}

DG_INLINE void dgDynamicBody::AddDampingAccelerationSimd()
{
	(dgSimd&)m_veloc = (dgSimd&)m_veloc - ((dgSimd&)m_veloc) * dgSimd(m_dampCoef.m_w);
	dgSimd omega (m_matrix.UnrotateVectorSimd(m_omega));
	omega = omega - omega * (dgSimd&)m_dampCoef;
	(dgSimd&)m_omega = m_matrix.RotateVectorSimd(omega);
}



DG_INLINE dgVector dgDynamicBody::PredictLinearVelocity(dgFloat32 timestep) const
{
	return 	m_veloc + m_accel.Scale (timestep * m_invMass.m_w);
}

DG_INLINE dgVector dgDynamicBody::PredictAngularVelocity(dgFloat32 timestep) const
{
	return m_omega + m_invWorldInertiaMatrix.RotateVector(m_alpha).Scale (timestep);
}


DG_INLINE OnApplyExtForceAndTorque dgDynamicBody::GetExtForceAndTorqueCallback () const
{
	return m_applyExtForces;
}

DG_INLINE void dgDynamicBody::SetExtForceAndTorqueCallback (OnApplyExtForceAndTorque callback)
{
	m_applyExtForces = callback;
}

#endif 


