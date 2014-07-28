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

#include "dStdAfxMath.h"
#include "dMathDefines.h"
#include "dLinearAlgebra.h"


#define COMPLEMENTARITY_VEL_DAMP			dFloat(100.0f)
#define COMPLEMENTARITY_POS_DAMP			dFloat(1500.0f)
#define COMPLEMENTARITY_PSD_DAMP_TOL		dFloat(1.0e-4f)
#define COMPLEMENTARITY_STACK_ENTRIES		64
#define COMPLEMENTARITY_MAX_FRICTION_BOUND	dFloat(1.0e15f)
#define COMPLEMENTARITY_MIN_FRICTION_BOUND	-COMPLEMENTARITY_MAX_FRICTION_BOUND


#ifdef _MSC_VER
#pragma warning (disable: 4100) //unreferenced formal parameter
#endif


void dSymmetricBiconjugateGradientSolve::ScaleAdd (int size, dFloat64* const a, const dFloat64* const b, dFloat64 scale, const dFloat64* const c) const
{
	for (int i = 0; i < size; i ++) {
		a[i] = b[i] + scale * c[i];
	}
}

void dSymmetricBiconjugateGradientSolve::Sub (int size, dFloat64* const a, const dFloat64* const b, const dFloat64* const c) const
{
	for (int i = 0; i < size; i ++) {
		a[i] = b[i] - c[i];
	}
}

dFloat64 dSymmetricBiconjugateGradientSolve::DotProduct (int size, const dFloat64* const b, const dFloat64* const c) const
{
	dFloat64 product = dFloat64 (0.0f);
	for (int i = 0; i < size; i ++) {
		product += b[i] * c[i];
	}
	return product;
}

dFloat64 dSymmetricBiconjugateGradientSolve::Solve (int size, dFloat64 tolerance, dFloat64* const x, const dFloat64* const b) const
{
//	dgStack<dFloat64> bufferR0(size);
//	dgStack<dFloat64> bufferP0(size);
//	dgStack<dFloat64> matrixTimesP0(size);
//	dgStack<dFloat64> bufferConditionerInverseTimesR0(size);
//	dFloat64* const r0 = &bufferR0[0];
//	dFloat64* const p0 = &bufferP0[0];
//	dFloat64* const MinvR0 = &bufferConditionerInverseTimesR0[0];
//	dFloat64* const matrixP0 = &matrixTimesP0[0];

	dFloat64* const r0 = new dFloat64 [size];
	dFloat64* const p0 = new dFloat64 [size];
	dFloat64* const MinvR0 = new dFloat64 [size];
	dFloat64* const matrixP0 = new dFloat64 [size];

	MatrixTimeVector (matrixP0, x);
	Sub(size, r0, b, matrixP0);
	bool continueExecution = InversePrecoditionerTimeVector (p0, r0);

	int iter = 0;
	dFloat64 num = DotProduct (size, r0, p0);
	dFloat64 error2 = num;
	for (int j = 0; (j < size) && (error2 > tolerance) && continueExecution; j ++) {

		MatrixTimeVector (matrixP0, p0);
		dFloat64 den = DotProduct (size, p0, matrixP0);

		dAssert (fabs(den) > dFloat64 (0.0f));
		dFloat64 alpha = num / den;

		ScaleAdd (size, x, x, alpha, p0);
        if ((j % 50) != 49) {
		    ScaleAdd (size, r0, r0, -alpha, matrixP0);
        } else {
            MatrixTimeVector (matrixP0, x);
            Sub(size, r0, b, matrixP0);
        }

		continueExecution = InversePrecoditionerTimeVector (MinvR0, r0);

		dFloat64 num1 = DotProduct (size, r0, MinvR0);
		dFloat64 beta = num1 / num;
		ScaleAdd (size, p0, MinvR0, beta, p0);
		num = DotProduct (size, r0, MinvR0);
		iter ++;
		error2 = num;
		if (j > 10) {
			error2 = dFloat64 (0.0f);
			for (int i = 0; i < size; i ++) {
				error2 = dMax (error2, r0[i] * r0[i]);
			}
		}
	}

	delete[] matrixP0;
	delete[] MinvR0;
	delete[] p0;
	delete[] r0;

	dAssert (iter < size);
	return num;
}


dComplemtaritySolver::dBodyState::dBodyState()
	:m_matrix(GetIdentityMatrix())
	,m_localFrame(GetZeroMatrix())
	,m_inertia(GetZeroMatrix())
	,m_invInertia(GetZeroMatrix())
	,m_localInertia (0.0f, 0.0f, 0.0f, 0.0f)
	,m_localInvInertia(0.0f, 0.0f, 0.0f, 0.0f)
	,m_veloc(0.0f, 0.0f, 0.0f, 0.0f)
	,m_omega(0.0f, 0.0f, 0.0f, 0.0f)
	,m_externalForce(0.0f, 0.0f, 0.0f, 0.0f)
	,m_externalTorque(0.0f, 0.0f, 0.0f, 0.0f)
	,m_globalCentreOfMass(0.0f, 0.0f, 0.0f, 0.0f)
	,m_mass(0.0f)
	,m_invMass(0.0f)
	,m_myIndex(0)
{
}


const dVector& dComplemtaritySolver::dBodyState::GetOmega() const
{
	return m_omega;
}

const dVector& dComplemtaritySolver::dBodyState::GetVelocity() const
{
	return m_veloc;
}


dFloat dComplemtaritySolver::dBodyState::GetMass () const
{
	return m_mass;
}

const dMatrix& dComplemtaritySolver::dBodyState::GetMatrix () const
{
	return m_matrix;
}

const dMatrix& dComplemtaritySolver::dBodyState::GetLocalMatrix () const
{
	return m_localFrame;
}

const dVector& dComplemtaritySolver::dBodyState::GetCenterOfMass () const
{
	return m_globalCentreOfMass;
}



void dComplemtaritySolver::dBodyState::UpdateInertia()
{
	dMatrix tmpMatrix (GetZeroMatrix());

	tmpMatrix[0] = m_localInertia.CompProduct (dVector (m_matrix[0][0], m_matrix[1][0], m_matrix[2][0], 0.0f));
	tmpMatrix[1] = m_localInertia.CompProduct (dVector (m_matrix[0][1], m_matrix[1][1], m_matrix[2][1], 0.0f));
	tmpMatrix[2] = m_localInertia.CompProduct (dVector (m_matrix[0][2], m_matrix[1][2], m_matrix[2][2], 0.0f));
	m_inertia = tmpMatrix * m_matrix;

	tmpMatrix[0] = m_localInvInertia.CompProduct (dVector (m_matrix[0][0], m_matrix[1][0], m_matrix[2][0], 0.0f));
	tmpMatrix[1] = m_localInvInertia.CompProduct (dVector (m_matrix[0][1], m_matrix[1][1], m_matrix[2][1], 0.0f));
	tmpMatrix[2] = m_localInvInertia.CompProduct (dVector (m_matrix[0][2], m_matrix[1][2], m_matrix[2][2], 0.0f));
	m_invInertia = tmpMatrix * m_matrix;
}

void dComplemtaritySolver::dBodyState::IntegrateForce (dFloat timestep, const dVector& force, const dVector& torque)
{
	dVector accel (force.Scale (m_invMass));
	dVector alpha (m_invInertia.RotateVector(torque));
	m_veloc += accel.Scale (timestep);
	m_omega += alpha.Scale (timestep);
}

void dComplemtaritySolver::dBodyState::ApplyNetForceAndTorque (dFloat invTimestep, const dVector& veloc, const dVector& omega)
{
	dVector accel = (m_veloc - veloc).Scale(invTimestep);
	dVector alpha = (m_omega - omega).Scale(invTimestep);

	m_externalForce = accel.Scale(m_mass);
	alpha = m_matrix.UnrotateVector(alpha);
	m_externalTorque = m_matrix.RotateVector(alpha.CompProduct(m_localInertia));
}



void dComplemtaritySolver::dBilateralJoint::Init(dBodyState* const state0, dBodyState* const state1)
{
	m_start = 0;
	m_count = 0;
	memset (m_rowIsMotor, 0, sizeof (m_rowIsMotor));
	memset (m_motorAcceleration, 0, sizeof (m_motorAcceleration));
	memset (m_jointFeebackForce, 0, sizeof (m_jointFeebackForce));

	m_state0 = state0;
	m_state1 = state1;
}


void dComplemtaritySolver::dBilateralJoint::InitPointParam (dPointDerivativeParam& param, const dVector& pivot) const
{
	dAssert (m_state0);
	dAssert (m_state1);

	param.m_posit0 = pivot;
	param.m_r0 = pivot - m_state0->m_globalCentreOfMass;
	param.m_veloc0 = m_state0->m_omega * param.m_r0;
	param.m_centripetal0 = m_state0->m_omega * param.m_veloc0;
	param.m_veloc0 += m_state0->m_veloc;

	param.m_posit1 = pivot;
	param.m_r1 = pivot - m_state1->m_globalCentreOfMass;
	param.m_veloc1 = m_state1->m_omega * param.m_r1;
	param.m_centripetal1 = m_state1->m_omega * param.m_veloc1;
	param.m_veloc1 += m_state1->m_veloc;
}


void dComplemtaritySolver::dBilateralJoint::CalculatePointDerivative (dParamInfo* const constraintParams, const dVector& dir, const dPointDerivativeParam& param)
{
	int index = constraintParams->m_count;

	dJacobian &jacobian0 = constraintParams->m_jacobians[index].m_jacobian_IM0; 
	dVector r0CrossDir (param.m_r0 * dir);
	jacobian0.m_linear[0] = dir.m_x;
	jacobian0.m_linear[1] = dir.m_y;
	jacobian0.m_linear[2] = dir.m_z;
	jacobian0.m_linear[3] = dFloat (0.0f);
	jacobian0.m_angular[0] = r0CrossDir.m_x;
	jacobian0.m_angular[1] = r0CrossDir.m_y;
	jacobian0.m_angular[2] = r0CrossDir.m_z;
	jacobian0.m_angular[3] = 0.0f;

	dJacobian &jacobian1 = constraintParams->m_jacobians[index].m_jacobian_IM1; 
	dVector r1CrossDir (dir * param.m_r1);
	jacobian1.m_linear[0] = -dir.m_x;
	jacobian1.m_linear[1] = -dir.m_y;
	jacobian1.m_linear[2] = -dir.m_z;
	jacobian1.m_linear[3] = dFloat (0.0f);
	jacobian1.m_angular[0] = r1CrossDir.m_x;
	jacobian1.m_angular[1] = r1CrossDir.m_y;
	jacobian1.m_angular[2] = r1CrossDir.m_z;
	jacobian1.m_angular[3] = 0.0f;

	dVector velocError (param.m_veloc1 - param.m_veloc0);
	dVector positError (param.m_posit1 - param.m_posit0);
	dVector centrError (param.m_centripetal1 - param.m_centripetal0);

	dFloat relPosit = positError % dir;
	dFloat relVeloc = velocError % dir;
	dFloat relCentr = centrError % dir; 

	dFloat dt = constraintParams->m_timestep;
	dFloat ks = COMPLEMENTARITY_POS_DAMP;
	dFloat kd = COMPLEMENTARITY_VEL_DAMP;
	dFloat ksd = dt * ks;
	dFloat num = ks * relPosit + kd * relVeloc + ksd * relVeloc;
	dFloat den = dFloat (1.0f) + dt * kd + dt * ksd;
	dFloat accelError = num / den;

	m_rowIsMotor[index] = false;
	m_motorAcceleration[index] = 0.0f;
	constraintParams->m_jointAccel[index] = accelError + relCentr;
	constraintParams->m_jointLowFriction[index] = COMPLEMENTARITY_MIN_FRICTION_BOUND;
	constraintParams->m_jointHighFriction[index] = COMPLEMENTARITY_MAX_FRICTION_BOUND;
	constraintParams->m_count = index + 1;
}


void dComplemtaritySolver::dBilateralJoint::AddAngularRowJacobian (dParamInfo* const constraintParams, const dVector& dir, dFloat jointAngle)
{
	int index = constraintParams->m_count;
	dJacobian &jacobian0 = constraintParams->m_jacobians[index].m_jacobian_IM0; 

	jacobian0.m_linear[0] = 0.0f;
	jacobian0.m_linear[1] = 0.0f;
	jacobian0.m_linear[2] = 0.0f;
	jacobian0.m_linear[3] = 0.0f;
	jacobian0.m_angular[0] = dir.m_x;
	jacobian0.m_angular[1] = dir.m_y;
	jacobian0.m_angular[2] = dir.m_z;
	jacobian0.m_angular[3] = 0.0f;

	dJacobian &jacobian1 = constraintParams->m_jacobians[index].m_jacobian_IM1; 
	jacobian1.m_linear[0] = 0.0f;
	jacobian1.m_linear[1] = 0.0f;
	jacobian1.m_linear[2] = 0.0f;
	jacobian1.m_linear[3] = 0.0f;
	jacobian1.m_angular[0] = -dir.m_x;
	jacobian1.m_angular[1] = -dir.m_y;
	jacobian1.m_angular[2] = -dir.m_z;
	jacobian1.m_angular[3] = 0.0f;

	const dVector& omega0 = m_state0->m_omega;
	const dVector& omega1 = m_state1->m_omega;
	dFloat omegaError = (omega1 - omega0) % dir;


	//at =  [- ks (x2 - x1) - kd * (v2 - v1) - dt * ks * (v2 - v1)] / [1 + dt * kd + dt * dt * ks] 
	dFloat dt = constraintParams->m_timestep;
	dFloat ks = COMPLEMENTARITY_POS_DAMP;
	dFloat kd = COMPLEMENTARITY_VEL_DAMP;
	dFloat ksd = dt * ks;
	dFloat num = ks * jointAngle + kd * omegaError + ksd * omegaError;
	dFloat den = dFloat (1.0f) + dt * kd + dt * ksd;
	dFloat alphaError = num / den;

	m_rowIsMotor[index] = false;
	m_motorAcceleration[index] = 0.0f;
	constraintParams->m_jointAccel[index] = alphaError;
	constraintParams->m_jointLowFriction[index] = COMPLEMENTARITY_MIN_FRICTION_BOUND;
	constraintParams->m_jointHighFriction[index] = COMPLEMENTARITY_MAX_FRICTION_BOUND;
	constraintParams->m_count = index + 1;
}


void dComplemtaritySolver::dBilateralJoint::AddAngularRowJacobian (dParamInfo* const constraintParams, const dVector& dir0, const dVector& dir1, dFloat accelerationRatio)
{
	int index = constraintParams->m_count;
	dJacobian &jacobian0 = constraintParams->m_jacobians[index].m_jacobian_IM0; 

	jacobian0.m_linear[0] = 0.0f;
	jacobian0.m_linear[1] = 0.0f;
	jacobian0.m_linear[2] = 0.0f;
	jacobian0.m_linear[3] = 0.0f;
	jacobian0.m_angular[0] = dir0.m_x;
	jacobian0.m_angular[1] = dir0.m_y;
	jacobian0.m_angular[2] = dir0.m_z;
	jacobian0.m_angular[3] = 0.0f;

	dJacobian &jacobian1 = constraintParams->m_jacobians[index].m_jacobian_IM1; 
	jacobian1.m_linear[0] = 0.0f;
	jacobian1.m_linear[1] = 0.0f;
	jacobian1.m_linear[2] = 0.0f;
	jacobian1.m_linear[3] = 0.0f;
	jacobian1.m_angular[0] = dir1.m_x;
	jacobian1.m_angular[1] = dir1.m_y;
	jacobian1.m_angular[2] = dir1.m_z;
	jacobian1.m_angular[3] = 0.0f;

	m_rowIsMotor[index] = true;
	m_motorAcceleration[index] = accelerationRatio;
	constraintParams->m_jointAccel[index] = 0.0f;
	constraintParams->m_jointLowFriction[index] = COMPLEMENTARITY_MIN_FRICTION_BOUND;
	constraintParams->m_jointHighFriction[index] = COMPLEMENTARITY_MAX_FRICTION_BOUND;
	constraintParams->m_count = index + 1;
}

void dComplemtaritySolver::dBilateralJoint::AddLinearRowJacobian (dParamInfo* const constraintParams, const dVector& pivot, const dVector& dir)
{
	dPointDerivativeParam pointData;
	InitPointParam (pointData, pivot);
	CalculatePointDerivative (constraintParams, dir, pointData); 
}


void dComplemtaritySolver::dBilateralJoint::JointAccelerations (dJointAccelerationDecriptor* const params)
{
	dJacobianColum* const jacobianColElements = params->m_colMatrix;
	dJacobianPair* const jacobianRowElements = params->m_rowMatrix;

	const dVector& bodyVeloc0 = m_state0->m_veloc;
	const dVector& bodyOmega0 = m_state0->m_omega;
	const dVector& bodyVeloc1 = m_state1->m_veloc;
	const dVector& bodyOmega1 = m_state1->m_omega;

	dFloat timestep = params->m_timeStep;
	dFloat kd = COMPLEMENTARITY_VEL_DAMP * dFloat (4.0f);
	dFloat ks = COMPLEMENTARITY_POS_DAMP * dFloat (0.25f);
	for (int k = 0; k < params->m_rowsCount; k ++) {
		if (m_rowIsMotor[k]) {
			jacobianColElements[k].m_coordenateAccel = m_motorAcceleration[k] + jacobianColElements[k].m_deltaAccel;
		} else {
			const dJacobianPair& Jt = jacobianRowElements[k];
			dVector relVeloc (Jt.m_jacobian_IM0.m_linear.CompProduct(bodyVeloc0) +
				Jt.m_jacobian_IM0.m_angular.CompProduct(bodyOmega0) + 
				Jt.m_jacobian_IM1.m_linear.CompProduct(bodyVeloc1) +
				Jt.m_jacobian_IM1.m_angular.CompProduct(bodyOmega1));

			dFloat vRel = relVeloc.m_x + relVeloc.m_y + relVeloc.m_z;
			dFloat aRel = jacobianColElements[k].m_deltaAccel;
			dFloat ksd = timestep * ks;
			dFloat relPosit = 0.0f - vRel * timestep * params->m_firstPassCoefFlag;

			dFloat num = ks * relPosit - kd * vRel - ksd * vRel;
			dFloat den = dFloat (1.0f) + timestep * kd + timestep * ksd;
			dFloat aRelErr = num / den;
			jacobianColElements[k].m_coordenateAccel = aRelErr + aRel;
		}
	}
}

int dComplemtaritySolver::GetActiveJoints (dBilateralJoint** const jointArray, int bufferSize) 
{
	return 0;
}


int dComplemtaritySolver::BuildJacobianMatrix (int jointCount, dBilateralJoint** const jointArray, dFloat timestep, dJacobianPair* const jacobianArray, dJacobianColum* const jacobianColumnArray, int maxRowCount)
{
	int rowCount = 0;

	dParamInfo constraintParams;
	constraintParams.m_timestep = timestep;
	constraintParams.m_timestepInv = 1.0f / timestep;

	// calculate Jacobian derivative for each active joint	
	for (int j = 0; j < jointCount; j ++) {
		dBilateralJoint* const joint = jointArray[j];
		constraintParams.m_count = 0;
		joint->JacobianDerivative (&constraintParams); 

		int dofCount = constraintParams.m_count;
		joint->m_count = dofCount;
		joint->m_start = rowCount;

		// complete the derivative matrix for this joint
		int index = joint->m_start;
		dBodyState* const state0 = joint->m_state0;
		dBodyState* const state1 = joint->m_state1;

		const dMatrix& invInertia0 = state0->m_invInertia;
		const dMatrix& invInertia1 = state1->m_invInertia;

		dFloat invMass0 = state0->m_invMass;
		dFloat invMass1 = state1->m_invMass;
		dFloat weight = 0.9f;

		for (int i = 0; i < dofCount; i ++) {
			dJacobianPair* const row = &jacobianArray[index];
			dJacobianColum* const col = &jacobianColumnArray[index];
			jacobianArray[rowCount] = constraintParams.m_jacobians[i]; 

			dVector JMinvIM0linear (row->m_jacobian_IM0.m_linear.Scale (invMass0));
			dVector JMinvIM1linear (row->m_jacobian_IM1.m_linear.Scale (invMass1));
			dVector JMinvIM0angular = invInertia0.UnrotateVector(row->m_jacobian_IM0.m_angular);
			dVector JMinvIM1angular = invInertia1.UnrotateVector(row->m_jacobian_IM1.m_angular);

			dVector tmpDiag (JMinvIM0linear.CompProduct(row->m_jacobian_IM0.m_linear) + JMinvIM0angular.CompProduct(row->m_jacobian_IM0.m_angular) + JMinvIM1linear.CompProduct(row->m_jacobian_IM1.m_linear) + JMinvIM1angular.CompProduct(row->m_jacobian_IM1.m_angular));
			dVector tmpAccel (JMinvIM0linear.CompProduct (state0->m_externalForce) + JMinvIM0angular.CompProduct(state0->m_externalTorque) + JMinvIM1linear.CompProduct (state1->m_externalForce) + JMinvIM1angular.CompProduct(state1->m_externalTorque));
			dFloat extenalAcceleration = -(tmpAccel[0] + tmpAccel[1] + tmpAccel[2]);

			col->m_diagDamp = 1.0f;
			col->m_coordenateAccel = constraintParams.m_jointAccel[i];
			col->m_jointLowFriction = constraintParams.m_jointLowFriction[i];
			col->m_jointHighFriction = constraintParams.m_jointHighFriction[i];

			col->m_deltaAccel = extenalAcceleration;
			col->m_coordenateAccel += extenalAcceleration;

			col->m_force = joint->m_jointFeebackForce[i] * weight;

			dFloat stiffness = COMPLEMENTARITY_PSD_DAMP_TOL * col->m_diagDamp;
			dFloat diag = (tmpDiag[0] + tmpDiag[1] + tmpDiag[2]);
			dAssert (diag > dFloat (0.0f));
			col->m_diagDamp = diag * stiffness;

			diag *= (dFloat(1.0f) + stiffness);
			col->m_invDJMinvJt = dFloat(1.0f) / diag;
			index ++;
			rowCount ++;
			dAssert (rowCount < maxRowCount);
		}
	}
	return rowCount;
}

void dComplemtaritySolver::CalculateReactionsForces (int bodyCount, dBodyState** const bodyArray, int jointCount, dBilateralJoint** const jointArray, dFloat timestepSrc, dJacobianPair* const jacobianArray, dJacobianColum* const jacobianColumnArray)
{
	dJacobian stateVeloc[COMPLEMENTARITY_STACK_ENTRIES];
	dJacobian internalForces [COMPLEMENTARITY_STACK_ENTRIES];

	int stateIndex = 0;
	dVector zero(dFloat (0.0f), dFloat (0.0f), dFloat (0.0f), dFloat (0.0f));
	for (int i = 0; i < bodyCount; i ++) {
		dBodyState* const state = bodyArray[i];
		stateVeloc[stateIndex].m_linear = state->m_veloc;
		stateVeloc[stateIndex].m_angular = state->m_omega;

		internalForces[stateIndex].m_linear = zero;
		internalForces[stateIndex].m_angular = zero;

		state->m_myIndex = stateIndex;
		stateIndex ++;
		dAssert (stateIndex < int (sizeof (stateVeloc)/sizeof (stateVeloc[0])));
	}

	for (int i = 0; i < jointCount; i ++) {
		dJacobian y0;
		dJacobian y1;
		y0.m_linear = zero;
		y0.m_angular = zero;
		y1.m_linear = zero;
		y1.m_angular = zero;
		dBilateralJoint* const constraint = jointArray[i];
		int first = constraint->m_start;
		int count = constraint->m_count;
		for (int j = 0; j < count; j ++) { 
			dJacobianPair* const row = &jacobianArray[j + first];
			const dJacobianColum* const col = &jacobianColumnArray[j + first];
			dFloat val = col->m_force; 
			y0.m_linear += row->m_jacobian_IM0.m_linear.Scale(val);
			y0.m_angular += row->m_jacobian_IM0.m_angular.Scale(val);
			y1.m_linear += row->m_jacobian_IM1.m_linear.Scale(val);
			y1.m_angular += row->m_jacobian_IM1.m_angular.Scale(val);
		}
		int m0 = constraint->m_state0->m_myIndex;
		int m1 = constraint->m_state1->m_myIndex;
		internalForces[m0].m_linear += y0.m_linear;
		internalForces[m0].m_angular += y0.m_angular;
		internalForces[m1].m_linear += y1.m_linear;
		internalForces[m1].m_angular += y1.m_angular;
	}


	dFloat invTimestepSrc = dFloat (1.0f) / timestepSrc;
	dFloat invStep = dFloat (0.25f);
	dFloat timestep = timestepSrc * invStep;
	dFloat invTimestep = invTimestepSrc * dFloat (4.0f);

	int maxPasses = 5;
	dFloat firstPassCoef = dFloat (0.0f);
	dFloat maxAccNorm = dFloat (1.0e-2f);

	for (int step = 0; step < 4; step ++) {
		dJointAccelerationDecriptor joindDesc;
		joindDesc.m_timeStep = timestep;
		joindDesc.m_invTimeStep = invTimestep;
		joindDesc.m_firstPassCoefFlag = firstPassCoef;

		for (int curJoint = 0; curJoint < jointCount; curJoint ++) {
			dBilateralJoint* const constraint = jointArray[curJoint];
			joindDesc.m_rowsCount = constraint->m_count;
			joindDesc.m_rowMatrix = &jacobianArray[constraint->m_start];
			joindDesc.m_colMatrix = &jacobianColumnArray[constraint->m_start];
			constraint->JointAccelerations (&joindDesc);
		}
		firstPassCoef = dFloat (1.0f);

		dFloat accNorm = dFloat (1.0e10f);
		for (int passes = 0; (passes < maxPasses) && (accNorm > maxAccNorm); passes ++) {
			accNorm = dFloat (0.0f);
			for (int curJoint = 0; curJoint < jointCount; curJoint ++) {

				dBilateralJoint* const constraint = jointArray[curJoint];
				int index = constraint->m_start;
				int rowsCount = constraint->m_count;
				int m0 = constraint->m_state0->m_myIndex;
				int m1 = constraint->m_state1->m_myIndex;

				dVector linearM0 (internalForces[m0].m_linear);
				dVector angularM0 (internalForces[m0].m_angular);
				dVector linearM1 (internalForces[m1].m_linear);
				dVector angularM1 (internalForces[m1].m_angular);

				dBodyState* const state0 = constraint->m_state0;
				dBodyState* const state1 = constraint->m_state1;
				const dMatrix& invInertia0 = state0->m_invInertia;
				const dMatrix& invInertia1 = state1->m_invInertia;
				dFloat invMass0 = state0->m_invMass;
				dFloat invMass1 = state1->m_invMass;

				for (int k = 0; k < rowsCount; k ++) {
					dJacobianPair* const row = &jacobianArray[index];
					dJacobianColum* const col = &jacobianColumnArray[index];

					dVector JMinvIM0linear (row->m_jacobian_IM0.m_linear.Scale (invMass0));
					dVector JMinvIM1linear (row->m_jacobian_IM1.m_linear.Scale (invMass1));
					dVector JMinvIM0angular = invInertia0.UnrotateVector(row->m_jacobian_IM0.m_angular);
					dVector JMinvIM1angular = invInertia1.UnrotateVector(row->m_jacobian_IM1.m_angular);
					dVector acc (JMinvIM0linear.CompProduct(linearM0) + JMinvIM0angular.CompProduct(angularM0) + JMinvIM1linear.CompProduct(linearM1) + JMinvIM1angular.CompProduct(angularM1));

					dFloat a = col->m_coordenateAccel - acc.m_x - acc.m_y - acc.m_z - col->m_force * col->m_diagDamp;
					dFloat f = col->m_force + col->m_invDJMinvJt * a;

					dFloat lowerFrictionForce = col->m_jointLowFriction;
					dFloat upperFrictionForce = col->m_jointHighFriction;

					if (f > upperFrictionForce) {
						a = dFloat (0.0f);
						f = upperFrictionForce;
					} else if (f < lowerFrictionForce) {
						a = dFloat (0.0f);
						f = lowerFrictionForce;
					}

					accNorm = dMax (accNorm, dAbs (a));
					dFloat prevValue = f - col->m_force;
					col->m_force = f;

					linearM0 += row->m_jacobian_IM0.m_linear.Scale (prevValue);
					angularM0 += row->m_jacobian_IM0.m_angular.Scale (prevValue);
					linearM1 += row->m_jacobian_IM1.m_linear.Scale (prevValue);
					angularM1 += row->m_jacobian_IM1.m_angular.Scale (prevValue);
					index ++;
				}
				internalForces[m0].m_linear = linearM0;
				internalForces[m0].m_angular = angularM0;
				internalForces[m1].m_linear = linearM1;
				internalForces[m1].m_angular = angularM1;
			}
		}

		for (int i = 0; i < bodyCount; i ++) {
			dBodyState* const state = bodyArray[i];
			//int index = state->m_myIndex;
			dAssert (state->m_myIndex == i);
			dVector force (state->m_externalForce + internalForces[i].m_linear);
			dVector torque (state->m_externalTorque + internalForces[i].m_angular);
			state->IntegrateForce(timestep, force, torque);
		}
	}

	for (int i = 0; i < jointCount; i ++) {
		dBilateralJoint* const constraint = jointArray[i];
		int first = constraint->m_start;
		int count = constraint->m_count;
		for (int j = 0; j < count; j ++) { 
			const dJacobianColum* const col = &jacobianColumnArray[j + first];
			dFloat val = col->m_force; 
			constraint->m_jointFeebackForce[j] = val;
		}
	}

	for (int i = 0; i < bodyCount; i ++) {
		dBodyState* const state = bodyArray[i];
		//int index = state->m_myIndex;
		dAssert (state->m_myIndex == i);
		state->ApplyNetForceAndTorque (invTimestepSrc, stateVeloc[i].m_linear, stateVeloc[i].m_angular);
	}

	for (int i = 0; i < jointCount; i ++) {
		dBilateralJoint* const constraint = jointArray[i];
		constraint->UpdateSolverForces (jacobianArray);
	}
}

