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


#ifndef _D_NEWTON_VEHICLE_MANAGER_H_
#define _D_NEWTON_VEHICLE_MANAGER_H_

#include "dStdAfxNewton.h"
#include "dNewtonDynamicBody.h"

class dNewtonVehicleManager: public CustomVehicleControllerManager
{
	public:
	class dNewtonVehicle: public dNewtonDynamicBody
	{
		public:
		CNEWTON_API dNewtonVehicle (dNewtonVehicleManager* const manager, void* const userData, dFloat mass, dFloat outerRadius, dFloat innerRadius, dFloat height, dFloat stairStep, const dFloat* const upDir, const dFloat* const frontDir, dLong collisionMask);
		CNEWTON_API ~dNewtonVehicle ();

//		virtual void OnVehicleMove (dFloat timestep) = 0;
//		CNEWTON_API dFloat GetVehicleHigh() const;
//		CNEWTON_API void SetVehicleVelocity (dFloat forwardSpeed, dFloat lateralSpeed, dFloat verticalSpeed, dFloat headingAngle, const dFloat* const gravity, dFloat timestep);

		private:
		CustomVehicleController* m_controller;
		friend class dNewtonVehicleManager;
	};

	CNEWTON_API dNewtonVehicleManager (dNewton* const world);
	CNEWTON_API virtual ~dNewtonVehicleManager ();

	CNEWTON_API dNewtonVehicle* GetFirstVehicle() const;
	CNEWTON_API dNewtonVehicle* GetNextVehicle(const dNewtonVehicle* const Vehicle) const;
//	CNEWTON_API virtual void ApplyVehicleMove (CustomVehicleController* const controller, dFloat timestep);
};




#endif
