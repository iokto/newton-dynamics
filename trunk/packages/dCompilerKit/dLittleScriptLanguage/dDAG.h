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

#ifndef __dDAG_h_
#define __dDAG_h_

#include "dLSCstdafx.h"

#ifdef _MSC_VER
#pragma warning (disable: 4100) // warning C4100: unreferenced formal parameter
#endif




#define D_SCOPE_PREFIX		"scope" 


class dDAGClassNode;
class dDAGFunctionNode;
class dDAGScopeBlockNode;

class dDAG
{
	public:
	dDAG(dList<dDAG*>& allNodes);
	virtual ~dDAG(void);
	
	virtual void CompileCIL(dCIL& cil)  {_ASSERTE (0);}
	virtual void ConnectParent(dDAG* const parent) {_ASSERTE (0);}
	virtual dDAG* Clone (dList<dDAG*>& allNodes) const {_ASSERTE (0); return NULL;}

	dDAGClassNode* GetClass() const;
	dDAGScopeBlockNode* GetScope() const;
	dDAGFunctionNode* GetFunction() const;
	bool RenameLocalVariable(dCIL& cil, dString& variable) const;
	
	dString m_name;
	dTreeAdressStmt::dArg m_result;
	dDAG* m_next;
	dDAG* m_parent;
	dList<dDAG*>::dListNode* m_myListNode;

	dRttiRootClassSupportDeclare(dDAG);
	
};


#endif