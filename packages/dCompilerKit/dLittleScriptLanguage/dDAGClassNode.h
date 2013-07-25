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

#ifndef __dDAGClassNode_H__
#define __dDAGClassNode_H__

#include "dDAG.h"
#include "dLSCstdafx.h"

class dDAGTypeNode;
class dDAGFunctionNode;
class dDAGParameterNode;
class dDAGExpressionNode;

class dDAGClassNode: public dDAG
{
	public:
	dDAGClassNode(dList<dDAG*>& allNodes);
	~dDAGClassNode(void);

	void FinalizeImplementation (const char* const visibility, const char* const name, dDAGClassNode* const baseClass);

	void AddVariable (dDAGParameterNode* const variable);
	void AddFunction (dDAGFunctionNode* const functionNode);

	dString GetFunctionName (const char* const functionName, dList<dDAGParameterNode*>& parameterNodeList) const;
	dString GetFunctionName (const char* const functionName, dList<dDAGExpressionNode*>& argumentList) const;
	dDAGTypeNode* GetFunctionReturnType(const char* const functionName, dList<dDAGExpressionNode*>& argumentList) const;

	dDAGFunctionNode* GetCurrentFunction ();
	dDAGParameterNode* FindVariable(const char* name) const;

	virtual void CompileCIL(dCIL& cil);
	virtual void ConnectParent(dDAG* const parent);

	bool m_isPublic;
	const dDAGClassNode* m_baseClass;
	dList<dDAGParameterNode*> m_variables;
	dList<dDAGFunctionNode*> m_functionList;
	dList<dDAGFunctionNode*> m_constructors;
	dList<dCIL::dListNode*> m_cilCodeList;
	dAddRtti(dDAG);
};


#endif