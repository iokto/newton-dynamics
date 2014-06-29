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

#include "dLSCstdafx.h"
#include "dDAG.h"
#include "dDAGTypeNode.h"
#include "dDAGDimensionNode.h"
#include "dDAGExpressionNodeAssigment.h"

dInitRtti(dDAGExpressionNodeAssigment);

dDAGExpressionNodeAssigment::dDAGExpressionNodeAssigment(dList<dDAG*>& allNodes, dDAGExpressionNodeVariable* const leftVariable, dDAGExpressionNode* const expression)
	:dDAGExpressionNode(allNodes)
	,m_expression(expression)
	,m_leftVariable(leftVariable)
{
}


dDAGExpressionNodeAssigment::~dDAGExpressionNodeAssigment()
{
}

void dDAGExpressionNodeAssigment::ConnectParent(dDAG* const parent)  
{
	m_parent = parent;
	m_expression->ConnectParent(this);
	m_leftVariable->ConnectParent(this);
}


dCIL::dReturnValue dDAGExpressionNodeAssigment::Evalue(dCIL& cil)
{
	dCIL::dReturnValue val (m_expression->Evalue(cil));
	dAssert (!m_leftVariable->m_dimExpressions.GetCount());

	if (m_leftVariable->m_type->m_intrinsicType != val.m_type) {
		// see if is ok to do cast promotions
		dAssert (0);
	}
	return val;
}

void dDAGExpressionNodeAssigment::CompileCIL(dCIL& cil)  
{
	m_expression->CompileCIL(cil); 
	if (m_leftVariable->m_dimExpressions.GetCount()) {
		dAssert(0);
/*
		dString& variable (m_leftVariable->m_name);
		int pos = variable.Find (D_SCOPE_PREFIX, 0, int (strlen (D_SCOPE_PREFIX)));
		if (pos != 0) {
			bool state = RenameLocalVariable(cil, variable);
			if (!state) {
				dTrace (("undefined local variable\n"));
				dAssert (0);
			}
		}

		dDAGDimensionNode* const dim = m_leftVariable->m_dimExpressions.GetFirst()->GetInfo();
		dim->CompileCIL(cil);
		dCIL::dListNode* const dimInstruction = cil.NewStatement();
		dTreeAdressStmt& addressIndex = dimInstruction->GetInfo();
		addressIndex.m_instruction = dTreeAdressStmt::m_assigment;
		addressIndex.m_arg0.m_label = cil.NewTemp();
		addressIndex.m_arg1 = dim->m_result; 

		dString result = addressIndex.m_arg0.m_label;
		DTRACE_INTRUCTION (&addressIndex);

		for (dList<dDAGDimensionNode*>::dListNode* node = m_leftVariable->m_dimExpressions.GetFirst()->GetNext(); node; node = node->GetNext()) {
			dAssert (0);
			dDAGDimensionNode* const dim = node->GetInfo();
			dim->CompileCIL(cil);
			
			dTreeAdressStmt& stmtMul = cil.NewStatement()->GetInfo();
			stmtMul.m_instruction = dTreeAdressStmt::m_assigment;
			stmtMul.m_operator = dTreeAdressStmt::m_mul;
			stmtMul.m_arg0.m_label = cil.NewTemp();
			stmtMul.m_arg1.m_label = result;
			stmtMul.m_arg2.m_label = dim->m_arraySize;

			DTRACE_INTRUCTION (&stmtMul);

			dTreeAdressStmt& stmtAdd = cil.NewStatement()->GetInfo();
			stmtAdd.m_instruction = dTreeAdressStmt::m_assigment;
			stmtAdd.m_operator = dTreeAdressStmt::m_add;
			stmtAdd.m_arg0.m_label = cil.NewTemp();
			stmtAdd.m_arg1.m_label = stmtMul.m_arg0.m_label;
			stmtAdd.m_arg2 = dim->m_result;

			result = stmtAdd.m_arg0.m_label;

			DTRACE_INTRUCTION (&stmtAdd);
		}

		dAssert (m_parent);
		#ifdef D_USE_COMPLEX_ADRESSING_MODE
			// emit an indirect addressing mode
			dTreeAdressStmt& tmp = cil.NewStatement()->GetInfo();
			tmp.m_instruction = dTreeAdressStmt::m_store;
			tmp.m_arg0 = m_expression->m_result;
			tmp.m_arg1.m_label = variable;
			tmp.m_arg2.m_label = result;
			// the size of the integer register in power of two
			tmp.m_extraInformation = 2;
			DTRACE_INTRUCTION (&tmp);
			m_result.m_label = tmp.m_arg0.m_label; 

		#else
			dTreeAdressStmt& dimSize = cil.NewStatement()->GetInfo();
			dimSize.m_instruction = dTreeAdressStmt::m_assigment;
			dimSize.m_operator = dTreeAdressStmt::m_mul;
			dimSize.m_arg0.m_label = cil.NewTemp();
			dimSize.m_arg1.m_label = result; 
			dimSize.m_arg2.m_type = dTreeAdressStmt::m_intConst;
			dimSize.m_arg2.m_label = "4"; 
			DTRACE_INTRUCTION (&dimSize);

			
			// emit an indirect addressing mode
			dTreeAdressStmt& tmp = cil.NewStatement()->GetInfo();
			tmp.m_instruction = dTreeAdressStmt::m_store;
			tmp.m_arg0 = m_expression->m_result;
			tmp.m_arg1.m_label = variable;
			tmp.m_arg2.m_label = dimSize.m_arg0.m_label;
			DTRACE_INTRUCTION (&tmp);
			m_result.m_label = tmp.m_arg0.m_label; 
		#endif
*/
	} else {
		m_leftVariable->CompileCIL(cil); 
		dTreeAdressStmt::dArg arg1 (LoadLocalVariable(cil, m_expression->m_result));
		dTreeAdressStmt& stmt = cil.NewStatement()->GetInfo();
		stmt.m_instruction = (m_leftVariable->m_result.m_label.Find(dScopePrefix) == 0) ? dTreeAdressStmt::m_storeBase : stmt.m_instruction = dTreeAdressStmt::m_assigment;
		stmt.m_arg0 = m_leftVariable->m_result;
		stmt.m_arg1 = arg1;
		DTRACE_INTRUCTION (&stmt);
	}
}


dDAGExpressionNodeVariable* dDAGExpressionNodeAssigment::FindLeftVariable()
{
	return m_leftVariable->FindLeftVariable();
}