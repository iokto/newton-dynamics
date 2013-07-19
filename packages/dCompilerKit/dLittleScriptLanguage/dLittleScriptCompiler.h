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


#ifndef __dScriptCompiler_h__
#define __dScriptCompiler_h__



#include <dCRC.h>
#include <dTree.h>
#include <dList.h>
#include <dRtti.h>
#include <dRefCounter.h>
#include <dContainersStdAfx.h>

#include "dDAG.h"
#include "dLittleScriptParser.h"

class dDAG;
class dDAGClassNode;
class dScriptPackage;
class dSyntaxTreeCode;
class dDAGFunctionNode;
class dDAGScopeBlockNode;
class dLittleScriptLexical;


class dScriptCompiler: public dLittleScriptParser
{
	public:
	dScriptCompiler(const char* const pakacgesRootNameDirectory);
	virtual ~dScriptCompiler();

	int CompileSource (const char* const sourceFileName);

	protected:
	virtual bool Parse(dLittleScriptLexical& scanner);


//	void DisplayError (const char* format, ...) const;
//	void SyntaxError (const dLittleScriptLexical& scanner, const dUserVariable& errorToken, const dUserVariable& errorTokenMarker);
//	dUserVariable NewConstructorPrototype (const dUserVariable& functionName, const dUserVariable& parameterList, const dUserVariable& initBaseClass);
//	dUserVariable NewFunctionPrototype (const dUserVariable& returnType, const dUserVariable& functionName, const dUserVariable& parameterList, const dUserVariable& isConst);
//	dUserVariable NewFunctionModifier (const dUserVariable& modifirer, const dUserVariable& appendModifirer);
//	dUserVariable TypeAddDimension(const dUserVariable& type, const dUserVariable& dimension);

	// called from parcel
	void ImportClass (const dString& className);
	void ImportAllClasses (const dString& className);
	void OpenPackage (const dString& packageName);


	dDAGClassNode* GetCurrentClass() const;
	dDAGScopeBlockNode* GetCurrentScope() const;

	dUserVariable CreateClass (const dString& visibility, const dString& classType, const dString& className, const dString& superClassName, const dString& interfaces);

	dUserVariable AddClassContructor (const dString& name, const dString& visibility);
	dUserVariable AddClassFunction (const dUserVariable& returnType, const dString& name, const dString& visibility);
	
	dUserVariable EmitTypeNode (const dUserVariable& type, const dUserVariable& dim = dUserVariable());

	dUserVariable NewParameterNode (const dUserVariable& primitiveType, const dString& name);
	dUserVariable FunctionAddParameterNode (const dUserVariable& parameter);
	dUserVariable FunctionAddBodyBlock (const dUserVariable& functionBody);

	dUserVariable BeginScopeBlock ();
	dUserVariable AddStatementToCurrentBlock(const dUserVariable& statement);
	dUserVariable NewVariableToCurrentBlock (const dString& modifiers, const dUserVariable& type, const dString& name);
	dUserVariable EndScopeBlock ();
	
//	dUserVariable AddClassVariable (const dString& modifiers, const dUserVariable& type, const dUserVariable& variableList);
//	dUserVariable ConcatenateVariables(const dUserVariable& variableA, const dUserVariable& variableB);
	dUserVariable ConcatenateExpressions(const dUserVariable& expressionA, const dUserVariable& expressionB);

	dUserVariable NewVariableStatement(const dString& name);
//	dUserVariable SetLocalVariablesType (const dUserVariable& type, const dUserVariable& variableList);
	dUserVariable NewExpressionNodeConstant (const dUserVariable& identifier);

	dUserVariable NewExpresionNodePrefixPostfixOperator (const dUserVariable& expression, bool isPrefix, bool isIncrement);
	dUserVariable NewExpressionNodeBinaryOperator (const dUserVariable& expressionA, const dUserVariable& binaryOperator, const dUserVariable& expressionB);

	dUserVariable NewExpressionNodeLogiOperator (const dUserVariable& expressionA, const dUserVariable& logigOperator, const dUserVariable& expressionB);

	dUserVariable NewDimensionNode(const dUserVariable& expression);
	dUserVariable ConcatenateDimensionNode(const dUserVariable& dim0, const dUserVariable& dim1);
	dUserVariable NewExpressionOperatorNew (const dString& typeName, const dUserVariable& dimension);

	dUserVariable NewExpressionNodeVariable (const dString& name, const dUserVariable& dimArray = dUserVariable());
	dUserVariable NewExpresionNodeAssigment (const dUserVariable& leftVariable, const dUserVariable& expression);
	dUserVariable NewExpresionNodeAssigment (const dUserVariable& leftVariable, const dUserVariable& assigmentOperator, const dUserVariable& expression);

	dUserVariable NewReturnStatement(const dUserVariable& expression);
	dUserVariable NewIFStatement(const dUserVariable& expression, const dUserVariable& thenExpression, const dUserVariable& elseExpression);

	
	dUserVariable NewExpressionFunctionCall (const dString& name, const dUserVariable& argumnetList);
	dUserVariable NewCaseStatement(const dString& constID, const dUserVariable& statement);

	dUserVariable ConcatenateCaseBlocks (const dUserVariable& firstStatement, const dUserVariable& lastStatement);
	dUserVariable ConcatenateCaseStatement (const dUserVariable& firstStatement, const dUserVariable& lastStatement);
	dUserVariable ConcatenateParametersExpressions(const dUserVariable& ExpressionA, const dUserVariable& ExpressionB);

	dUserVariable NewEmptyStatement();
	dUserVariable NewBreakStatement();
	dUserVariable NewContinueStatement();
	dUserVariable NewDoStatement(const dUserVariable& expression, const dUserVariable& statement);
	dUserVariable NewWhileStatement(const dUserVariable& expression, const dUserVariable& statement);
	dUserVariable NewSwitchStatement(const dUserVariable& expression, const dUserVariable& caseStatementList);
	dUserVariable NewForStatement(const dUserVariable& init_exp, const dUserVariable& conditional, const dUserVariable& step_Exp, const dUserVariable& statement);

	

	dString m_packageFileName;
	dString m_packageRootDirectory;
	dScriptPackage* m_currentPackage;
	dDAGFunctionNode* m_currentFunction;
	dList<dDAGClassNode*> m_classList;
	dList<dDAGScopeBlockNode*> m_scopeStack;
	dList<dDAG*> m_allNodes;


	friend class dLittleScriptParser;
};

#endif