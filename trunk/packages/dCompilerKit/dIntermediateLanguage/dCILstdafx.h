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

#ifndef __dCILstdafx_H_
#define __dCILstdafx_H_


#include <dCRC.h>
#include <dTree.h>
#include <dList.h>
#include <dHeap.h>
#include <dString.h>
#include <dVirtualMachine.h>
#include <dContainersStdAfx.h>


#ifdef _MSC_VER
	#ifdef max
		#undef max
	#endif

	#ifdef min
		#undef min
	#endif

#pragma warning (disable: 4100)		// 'O' : unreferenced formal parameter
#pragma warning (disable: 4244)		//'argument' : conversion from 'unsigned int' to 'unsigned short', possible loss of data
#pragma warning (disable: 4480)		// nonstandard extension used: specifying underlying type for enum ''
#pragma warning (disable: 4355)		//'this' : used in base member initializer list
#pragma warning (disable: 4800)		//'unsigned int' : forcing value to bool 'true' or 'false' (performance warning)
#pragma warning (disable: 4512)		//'llvm::Type' : assignment operator could not be generated

//4275
//4146
//4180
//4267
//4345
//4351
//4503
//4624
//4291
#endif

#include <llvm\ADT\Triple.h>

#include <llvm/IR/Verifier.h>
#include <llvm/ExecutionEngine/GenericValue.h>
#include <llvm/ExecutionEngine/Interpreter.h>
#include <llvm/ExecutionEngine/JIT.h>
#include <llvm/IR/Constants.h>
#include <llvm/IR/DerivedTypes.h>
#include <llvm/IR/Instructions.h>
#include <llvm/IR/LLVMContext.h>
#include <llvm/IR/Module.h>
#include <llvm/Target/TargetMachine.h>
#include <llvm\Target\TargetOptions.h>

//#include <llvm/Support/TargetSelect.h>
#include <llvm/Support/CodeGen.h>
#include <llvm/Support/raw_ostream.h>
#include <llvm\Support/TargetRegistry.h>




#endif