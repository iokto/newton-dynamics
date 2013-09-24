/////////////////////////////////////////////////////////////////////////////
// Name:        dNodeInfo.h
// Purpose:     
// Author:      Julio Jerez
// Modified by: 
// Created:     22/05/2010 08:02:08
// RCS-ID:      
// Copyright:   Copyright (c) <2010> <Newton Game Dynamics>
// License:     
// This software is provided 'as-is', without any express or implied
// warranty. In no event will the authors be held liable for any damages
// arising from the use of this software.
// 
// Permission is granted to anyone to use this software for any purpose,
// including commercial applications, and to alter it and redistribute it
// freely
/////////////////////////////////////////////////////////////////////////////

#ifndef _D_NODEINFO_H_
#define _D_NODEINFO_H_

#include "dScene.h"
#include "dVariable.h"
#include <dString.h>


class dNodeInfo;
class dSceneRender;


#define D_DEFINE_CLASS_NODE_ESSENCIALS(className,baseClass,exportType)		\
	dAddRtti(baseClass,exportType);											\
	virtual exportType dNodeInfo* MakeCopy () const							\
	{																		\
		return new className(*this);										\
	}																		\
	virtual exportType dNodeInfo* MetaFunction(dScene* const world) const	\
	{																		\
		return new className(world);										\
	}																		\
	static exportType const char* BaseClassName ()							\
	{																		\
		return #baseClass;													\
	}																		\
	static exportType const className& GetSingleton()						\
	{																		\
		return m_singletonClass;											\
	}																		\
	static className m_singletonClass;


#define D_DEFINE_CLASS_NODE(className,baseClass,exportType)			\
	virtual exportType const char* GetClassName () const			\
	{																\
		return #className;											\
	}																\
	D_DEFINE_CLASS_NODE_ESSENCIALS(className,baseClass,exportType)		




#define D_IMPLEMENT_CLASS_NODE(className)						\
	dInitRtti(className);										\
	className className::m_singletonClass;						\
	static className::dRegisterSingleton m_registerSingletonAgent (#className, &className::m_singletonClass);


#define SerialiseBase(baseClass,rootNode)								\
	TiXmlElement* const baseClassNode = new TiXmlElement (#baseClass);	\
	rootNode->LinkEndChild(baseClassNode);								\
	baseClass::Serialize(baseClassNode);							

#define DeserialiseBase(scene,baseClass,rootNode)															\
	TiXmlElement* const baseClassNode = (TiXmlElement*) rootNode->FirstChild (baseClass::GetClassName());	\
	baseClass::Deserialize (scene, baseClassNode);



class dNodeInfo: public dClassInfo, public dVariableList
{
	public:
	class dRegisterSingleton
	{	
		public:
		dRegisterSingleton (const char* const className, const dNodeInfo* const singleton);
	};

	
	dNodeInfo();
	dNodeInfo(const dNodeInfo& me);
	virtual ~dNodeInfo(void);
	virtual dNodeInfo* MakeCopy () const;
	virtual const char* GetClassName () const;		
	virtual const char* GetBaseClassName ()	const;
	virtual dNodeInfo* MetaFunction(dScene* const world) const;

	virtual const char* GetName () const;
	virtual void SetName (const char* const name);
	
	virtual void Serialize (TiXmlElement* const rootNode) const; 
	virtual bool Deserialize (const dScene* const scene, TiXmlElement* const rootNode);

	// draw scene in wire frame mode
	virtual void DrawWireFrame(dSceneRender* const render, dScene* const scene, dScene::dTreeNode* const myNode) const{dAssert (0);}
	virtual void DrawFlatShaded(dSceneRender* const render, dScene* const scene, dScene::dTreeNode* const myNode) const{dAssert (0);}

	virtual void BakeTransform (const dMatrix& transform){};
	virtual unsigned GetUniqueID() const {return m_uniqueID;}

	static dNodeInfo* CreateFromClassName (const char* const className, dScene* const world);
	static dTree<const dNodeInfo*, dCRCTYPE>& GetSingletonDictionary();
	static void ReplaceSingletonClass (const char* const className, const dNodeInfo* const singleton);

	dAddRtti(dClassInfo,);

	private:
	dString m_name;

	unsigned m_uniqueID;
	static unsigned m_uniqueIDCounter;
};





#endif
