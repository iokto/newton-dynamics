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


#ifndef __DSTRING_H_
#define __DSTRING_H_

#include "dContainersStdAfx.h"
#include "dContainersAlloc.h"


class dString: public dContainersAlloc
{
	class dStringAllocator;
	public:
	DCONTAINER_API dString ();
	DCONTAINER_API dString (char chr);
	DCONTAINER_API dString (const dString& src);
	DCONTAINER_API dString (const char* const data);
	DCONTAINER_API dString (const char* const data, int maxSize);
	DCONTAINER_API dString (int val);
	DCONTAINER_API dString (long long val);
	DCONTAINER_API ~dString ();

	char& operator[] (int index);
	char operator[] (int index) const;
	
	DCONTAINER_API dString& operator= (const dString& src);
	bool operator== (const dString& src) const;
	bool operator!= (const dString& src) const;
	bool operator< (const dString& src) const;
	bool operator> (const dString& src) const;
	bool operator<= (const dString& src) const;
	bool operator>= (const dString& src) const;

	DCONTAINER_API void operator+= (const char* const src);
	void operator+= (const dString& src);

	dString operator+ (const char* const src) const;
	dString operator+ (const dString& src) const;

	DCONTAINER_API int Find (char ch, int from = 0) const;
	int Find (const dString& subString, int from = 0) const;
	DCONTAINER_API int Find (const char* const subString, int from = 0, int lenght = 0x7ffffff) const;

	DCONTAINER_API void Replace (int start, int size, const char* const str, int strSize);
	void Replace (int start, int size, const dString& str);
	void Empty();

	DCONTAINER_API void ToUpper();
	DCONTAINER_API void ToLower();
	DCONTAINER_API int ToInteger() const;
	DCONTAINER_API long long ToInteger64() const;

	int Size() const;
	int Capacity() const;
	DCONTAINER_API void Expand (int size);

	DCONTAINER_API void LoadFile (FILE* const file);
	dString SubString(int start = 0, int size = 0x7fffffff) const;

	const char* GetStr () const;

	private:
	int CalculateSize (const char* const data) const;
	int Compare (const char* const str0, const char* const str1) const;
	void CopyData (char* const dst, const char* const src, int size) const;

	int Find (const char* const subString, int stringSize, int from, int lenght) const;


	protected:
	char* AllocMem(int size);
	void FreeMem (char* const ptr);
	dString (const dString& src, const char* const concatenate, int maxSize);
	
	char* m_string;
	int m_size;
	int m_capacity;

	private:
	dStringAllocator& GetAllocator() const;
};


inline char& dString::operator[] (int index)
{
	dAssert (m_string);
	dAssert (index >= 0);
	dAssert (index < m_size);
	return m_string[index];
}

inline char dString::operator[] (int index) const
{
	dAssert (m_string);
	dAssert (index >= 0);
	dAssert (index < m_size);
	return m_string[index];
}

inline const char* dString::GetStr () const
{
	return m_string;
}

inline int dString::Size() const
{
	return m_size;
}


inline int dString::Find (const char* const subString, int from, int lenght) const
{
	return Find (subString, CalculateSize(subString), from, lenght);
}

inline int dString::Find (const dString& subStream, int from) const
{
	dAssert (subStream.m_string);
	return Find (subStream.m_string, subStream.m_size, from, subStream.m_size);
}

inline void dString::Replace (int start, int size, const dString& str)
{
	Replace(start, size, str.m_string, str.m_size);
}

inline void dString::operator+= (const dString& src)
{
	*this += src.m_string;
}

inline dString dString::operator+ (const dString& src) const
{
	return dString (*this, src.m_string, src.m_size);
}

inline dString dString::operator+ (const char* const copy) const
{
	return dString (*this, copy, CalculateSize (copy));
}


inline int dString::Capacity() const
{
	return m_capacity;
}

inline void dString::CopyData (char* const dst, const char* const src, int size) const
{
	dAssert (dst);
	dAssert (src);
	memcpy (dst, src, size);
}

inline int dString::Compare (const char* const str0, const char* const str1) const
{
	dAssert (str0);
	dAssert (str1);
	return strcmp (str0, str1);
}


inline bool dString::operator== (const dString& src) const
{
	return Compare (m_string, src.m_string) == 0;
}

inline bool dString::operator!= (const dString& src) const
{
	return Compare (m_string, src.m_string) != 0;
}


inline bool dString::operator< (const dString& src) const
{
	return Compare (m_string, src.m_string) < 0;
}

inline bool dString::operator> (const dString& src) const
{
	return Compare (m_string, src.m_string) > 0;
}

inline bool dString::operator<= (const dString& src) const
{
	return Compare (m_string, src.m_string) <= 0;
}

inline bool dString::operator>= (const dString& src) const
{
	return Compare (m_string, src.m_string) >= 0;
}

inline dString dString::SubString(int start, int size) const
{
	dAssert (m_string);
	return dString (&m_string[start], size);
}


#endif


