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


#include "dLexCompiler.h"


int main(int argc, char* argv[])
{
#ifdef _MSC_VER
	_CrtSetDbgFlag( _CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF );
#endif

	if (argc < 3) {
		fprintf (stdout, "usage: dLexCompiler [inputRulesFileName] [outputFileName]\n");
		fprintf (stdout, "[inputRulesFileName] name of the file containing a Lex or Flex rules\n");
		fprintf (stdout, "[outputFileName] name of the file cpp output file\n");
		exit (0);
	}


	const char* const inputRulesFileName = argv[1];
	const char* const outputFileName = argv[2];

	FILE* const rules = fopen (inputRulesFileName, "rb");
	if (!rules) {
		fprintf (stdout, "Rule file \"%s\" not found\n",  inputRulesFileName);
		exit (0);
	}

	dString buffer;
	buffer.LoadFile(rules);
	fclose (rules);
	dLexCompiler lexical (buffer, outputFileName);

	return 0;
}
