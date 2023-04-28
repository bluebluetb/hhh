#ifndef DATAGENERATOR_H
#define DATAGENERATOR_H
#include "structures.h"

namespace Generator {
	void Triangulate(ProblemSet &p, bool addEdges=false);
	typedef ProblemSet (*generator_function)(GenParams &genParams);
	struct type {
		const char * name;
		generator_function function;
	};
	extern type generators[];
}

#endif // DATAGENERATOR_H
