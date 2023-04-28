#ifndef SOLVER_H
#define SOLVER_H

#include <fstream>
#include <iostream>
#include "structures.h"

class BaseTask;
namespace Solver {
	typedef int (*solver_function)(BaseTask * t);
	struct type {
		const char * name;
		solver_function function;
	};
	extern type solvers[];
	static std::ofstream logFile;
};

#endif // SOLVER_H
