// Copyright 2009-2021 NTESS. Under the terms
// of Contract DE-NA0003525 with NTESS, the U.S.
// Government retains certain rights in this software.
//
// Copyright (c) 2009-2021, NTESS
// All rights reserved.
//
// Portions are copyright of other developers:
// See the file CONTRIBUTORS.TXT in the top level directory
// the distribution for more information.
//
// This file is part of the SST software package. For license
// information, see the LICENSE file in the top level directory of the
// distribution.


#ifndef _H_EMBER_COMPUTE
#define _H_EMBER_COMPUTE

#include <sst/core/rng/gaussian.h>
#include "mpi/embermpigen.h"
#include <string.h>
#include <map>
#include <Python.h>
namespace SST {
namespace Ember {

class EmberCompute {

public:
	EmberCompute(std::string modelfile);
//	~EmberCMT3DGenerator();
	std::map<std::string,uint32_t> staticcompute(std::map<std::string,uint32_t>);
	std::map<std::string,SSTGaussianDistribution*> stochasticcompute(std::map<std::string,uint32_t>);
	std::string modelinput;

	uint32_t my_rank;

	std::map<std::string,uint32_t> outputvalues;

	std::map<std::string,SSTGaussianDistribution*> outputvaluesdist;




};

}
}

#endif
