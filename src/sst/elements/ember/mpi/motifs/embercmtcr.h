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


#ifndef _H_EMBER_CMT_CR
#define _H_EMBER_CMT_CR

#include <sst/core/rng/gaussian.h>
#include "mpi/embermpigen.h"
#include <string.h>
#include <map>
#include <Python.h>


namespace SST {
namespace Ember {

class EmberCMTCRGenerator : public EmberMessagePassingGenerator {

public:
    SST_ELI_REGISTER_SUBCOMPONENT_DERIVED(
        EmberCMTCRGenerator,
        "ember",
        "CMTCRMotif",
        SST_ELI_ELEMENT_VERSION(1,0,0),
        "Performs all-to-all communication using Crystal Router",
        SST::Ember::EmberGenerator
    )

    SST_ELI_DOCUMENT_PARAMS(
        {   "arg.iterations",   "Sets the number of data exchanges to perform", "1"},
        {   "arg.elementsize",  "Sets the number of gridpoints per element", "10"},
        {   "arg.variables",    "Sets the number of physical quantities for which derivatives are calculated", "1"},
        {   "arg.px",           "Sets the size of the processors in the machine in x_dim", "4"},
        {   "arg.py",           "Sets the size of the processors in the machine in y_dim", "4"},
        {   "arg.pz",           "Sets the size of the processors in the machine in z_dim", "4"},
        {   "arg.threads",          "Sets the number of MPI threads per processor", "1"},
        {   "arg.mx",           "Sets the number of elements per processor in x dim",   "10"},
        {   "arg.my",           "Sets the number of elements per processor in y dim",   "10"},
        {   "arg.mz",           "Sets the number of elements per processor in z dim",   "10"},
        {   "arg.processorflops",   "Sets the processor flops for compute time estimation", "4"},
        {    "arg.processorfreq",    "Sets the processor frequency for compute time estimation", "2.5"},
        {   "arg.nsComputeMean",    "Sets the mean compute time per processor", "1000"},
        {   "arg.nsComputeStddev",  "Sets the stddev in compute time per processor", "50"},
    )

    SST_ELI_DOCUMENT_STATISTICS(
        { "time-Init", "Time spent in Init event",          "ns",  0},
        { "time-Finalize", "Time spent in Finalize event",  "ns", 0},
        { "time-Rank", "Time spent in Rank event",          "ns", 0},
        { "time-Size", "Time spent in Size event",          "ns", 0},
        { "time-Send", "Time spent in Recv event",          "ns", 0},
        { "time-Recv", "Time spent in Recv event",          "ns", 0},
        { "time-Irecv", "Time spent in Irecv event",        "ns", 0},
        { "time-Isend", "Time spent in Isend event",        "ns", 0},
        { "time-Wait", "Time spent in Wait event",          "ns", 0},
        { "time-Waitall", "Time spent in Waitall event",    "ns", 0},
        { "time-Waitany", "Time spent in Waitany event",    "ns", 0},
        { "time-Compute", "Time spent in Compute event",    "ns", 0},
        { "time-Barrier", "Time spent in Barrier event",    "ns", 0},
        { "time-Alltoallv", "Time spent in Alltoallv event", "ns", 0},
        { "time-Alltoall", "Time spent in Alltoall event",  "ns", 0},
        { "time-Allreduce", "Time spent in Allreduce event", "ns", 0},
        { "time-Reduce", "Time spent in Reduce event",      "ns", 0},
        { "time-Bcast", "Time spent in Bcast event",        "ns", 0},
        { "time-Gettime", "Time spent in Gettime event",    "ns", 0},
        { "time-Commsplit", "Time spent in Commsplit event", "ns", 0},
        { "time-Commcreate", "Time spent in Commcreate event", "ns", 0},
    )

public:
	EmberCMTCRGenerator(SST::ComponentId_t, Params& params);
//	~EmberCMT3DGenerator();
	void configure();
	bool generate( std::queue<EmberEvent*>& evQ);

private:
// User parameters - application
	uint32_t iterations;	// Total no. of timesteps being simulated
	uint32_t eltSize;		// Size of element (5-20)
	uint32_t variables;     // No. of physical quantities
	uint32_t rkstages;      // No. of range-kutta stages (3)
	uint32_t nelt;			// Total no. of elements per process (100-10,000)
	std::string equationsfile; //filename containing analytic models for various compute kernels
	std::string tracefile; //filename of the main trace file containing details for trace-based simulation
	uint32_t npart;      // No. of particles per process
	uint32_t nwallgpart;			// Total no. of ghost particles per process
	uint32_t lr; // no. of real attributes associated with each particle
	uint32_t li; // no. of integer attributes associated with each particle


// User parameters - machine
	int32_t px;				// Machine size (no. of nodes in 3d dimensions)
	int32_t py;
	int32_t pz;
    uint32_t threads;

// User parameters - mpi rank
	uint32_t mx;			// Local distribution of the elements on a MPI rank
	uint32_t my;
	uint32_t mz;
//	uint64_t nsCompute;		// computation time (or delay)

// User parameters - processor
	uint64_t procFlops;		// no. of FLOPS/cycle for the processor
	uint64_t procFreq;		// operating frequency of the processor
	double m_mean;
	double m_stddev;

// Model parameters
	uint32_t m_loopIndex;   // Loop over 'iterations'
    uint32_t m_phyIndex;    // Loop over 'variables'
    uint32_t m_rkstage;     // Loop over 'rkstages'
	int32_t myX;			// Local (x,y,z) coordinates and rank
	int32_t myY;
	int32_t myZ;
	uint64_t myID;          // my MPI rank
	uint64_t xferSize;      // Average transfer size
	uint32_t stages;        // No. of Stages in cyrstal router reduction
	SSTGaussianDistribution* m_random;
	bool docompute;       // Set this flag to true to include computation blocks during simulation
	bool istrace;       // Set this flag to true to include trace based simulation
	std::map<std::vector<int>,int> partcounthashmap;
	std::map<std::vector<int>,int> wallpartcounthashmap;
	std::map<std::vector<int>,int> partmvmthashmap;
	std::map<std::vector<int>,int> wallpartmvmthashmap;
	std::map<std::string,std::string> equationshashmap;
	std::map<std::string,int> parameterhashmap;
	//std::map<uint32_t,uint32_t> recvbuffer;



// Helper methods
	uint32_t getparticles(uint32_t);
	uint32_t getwallghostparticles(uint32_t);
	double getUPL();
	double getcompute1();
	double getCGP();
	double getIPPL();
	double getUPF();
	double getRK3();
	//uint32_t getgoingparticles(uint64_t,uint32_t,uint32_t);
	//uint32_t getcomingparticles(uint64_t,uint32_t,uint32_t);
	//uint32_t getgoingwallparticles(uint64_t,uint32_t,uint32_t);
	//uint32_t getcomingwallparticles(uint64_t,uint32_t,uint32_t);
	std::vector<uint32_t> getmovingparticles(uint64_t,uint32_t,uint32_t);
	void setuptracefile();
	void setparameterhashmap();
	std::map<std::vector<int>,int> workloaddictTohashmap(PyObject*);
	std::map<std::string,std::string> equationsdictTohashmap(PyObject*);
	PyObject* equationshashMapTodict(std::map<std::string,std::string>);
	std::vector<std::string> listToVector(PyObject*);
	PyObject* paramshashmapTodict(std::map<std::string,int>);
	double getcomputetime(std::string);

};

}
}

#endif
