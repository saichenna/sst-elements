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


#include <math.h>
#include <sst_config.h>
#include <sst/core/rng/marsaglia.h>

#include "embercmtcr.h"

/*
	CMT communication motif for Crystal Router algorithm using a 3d decom-
	position of elements over any network topology.

	todo: replace with isend and irecv
*/

using namespace SST::Ember;

EmberCMTCRGenerator::EmberCMTCRGenerator(SST::ComponentId_t id, Params& params) :
	EmberMessagePassingGenerator(id, params, "CMTCR"),
	m_loopIndex(0),
	m_phyIndex(0),
	m_rkstage(0)
{
	iterations = (uint32_t) params.find("arg.iterations", 1);
	rkstages = (uint32_t) params.find("arg.rkstages", 3);   //no. of rk stages default 3
	eltSize = (uint32_t) params.find("arg.elementsize", 5);    //polynomial degree = eltSize-1
	variables = (uint32_t) params.find("arg.variables", 1);    //no. of physical quantities
	nelt = (uint32_t) params.find("arg.nelt",1);               // no. of elements per processor
	npart = (uint32_t) params.find("arg.npart", 100);          // no of particles per processor
	nwallgpart = (uint32_t) params.find("arg.wallgpart", 50);  // no of ghost particles per processor
	//tracefile = (std::string) params.find("arg.tracefile", " "); // file containing all the trace file locations
	//tracefile = "heleshawtraceinput.txt";
	tracefile = "testtraceinput.txt";
	docompute = (uint32_t) params.find("arg.docompute", true);; // check if compute block delays are included in simulation
	istrace = (uint32_t) params.find("arg.istrace", false);;  // Check if trace-based simulation is enabled
	lr = (uint32_t) params.find("arg.lr", 77); // no. of real attributes associated with each particle
	li = (uint32_t) params.find("arg.li", 40); // no. of integer attributes associated with each particle
	equationsfile = "equations.txt";
	//equationsfile = "equationstest.txt";



	//Distribution of processors in the machine
	px  = (uint32_t) params.find("arg.px", 8);
	py  = (uint32_t) params.find("arg.py", 8);
	pz  = (uint32_t) params.find("arg.pz", 8);
	threads = (uint32_t) params.find("arg.threads", 1);

	// Local distribution of elements on each MPI rank
	/*
	mx = (uint32_t) params.find("arg.mx", 10);
	my = (uint32_t) params.find("arg.my", 10);
	mz = (uint32_t) params.find("arg.mz", 10);
	nelt  = mx * my * mz; //Default : 1000
	*/
	// Calculate computation time in nanoseconds
	/*
	procFlops = (uint64_t) params.find("arg.processorflops", 4); //4 DP FLOPS/cycle
	procFreq = (uint64_t) params.find("arg.processorfreq", 2);
	const int64_t flopCount		= pow(eltSize, 3) * (nelt) * (2*eltSize-1); // Not assuming FMAs
	const double time = (double)flopCount / ( (double)procFlops * (double)procFreq ); //Time in ns

    m_mean = params.find("arg.nsComputeMean", time);
    m_stddev = params.find("arg.nsComputeStddev", (m_mean*0.05));
	*/

	xferSize = lr*sizeof(double)*li*sizeof(int);

	configure();
}



void EmberCMTCRGenerator::configure()
{

	/*
	if( (px * py *pz) != (signed)size() ) {
		fatal(CALL_INFO, -1, "Error: CMTCR motif checked processor decomposition: %" \
			PRIu32 "x%" PRIu32 "x%" PRIu32 " != MPI World %" PRIu32 "\n",
			px, py, pz, size());
	}
	*/

	// Get our (x,y) position in a 3D decomposition
	myX = -1; myY = -1; myZ = -1;
	myID = rank();
	
	//getPosition(myID, px, py, pz, &myX, &myY, &myZ);

    stages = (uint32_t) log2( size() );

    //m_phyIndex = 0;

    m_random = new SSTGaussianDistribution( m_mean, m_stddev,
                                new RNG::MarsagliaRNG( 7+myID, getJobId() ) );

    /*
	if(rank() == 0) {
		output( "CMTCR (Crystal Router Reduction) Motif \n" \
		    "element_size = %" PRIu32 ", elements_per_proc = %" PRIu32 ", total processes = %" PRIu32 \
			"\n npx: %" PRIu32 " ,py: %" PRIu32 " ,pz: %" PRIu32 "\n",
			eltSize, nelt, size(),
			px, py, pz );
	}

	verbose(CALL_INFO, 2, 0, "Rank: %" PRIu64 " is located at coordinates \
		(%" PRId32 ", %" PRId32 ", %" PRId32") in the 3D grid, \n",
		myID, myX, myY, myZ );
	*/
	// 1. Read the tracefile in case trace-based simulations is ON and setup the tracefiles

	setuptracefile();

	//std::cout << "Rank " << myID << " finished setuptracefile method!!" << std::endl;




}



bool EmberCMTCRGenerator::generate( std::queue<EmberEvent*>& evQ)
{
        if (m_loopIndex == 0) {
            verbose(CALL_INFO, 2, 0, "rank=%" PRIu64 ", size=%d\n", myID, size());
        }

        //Do this for stage 1
        if (m_rkstage == 0)
	{

        //Get real and wall ghost particles for each rank at given iteration
        npart = getparticles(m_loopIndex+1);
        nwallgpart = getwallghostparticles(m_loopIndex+1);

        setparameterhashmap();


        if (docompute) {
        double nsCompute = getUPL();
	//std::cout << "Mu rank: " << myID << " compute time (ns) for UPL: " << nsCompute << std::endl;
    	enQ_compute( evQ, nsCompute );  	// Delay block for compute
        }

    	if (docompute) {
    	double nsCompute = getcompute1();
	//std::cout << "Mu rank: " << myID << " compute time (ns) for compute1: " << nsCompute << std::endl;
    	enQ_compute( evQ, nsCompute );  	// Delay block for compute1
    	}



    	//Start crystal router communication for moving particles

        uint32_t i = 0;
        uint64_t MASK = 0;
        uint64_t myDest = 0;

	//std::cout << "Rank = " << myID << " beginning particle CR communication!!!" << std::endl;

        for (i=1; i<=stages; i++) {
            MASK = (uint64_t) exp2( stages - i );
            myDest = myID ^ MASK;


            std::vector<uint32_t> movingparticles = getmovingparticles(myDest,i,m_loopIndex+1);

            uint32_t outparts = movingparticles[0];
            uint32_t inparts = movingparticles[2];

	    //std::cout << " My rank : " << myID << " Incoming particles CR message size: " << inparts*xferSize << " timestep : " << m_loopIndex << " commstage : " << i << std::endl;
	    //std::cout << " My rank : " << myID << " Outgoing particles CR message size: " << outparts*xferSize << " timestep : " << m_loopIndex << " commstage : " << i << std::endl;
	    //std::cout << " My rank : " << myID << " Incoming ghost particles CR message size: " << movingparticles[1]*xferSize << " timestep : " << m_loopIndex << " commstage : " << i << std::endl;
	    //std::cout << " My rank : " << myID << " OUtgoing ghost particles CR message size: " << movingparticles[3]*xferSize << " timestep : " << m_loopIndex << " commstage : " << i << std::endl;



        //verbose(CALL_INFO, 2, 0, "stage: %" PRIu32 " \trank: %" PRIu64 " \trecv/send from %" PRIu64 "\n",
        //		    myID, i, myDest);

            if (myID < myDest) {
                enQ_send( evQ, myDest, outparts*xferSize, 0, GroupWorld);
                enQ_recv( evQ, myDest, inparts*xferSize, 0, GroupWorld);
            } else {
                enQ_recv( evQ, myDest, inparts*xferSize, 0, GroupWorld);
                enQ_send( evQ, myDest, outparts*xferSize, 0, GroupWorld);
            }
        }

        enQ_barrier( evQ, GroupWorld);
	//
	//std::cout << "Rank = " << myID << " end of particle CR communication!!!"<< std::endl;

        //Compute event for create ghost particles kernel
        if (docompute) {
        double nsCompute = getCGP();
            	enQ_compute( evQ, nsCompute );  	// Delay block for compute
        }
        //moving wall ghost particles using crystal router
	
	//std::cout << "Rank = " << myID << " beginning ghost particle CR communication!!!" << std::endl;

        	i = 0;
		MASK = 0;
		myDest = 0;

		for (i=1; i<=stages; i++) {
			MASK = (uint64_t) exp2( stages - i );
			myDest = myID ^ MASK;


			std::vector<uint32_t> movingparticles = getmovingparticles(myDest,i,m_loopIndex+1);

			uint32_t outparts = movingparticles[1];
			uint32_t inparts = movingparticles[3];

		//verbose(CALL_INFO, 2, 0, "stage: %" PRIu32 " \trank: %" PRIu64 " \trecv/send from %" PRIu64 "\n",
		//		    myID, i, myDest);

			if (myID < myDest) {
				enQ_send( evQ, myDest, outparts*xferSize, 0, GroupWorld);
				enQ_recv( evQ, myDest, inparts*xferSize, 0, GroupWorld);
			}
			else {
				enQ_recv( evQ, myDest, inparts*xferSize, 0, GroupWorld);
				enQ_send( evQ, myDest, outparts*xferSize, 0, GroupWorld);
			}
		}

		//std::cout << "Rank = " << myID << " end ghost particle CR communication!!!" << std::endl;

        }

        enQ_barrier( evQ, GroupWorld);
        //Compute event for interpolation kernel
        if (docompute) {
		double nsCompute = getIPPL();
	//std::cout << "Mu rank: " << myID << " compute time (ns) for IPPL: " << nsCompute << std::endl;
		enQ_compute( evQ, nsCompute );  	// Delay block for compute
        }

		//Compute event for particle forces (no collision) kernel
		if (docompute) {
	
		double nsCompute = getUPF();
		//std::cout << "Mu rank: " << myID << " compute time (ns) for UPF: " << nsCompute << std::endl;
		enQ_compute( evQ, nsCompute );  	// Delay block for compute
		}

		//Compute event for rk stage computation kernel
		if (docompute) {
		double nsCompute = getRK3();
		//std::cout << "Mu rank: " << myID << " compute time (ns) for getRK: " << nsCompute << std::endl;
		enQ_compute( evQ, nsCompute );  	// Delay block for compute
		}

        //enQ_barrier( evQ, GroupWorld);


        //verbose(CALL_INFO, 2, 0, "Completed CR simulation\n");
        if ( ++m_rkstage == rkstages)
        {
            if ( ++m_loopIndex == iterations)
            {
		    
		if (myID == 0) {
		std::cout << "Mu rank: " << myID << " m_loopIndex: " << m_loopIndex << " m_rkstage: " << m_rkstage << std::endl;
		}
               return true;
            }
            else
            {
        	m_rkstage = 0;
		if (myID == 0){
		std::cout << "Mu rank: " << myID << " m_loopIndex: " << m_loopIndex << " m_rkstage: " << m_rkstage << std::endl;
		}
                return false;
            }
        }
        else
        {

		if (myID == 0) {
		std::cout << "Mu rank: " << myID << " m_loopIndex: " << m_loopIndex << " m_rkstage: " << m_rkstage << std::endl;
		}
        	return false;
      	}




}

//Generate a hashmap containing the key as the application parameters and values as the actual parameter values
void EmberCMTCRGenerator::setparameterhashmap(){

	parameterhashmap["nelt"] = nelt;
	parameterhashmap["eltSize"] = eltSize;
	parameterhashmap["npart"] = npart;
	parameterhashmap["nwallgpart"] = nwallgpart;

}

// Return a pair containing incoming and outgoing real particles for the given rank at a specific crystal router stage
std::vector<uint32_t> EmberCMTCRGenerator::getmovingparticles(uint64_t myDest, uint32_t commstage, uint32_t loopindex) {

	//1. Generate two bins of ranks to define the source and destination group at this CR stage
	uint32_t binsize = exp2( stages - commstage );

	std::vector<uint32_t> low;
	std::vector<uint32_t> high;
	uint32_t accumulator1,accumulator2,accumulator3,accumulator4;

	uint32_t index;

	if (myID < myDest) {
		index = myID;
	}

	else {

		index = myDest;

	}

	int rem = (int) index / binsize;

	for(int i=0; i < binsize; i++) {

		low.push_back(rem*binsize+i);
		high.push_back((rem+1)*binsize+i);
	}


// Calculate the send particles

	if (myID < myDest)
	{

		//1. Send to the myDest rank all the particles which myID is supposed to send to the ranks in high group

		accumulator1 = 0;
		accumulator2 = 0;



		for(int i=0; i < high.size(); i++)
		{

			std::vector<int> key{myID,high[i],loopindex};



			uint32_t value,value2;

			std::map<std::vector<int>,int>::const_iterator it = partmvmthashmap.find(key);

			if(it == partmvmthashmap.end())
			{

				value = 0;

			}

			else
			{

				value = (uint32_t) it->second;

			}

			accumulator1 += value;

			std::map<std::vector<int>,int>::const_iterator it2 = wallpartmvmthashmap.find(key);

			if(it2 == wallpartmvmthashmap.end())
			{

				value2 = 0;

			}

			else
			{

				value2 = (uint32_t) it2->second;

			}

			accumulator2 += value2;


		}

	}

	else
	{

			//1. Send to the myDest rank all the particles which myID is supposed to send to the ranks in low group

			accumulator1 = 0;
			accumulator2 = 0;


			for(int i=0; i < low.size(); i++)
			{
				std::vector<int> key{myID,low[i],loopindex};

				uint32_t value,value2;

				std::map<std::vector<int>,int>::const_iterator it = partmvmthashmap.find(key);

				if(it == partmvmthashmap.end())
				{

					value = 0;

				}

				else
				{

					value = (uint32_t) it->second;

				}

				accumulator1 += value;

				std::map<std::vector<int>,int>::const_iterator it2 = wallpartmvmthashmap.find(key);

				if(it2 == wallpartmvmthashmap.end())
				{

					value2 = 0;

				}

				else
				{

					value2 = (uint32_t) it2->second;

				}

				accumulator2 += value2;
			}

	}


//Calculate the receive particles

	if (myID < myDest)
	{

		//1. Calculate how many particles are being sent to you from your myDest rank, it follows the same logic myRank used to send

		accumulator3 = 0;
		accumulator4 = 0;



		for(int i=0; i < low.size(); i++)
		{

			std::vector<int> key{myDest,low[i],loopindex};



			uint32_t value3,value4;

			std::map<std::vector<int>,int>::const_iterator it3 = partmvmthashmap.find(key);

			if(it3 == partmvmthashmap.end())
			{

				value3 = 0;

			}

			else
			{

				value3 = (uint32_t) it3->second;

			}

			accumulator3 += value3;

			std::map<std::vector<int>,int>::const_iterator it4 = wallpartmvmthashmap.find(key);

			if(it4 == wallpartmvmthashmap.end())
			{

				value4 = 0;

			}

			else
			{

				value4 = (uint32_t) it4->second;

			}

			accumulator4 += value4;


		}

	}

	else
	{

			accumulator3 = 0;
			accumulator4 = 0;


			for(int i=0; i < high.size(); i++)
			{
				std::vector<int> key{myDest,high[i],loopindex};

				uint32_t value3,value4;

				std::map<std::vector<int>,int>::const_iterator it3 = partmvmthashmap.find(key);

				if(it3 == partmvmthashmap.end())
				{

					value3 = 0;

				}

				else
				{

					value3 = (uint32_t) it3->second;

				}

				accumulator3 += value3;

				std::map<std::vector<int>,int>::const_iterator it4 = wallpartmvmthashmap.find(key);

				if(it4 == wallpartmvmthashmap.end())
				{

					value4 = 0;

				}

				else
				{

					value4 = (uint32_t) it4->second;

				}

				accumulator4 += value4;
			}

	}




	std::vector<uint32_t> movparts{accumulator1,accumulator2,accumulator3,accumulator4};

	return movparts;


}


void EmberCMTCRGenerator::setuptracefile(){

	//1. Read the main trace file. This should contain some basic information like # of timesteps, followed by the names of other trace files
	//2. Setup a hash-table for each of the individual workload traces, such as:
	//    a. particle-count trace file
	//    b. particle-movement trace file
	//    c. wall ghost particle-count trace file
	//    d. wall ghost particle-movement trace file
	//3. Read the equations file and setup the SR analytic models
	//4. Generate the parameter dictionary with application parameters as values
	// Hint: Maybe use call a python function which makes it easier to parse the trace files and generate corresponding dictionaries

	Py_Initialize();

		/* Python evaluate module setup */
		PyRun_SimpleString("import sys; import os");

		//TO-DO: Currently hardcoding the tests folder path. Change it in the future as $SST_EMBER/tests

		PyRun_SimpleString("sys.path.insert(0,'/home/saichenna/scratch/src/sst-elements-library-11.0.0/src/sst/elements/ember/mpi/motifs/equations/')");

		PyRun_SimpleString("sys.path.append(os.getcwd())");


	   	PyObject* myModule = PyImport_ImportModule((char*)"evaluate");

		if (myModule == NULL) {

			std::cout << "ERROR importing module" << std::endl;
			PyErr_Print();
	    		exit(-1);
	    	}

		PyObject* initTraceFunction = PyObject_GetAttrString(myModule,(char*)"setuptrace");

		if (initTraceFunction == NULL) {

			std::cout << "ERROR getting setuptrace attribute" << std::endl;
			exit(-1);

		}

		PyObject* filename = PyUnicode_FromString(tracefile.c_str());

		PyObject* outputdict = PyObject_CallFunctionObjArgs(initTraceFunction, filename, NULL);

		std::vector<std::string> outputlist = listToVector(outputdict);

		PyObject* partcountTraceFunction = PyObject_GetAttrString(myModule,(char*)"partcounttrace");

		PyObject* partmvmtTraceFunction = PyObject_GetAttrString(myModule,(char*)"partmvmttrace");

		//PyObject* wallpartcountTraceFunction = PyObject_GetAttrString(myModule,(char*)"wallpartcounttrace");

		//PyObject* wallpartmvmtTraceFunction = PyObject_GetAttrString(myModule,(char*)"wallpartmvmttrace");

		PyObject* equationsFunction = PyObject_GetAttrString(myModule,(char*)"setupequations");


		PyObject* parcountfilename = PyUnicode_FromString(outputlist[0].c_str());

		PyObject* parmvmtfilename = PyUnicode_FromString(outputlist[1].c_str());

		PyObject* wallparcountfilename = PyUnicode_FromString(outputlist[2].c_str());

		PyObject* wallparmvmtfilename = PyUnicode_FromString(outputlist[3].c_str());

		PyObject* equationfilename = PyUnicode_FromString(outputlist[4].c_str());



		//Compute the dynamic workload on each PE
		PyObject* partcountdict = PyObject_CallFunctionObjArgs(partcountTraceFunction, parcountfilename, NULL);

		PyObject* partmvmtdict = PyObject_CallFunctionObjArgs(partmvmtTraceFunction, parmvmtfilename, NULL);

		PyObject* wallpartcountdict = PyObject_CallFunctionObjArgs(partcountTraceFunction, wallparcountfilename, NULL);

		PyObject* wallpartmvmtdict = PyObject_CallFunctionObjArgs(partmvmtTraceFunction, wallparmvmtfilename, NULL);

		//PyObject* wallpartcountdict = PyObject_CallFunctionObjArgs(wallpartcountTraceFunction, wallparcountfilename, NULL);

		//PyObject* wallpartmvmttdict = PyObject_CallFunctionObjArgs(wallpartmvmtTraceFunction, wallparmvmtfilename, NULL);

		PyObject* equationsdict = PyObject_CallFunctionObjArgs(equationsFunction, equationfilename, NULL);

		//Convert the python dictionaries to the corresponding hashmaps

		partcounthashmap = workloaddictTohashmap(partcountdict);
		wallpartcounthashmap = workloaddictTohashmap(wallpartcountdict);

		partmvmthashmap = workloaddictTohashmap(partmvmtdict);
		wallpartmvmthashmap = workloaddictTohashmap(wallpartmvmtdict);

		equationshashmap = equationsdictTohashmap(equationsdict);


		Py_Finalize();



}


//Convert the python dictionary int C++ hashmap
std::map<std::vector<int>,int> EmberCMTCRGenerator::workloaddictTohashmap(PyObject* outputdict){


	std::map<std::vector<int>,int> retmap;

	int tmp = 0;

	//std::cout << "Entered dictTohashMap routine" << std::endl;

	if (PyDict_Check(outputdict)){

		//Generate a list of keys

		PyObject* list = PyDict_Keys(outputdict);

		if (list == NULL) {

			std::cout << "Error generating list of keys from dict" << std::endl;

			PyErr_Print();
			std::exit(1);


		}

		//Iterate over the keys and insert the key-value pair in the hashMap

		for(Py_ssize_t i = 0; i < PyList_Size(list); i++){

				//std::cout << "Iterating over the keys" << std::endl;
				PyObject* key = PyList_GetItem(list, i);
				PyObject* value = PyDict_GetItem(outputdict,key);

				if (PyTuple_Check(key)){

					std::vector<int> hashmapkey;

					for (Py_ssize_t i = 0;  i < PyTuple_Size(key); i++ ) {

						PyObject* val = PyTuple_GetItem(key,i);

						tmp = (int) PyLong_AsLong(val);

						hashmapkey.push_back(tmp);


					}

					retmap[hashmapkey] = (int) PyLong_AsLong(value);


				}



				if (key == NULL || value == NULL) {

					std::cout << "Error extracting key value from python dictionary" << std::endl;

		                        PyErr_Print();
                		        std::exit(1);


				}

		}


	}

	else
	{
		throw std::runtime_error("Passed PyObject pointer was not a dictionary!");
		std::exit(1);
	}


	return retmap;

}


//Convert the python dictionary int C++ hashmap
std::map<std::string,std::string> EmberCMTCRGenerator::equationsdictTohashmap(PyObject* dict){


	std::map<std::string,std::string> retmap;

	//std::cout << "Entered dictTohashMap routine" << std::endl;

	if (PyDict_Check(dict)){

		//Generate a list of keys

		PyObject* list = PyDict_Keys(dict);

		if (list == NULL) {

			std::cout << "Error generating list of keys from outputdict" << std::endl;

			PyErr_Print();
			std::exit(1);


		}

		//Iterate over the keys and insert the key-value pair in the hashMap

		for(Py_ssize_t i = 0; i < PyList_Size(list); i++){

				//std::cout << "Iterating over the keys" << std::endl;
				PyObject* key = PyList_GetItem(list, i);
				PyObject* value = PyDict_GetItem(dict,key);

				if (key == NULL || value == NULL) {

					std::cout << "Error extracting key value from python dictionary" << std::endl;

		                        PyErr_Print();
                		        std::exit(1);


				}

				retmap[PyUnicode_AsUTF8(key)] = PyUnicode_AsUTF8(value);



		}


	}

	else
	{
		throw std::runtime_error("Passed PyObject pointer was not a dictionary!");
		retmap.insert(std::pair<std::string,std::string>("error","0.0"));
	}


	return retmap;

}


PyObject* EmberCMTCRGenerator::equationshashMapTodict(std::map<std::string,std::string> parammap){


	PyObject* dictObj = PyDict_New();

	for (auto iter = parammap.begin(); iter != parammap.end(); ++iter){

		PyObject* parameter = PyUnicode_FromString(iter->first.c_str());
		PyObject* value = PyUnicode_FromString(iter->second.c_str());

		if (!parameter || !value) {

			Py_DECREF(dictObj);

			throw std::runtime_error("Unable to allocate memory for Python dict");
		}

		PyDict_SetItem(dictObj, parameter, value);


    	}


	//PyRun_SimpleString("print(dictObj)");



	return dictObj;


}




std::vector<std::string> EmberCMTCRGenerator::listToVector(PyObject* incoming){

	std::vector<std::string> data;

	if(PyList_Check(incoming))
	{
	    for(Py_ssize_t i = 0; i < PyList_Size(incoming); i++){
				PyObject *value = PyList_GetItem(incoming, i);
				data.push_back(PyUnicode_AsUTF8(value));
		}
	}
	else
	{
		throw std::runtime_error("Passed PyObject pointer was not a list!");
		data.push_back("error");
	}

	return data;

}


uint32_t EmberCMTCRGenerator::getparticles(uint32_t timestep){


	std::vector<int> key;

	key.push_back((int)myID);
	key.push_back((int)timestep);
	uint32_t val = partcounthashmap[key];

	return val;

	//return 100;
}

uint32_t EmberCMTCRGenerator::getwallghostparticles(uint32_t timestep){

	std::vector<int> key;

	key.push_back((int)myID);
	key.push_back((int)timestep);
	uint32_t val = wallpartcounthashmap[key];

	return val;


	//return 50;
}

PyObject* EmberCMTCRGenerator::paramshashmapTodict(std::map<std::string,int> parammap) {

	PyObject* dictObj = PyDict_New();

		for (auto iter = parammap.begin(); iter != parammap.end(); ++iter){

			PyObject* parameter = PyUnicode_FromString(iter->first.c_str());
			PyObject* value = PyFloat_FromDouble( (double) iter->second);

			if (!parameter || !value) {

				Py_DECREF(dictObj);

				throw std::runtime_error("Unable to allocate memory for Python dict");
			}

			PyDict_SetItem(dictObj, parameter, value);


	    }

	return dictObj;


}

double EmberCMTCRGenerator::getcomputetime(std::string eq_string){

	std::cout << "My rank " << myID << " entered getcomputetime method!!!" << std::endl;


	//1. Setup the python environment
	//2. Send the equation string along with the python dictionary containing parameter-name and values

				Py_Initialize();


				PyErr_Print();

				/* Python evaluate module setup */
				PyRun_SimpleString("import sys; import os");

				//TO-DO: Currently hardcoding the tests folder path. Change it in the future as $SST_EMBER/tests

				PyRun_SimpleString("sys.path.insert(0,'/home/saichenna/scratch/src/sst-elements-library-11.0.0/src/sst/elements/ember/mpi/motifs/equations/')");

				PyRun_SimpleString("sys.path.append(os.getcwd())");


			   	PyObject* myModule = PyImport_ImportModule((char*)"evaluate");

				if (myModule == NULL) {

					std::cout << "ERROR importing module" << std::endl;
					PyErr_Print();
			    		exit(-1);
			    	}

				//std::cout << "My rank " << myID << " is about to enter dynamiccomputetime method" << std::endl;

				PyObject* equationFunction = PyObject_GetAttrString(myModule,(char*)"dynamiccomputetime");

				if (equationFunction == NULL) {

					std::cout << "ERROR getting dynamiccomputetime attribute" << std::endl;
					exit(-1);

				}

				//Get the equation string


				PyObject* equation_string = PyUnicode_FromString(eq_string.c_str());


				PyObject* paramdict = paramshashmapTodict(parameterhashmap);

				std::cout << "My rank " << myID << " is about to enter dynamiccomputetime method" << std::endl;
				PyObject* computetime = PyObject_CallFunctionObjArgs(equationFunction, equation_string, paramdict, NULL);

				double val = PyFloat_AsDouble(computetime);

				Py_Finalize();

				return val;



}

double EmberCMTCRGenerator::getUPL(){

			//Get the equation string

			double value;

			std::string eq_string = equationshashmap["upl"];

			value = getcomputetime(eq_string);


			return value;


			//return 100;


}

double EmberCMTCRGenerator::getcompute1(){

	//Get the equation string

	double value;

	std::string eq_string = equationshashmap["compute1"];

	value = getcomputetime(eq_string);


	return value;


	//return 100;


	//return 100;
}

double EmberCMTCRGenerator::getCGP(){

	//Get the equation string

	double value;

	std::string eq_string = equationshashmap["cgp"];

	value = getcomputetime(eq_string);


	return value;


	//return 100;


}

double EmberCMTCRGenerator::getIPPL(){

	//Get the equation string

	double value;

	std::string eq_string = equationshashmap["ippl"];

	value = getcomputetime(eq_string);


	return value;


	//return 100;
}

double EmberCMTCRGenerator::getUPF(){

	//Get the equation string

	double value;

	std::string eq_string = equationshashmap["upf"];

	value = getcomputetime(eq_string);


	return value;


	//return 100;
}

double EmberCMTCRGenerator::getRK3(){

	//Get the equation string

	double value;

	std::string eq_string = equationshashmap["rk3"];

	value = getcomputetime(eq_string);


	return value;


	//return 100;
}

/*

uint32_t EmberCMTCRGenerator::getgoingparticles(uint64_t mydest,uint32_t stage){

	return 100;
}

uint32_t EmberCMTCRGenerator::getcomingparticles(uint64_t mydest,uint32_t stage){

	return 100;
}

uint32_t EmberCMTCRGenerator::getgoingwallparticles(uint64_t mydest,uint32_t stage){

	return 100;
}

uint32_t EmberCMTCRGenerator::getcomingwallparticles(uint64_t mydest,uint32_t stage){

	return 100;
}

*/












