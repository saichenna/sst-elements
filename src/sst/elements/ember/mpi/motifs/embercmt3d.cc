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

#include "embercmt3d.h"

/*
	CMT communication motif for a 3d decomposition of elements over any
	network topology. Each process communicates with its +-x, +-y & +-z
	neighbors. Note: The motif works only when all the processors in the
	grid are used.

	todo: replace with isend and irecv
*/

using namespace SST::Ember;

EmberCMT3DGenerator::EmberCMT3DGenerator(SST::ComponentId_t id, Params& params) :
	EmberMessagePassingGenerator(id, params, "CMT3D"),
        m_loopIndex(0),
	m_phyIndex(0),
	m_rkstage(0),
        x_pos(-1), x_neg(-1),
        y_pos(-1), y_neg(-1),
        z_pos(-1), z_neg(-1),
        sendx_pos(false), sendx_neg(false),
        sendy_pos(false), sendy_neg(false),
        sendz_pos(false), sendz_neg(false)
{
    	iterations = (uint32_t) params.find("arg.iterations", 1);
    	eltSize = (uint32_t) params.find("arg.elementsize", 5); //polynomial degree = eltSize-1
    	variables = (uint32_t) params.find("arg.variables", 1); //no. of physical quantities
	rkstages = (uint32_t) params.find("arg.rkstages", 3);   //no. of rk stages default 3
	nelt = (uint32_t) params.find("arg.nelt", 1); //no. of elements per processor
	sizeof_cell = params.find<uint32_t>("arg.datatype_width", 8);
	//equationsfile = (std::string) params.find("arg.eqfile", (char*)""); //filename containing analytic models for various compute kernels
	equationsfile = "test.txt";

    	//Distribution of MPI ranks
    	px  = (uint32_t) params.find("arg.px", 8);
    	py  = (uint32_t) params.find("arg.py", 8);
    	pz  = (uint32_t) params.find("arg.pz", 8);
    	threads = (uint32_t) params.find("arg.threads", 1);

    	//Local distribution of elements on each MPI rank
    	mx = (uint32_t) params.find("arg.mx", 10);
    	my = (uint32_t) params.find("arg.my", 10);
    	mz = (uint32_t) params.find("arg.mz", 10);


	//Check to enable computation in the simulation
	
	docompute = (uint32_t) params.find("arg.docompute", true);;
	
	//Set nelt to default (1000) if not specified at runtime
	if (nelt == 1){
    		nelt  = mx * my * mz; //Default : 1000

	}

	// if nelt is specified, setup mx,my,mz 
	else {

		mx = (uint32_t) cbrt(nelt);
		my = (uint32_t) cbrt(nelt);
		mz = (uint32_t) cbrt(nelt);

	}

    	// Calculate computation time in nanoseconds
    	procFlops = (uint64_t) params.find("arg.processorflops", 4); // Nehalem 4 DP FLOPS/cycle
    	procFreq = (uint64_t) params.find("arg.processorfreq", 2);
    	const int64_t flopCount	= pow(eltSize, 3) * (nelt) * (2*eltSize-1); // Not assuming FMAs
    	const double time = (double)flopCount / ( (double)procFlops * (double)procFreq ); //Time in ns

        m_mean = params.find("arg.nsComputeMean", time);
        m_stddev = params.find("arg.nsComputeStddev", (m_mean*0.05));

    	x_xferSize = eltSize*eltSize*my*mz*sizeof_cell;
    	y_xferSize = eltSize*eltSize*mx*mz*sizeof_cell;
    	z_xferSize = eltSize*eltSize*mx*my*sizeof_cell;

        configure();
}


void EmberCMT3DGenerator::configure()
{

    	// Check that we are using all the processors or else lock up will happen :(.
//    	if( (px * py *pz *threads) != (signed)size() ) {
    	if( (px * py *pz *threads) != (signed)size() ) {
    		fatal(CALL_INFO, -1, "Error: CMT3D motif checked processor decomposition: %" \
    			PRId32 "x%" PRId32 "x%" PRId32 "x%" PRIu32 " != MPI World %" PRIu32 "\n",
    			px, py, pz, threads, size());
    	}

    	if(rank() == 0) {
    		output( "CMT3D (Pairwise Exchange) Motif \n" \
    		    "nx1 (elt_size) = %" PRIu32 ", nelt (elts/proc) = %" PRIu32 ", np (total processes) = %" PRIu32 \
    		    ", elements (total) = %" PRIu32 \
     			"\npx: %" PRIu32 " ,py: %" PRIu32 " ,pz: %" PRIu32 " ,threads/proc: %" PRIu32 \
     			"\nmx: %" PRIu32 " ,my: %" PRIu32 " ,mz: %" PRIu32 \
    			"\ncompute time: mean = %fns ,stddev = %fns \
    			\nx_xfer: %" PRIu64 " ,y_xfer: %" PRIu64 " ,z_xfer: %" PRIu64 "\n",
    			eltSize, nelt, size(), nelt*size(),
    			px, py, pz, threads,
    			mx, my, mz,
                m_mean, m_stddev,
    			x_xferSize, y_xferSize, z_xferSize );
    	}

    	// Get our (x,y,z) position and neighboring ranks in a 3D decomposition
	//
	
	//Sai Chenna: bug-fix. initialize m_phyIndex to 0 
	//m_phyIndex = 0;
    	myX=-1; myY=-1; myZ=-1;
    	myID = rank();
    	getPosition(myID, px, py, pz, &myX, &myY, &myZ);

	//Sai Chenna: Create stochastic compute block times for all the compute kernels in CMT3D
	//
	//
	//
	
	generate_computetimes();

	generate_distributions();
	
	/*
	Py_Initialize();
	PyRun_SimpleString("print('Hello world! Successfully able to use python within embercmt3d!!')");
	PyRun_SimpleString("import os");
	PyRun_SimpleString("print(os.getcwd())");
	Py_Finalize();

	*/
	/*
	m_computeconv = generate_computeconvdist();
	m_computedr = generate_computedrdist();
	m_computeds = generate_computedsdist();
	m_computedt = generate_computedtdist();
	m_comminit = generate_comminitdist();
	m_comminitaxis = generate_comminitaxisdist();
	m_computeprepfaces = generate_computeprepfacesdist();
	m_computecleanfaces = generate_computecleanfacesdist();
	m_computesum = generate_computesumdist();
	m_computerk = generate_computerkdist();
	
	*/


        // Initialize Marsaglia RNG for compute time
        m_random = new SSTGaussianDistribution( m_mean, m_stddev,
                                    new RNG::MarsagliaRNG( 7+myID, getJobId() ) );

    	// Set which direction to transfer and the neighbors in that direction
    	if( myX > 0 ) {
    		sendx_neg = true;
    		x_neg	= convertPositionToRank(px, py, pz, myX-1, myY, myZ);
    	}

    	if( myX < px-1 ) {
    		sendx_pos = true;
    		x_pos	= convertPositionToRank(px, py, pz, myX+1, myY, myZ);
    	}

    	if( myY > 0 ) {
    		sendy_neg = true;
    		y_neg 	= convertPositionToRank(px, py, pz, myX, myY-1, myZ);
    	}

    	if( myY < py-1 ) {
    		sendy_pos = true;
    		y_pos	= convertPositionToRank(px, py, pz, myX, myY+1, myZ);
    	}

    	if( myZ > 0 ) {
    		sendz_neg = true;
    		z_neg 	= convertPositionToRank(px, py, pz, myX, myY, myZ-1);
    	}

    	if( myZ < pz-1 ) {
    		sendz_pos = true;
    		z_pos	= convertPositionToRank(px, py, pz, myX, myY, myZ+1);
    	}

    	verbose(CALL_INFO, 2, 0, "Rank: %" PRIu64 " is located at coordinates \
    		(%" PRId32 ", %" PRId32 ", %" PRId32") in the 3D grid,\
    		X+:%s %" PRId32 ", X-:%s %" PRId32 ", \
    		Y+:%s %" PRId32 ", Y-:%s %" PRId32 ", \
    		Z+:%s %" PRId32 ", Z-:%s %" PRId32 "\n",
    		myID,
    		myX, myY, myZ,
    		(sendx_pos ? "Y" : "N"), x_pos,
    		(sendx_neg ? "Y" : "N"), x_neg,
    		(sendy_pos ? "Y" : "N"), y_pos,
    		(sendy_neg ? "Y" : "N"), y_neg,
    		(sendz_pos ? "Y" : "N"), z_pos,
    		(sendz_neg ? "Y" : "N"), z_neg);


	if (rank() == 0){

		output( "Debug: CMT3D Halo3D Generator configure method!! \n" );
	
	
	}

}



bool EmberCMT3DGenerator::generate( std::queue<EmberEvent*>& evQ)
{
        if (m_loopIndex == 0) {
    		verbose(CALL_INFO, 2,0, "rank=%d, size=%d\n", rank(), size());
        }

        //double nsCompute = m_random->getNextDouble();
    	//enQ_compute( evQ, nsCompute );  	// Delay block for compute
	//
	
	double nsCompute;

	if (docompute) {

	nsCompute = m_computeconv->getNextDouble();
	enQ_compute(evQ,nsCompute);

	nsCompute = m_computedr->getNextDouble();
        enQ_compute(evQ,nsCompute);

	nsCompute = m_computeds->getNextDouble();
        enQ_compute(evQ,nsCompute);

	nsCompute = m_computedt->getNextDouble();
        enQ_compute(evQ,nsCompute);

	nsCompute = m_computesum->getNextDouble();
        enQ_compute(evQ,nsCompute);



	nsCompute = m_comminit->getNextDouble();
        enQ_compute(evQ,nsCompute);

	}


    	// +x/-x transfers
	//
	
	if (docompute) {
	
	nsCompute = m_comminitaxis->getNextDouble();
        enQ_compute(evQ,nsCompute);

	}

    	if ( myX % 2 == 0){
    		if (sendx_pos) {

			if (docompute) {
			nsCompute = m_preparefaces->getNextDouble();
        		enQ_compute(evQ,nsCompute);
			}

    			enQ_recv( evQ, x_pos, x_xferSize, 0, GroupWorld );
    			enQ_send( evQ, x_pos, x_xferSize, 0, GroupWorld );

			if (docompute) {
			nsCompute = m_cleanfaces->getNextDouble();
                        enQ_compute(evQ,nsCompute);
			}
    		}
    		if (sendx_neg) {

			if (docompute) {
			nsCompute = m_preparefaces->getNextDouble();
                        enQ_compute(evQ,nsCompute);
			}

    			enQ_recv( evQ, x_neg, x_xferSize, 0, GroupWorld );
    			enQ_send( evQ, x_neg, x_xferSize, 0, GroupWorld );

			if (docompute){
			nsCompute = m_cleanfaces->getNextDouble();
                        enQ_compute(evQ,nsCompute);
			}

    		}
    	}
    	else {
    		if (sendx_pos) {

			if (docompute){
			nsCompute = m_preparefaces->getNextDouble();
                        enQ_compute(evQ,nsCompute);
			}

    			enQ_send( evQ, x_pos, x_xferSize, 0, GroupWorld );
    			enQ_recv( evQ, x_pos, x_xferSize, 0, GroupWorld );

			if (docompute){
			nsCompute = m_cleanfaces->getNextDouble();
                        enQ_compute(evQ,nsCompute);
			}


    		}
    		if (sendx_neg) {

			if (docompute){
			nsCompute = m_preparefaces->getNextDouble();
                        enQ_compute(evQ,nsCompute);
			}

    			enQ_send( evQ, x_neg, x_xferSize, 0, GroupWorld );
    			enQ_recv( evQ, x_neg, x_xferSize, 0, GroupWorld );

			if (docompute){
			nsCompute = m_cleanfaces->getNextDouble();
                        enQ_compute(evQ,nsCompute);
			}
    		}
    	}

    	// +y/-y transfers
	//
	if (docompute){
	nsCompute = m_comminitaxis->getNextDouble();
        enQ_compute(evQ,nsCompute);
	}
    	if ( myY % 2 == 0){
    		if (sendy_pos) {


			if (docompute){
			nsCompute = m_preparefaces->getNextDouble();
                        enQ_compute(evQ,nsCompute);
			}

    			enQ_recv( evQ, y_pos, y_xferSize, 0, GroupWorld );
    			enQ_send( evQ, y_pos, y_xferSize, 0, GroupWorld );

			if (docompute){
			nsCompute = m_cleanfaces->getNextDouble();
                        enQ_compute(evQ,nsCompute);
			}
    		}
    		if (sendy_neg) {


			if (docompute){
			nsCompute = m_preparefaces->getNextDouble();
                        enQ_compute(evQ,nsCompute);
			}

    			enQ_recv( evQ, y_neg, y_xferSize, 0, GroupWorld );
    			enQ_send( evQ, y_neg, y_xferSize, 0, GroupWorld );

			if (docompute){
			nsCompute = m_cleanfaces->getNextDouble();
                        enQ_compute(evQ,nsCompute);
			}
    		}
    	}
    	else {
    		if (sendy_pos) {


			if (docompute){
			nsCompute = m_preparefaces->getNextDouble();
                        enQ_compute(evQ,nsCompute);
			}

    			enQ_send( evQ, y_pos, y_xferSize, 0, GroupWorld );
    			enQ_recv( evQ, y_pos, y_xferSize, 0, GroupWorld );

			if (docompute){
			nsCompute = m_cleanfaces->getNextDouble();
                        enQ_compute(evQ,nsCompute);
			}


    		}
    		if (sendy_neg) {

			if (docompute){
			nsCompute = m_preparefaces->getNextDouble();
                        enQ_compute(evQ,nsCompute);
			}

    			enQ_send( evQ, y_neg, y_xferSize, 0, GroupWorld );
    			enQ_recv( evQ, y_neg, y_xferSize, 0, GroupWorld );

			if (docompute){
			nsCompute = m_cleanfaces->getNextDouble();
                        enQ_compute(evQ,nsCompute);
			}
    		}
    	}

    	// +z/-z transfers
	//
	if (docompute){
	nsCompute = m_comminitaxis->getNextDouble();
        enQ_compute(evQ,nsCompute);
	}
    	if ( myZ % 2 == 0){
    		if (sendz_pos) {

			if (docompute){
			nsCompute = m_preparefaces->getNextDouble();
                        enQ_compute(evQ,nsCompute);
			}


    			enQ_recv( evQ, z_pos, z_xferSize, 0, GroupWorld );
    			enQ_send( evQ, z_pos, z_xferSize, 0, GroupWorld );

			if (docompute){
			nsCompute = m_cleanfaces->getNextDouble();
                        enQ_compute(evQ,nsCompute);
			}
    		}
    		if (sendz_neg) {

			if (docompute){
			nsCompute = m_preparefaces->getNextDouble();
                        enQ_compute(evQ,nsCompute);
			}


    			enQ_recv( evQ, z_neg, z_xferSize, 0, GroupWorld );
    			enQ_send( evQ, z_neg, z_xferSize, 0, GroupWorld );

			if (docompute){
			nsCompute = m_cleanfaces->getNextDouble();
                        enQ_compute(evQ,nsCompute);
			}
    		}
    	}
    	else {
    		if (sendz_pos) {


			if (docompute){
			nsCompute = m_preparefaces->getNextDouble();
                        enQ_compute(evQ,nsCompute);
			}


    			enQ_send( evQ, z_pos, z_xferSize, 0, GroupWorld );
    			enQ_recv( evQ, z_pos, z_xferSize, 0, GroupWorld );

			if (docompute){
			nsCompute = m_cleanfaces->getNextDouble();
                        enQ_compute(evQ,nsCompute);
			}
    		}
    		if (sendz_neg) {

			if (docompute){
			nsCompute = m_preparefaces->getNextDouble();
                        enQ_compute(evQ,nsCompute);
			}


    			enQ_send( evQ, z_neg, z_xferSize, 0, GroupWorld );
    			enQ_recv( evQ, z_neg, z_xferSize, 0, GroupWorld );

			if (docompute){
			nsCompute = m_cleanfaces->getNextDouble();
                        enQ_compute(evQ,nsCompute);
			}
    		}
    	}

        //enQ_barrier( evQ, GroupWorld);

	if (docompute){
	nsCompute = m_computerk->getNextDouble();
        enQ_compute(evQ,nsCompute);
	}

	if ( ++m_rkstage == rkstages) {

    		if ( ++m_loopIndex == iterations) {
            		if ( ++m_phyIndex == variables) {
                		return true;
            		} else {
                		m_loopIndex = 0;
				m_rkstage = 0;
                		return false;
            		}
        	} else {
			m_rkstage = 0;
            		return false;
        	}
	}

	else {
		return false;
	}

	

	//Sai Chenna: Updated
	//
	/*
	if (++m_loopIndex == iterations) {

		return true;
	}

	else {
		return false;

	}
	*/


}



void EmberCMT3DGenerator::generate_computetimes(){

	//1. Setup the python instance.
	//2. In EmberCMT3D, fluid solver phase, workload is static across the ranks and does not change during the execution.
	//   Hence we generate the gaussian distributions for each function at the configure step. 
	//
	

	//if (rank() == 0){
	
	
		//std::cout << "Entered generate_computetimes" << std::endl;
	
	//}

	

	Py_Initialize();

	/* Python evaluate module setup */
	PyRun_SimpleString("import sys; import os");

	//TO-DO: Currently hardcoding the tests folder path. Change it in the future as $SST_EMBER/tests 
	
	PyRun_SimpleString("sys.path.insert(0,'/home/saichenna/scratch/src/sst-elements-library-11.0.0/src/sst/elements/ember/mpi/motifs/equations/')");

	PyRun_SimpleString("sys.path.append(os.getcwd())");

	//PyRun_SimpleString("print(sys.path)");



	//PyObject* myModuleString = PyBytes_FromString((char*)"evaluate");

	//std::string test = PyBytes_AsString(myModuleString);

	//std::cout << "myModuleString = " << test << std::endl;

	/*if (myModuleString == NULL) {

                std::cout << "ERROR generating python string" << std::endl;
                exit(-1);
        }*/

	//PyRun_SimpleString("print(myModuleString)");

   	//PyObject* myModule = PyImport_Import(myModuleString);
   	//PyObject* myModule = PyImport_Import(PyBytes_FromString((char*)"evaluate"));
   	PyObject* myModule = PyImport_ImportModule((char*)"evaluate");

	if (myModule == NULL) {

		std::cout << "ERROR importing module" << std::endl;
		PyErr_Print();
    		exit(-1);
    	}

	PyObject* myFunction = PyObject_GetAttrString(myModule,(char*)"staticcomputetime");

	if (myFunction == NULL) {

		std::cout << "ERROR getting computetime attribute" << std::endl;
		exit(-1);

	}



	//Generate the hashmap/dictionary of input parameter values that are used in the compute model parameters
	
	std::map<std::string,uint32_t> modelparams;

	modelparams.insert(std::pair<std::string,uint32_t>("nelt",nelt));
	modelparams.insert(std::pair<std::string,uint32_t>("eltSize",eltSize));

	// Convert the hashmap into dict (for python)
	PyObject* inputdict = hashMapTodict(modelparams);


	//Convert filename string to python format
        PyObject* filename = PyUnicode_FromString(equationsfile.c_str());

	//Call the corresponding python function which calculates the compute time for all the compute kernels and returns the dictionary
	//
	//std::cout << "File name = " << equationsfile << std::endl;
	
	PyObject* outputdict = PyObject_CallFunctionObjArgs(myFunction, filename, inputdict, NULL);

	if (outputdict == NULL) {


		std::cout << "Error calling the computetimes python function" << std::endl;
		PyErr_Print();
                exit(-1);
		
	
	
	}

	//std::cout << "Finished generating python output dictionary" << std::endl;




	outputvalues = dictTohashMap(outputdict);

	if (rank() == 0){


		std::cout << "Printing the output hashMap containing the compute times for various kernels" << std::endl;

		for (auto iter = outputvalues.begin(); iter != outputvalues.end(); ++iter){

			std::cout << "Kernel: " << iter->first << "\t" << "value: " << iter->second << std::endl;

		}

	}

	Py_Finalize();
	


}


//Convert C++ hashmap into Python dictionary

PyObject* EmberCMT3DGenerator::hashMapTodict(std::map<std::string,uint32_t> parammap){


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


	//PyRun_SimpleString("print(dictObj)");



	return dictObj;


}

//Convert Python dictionary into C++ hashmap

std::map<std::string,double> EmberCMT3DGenerator::dictTohashMap(PyObject* outputdict){


	std::map<std::string,double> retmap;

	//std::cout << "Entered dictTohashMap routine" << std::endl;

	if (PyDict_Check(outputdict)){

		//Generate a list of keys

		PyObject* list = PyDict_Keys(outputdict);

		if (list == NULL) {

			std::cout << "Error generating list of keys from outputdict" << std::endl;
		
			PyErr_Print();
			std::exit(1);
		
		
		}

		//Iterate over the keys and insert the key-value pair in the hashMap

		for(Py_ssize_t i = 0; i < PyList_Size(list); i++){

				//std::cout << "Iterating over the keys" << std::endl;
				PyObject* key = PyList_GetItem(list, i);
				PyObject* value = PyDict_GetItem(outputdict,key);
				if (key == NULL || value == NULL) {
				
					std::cout << "Error extracting key value from python dictionary" << std::endl;

		                        PyErr_Print();
                		        std::exit(1);
				
				
				}
				//std::string key_string = PyUnicode_AsString(key);
				std::string key_string = PyUnicode_AsUTF8(key);
				double val = PyFloat_AsDouble(value);
				retmap.insert(std::pair<std::string,double>(key_string,val));
		}


	}

	else
	{
		throw std::runtime_error("Passed PyObject pointer was not a dictionary!");
		retmap.insert(std::pair<std::string,double>("error",0.0));
	}


	return retmap;


}


void EmberCMT3DGenerator::generate_distributions() {


	//1. Generate the gaussian distributions for each compute kernel based on the time obtained by generate_computetimes routine
	//
	

	double val = outputvalues.find("computedr")->second;
	m_computedr = new SSTGaussianDistribution( val, (val*0.05),new RNG::MarsagliaRNG( 7+rank(), getJobId() ) );
	val = outputvalues.find("computeds")->second;
	m_computeds = new SSTGaussianDistribution( val, (val*0.05),new RNG::MarsagliaRNG( 7+rank(), getJobId() ) );
	val = outputvalues.find("computedt")->second;
	m_computedt = new SSTGaussianDistribution( val, (val*0.05),new RNG::MarsagliaRNG( 7+rank(), getJobId() ) );
	val = outputvalues.find("computeconv")->second;
	m_computeconv = new SSTGaussianDistribution( val, (val*0.05),new RNG::MarsagliaRNG( 7+rank(), getJobId() ) );
	val = outputvalues.find("comminit")->second;
	m_comminit = new SSTGaussianDistribution( val, (val*0.05),new RNG::MarsagliaRNG( 7+rank(), getJobId() ) );
	val = outputvalues.find("comminitaxis")->second;
	m_comminitaxis = new SSTGaussianDistribution( val, (val*0.05),new RNG::MarsagliaRNG( 7+rank(), getJobId() ) );
	val = outputvalues.find("cleanfaces")->second;
	m_cleanfaces = new SSTGaussianDistribution( val, (val*0.05),new RNG::MarsagliaRNG( 7+rank(), getJobId() ) );
	val = outputvalues.find("preparefaces")->second;
	m_preparefaces = new SSTGaussianDistribution( val, (val*0.05),new RNG::MarsagliaRNG( 7+rank(), getJobId() ) );
	val = outputvalues.find("computesum")->second;
	m_computesum = new SSTGaussianDistribution( val, (val*0.05),new RNG::MarsagliaRNG( 7+rank(), getJobId() ) );
	val = outputvalues.find("computerk")->second;
	m_computerk = new SSTGaussianDistribution( val, (val*0.05),new RNG::MarsagliaRNG( 7+rank(), getJobId() ) );
	

}


	


/*
SSTGaussianDistribution* EmberCMT3DGenerator::generate_computeconvdist()
{


	//Compute the time taken for the kernel using analytic model and convert it into nanoseconds
	double mean = ((pow(((1.20612888292*eltSize)-4.29444533027),3)*(nelt)*(3.22111942021e-07))+(9.17890674986e-05))*1e9;


	//Generate a gaussian distribution with S% standard deviation around mean
	SSTGaussianDistribution* random = new SSTGaussianDistribution( mean, (mean*0.05),
                                    new RNG::MarsagliaRNG( 7+rank(), getJobId() ) );


	return random;


}



SSTGaussianDistribution* EmberCMT3DGenerator::generate_computedrdist()
{


        //Compute the time taken for the kernel using analytic model and convert it into nanoseconds
        double mean = (((5.897172491e-09)*(((0.971939537485*nelt)+(2.46932067451))*(pow(eltSize,4)+1.92818093689-eltSize)))+(4.56245565708e-05))*1e9;


        //Generate a gaussian distribution with S% standard deviation around mean
        SSTGaussianDistribution* random = new SSTGaussianDistribution( mean, (mean*0.05),
                                    new RNG::MarsagliaRNG( 7+rank(), getJobId() ) );


        return random;


}

SSTGaussianDistribution* EmberCMT3DGenerator::generate_computeconvdist()
{


        //Compute the time taken for the kernel using analytic model and convert it into nanoseconds
        double mean = ((pow(((1.20612888292*eltSize)-4.29444533027),3)*(nelt)*(3.22111942021e-07))+(9.17890674986e-05))*1e9;


        //Generate a gaussian distribution with S% standard deviation around mean
        SSTGaussianDistribution* random = new SSTGaussianDistribution( mean, (mean*0.05),
                                    new RNG::MarsagliaRNG( 7+rank(), getJobId() ) );


        return random;


}

SSTGaussianDistribution* EmberCMT3DGenerator::generate_computeconvdist()
{


        //Compute the time taken for the kernel using analytic model and convert it into nanoseconds
        double mean = ((pow(((1.20612888292*eltSize)-4.29444533027),3)*(nelt)*(3.22111942021e-07))+(9.17890674986e-05))*1e9;


        //Generate a gaussian distribution with S% standard deviation around mean
        SSTGaussianDistribution* random = new SSTGaussianDistribution( mean, (mean*0.05),
                                    new RNG::MarsagliaRNG( 7+rank(), getJobId() ) );


        return random;


}

SSTGaussianDistribution* EmberCMT3DGenerator::generate_computeconvdist()
{


        //Compute the time taken for the kernel using analytic model and convert it into nanoseconds
        double mean = ((pow(((1.20612888292*eltSize)-4.29444533027),3)*(nelt)*(3.22111942021e-07))+(9.17890674986e-05))*1e9;


        //Generate a gaussian distribution with S% standard deviation around mean
        SSTGaussianDistribution* random = new SSTGaussianDistribution( mean, (mean*0.05),
                                    new RNG::MarsagliaRNG( 7+rank(), getJobId() ) );


        return random;


}

SSTGaussianDistribution* EmberCMT3DGenerator::generate_computeconvdist()
{


        //Compute the time taken for the kernel using analytic model and convert it into nanoseconds
        double mean = ((pow(((1.20612888292*eltSize)-4.29444533027),3)*(nelt)*(3.22111942021e-07))+(9.17890674986e-05))*1e9;


        //Generate a gaussian distribution with S% standard deviation around mean
        SSTGaussianDistribution* random = new SSTGaussianDistribution( mean, (mean*0.05),
                                    new RNG::MarsagliaRNG( 7+rank(), getJobId() ) );


        return random;


}

SSTGaussianDistribution* EmberCMT3DGenerator::generate_computeconvdist()
{


        //Compute the time taken for the kernel using analytic model and convert it into nanoseconds
        double mean = ((pow(((1.20612888292*eltSize)-4.29444533027),3)*(nelt)*(3.22111942021e-07))+(9.17890674986e-05))*1e9;


        //Generate a gaussian distribution with S% standard deviation around mean
        SSTGaussianDistribution* random = new SSTGaussianDistribution( mean, (mean*0.05),
                                    new RNG::MarsagliaRNG( 7+rank(), getJobId() ) );


        return random;


}

SSTGaussianDistribution* EmberCMT3DGenerator::generate_computeconvdist()
{


        //Compute the time taken for the kernel using analytic model and convert it into nanoseconds
        double mean = ((pow(((1.20612888292*eltSize)-4.29444533027),3)*(nelt)*(3.22111942021e-07))+(9.17890674986e-05))*1e9;


        //Generate a gaussian distribution with S% standard deviation around mean
        SSTGaussianDistribution* random = new SSTGaussianDistribution( mean, (mean*0.05),
                                    new RNG::MarsagliaRNG( 7+rank(), getJobId() ) );


        return random;


}

SSTGaussianDistribution* EmberCMT3DGenerator::generate_computeconvdist()
{


        //Compute the time taken for the kernel using analytic model and convert it into nanoseconds
        double mean = ((pow(((1.20612888292*eltSize)-4.29444533027),3)*(nelt)*(3.22111942021e-07))+(9.17890674986e-05))*1e9;


        //Generate a gaussian distribution with S% standard deviation around mean
        SSTGaussianDistribution* random = new SSTGaussianDistribution( mean, (mean*0.05),
                                    new RNG::MarsagliaRNG( 7+rank(), getJobId() ) );


        return random;


}

SSTGaussianDistribution* EmberCMT3DGenerator::generate_computeconvdist()
{


        //Compute the time taken for the kernel using analytic model and convert it into nanoseconds
        double mean = ((pow(((1.20612888292*eltSize)-4.29444533027),3)*(nelt)*(3.22111942021e-07))+(9.17890674986e-05))*1e9;


        //Generate a gaussian distribution with S% standard deviation around mean
        SSTGaussianDistribution* random = new SSTGaussianDistribution( mean, (mean*0.05),
                                    new RNG::MarsagliaRNG( 7+rank(), getJobId() ) );


        return random;


}

*/


