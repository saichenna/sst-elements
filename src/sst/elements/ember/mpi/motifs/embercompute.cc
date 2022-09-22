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

#include "embercompute.h"

/*
	Ember compute model API to provide plug-in analytic compute models for Ember motifs.
	
*/

using namespace SST::Ember;

EmberCompute::EmberCompute(std::string modelfile,uint32_t id) 
{
	modelinput = modelfile;

	my_rank = id;


        //configure();
}


std::map<std::string,uint32_t> EmberCompute::staticcompute(std::map<std::string,uint32_t> mparams)
{

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

	modelparams = mparams;

	// Convert the hashmap into dict (for python)
        PyObject* inputdict = hashMapTodict(modelparams);

	//Convert filename string to python format
        PyObject* filename = PyUnicode_FromString(modelinput.c_str());

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

        if (my_rank == 0){


                std::cout << "Printing the output hashMap containing the compute times for various kernels" << std::endl;

                for (auto iter = outputvalues.begin(); iter != outputvalues.end(); ++iter){

                        std::cout << "Kernel: " << iter->first << "\t" << "value: " << iter->second << std::endl;

                }

        }

        Py_Finalize();


	return outputvalues;


}


std::map<std::string,SSTGaussianDistribution*> EmberCompute::stochasticcompute()
{

	// Go over the outputvalues hashmap and generate a gaussian distribution for each key-value pair
	//
	for (auto iter = outputvalues.begin(); iter != outputvalues.end(); ++iter){

		double val = iter->second;

		SSTGaussianDistribution* valdist = new SSTGaussianDistribution(val,(val*0.05),new RNG::MarsagliaRNG(7 + my_rank));


		outputvaluesdist.insert(std::pair<std::string,SSTGaussianDistribution*>(iter->first,valdist));


        }

	return outputvaluesdist;





}






//Convert C++ hashmap into Python dictionary

PyObject* EmberCompute::hashMapTodict(std::map<std::string,uint32_t> parammap){


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

std::map<std::string,double> EmberCompute::dictTohashMap(PyObject* outputdict){


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

