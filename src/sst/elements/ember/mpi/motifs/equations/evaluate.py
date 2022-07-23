import math
import os
import sys

		 
def dynamiccomputetime(equation_string,input_dict):

    def add(a,b):
        return (a+b)
    def sub(a,b):
        return (a-b)
    def mul(a,b):
        return (a*b)
    def div(a,b):
        return (a/b)
    def square(a):
        return (a*a)
    def cube(a):
        return (a*a*a)
    def ln(a):
        return (math.log(a))
        
        
    equation_expression = equation_string

    print("Inside python dynamiccomputetime method!!")

    print("Inside python input_dict")

    print(input_dict)
    
    print("Equation string")
    
    print(equation_string)
    
    for key in input_dict.keys():
    
    	#print("key = "+key+" value = "+str(inputdict[key]))
    
    	equation_expression = equation_expression.replace(key,str(input_dict[key]))
    #print("equation_expression")
    #print(equation_expression)
    	
    	
    predtime = eval(equation_expression)

    #Assuming that the analytical model output is in seconds, convert it into nanoseconds for our ember compute block
    #predtime_ns = predtime*10e9
    predtime_ns = predtime*10e6

    return float(predtime_ns)


def staticcomputetime(equationfile, inputdict):
#def lookupEquation(equation, input1, i_scheme):
    #print ("[lookupEquation] ---  ", equation)
    #print ("inputs")
    #print (input1)              
    #print (i_scheme)
    def add(a,b):
        return (a+b)
    def sub(a,b):
        return (a-b)
    def mul(a,b):
        return (a*b)
    def div(a,b):
        return (a/b)
    def square (a):
        return (a*a)
    def cube (a):
        return (a*a*a)
    def ln(a):
        return (math.log(a))

    #print("File name is "+str(equationfile))

    #print(os.getcwd())
    #print(sys.path)
    #print("inputdict = ")
    #print(inputdict)
    
    fp = open(equationfile,'r')
    
    length = len(inputdict.keys())

    lines = fp.read().splitlines()

    equations = dict()

    outputdict = dict()

    for l in lines:
        if len(l) != 0:
            if l[0] != '#':
                variable=l.split("=")
                if variable[0] == "Parameters" or variable[0] == "parameters" or variable[0] == "PARAMETERS":
                    params = variable[1].split(",")
                else:
                    #Store the equation models in a dictionary
                    equations[variable[0]] = variable[1]

        
    # print(params)
    # print(equation_expression)

    fp.close()

    if len(params) != length:
        print ("Number of app parameters in " +equation+ " file does not match the one specified in Ember motif")
        exit()

    for key in equations.keys():

        equation_expression = equations[key]

        for x in range(length):
            equation_expression = equation_expression.replace(params[x],str(inputdict[params[x]]))

        predtime = eval(equation_expression)

        #Assuming that the analytical model output is in seconds, convert it into nanoseconds for our ember compute block
        predtime_ns = predtime*10e9

        outputdict[key] = float(predtime_ns)


    #print(outputdict)

    return outputdict


#Method to read the main trace files and return the corresponding filename strings back
def setuptrace(filename):

    fp = open(filename,'r')

    lines = fp.read().splitlines()

    for l in lines:
    
    	if len(l) != 0:
    	
    		if l[0] != '#':
    		
    			variable=l.split("=")
    			
    			if variable[0] == "particlecount" :
    				parcountfile = variable[1]
    				
    			elif variable[0] == "particlemovement" :
    				parmvmtfile = variable[1]
    				
    			elif variable[0] == "wallparticlecount" :
    				wallparcountfile = variable[1]
    				
    			elif variable[0] == "wallparticlemovement" :
    				wallparmvmtfile = variable[1]
    				
    			elif variable[0] == "equations" :
    				equationsfile = variable[1]
    				
    outputlist = [parcountfile,parmvmtfile,wallparcountfile,wallparmvmtfile,equationsfile]
    
    fp.close()
    
    return outputlist
    

#Read the trace file and generate a python dictionary in the following format: key=(rank,timestep), value = # of particles
#input trace format: <rank>,<time-step>,<particle-count>
def partcounttrace(filename):

	fp = open(filename,'r')
	
	countdict = dict()
	
	#Divide time-step value by delta, as the trace is collected for every delta time-steps. Most of the trace files have a delta value of 100
	delta = 100
	
	lines = fp.read().splitlines()
	
	for l in lines:
	
		l = l.rstrip()
	
		ls = l.split(",")
	
		countdict[(int(ls[0]),int(int(ls[1])/delta))] = int(ls[2])
		
		
	fp.close()
	
	return countdict
	

#Read the particle movement trace file and generate a python dictionary in the following format: key=(source-rank,destination-rank,timestep); value = # of particles moved	
#input trace format: <source-rank>,<destination-rank>,<time-step>,<# of particles moving across processors>
def partmvmttrace(filename):

	fp = open(filename,'r')
	
	mvmtdict = dict()
	
	#Divide time-step value by delta, as the trace is collected for every delta time-steps. Most of the trace files have a delta value of 100
	delta = 100
	
	lines = fp.read().splitlines()
	
	for l in lines:
	
		l = l.rstrip()
		
		ls = l.split(",")
		
		mvmtdict[(int(ls[0]),int(ls[1]),int(int(ls[2])/delta))] = int(ls[3])
		
		
	fp.close()
	
	return mvmtdict
	

#Read the equations file and generate a python dictionary in the following format: key=(name of the kernel); value= string containing the SR model
# input trace format: <kernel-name>=<equation>
def setupequations(filename):

	fp = open(filename,'r')
	
	eqdict = dict()
	
	lines = fp.read().splitlines()
	
	for l in lines:
	
		if len(l) != 0:
    	
    		        if l[0] != '#':
    		
    			        variable=l.split("=")
    			
    			        eqdict[variable[0]] = variable[1]
    			
    			
	fp.close()
	
	return eqdict
		

		
		

	




