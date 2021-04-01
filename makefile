de_ePX_bc : aux_functions.o cec17_test_func.o de.o file_man.o global.o statistics.o transformation.o 
	g++ -Wall aux_functions.o cec17_test_func.o de.o file_man.o global.o statistics.o transformation.o -o de_ePX_bc

aux_functions.o : aux_functions.cpp	
	g++ -Wall -o aux_functions.o -c aux_functions.cpp

cec17_test_func.o : cec17_test_func.cpp	
	g++ -Wall -o cec17_test_func.o -c cec17_test_func.cpp

de.o : de.cpp	
	g++ -Wall -o de.o -c de.cpp

file_man.o : file_man.cpp	
	g++ -Wall -o file_man.o -c file_man.cpp

global.o : global.cpp	
	g++ -Wall -o global.o -c global.cpp

statistics.o : statistics.cpp	
	g++ -Wall -o statistics.o -c statistics.cpp

transformation.o : transformation.cpp	
	g++ -Wall -o transformation.o -c transformation.cpp



