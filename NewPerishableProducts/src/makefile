#INCCPLEX= -I/opt/ibm/ILOG/CPLEX_Studio124/cplex/include/
#LIBCPLEX= -L/opt/ibm/ILOG/CPLEX_Studio124/cplex/lib/x86_sles10_4.1/static_pic/ -lcplex -m32 -lm -pthread
#CPLEXBIN=opt/ibm/ILOG/CPLEX_Studio124/cplex/bin/x86_sles10_4.1/

#folder include
INCCPLEX= -I/opt/ibm/ILOG/CPLEX_Studio1251/cplex/include/
#folder lib
LIBCPLEX= -L/opt/ibm/ILOG/CPLEX_Studio1251/cplex/lib/x86-64_sles10_4.1/static_pic/ -m64 -lcplex -lm -lpthread

#INCCPLEX= -I/home/lab/santosma/cplex/cplex/include
#LIBCPLEX= -L/home/lab/santosma/cplex/cplex/lib/x86-64_sles10_4.1/static_pic/ -lcplex -lm -pthread
#CPLEXBIN=opt/ibm/ILOG/CPLEX_Studio124/cplex/bin/x86_sles10_4.1/

comp: adv lbs main FORCE

adv: FORCE
	g++ -c -Wextra -ansi -O3 -funroll-all-loops $(INCCPLEX) adversarial/Master.cpp -o ../lib/master.o
	g++ -c -Wextra -ansi -O3 -funroll-all-loops $(INCCPLEX) adversarial/Adversarial.cpp -o ../lib/adv.o
	g++ -c -Wextra -ansi -O3 -funroll-all-loops $(INCCPLEX) adversarial/MILP.cpp -o ../lib/mip.o
	g++ -c -Wextra -ansi -O3 -funroll-all-loops $(INCCPLEX) adversarial/HeurMILP1.cpp -o ../lib/mip1.o
	g++ -c -Wextra -ansi -O3 -funroll-all-loops $(INCCPLEX) adversarial/HeurMILP2.cpp -o ../lib/mip2.o
	g++ -c -Wextra -ansi -O3 -funroll-all-loops $(INCCPLEX) adversarial/DPA.cpp -o ../lib/dpa.o

lbs: FORCE
	g++ -c -Wextra -ansi -O3 -funroll-all-loops $(INCCPLEX) LBs/AffDecRules.cpp -o ../lib/aff.o	
	
main: FORCE
	g++ -c -Wextra -ansi -O3 -funroll-all-loops main/Timer.cpp -o ../lib/time.o
	g++ -c -Wextra -ansi -O3 -funroll-all-loops $(INCCPLEX) main/Solver.cpp -o ../lib/solver.o
	g++ -Wextra -ansi -O3 -funroll-all-loops $(INCCPLEX) main/instMaker.cpp -o ../bin/INSTANCEMAKER
	g++ -Wextra -ansi -O3 -funroll-all-loops $(INCCPLEX) ../lib/time.o ../lib/adv.o ../lib/mip.o ../lib/mip1.o ../lib/mip2.o ../lib/dpa.o ../lib/master.o ../lib/solver.o main/main.cpp -o ../bin/main $(LIBCPLEX)
	g++ -Wextra -ansi -O3 -funroll-all-loops $(INCCPLEX) ../lib/time.o ../lib/aff.o main/mainAff.cpp -o ../bin/mainaff $(LIBCPLEX)
FORCE:
