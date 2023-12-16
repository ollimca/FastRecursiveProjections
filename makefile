# Makefile for the FRP shared library. 
# main.pricing is an example of an executable that uses the FRP library

CXX=g++ 
CXXFLAGS =  -I./include/    
CXXLFLAGS = -L. -L./lib/ 

OBJS += FRP_Pricing.o FRP_Progs.o Seq.o FRP_Merton_class.o FRP_Heston_class.o FRP_BS_class.o Matrutils.o Model.o haar.o
LIB = FRP
EXEC = main_pricing

%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c -g -O2 -Wall -o $@ $^

$(OBJS):
exec : lib
	$(CXX) -g $(CXXFLAGS) -o $(EXEC) $(EXEC).cpp -L. -l$(LIB)

lib : $(LIB) 
$(LIB): $(OBJS)
	$(CXX) -shared -o lib$(LIB).so $(OBJS) $(CXXLFLAGS) -lfftw3 -lnewmat -lrecipes 

clean : 
	rm -f $(OBJS)
realclean: clean
	rm -f $(OBJS) $(EXEC) $(LIB)

run: exec
	./$(EXEC)
	
install: all
