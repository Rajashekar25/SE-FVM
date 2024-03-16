# Makefile
CXX = g++
CXXFLAGS = -O3

objects = main.o pressure.o uvnew.o uvquick.o uvupwind.o

lidDrivenCavity: $(objects)
	$(CXX) -o lidDrivenCavity $(objects)

main.o: main.cpp uvupwind.h uvquick.h pressure.h uvnew.h
	$(CXX) $(CXXFLAGS) -c main.cpp

uvupwind.o: uvupwind.cpp 
	$(CXX) $(CXXFLAGS) -c uvupwind.cpp

uvquick.o: uvquick.cpp 
	$(CXX) $(CXXFLAGS) -c uvquick.cpp

pressure.o: pressure.cpp 
	$(CXX) $(CXXFLAGS) -c pressure.cpp

uvnew.o: uvnew.cpp 
	$(CXX) $(CXXFLAGS) -c uvnew.cpp

clean:
	rm lidDrivenCavity $(objects)




