ROOT = /home/chris/root
SRC = src/json.c src/main.cpp src/PJMCoords.cc
INC = -I/usr/include/postgresql -I/home/chris/Healpix_3.50/src/cxx/optimized_gcc/include `$(ROOT)/bin/root-config --cflags`
DEF = -DLIBUS_NO_SSL -DHAVE_INLINE
# -D_GLIBCXX_PARALLEL
LIBS = -lstdc++fs -lpq -luWS -lcurl -lcrypto -lssl -lz -lnuma -lpthread -luuid -L/home/chris/Healpix_3.50/src/cxx/optimized_gcc/lib -lhealpix_cxx -lcxxsupport `$(ROOT)/bin/root-config --libs`
JEMALLOC = -L`jemalloc-config --libdir` -Wl,-rpath,`jemalloc-config --libdir` -ljemalloc `jemalloc-config --libs`
TARGET = gaiawebql

dev:
	icpc -g -O3 -xCORE-AVX2 -mcmodel large -qopenmp -qopenmp-simd -qopt-streaming-stores auto -funroll-loops -ipo -std=c++17 -fp-model fast -qopt-report=5 -qopt-report-phase=vec $(DEF) $(INC) $(SRC) -o $(TARGET) $(LIBS) $(JEMALLOC) 

#	how to run:
#	enable OpenMP for loop cancellation
#	OMP_CANCELLATION=true ./gaiawebql

gcc:
	g++ -march=native -g -O3 -std=c++17 -fopenmp -fopenmp-simd -funroll-loops -ftree-vectorize -Wno-register $(DEF) $(INC) $(SRC) -o $(TARGET) $(LIBS) $(JEMALLOC)

test:
	icpc -g -O3 -xCORE-AVX2 -I/usr/include/python2.7 src/PJMCoords.cc src/testAstroPy.cpp -o astropy