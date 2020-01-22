ROOT = /home/chris/root
SRC = src/json.c src/main.cpp src/SeedHist2D.cpp src/SeedHist3D.cpp src/PJMCoords.cc
INC = -I/usr/include/postgresql `$(ROOT)/bin/root-config --cflags`
DEF = -DLIBUS_NO_SSL -DHAVE_INLINE
# -D_GLIBCXX_PARALLEL
LIBS = -lstdc++fs -lpq -lcurl -lcrypto -lssl -lz -l:libnuma.so.1 -lpthread `$(ROOT)/bin/root-config --libs` -lnghttp2_asio -lboost_system -lboost_iostreams
JEMALLOC = -L`jemalloc-config --libdir` -Wl,-rpath,`jemalloc-config --libdir` -l:libjemalloc.so.2 `jemalloc-config --libs`
TARGET = gaiawebql

dev:
	icpc -g -O3 -xCORE-AVX2 -mcmodel=large -qopenmp -qopenmp-simd -qopt-streaming-stores auto -funroll-loops -ipo -std=c++17 -fp-model fast -qopt-report=5 -qopt-report-phase=vec $(DEF) $(INC) $(SRC) -o $(TARGET) $(LIBS) $(JEMALLOC) 

#	how to run:
#	enable OpenMP for loop cancellation
#	OMP_CANCELLATION=true ./gaiawebql

gcc:
	g++ -march=native -g -O3 -std=c++17 -fopenmp -fopenmp-simd -funroll-loops -ftree-vectorize -Wno-register $(DEF) $(INC) $(SRC) -o $(TARGET) $(LIBS) $(JEMALLOC)

test:
	icpc -g -O3 -xCORE-AVX2 -I/usr/include/python2.7 src/PJMCoords.cc src/testAstroPy.cpp -o astropy
