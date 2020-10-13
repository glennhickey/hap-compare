LIB_DIR:=lib
INC_DIR:=include
CWD:=$(shell pwd)
CXX ?= g++

# Multithreading with OpenMP.
PARALLEL_FLAGS=-fopenmp -pthread
LIBS=-L$(LIB_DIR) -lsdsl -ldivsufsort -ldivsufsort64

ifeq ($(shell uname -s),Darwin)
    # Our compiler might be clang that lacks -fopenmp support.
    # Sniff that
    ifeq ($(strip $(shell $(MY_CXX) -fopenmp /dev/null -o/dev/null 2>&1 | grep fopenmp | wc -l)), 1)
        # The compiler complained about fopenmp instead of its nonsense input file.
        # We need to use the hard way of getting OpenMP not bundled with the compiler.

        # The compiler only needs to do the preprocessing
        PARALLEL_FLAGS = -Xpreprocessor -fopenmp -pthread

        ifeq ($(shell if [ -d /opt/local/lib/libomp ];then echo 1;else echo 0;fi), 1)
            # Use /opt/local/lib/libomp if present, because Macports installs libomp there.
            # Brew is supposed to put it somewhere the compiler can find it by default.
            LIBS += -L/opt/local/lib/libomp
            # And we need to find the includes. Homebrew puts them in the normal place
            # but Macports hides them in "libomp"
            PARALLEL_FLAGS += -I/opt/local/include/libomp
        endif

        # We also need to link it
        LIBS += -lomp

    endif
endif

CXXFLAGS := -O0 -Werror=return-type -std=c++14 -ggdb -g -MMD -MP $(PARALLEL_FLAGS) $(CXXFLAGS)

XG_DIR = $(CWD)/deps/xg
LIBVGIO_DIR = $(CWD)/deps/libvgio

LIB_DEPS = $(LIB_DIR)/libxg.a $(LIB_DIR)/libvgio.a

LIB_FLAGS = -lxg -lvgio -lsdsl -lhandlegraph -ldivsufsort -ldivsufsort64 -latomic -lhts -lprotoc -lprotobuf -ljansson -L$(CWD)/$(LIB_DIR) 
INC_FLAGS = -I$(CWD)/$(INC_DIR)

all : bin/pg-pathcomp venv/bin/activate

$(LIB_DIR)/libxg.a: $(XG_DIR)/src/*.hpp $(XG_DIR)/src/*.cpp
	+mkdir -p include lib
	+cd $(XG_DIR) && cmake -H. -Bbuild && cmake --build build -- -j4
	+cp -r $(XG_DIR)/src/*.hpp $(CWD)/$(INC_DIR)
	+cp $(XG_DIR)/lib/libxg.a $(LIB_DIR)
   #todo use cmake properly!
	+cp $(XG_DIR)/build/*prefix/lib/lib*.a $(LIB_DIR)
	+cp -r $(XG_DIR)/build/*prefix/include/* $(INC_DIR)
	+cp $(XG_DIR)/build/*prefix/src/*build/lib/lib*.a $(LIB_DIR)
	+cp -r $(XG_DIR)/build/*prefix/src/*build/include/* $(INC_DIR)
	+cp $(XG_DIR)/build/sdsl-lite-prefix/src/sdsl-lite-build/external/*/lib/*.a $(LIB_DIR)
	+cp deps/xg/deps/*/src/*.hpp deps/xg/deps/args/*.hxx $(INC_DIR)

$(LIB_DIR)/libvgio.a: $(LIBVGIO_DIR)/src/*.cpp
	+mkdir -p include lib
	+cd $(LIBVGIO_DIR) && cmake -H. -Bbuild && cmake --build build -- -j4
	+cp -r $(LIBVGIO_DIR)/include/* $(CWD)/$(INC_DIR)
	+cp -r $(LIBVGIO_DIR)/build/libvgio.a $(CWD)/$(LIB_DIR)

pg-pathcomp.o:$(LIB_DEPS) pg-pathcomp.hpp pg-pathcomp.cpp
	$(CXX) $(INCLUDE_FLAGS) $(CXXFLAGS) $(CPPFLAGS) -c pg-pathcomp.cpp $(INC_FLAGS)

main.o:$(LIB_DEPS) main.cpp pg-pathcomp.hpp
	$(CXX) $(INCLUDE_FLAGS) $(CXXFLAGS) $(CPPFLAGS) -c main.cpp $(INC_FLAGS)

bin/pg-pathcomp:$(LIB_DEPS) pg-pathcomp.o main.o
	$(CXX) $(INCLUDE_FLAGS) $(CXXFLAGS) $(CPPFLAGS) -o bin/pg-pathcomp main.o pg-pathcomp.o $(LIB_FLAGS)

venv/bin/activate:
	virtualenv -p python2 venv
	. venv/bin/activate && pip install numpy scipy matplotlib

clean:
	rm -rf pg-pathcomp pg-pathcomp.o main.o $(LIB_DIR)/* $(INC_DIR)/*
