CXX = g++-10
CXXFLAGS = -g -Wall -Wextra --std=c++17 -fno-finite-math-only -ffast-math -march=native -Wfatal-errors -MMD -fno-omit-frame-pointer
INCFLAGS = -I./include -I./lodepng -I/usr/local/include
LINKFLAGS = -L/usr/local/lib

GCCVERSION = $(shell $(CXX) -dumpversion)

# clang on Mac is aliased to g++.
# It also doesn't support OpenMP.
# This ugly hack enables OpenMP with my version of gcc:
ifeq "$(GCCVERSION)" "10.2.0"
	CXXFLAGS += -fopenmp
endif

ifdef DEBUG
	CXXFLAGS += -O2 -fsanitize=address -fsanitize=undefined -DDEBUG
else
	CXXFLAGS += -O3
endif

SRCS := $(wildcard src/*.cpp)
EXECS := $(patsubst src/%.cpp,%.x,$(SRCS))
DEPS :=  $(patsubst src/%.cpp,%.d,$(SRCS))
DEPS += tests.d

all: $(EXECS)

#This builds all .cpp files into separate executables:
%.x: src/%.cpp
	$(CXX) $(CXXFLAGS) $(INCFLAGS) -o $@ lodepng/lodepng.cpp $< $(LINKFLAGS)

tests.x: tests/tests.cpp
	$(CXX) $(CXXFLAGS) $(INCFLAGS) -o tests.x lodepng/lodepng.cpp $< $(LINKFLAGS) -lgtest -lgtest_main

test: tests.x
	./tests.x

-include $(DEPS)

.PHONY: clean
clean:
	rm -rf *.x *.o *.d *.x.dSYM *.png
