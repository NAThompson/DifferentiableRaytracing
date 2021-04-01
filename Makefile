CXX = g++
CXXFLAGS = -g -Wall -Wextra --std=c++17 -fno-finite-math-only -ffast-math -march=native -Wfatal-errors -MMD -fno-omit-frame-pointer -fopenmp
INCFLAGS = -I./include -I./lodepng

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
