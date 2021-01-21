CXX = g++
CXXFLAGS = -g -Wall -Wextra --std=gnu++17 -fno-finite-math-only -march=native -Wfatal-errors -MMD -fno-omit-frame-pointer
INCFLAGS = -I./include -I./lodepng

ifdef DEBUG
	CXXFLAGS += -O2 -fsanitize=address -fsanitize=undefined
else
	CXXFLAGS += -O3
endif

SRCS := $(wildcard src/*.cpp)
EXECS := $(patsubst src/%.cpp,%.x,$(SRCS))

all: $(EXECS)

#This builds all .cpp files into separate executables:
%.x: src/%.cpp
	$(CXX) $(CXXFLAGS) $(INCFLAGS) -o $@ $< lodepng/lodepng.cpp $(LINKFLAGS)

-include $(SRCS:.cpp=.d)

.PHONY: clean
clean:
	rm -rf *.x *.o *.d *.x.dSYM
