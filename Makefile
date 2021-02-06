CXX = g++
CXXFLAGS = -g -Wall -Wextra --std=gnu++17 -fno-finite-math-only -ffast-math -march=native -Wfatal-errors -MMD -fno-omit-frame-pointer
INCFLAGS = -I./include -I./lodepng

ifdef DEBUG
	CXXFLAGS += -O2 -fsanitize=address -fsanitize=undefined -DDEBUG
else
	CXXFLAGS += -O3
endif

SRCS := $(wildcard src/*.cpp)
EXECS := $(patsubst src/%.cpp,%.x,$(SRCS))

all: $(EXECS)

#This builds all .cpp files into separate executables:
%.x: src/%.cpp
	$(CXX) $(CXXFLAGS) $(INCFLAGS) -o $@ lodepng/lodepng.cpp $< $(LINKFLAGS)

-include $(SRCS:.cpp=.d)

test: tests/tests.cpp
	$(CXX) $(CXXFLAGS) $(INCFLAGS) -o tests.x lodepng/lodepng.cpp $< $(LINKFLAGS) -lgtest -lgtest_main
	./tests.x

.PHONY: clean
clean:
	rm -rf *.x *.o *.d *.x.dSYM *.png
