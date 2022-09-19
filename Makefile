.POSIX:

BIN = fit

CXX = g++
CXXFLAGS = -std=c++17 -O3 -g -Wall -pedantic -DBINARY=\"$(BIN)\" $(shell root-config --cflags)
LDFLAGS = -L$(shell root-config --libdir)
LDLIBS = $(shell root-config --libs --noldflags) -lboost_program_options

$(BIN): main.o
	$(CXX) $(LDFLAGS) $^ -o $@ $(LDLIBS)

clean:
	$(RM) *.o $(BIN)

.PHONY: clean
