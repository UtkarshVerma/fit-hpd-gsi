.POSIX:

CXX = g++
CXXFLAGS = -O3 -g -Wall -pedantic $(shell root-config --cflags) -std=c++17
LDFLAGS = -L$(shell root-config --libdir)
LDLIBS = $(shell root-config --libs --noldflags)
BIN = fit

$(BIN): main.o
	$(CXX) $(LDFLAGS) $^ -o $@ $(LDLIBS)

clean:
	$(RM) *.o $(BIN)

.PHONY: clean
