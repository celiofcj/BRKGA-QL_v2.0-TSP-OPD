UNAME_S := $(shell uname -s)
CCC = g++
TARGETDIR=.
CCFLAGS = -O3 -Wall -std=c++11
ifeq ($(UNAME_S),Linux)
    CCFLAGS += -fopenmp
endif
ifeq ($(UNAME_S),Darwin)
    CCFLAGS += -openmp
endif

debug: CCFLAGS += -g
release: CCFLAGS += -O3 -march=native

debug: runTest
release: runTest
 
OBJS2 = \
		$(TARGETDIR)/MCMF.o \
        $(TARGETDIR)/BRKGA_QL.o \
        $(TARGETDIR)/Problem.o  \
        $(TARGETDIR)/LazyPropagationSegmentTree.o  \
        $(TARGETDIR)/Searcher.o \

$(TARGETDIR)/runTest: $(OBJS2)
	$(CCC) $(CCFLAGS) -o $(TARGETDIR)/runTest $(OBJS2)

$(TARGETDIR)/Problem.o: MCMF.h MCMF.cpp
	$(CCC) $(CCFLAGS) -c MCMF.cpp -o $(TARGETDIR)/MCMF.o

$(TARGETDIR)/Problem.o: Searcher.h Searcher.cpp
	$(CCC) $(CCFLAGS) -c Searcher.cpp -o $(TARGETDIR)/Searcher.o

$(TARGETDIR)/Problem.o: LazyPropagationSegmentTree.h LazyPropagationSegmentTree.cpp
	$(CCC) $(CCFLAGS) -c LazyPropagationSegmentTree.cpp -o $(TARGETDIR)/LazyPropagationSegmentTree.o

$(TARGETDIR)/Problem.o: Problem.h Problem.cpp
	$(CCC) $(CCFLAGS) -c Problem.cpp -o $(TARGETDIR)/Problem.o

$(TARGETDIR)/BRKGA_QL.o: BRKGA_QL.h BRKGA_QL.cpp
	$(CCC) $(CCFLAGS) -c BRKGA_QL.cpp -o $(TARGETDIR)/BRKGA_QL.o

test: runTest
	./runTest

clean:
	$(RM) \
	$(TARGETDIR)/BRKGA_QL.o \
	$(TARGETDIR)/Problem.o  \
	$(TARGETDIR)/LazyPropagationSegmentTree.o
	$(TARGETDIR)/Searcher.o
    
rebuild:
	$(MAKE) clean
	$(MAKE) debug