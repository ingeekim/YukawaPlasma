# Project: YBIM
# Yukawa binary ionic mixture
# In-Gee Kim, July, 2015

CXX = clang++ 

CFLAGS = -O2 -g -Wall -c
LDFLAGS =
LIBS = 

SOURCES = main.cpp YukawaPlasma.cpp parseInput.cpp evolve.cpp\
	initialize.cpp verlet.cpp thermostats.cpp yukawa.cpp\
	postEvolve.cpp vel_dist.cpp savFinal.cpp rad_dist.cpp\
	calcDiffusion.cpp AutoCorrelations.cpp\
	accAutoCorr.cpp avgAutoCorr.cpp\
	readVelocities.cpp setCurrents.cpp calcAutoCorr.cpp\
	savAutoCorr.cpp chkStatFit.cpp getDarken.cpp getFick.cpp

OBJECTS = $(SOURCES:.cpp=.o)

TARGET = ybim

all: $(TARGET)
    
$(TARGET): $(OBJECTS) 
	$(CXX) $(LDFLAGS) $(OBJECTS) -o $(TARGET) $(LIBS)

.cpp.o:
	$(CXX) $(CFLAGS) -c $< 

.SUFFIXES: .cpp

clean:
	rm -f *.o
