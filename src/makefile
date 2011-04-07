# Possible compilers: g++, icpc
CXX    = icpc

#
# DONT CHANGE ANYTHING BELOW
# UNLESS YOU KNOW WHAT YOU 
# ARE DOING !!!
#

# GNU compiler
ifeq ($(CXX),g++)
	CFLAGS = -O3 -Wall
endif

# LLVM compiler
ifeq ($(CXX),llvm-g++)
	CFLAGS = -O3 -Wall
endif

# Intel compiler
ifeq ($(CXX),icpc)
	CFLAGS = -fast
endif

HDR = $(wildcard *.h)
SRC = $(wildcard *.cc)
OBJ = $(patsubst %.cc,%.o,$(SRC))

TARGET = flo3d

all: $(TARGET)

flo3d: $(OBJ)
	$(CXX) -o $(TARGET) $(OBJ)

$(OBJ): $(HDR)

.cc.o: $(HDR)
	$(CXX) $(CFLAGS) -c $*.cc

clean:
	rm -f *.o $(TARGET)

allclean:
	make clean
	rm -f *.vtk
