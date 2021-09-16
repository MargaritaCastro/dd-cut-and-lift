# --- SYSTEM ---

#linux
#SYSTEM = x86-64_linux
#mac
SYSTEM = x86-64_osx
LIBFORMAT  = static_pic

# --- DIRECTORIES ---

CCC = g++ -std=c++14

BASISILOG = $(HOME)/Applications/IBM/ILOG/CPLEX_Studio129
CPOPTDIR   = $(BASISILOG)/cpoptimizer
CONCERTDIR = $(BASISILOG)/concert
CPLEXDIR   = $(BASISILOG)/cplex
BOOSTDIR   = /usr/local/Cellar/boost/1.67.0_1/include


# --- FLAGS ---

CCOPT = -fPIC -fstrict-aliasing -pedantic -Wall -fexceptions -Wno-long-long -ffloat-store -m64 -DILOUSEMT -D_REENTRANT -DILM_REENTRANT
CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CPOPTLIBDIR = $(CPOPTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
#GUROBILIB   = $(GUROBIDIR)/lib/

CCLNFLAGS = -L$(CPLEXLIBDIR) -lilocplex -lcplex -L$(CONCERTLIBDIR) -lconcert -lm -pthread
CLNFLAGS  = -L$(CPLEXLIBDIR) -lcplex -lm -pthread

# --- OPTIMIZATION FLAGS ---

##Debug
DEBUG_OPT = -g -O0 -DDEBUG

##Release
#DEBUG_OPT = -O -DNDEBUG

#DEBUG_OPT = -g3 -O0
#PROF = -pg
#PROF =

CFLAGS = -DIL_STD $(DEBUG_OPT) -I$(BOOSTDIR) -c -I$(CPLEXDIR)/include  -I$(CPOPTDIR)/include -I$(CONCERTDIR)/include  -I$(GUROBIDIR)/include/ -fPIC -pedantic -Wall -Wno-long-long -fexceptions -m64 -fno-strict-aliasing -DILOUSEMT -D_REENTRANT -DILM_REENTRANT -Wno-c99-extensions
LDFLAGS = -L$(CPOPTLIBDIR) -lcp -L$(CPLEXLIBDIR) -lcplex -L$(CONCERTLIBDIR) -lconcert -lilocplex -lpthread -framework CoreFoundation -framework IOKit -stdlib=libc++ 

# ---- COMPILE  ----
SRC_DIR   := src
OBJ_DIR   := obj

SRC_DIRS  := $(shell find $(SRC_DIR) -type d)
OBJ_DIRS  := $(addprefix $(OBJ_DIR)/,$(SRC_DIRS))

SOURCES   := $(shell find $(SRC_DIR) -name '*.cpp')
OBJ_FILES := $(addprefix $(OBJ_DIR)/, $(SOURCES:.cpp=.o))

vpath %.cpp $(SRC_DIRS)

# ---- TARGETS ----

EXECUTABLE = bdd_conic

all: $(EXECUTABLE)

$(EXECUTABLE): makedir $(SOURCES) $(OBJ_FILES) 
	$(CCC) $(OBJ_FILES) $(LDFLAGS) $(PROF) -o $@

$(OBJ_DIR)/%.o: %.cpp
	$(CCC) $(CFLAGS) $< -o $@

makedir: $(OBJ_DIRS)

$(OBJ_DIRS):
	@mkdir -p $@

clean:
	@rm -rf obj 
	@rm -rf $(EXECUTABLE)
