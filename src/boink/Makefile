# Profile?
# Set this variable to true if you wish to profile the codes.
WANT_PROFILING=false

# Which profiling tool to use?
# Assuming you have TAU installed and setup properly,
# you can instrument codes with it to get detailed multi-threaded profiling.
# Otherwise, gprof is able to give you some information without threading info.
# Choose one of: gprof, TAU
PROFILER_OF_CHOICE=gprof

# Perform extra sanity checking?
# Set this variable to true
# if you wish the codes to perform extra sanity checking
# (to the possible detriment of performance).
WANT_EXTRA_SANITY_CHECKING=false

# Compile with debugging symbols?
# Set this variable to true
# if you wish the codes to be built with debugging symbols
# (increases code size
# and does not always produce accurate stepping in a debugger
# when optimization is turned on).
WANT_DEBUGGING=false

PREFIX=/usr/local

### NOTE: No user-serviceable parts below this line! ###

INCLUDES= -I ../../include/ 


# Warnings in common to C and C++
WARNINGS=-Wall

# Flags in common to C and C++
COMMON_FLAGS=-O3 -fPIC -fno-omit-frame-pointer

# Base C/CXXFLAGS
CPPFLAGS ?=

CXXFLAGS ?=
CXXFLAGS += $(COMMON_FLAGS) $(WARNINGS)
CXXFLAGS += -Wstrict-null-sentinel -std=c++11
CXXFLAGS += $(INCLUDES) $(CPPFLAGS)

CFLAGS	 ?=
CXXFLAGS += $(COMMON_FLAGS) $(WARNINGS)
CFLAGS   += -Wshadow -Wcast-align -Wstrict-prototypes
CFLAGS   += $(INCLUDES) $(CPPFLAGS)

LDFLAGS  ?=

ifeq ($(WANT_DEBUGGING), true)
DEBUG_FLAGS=-g
CXXFLAGS += $(DEBUG_FLAGS)
CFLAGS   += $(DEBUG_FLAGS)
endif


ifeq ($(WANT_PROFILING), true)
ifeq ($(PROFILER_OF_CHOICE), TAU)
CXX=tau_cxx.sh
endif
ifeq ($(PROFILER_OF_CHOICE), gprof)
CXXFLAGS += -pg
CFLAGS   += -pg
LDFLAGS  += -pg
endif
endif

# Place POSIX threads last in linking order, if needed.
ifneq ($(shell uname), Linux)
LDFLAGS  += -pthread
endif


HAVE_OPENMP=$(shell 						\
	      $(CXX) -fopenmp -o chkomp .check_openmp.cc 	\
	      2>/dev/null && echo true || echo false; 		\
	      rm -f chkomp)

ifeq ($(HAVE_OPENMP), true)
CXXFLAGS +=-fopenmp
CFLAGS   +=-fopenmp
endif

ifneq ($(PACKAGE_VERSION),)
VERSION = $(PACKAGE_VERSION)
else
VERSION = $(pushd ../../ && python setup.py --version && popd)
endif

# The ABI version of libboink
LIB_VERSION = 1

ifeq ($(shell uname), Darwin)
SHARED_EXT   = dylib
SONAME       = libboink.$(SHARED_EXT).$(LIB_VERSION)
SONAME_FLAGS = -install_name $(PREFIX)/lib/$(SONAME) \
	       -compatibility_version $(LIB_VERSION) \
	       -current_version $(LIB_VERSION)
else
SHARED_EXT   = so
SONAME       = libboink.$(SHARED_EXT).$(LIB_VERSION)
SONAME_FLAGS = -Wl,-soname=$(SONAME)
endif

LIBBOINKSO=$(SONAME)

#CXXFLAGS += -DVERSION=$(VERSION)

export CXX
export CFLAGS
export CXXFLAGS
export LDFLAGS
export VERSION


LIBBOINK_OBJS= \
	dbg.o \
	hashing.o \
	assembly.o

TEST_OBJS= \
	test_assembly.o \
	benchmark_boink_consume.o \
	benchmark_oxli_consume.o

PRECOMILE_OBJS ?=
PRECLEAN_TARGS ?=

HEADERS= \
	hashing.hh \
	dbg.hh \
	boink.hh \
	assembly.hh

BOINK_HEADERS = $(addprefix ../../include/boink/,$(HEADERS))

# START OF RULES #

# The all rule comes first!
all: $(LIBBOINKSO) libboink.a boink.pc

clean: $(PRECLEAN_TARGS)
	rm -f *.o *.a *.$(SHARED_EXT)* boink.pc $(TEST_PROGS)

install: $(LIBBOINKSO) libboink.a boink.pc $(BOINK_HEADERS)
	rm -rf $(PREFIX)/include/boink
	mkdir -p $(PREFIX)/lib $(PREFIX)/lib/pkgconfig $(PREFIX)/include/boink
	cp -r $(BOINK_HEADERS) 				\
		$(PREFIX)/include/boink/
	cp boink.pc $(PREFIX)/lib/pkgconfig/
	cp $(LIBBOINKSO) libboink.a $(PREFIX)/lib/
	ln -sf $(PREFIX)/lib/$(LIBBOINKSO) $(PREFIX)/lib/libboink.$(SHARED_EXT)

boink.pc: boink.pc.in
	sed -e 's,@prefix@,$(PREFIX),'  -e 's,@VERSION@,$(VERSION),' $< >$@


%.o: %.cc $(PRECOMILE_OBJS) $(BOINK_HEADERS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -c -o $@ $< -loxli


$(LIBBOINKSO): $(LIBBOINK_OBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(SONAME_FLAGS) -shared -o $@ $^
	ln -sf $(SONAME) libboink.$(SHARED_EXT)

test_assembly: $(BOINK_HEADERS) $(LIBBOINK_OBJS) test_assembly.o
	$(CXX) $(LDFLAGS) -o $@ test_assembly.o $(LIBBOINK_OBJS) -loxli -L.

benchmark_boink_consume: $(BOINK_HEADERS) $(LIBBOINK_OBJS) benchmark_boink_consume.o
	$(CXX) $(LDFLAGS) -o $@ benchmark_boink_consume.o $(LIBBOINK_OBJS) -loxli -L.

benchmark_oxli_consume: $(BOINK_HEADERS) $(LIBBOINK_OBJS) benchmark_oxli_consume.o
	$(CXX) $(LDFLAGS) -o $@ benchmark_oxli_consume.o $(LIBBOINK_OBJS) -loxli -L.

test: test_assembly 
	./$<

BM_FILE_1 = /store/biodb/genomes/fugu/Takifugu_rubripes.FUGU5.dna_rm.toplevel.fa
BM_FILE_SM = /store/biodb/genomes/sacPom.pombase.fa
BM_FILE_LG = /store/biodb/MMETSP/data/Asterionellopsis_glacialis/SRR1300451/SRR1300451_1.fastq

benchmark: benchmark_boink_consume benchmark_oxli_consume
	echo "Boink Consume:"
	time ./benchmark_boink_consume $(BM_FILE_LG)
	echo "Oxli Consume:"
	time ./benchmark_oxli_consume $(BM_FILE_LG)

libboink.a: $(LIBBOINK_OBJS)
	ar rcs $@ $^
	ranlib $@
