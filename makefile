CXX ?= c++
PYTHON ?= python3
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
PREFIX ?= $(shell brew --prefix 2>/dev/null || printf '%s' /usr/local)
CPLEX_ARCH ?= $(if $(filter arm64,$(shell uname -m)),arm64_osx,x86-64_osx)
CPLEX_SYSTEM_LIBS = -lpthread -framework CoreFoundation -framework IOKit -framework Accelerate
else
PREFIX ?= /usr/local
CPLEX_ARCH ?= x86-64_linux
CPLEX_SYSTEM_LIBS = -lpthread -ldl
endif
HIGHS_ROOT ?= $(PREFIX)
GSL_ROOT ?= $(PREFIX)
HIGHS_LIBDIR ?= $(HIGHS_ROOT)/lib
GSL_LIBDIR ?= $(GSL_ROOT)/lib
CPLEX_ROOT ?=
CPLEX_LIBDIR ?= $(CPLEX_ROOT)/cplex/lib/$(CPLEX_ARCH)/static_pic
CPPFLAGS += -Iinclude -I$(HIGHS_ROOT)/include -I$(HIGHS_ROOT)/include/highs -I$(GSL_ROOT)/include
CXXFLAGS ?= -O2 -g -std=c++11 -Wall -Wextra -Wno-unused-parameter -Wno-sign-compare
LDLIBS += -L$(HIGHS_LIBDIR) -Wl,-rpath,$(HIGHS_LIBDIR) -lhighs -L$(GSL_LIBDIR) -Wl,-rpath,$(GSL_LIBDIR) -lgsl -lgslcblas -lm -pthread

SOURCES = utilities DELTA_solver global_phase cpm reverse_search cpm_controller warp lp_backend input resource_scope decision_path
OBJECTS = $(addprefix build/,$(addsuffix .o,$(SOURCES)))
CORE_OBJECTS = $(filter-out build/DELTA_solver.o,$(OBJECTS))
TESTS = fixed_reduction_test reverse_search_test conservative_cut_test oracle_driver

.PHONY: all test cplex rlt-replay clean
all: bin/DELTA_solver bin/generate_dblp

build/%.o: src_cpp/%.cpp
	@mkdir -p build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -c $< -o $@

bin/DELTA_solver: $(OBJECTS)
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)

bin/generate_dblp: tools/generate_dblp.cpp
	@mkdir -p bin
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $<

bin/rlt_replay: tools/rlt_replay.cpp $(CORE_OBJECTS)
	@mkdir -p bin
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)
rlt-replay: bin/rlt_replay

bin/%: tests/%.cpp $(CORE_OBJECTS)
	@mkdir -p bin
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)

test: all $(addprefix bin/,$(TESTS))
	bin/fixed_reduction_test
	bin/reverse_search_test
	bin/conservative_cut_test
	$(PYTHON) tests/oracle_regression.py --bin bin/oracle_driver --allow-controlled-case 13
	DECIDEIT_SEQUENTIAL_DIRECTIONS=1 $(PYTHON) tests/converter_regression.py --bin bin/DELTA_solver --converter tools/convert_group_b.py
	DECIDEIT_SEQUENTIAL_DIRECTIONS=1 $(PYTHON) tests/generator_regression.py --bin bin/DELTA_solver --generator bin/generate_dblp
	$(PYTHON) tests/decision_path_test.py --bin bin/DELTA_solver
	$(PYTHON) tests/check_illustrative.py --solver bin/DELTA_solver

CPLEX_OBJECTS = build/with_cplex/DELTA_solver.o build/with_cplex/cplex_reference.o
CPLEX_CPPFLAGS = $(CPPFLAGS) -I$(CPLEX_ROOT)/cplex/include -DDECIDEIT_WITH_CPLEX
build/with_cplex/%.o: src_cpp/%.cpp
	@test -f "$(CPLEX_ROOT)/cplex/include/ilcplex/cplex.h" || { echo "Set CPLEX_ROOT to an installed CPLEX Studio directory." >&2; exit 2; }
	@mkdir -p build/with_cplex
	$(CXX) $(CPLEX_CPPFLAGS) $(CXXFLAGS) -MMD -MP -c $< -o $@

bin/DELTA_solver_cplex: $(CORE_OBJECTS) $(CPLEX_OBJECTS)
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS) -L$(CPLEX_LIBDIR) -lcplex $(CPLEX_SYSTEM_LIBS)

cplex: bin/DELTA_solver_cplex

-include $(OBJECTS:.o=.d) $(CPLEX_OBJECTS:.o=.d)
clean:
	rm -rf build bin
