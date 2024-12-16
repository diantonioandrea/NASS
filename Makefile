.PHONY: all distclean

CXXFLAGS += -Wall -std=c++23 -pedantic -I./include -march=native -Ofast

# # Verbosity.
# CXXFLAGS += -DNVERBOSE

# # Debugging.
# CXXFLAGS += -DNDEBUG

ifeq ($(shell uname),Darwin)
ifneq ($(OpenMP),) # Apple's Clang, custom OpenMP installation under $OpenMP.
CXXFLAGS += -Xclang -fopenmp
CPPFLAGS += -I$(OpenMP)/include
LDFLAGS += -L$(OpenMP)/lib
LDLIBS += -lomp
endif
else # GCC.
CXXFLAGS += -fopenmp
LDLIBS += -lgomp
endif

# Headers.
HEADERS = ./include/*.hpp

# Executables.
TESTS = $(subst src/,executables/,$(subst .cpp,.out,$(shell find src -name "Test_*.cpp")))

# Objects.
OBJECTS = $(subst src/,objects/,$(subst .cpp,.o,$(shell find src -name "NASS_*.cpp")))

# Directories.
DIRECTORIES = ./objects ./executables

# All.
all: $(DIRECTORIES) $(TESTS)
	@echo "Compiled everything!"

# Tests.
$(TESTS): executables/Test_%.out: objects/Test_%.o $(OBJECTS) 
	@echo "Linking to $@"
	@$(CXX) $(LDFLAGS) $(LDLIBS) $^ -o $@

# Objects.
$(OBJECTS): objects/%.o: src/%.cpp $(HEADERS)
	@echo "Compiling $<"
	@$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

$(subst src/,objects/,$(subst .cpp,.o,$(shell find src -name "Test_*.cpp"))): objects/%.o: src/%.cpp $(HEADERS)
	@echo "Compiling $<"
	@$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

# Directories.
$(DIRECTORIES):
	@mkdir -p $(DIRECTORIES)

# Clean.
distclean:
	@echo "Cleaning the repo."
	@$(RM) -r $(DIRECTORIES)