.PHONY: all distclean

ifeq ($(shell uname),Darwin) # Forces GCC (g++-14, homebrew) instead of Clang.
CXX = g++-14
endif

CXXFLAGS += -Wall -std=c++23 -pedantic -I./include -march=native -Ofast -fopenmp
LDLIBS += -lgomp

# # Verbosity.
# CXXFLAGS += -DNVERBOSE

# # Debugging.
# CXXFLAGS += -DNDEBUG

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