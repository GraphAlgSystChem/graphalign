CXX = g++
#O2: Optimisation
#Wall, Wextra: Address, comments, formatting, bool-comparison, type-limits, uninitialised...
CXXFLAGS = -std=c++20 -g -O2

LIBS =  # insert libraries here
OBJS = main.o
EXE = align

$(EXE): $(OBJS)
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS)

# $@ refers to target
# $^ list of prerequisites (that must be up-to-date with the lates compilation (or exist at all))
# implicit rules are created for the prerequisite objects in OBJS as they should all be o files.

.PHONY : clean 
clean :
	rm -f $(EXE) *.o a.out
