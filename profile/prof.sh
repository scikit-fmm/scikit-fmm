# make a temporary environment for profiling only the extension module

cp ../skfmm/base_marcher.cpp .
cp ../skfmm/base_marcher.h .
cp ../skfmm/distance_marcher.cpp .
cp ../skfmm/distance_marcher.h .
cp ../skfmm/heap.cpp .
cp ../skfmm/heap.h .

g++ -O2 base_marcher.cpp distance_marcher.cpp heap.cpp prof.cpp -g -pg  -o a.out
g++ -O2 base_marcher.cpp distance_marcher.cpp heap.cpp prof.cpp  -o clean.out

time ./clean.out
time ./a.out
gprof ./a.out > output.txt
shark -i -1 -G ./clean.out

#########################################################################
# notes on profiling:                                                   #
# profile.py corresponds to the case in prof.cpp                        #
# the python version and the c version are more or less the same speed  #
# (when compiled with the same level of optimization!)                  #
# the call to gprof here lets us to profile the c++ code                #
# the python version only sees the cfmm module as a black box.          #
# beware that compiling c code with the profiler flags makes it slower  #
#########################################################################

# on os x: shark -i -1 -G my_program <my_program options>

# to check for memory leaks
#g++ -O1 prof.cpp fast_marching.cpp heap.cpp -g -o a2.out
#valgrind --leak-check=yes ./a2.out

