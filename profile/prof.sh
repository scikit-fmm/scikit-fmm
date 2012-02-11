# make a temporary environment for profiling only the extension module

cp ../skfmm/fast_marching.cpp .
cp ../skfmm/fast_marching.h .
cp ../skfmm/heap.h .
cp ../skfmm/heap.cpp .


g++ -O2 prof.cpp fast_marching.cpp heap.cpp -g -pg  -o a.out
g++ -O2 prof.cpp fast_marching.cpp heap.cpp -o clean.out

time ./clean.out
time ./a.out
gprof ./a.out > output.txt
shark -i -1 -G ./a.out

#########################################################################
# notes on profiling:						        #
# profile.py corresponds to the case in prof.cpp                        #
# the python version and the c version are more or less the same speed  #
# (when compiled with the same level of optimization!)                  #
# the call to gprof here lets us to profile the c++ code	        #
# the python version only sees the cfmm module as a black box.	        #
# beware that compiling c code with the profiler flags makes it slower  #
#########################################################################

# on os x: shark -i -1 -G my_program <my_program options>

# to check for memory leaks
#g++ -O1 prof.cpp fast_marching.cpp heap.cpp -g -o a2.out
#valgrind --leak-check=yes ./a2.out

