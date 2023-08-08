g++ -c radcor_wrap.cpp `root-config --cflags --glibs` -std=c++11 -fpic
echo here 1
gfortran -lowercase -std=legacy -ffixed-line-length-0 -extend_source -c sigrad_sim.f 
echo here 2
g++ -o radsbr radcor_wrap.o sigrad_sim.o -lgfortran `root-config --cflags --glibs` -std=c++11 -fpic
