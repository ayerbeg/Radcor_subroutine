# Radcor_subroutine

This is a C++ wrapper to a subroutine, based on Karl Slifer nradcor.f, to radiate cross-section. The subroutine sigrad_sim.f was coded by Ryan Zielinski
and David Ruth <David.Ruth@unh.edu>. 

There are two ways to compile: 
1. just run the script compile.sh typing ,/compile.sh . If everything is ok, it will produce an executable called radsbr. Run it typing ./radsbr. That's all.
2. The second way makes use of cmake. Create another folder, i.e. build. Change to the new folder. Type cmake <folder where the sources are located>. It will configure the system. If there are no errors in the configuration, type make. If no errors, it will produce the executbale radcor_sbr. Run it typing ./radcor_sbr
