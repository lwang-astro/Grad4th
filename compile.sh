#!/bin/sh
g++ -Wall -O3 -std=c++11 -fopenmp -DINTGRT=DKD_2nd   -DENERGY=energy main.cpp -o dkd
g++ -Wall -O3 -std=c++11 -fopenmp -DINTGRT=KDK_2nd   -DENERGY=energy main.cpp -o kdk
g++ -Wall -O3 -std=c++11 -fopenmp -DINTGRT=KDKDK_2nd -DENERGY=energy main.cpp -o kdkdk2
g++ -Wall -O3 -std=c++11 -fopenmp -DINTGRT=KDKDK_4th -DENERGY=energy main.cpp -o kdkdk4

g++ -Wall -O3 -std=c++11 -fopenmp -DINTGRT=DKD_2nd_changeover   -DENERGY=energy_changeover -DCHANGEOVER main.cpp -o dkdc
g++ -Wall -O3 -std=c++11 -fopenmp -DINTGRT=KDK_2nd_changeover   -DENERGY=energy_changeover -DCHANGEOVER main.cpp -o kdkc
g++ -Wall -O3 -std=c++11 -fopenmp -DINTGRT=KDKDK_2nd_changeover -DENERGY=energy_changeover -DCHANGEOVER main.cpp -o kdkdk2c
g++ -Wall -O3 -std=c++11 -fopenmp -DINTGRT=KDKDK_4th_changeover -DENERGY=energy_changeover -DCHANGEOVER main.cpp -o kdkdk4c
