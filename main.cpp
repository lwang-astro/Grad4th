#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <algorithm>

#include "Nbody.h"

#ifdef CHANGEOVER

#ifndef INTGRT
#define INTGRT KDKDK_4th_changeover
#endif

#ifndef ENERGY
#define ENERGY energy_changeover
#endif

#else

#ifndef INTGRT
#define INTGRT KDKDK_4th
#endif

#ifndef ENERGY
#define ENERGY energy
#endif

#endif

int main(int argc, char **argv){
	Nbody sys;

	assert(argc > 2);
	FILE *fp = fopen(argv[1], "r");
	assert(fp);
	sys.read(fp);
	fclose(fp);

	sys.set_eps(1./64.);
    sys.changeover.setR(1./64., 1./8.);

	double en0 = sys.ENERGY(stderr);

#ifdef CHANGEOVER
    sys.calc_acc_changeover();
#else
	sys.calc_acc();
#endif

	const int invtick = atoi(argv[2]);
	const double tick = 1.0 / invtick;
	const int nloop = invtick / 64;


	double err_max = 0.0;

	for(int i=0; i<nloop; i++){
		sys.INTGRT(tick);
		// sys.DKD_2nd(tick);
		// sys.KDK_2nd(tick);
		// sys.KDKDK_2nd(tick);
		// sys.KDKDK_4th(tick);
		double en1 = sys.ENERGY();
		double de = (en1 - en0) / en0;
		printf("%e %e\n", sys.tsys, de);

		err_max = std::max(err_max, fabs(de));
	}
	printf("%e %e (dt,error_max)\n", 12.*tick, err_max);

	return 0;
}


