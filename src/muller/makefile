all:
	make -f test11c4_omp.make HOST=linux-gfortran-openmp -j16
	make -f test12c4_omp.make HOST=linux-gfortran-openmp -j16

run_muller:
	run_omp.sh 15 muller_lnx64_omp config_muller.txt

run_mfie:
	run_omp.sh 15 mfie_lnx64_omp config_mfie.txt

clean:
	rm -f *.o fort.*

distclean: 
	rm -f *.o fort.*
	rm -f mfie_lnx64_omp muller_lnx64_omp
