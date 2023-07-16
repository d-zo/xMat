output:
	python preparation.py;

run:
	cd output/C; make;
	cd output/Fortran; make;
	cd output/Fortran-UMAT; make;
	cd output/Fortran-USER_MOD; make;
	cd output/Fortran-VUMAT; make;
	cd output/Octave; octave compile_lnx.m; octave tester.m;
	cd output/Python; make;
