 FC=gfortran -O3 -fopenmp
# POSSIBLE FLAGS
# -DRADTRAN 
# -DTHERMALHIST 
# -DWITHDRAG 
# -DPARTICLE 
# -DUSEPERT 
# -DSELFGRAVITY 
# -DVERBOSE
# -DFASTGRAVITY 
# -DEXTRAANCHORS 
# -DSUPPRESSDRIFT  
# -DROTATE
# -DEXTERNALPHI
# -DRUN_TEST_PHI
# -DPOLYEOS
# -DSLOPE_THETA
# -DBACKREACTION_DIRECT
# -DFLUX_CYL_Y --> allows 3D sim in 2D by using the y coord as a radial coord and assuming axisymmetry about x axis.
 FLAGS=-frecord-marker=4 -DSLOPE_THETA=1.25 -DVERBOSE -x f95-cpp-input  -Wall #-ffpe-trap=zero,overflow,invalid

 OBJ = parameters.o derived_types.o \
       grid_commons.o eos.o input.o \
       utils.o mcrtfld.o gravity.o \
       pdrag.o particle.o init_grid.o \
       units.o flux.o flux_ang.o source.o \
       velocity.o state.o init_conditions.o \
       read_hydro.o write_files.o courant.o avisc.o cleanup.o main.o

hexcake: $(OBJ)
	$(FC) $(LFLAGS) -o boxzy $(OBJ)

clean:
	rm -f *.o *.mod boxzy

%.o:%.f90
	$(FC) $(FLAGS) -c $^ -o $@
