 FC=gfortran -O3 -g -fopenmp
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
 FLAGS=-frecord-marker=4 -DTHERMALHIST -DWITHDRAG -DPARTICLE -DPOLYEOS  -DSLOPE_THETA=1.25 -DVERBOSE -x f95-cpp-input  -Wall -ffpe-trap=zero,overflow,invalid
# FLAGS=-frecord-marker=4 -DPOLYEOS  -DSLOPE_THETA=1.25 -DVERBOSE -x f95-cpp-input  -Wall #-DPOLYEOS #-ffpe-trap=zero,overflow,invalid

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
