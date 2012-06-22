 FC=gfortran -O3 -fopenmp 
# POSSIBLE FLAGS
# -DRADTRAN -DTHERMALHIST -DWITHDRAG -DPARTICLE -DFREEFLOW -DUSEPERT -DSELFGRAVITY -DVERBOSE
# -DSELFGRAVITY -DFASTGRAVITY -DEXTRAANCHORS -DSUPPRESSDRIFT  -DTCDIFFERENCE -DUPWIND
 FLAGS=-frecord-marker=4 -DVERBOSE -DTCDIFFERENCE -x f95-cpp-input  -Wall 

 OBJ = parameters.o derived_types.o \
       grid_commons.o eos.o input.o \
       utils.o mcrtfld.o gravity.o \
       pdrag.o particle.o init_grid.o \
       units.o flux.o source.o \
       velocity.o state.o init_conditions.o \
       read_hydro.o courant.o avisc.o cleanup.o main.o

hexcake: $(OBJ)
	$(FC) $(LFLAGS) -o boxzy $(OBJ)

clean:
	rm -f *.o *.mod boxzy

%.o:%.f90
	$(FC) $(FLAGS) -c $^ -o $@


