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
# -DTCDIFFERENCE 
# -DUPWIND 
# -DEXPANSION_LIMITED
 FLAGS=-frecord-marker=4 -DVERBOSE -DTCDIFFERENCE -x f95-cpp-input  -Wall  #-ffpe-trap=zero,overflow,invalid

 OBJ = parameters.o derived_types.o \
       grid_commons.o eos.o input.o \
       utils.o mcrtfld.o gravity.o \
       pdrag.o particle.o init_grid.o \
       units.o flux.o source.o \
       velocity.o state.o init_conditions.o \
       read_hydro.o write_files.o courant.o avisc.o cleanup.o main.o

hexcake: $(OBJ)
	$(FC) $(LFLAGS) -o boxzy $(OBJ)

clean:
	rm -f *.o *.mod boxzy

%.o:%.f90
	$(FC) $(FLAGS) -c $^ -o $@
