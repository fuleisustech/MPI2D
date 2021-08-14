include ../Makefile.config
include ../OBJS.config

all: ${OBJS} 

$(DIR_OBJ)/%.o: $(DIR_CORE)/%.f90 
	$(FCC) -o $@ $< $(F90FLAGS)

clean:
	rm -f $(DIR_OBJ)/*.o $(DIR_MOD)/*.mod

