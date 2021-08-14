include ../Makefile.config
include ../OBJS.config


a2ds: a2d_modeling.x a2d_rtm.x  a2d_rtm_h.x  PhaInv2d_Field.x  PhaInv2d_Field_2dfft.x fwi2d_ms_Volve.x 
inv2d: PhaInv2d_Field_2dfft.x
inv1d: PhaInv2d_Field.x


# a2ds -- acoustic 2d wave equation modeling and migration
a2d_modeling.x: $(OBJS) a2d_modeling.f90
	$(FCC) -o $@ $^ $(LDFLAGS) -I$(DIR_MOD)
	mv a2d_modeling.x $(DIR_BIN)/
e2d_modeling.x: $(OBJS) e2d_modeling.f90
	$(FCC) -o $@ $^ $(LDFLAGS) -I$(DIR_MOD)
	mv e2d_modeling.x $(DIR_BIN)/
a2d_rtm.x: $(OBJS) a2d_rtm.f90
	$(FCC) -o $@ $^ $(LDFLAGS) -I$(DIR_MOD)
	mv a2d_rtm.x $(DIR_BIN)/
a2d_rtm_h.x: $(OBJS) a2d_rtm_h.f90
	$(FCC) -o $@ $^ $(LDFLAGS) -I$(DIR_MOD)
	mv a2d_rtm_h.x $(DIR_BIN)/
wt2d_Field.x: $(OBJS) wt2d_Field.f90
	$(FCC) -o $@ $^ $(LDFLAGS) -I$(DIR_MOD)
	mv wt2d_Field.x $(DIR_BIN)/

PhaInv2d_Field.x: $(OBJS) PhaInv2d_Field.f90
	$(FCC) -o $@ $^ $(LDFLAGS) -I$(DIR_MOD)
	mv PhaInv2d_Field.x $(DIR_BIN)/
PhaInv2d_Field_2dfft.x: $(OBJS) PhaInv2d_Field_2dfft.f90
	$(FCC) -o $@ $^ $(LDFLAGS) -I$(DIR_MOD)
	mv PhaInv2d_Field_2dfft.x $(DIR_BIN)/
fwi2d_ms_Volve.x: $(OBJS) fwi2d_ms_Volve.f90
	$(FCC) -o $@ $^ $(LDFLAGS) -I$(DIR_MOD)
	mv fwi2d_ms_Volve.x $(DIR_BIN)/

testfft.x: $(OBJS) testfft.f90
	$(FCC) -o $@ $^ $(LDFLAGS) -I$(DIR_MOD)
	mv testfft.x $(DIR_BIN)/
clean:
	rm -f $(DIR_BIN)/*.x 
