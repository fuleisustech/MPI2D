include ./Makefile.config

# Acoustic WE Modeling and Migration
e2dmod:
	${RUN} ${DIR_BIN}/e2d_modeling.x par=${DIR_CURRENT}/parfile_e2dmod.sh
a2dmod:
	${RUN} ${DIR_BIN}/a2d_modeling.x par=${DIR_CURRENT}/parfile_a2dmod.sh
a2drtm :
	${RUN} ${DIR_BIN}/a2d_rtm.x par=${DIR_CURRENT}/parfile_rtm.sh
a2drtm_h :
	${RUN} ${DIR_BIN}/a2d_rtm_h.x par=${DIR_CURRENT}/parfile_rtm.sh
a2dfwi :
	${RUN} ${DIR_BIN}/fwi2d.x par=${DIR_CURRENT}/parfile_fwi.sh
a2dfwi_hil :
	${RUN} ${DIR_BIN}/fwi2d_hil.x par=${DIR_CURRENT}/parfile_fwi_hil.sh
a2dwt_Field :
	${RUN} ${DIR_BIN}/wt2d_Field.x par=${DIR_CURRENT}/parfile_wt.sh
a2dfwi_ms :
	${RUN} ${DIR_BIN}/fwi2d_ms.x par=${DIR_CURRENT}/parfile_fwi_ms.sh
a2dfwi_ms_field :
	${RUN} ${DIR_BIN}/fwi2d_ms_Volve.x par=${DIR_CURRENT}/parfile_fwi_ms.sh
a2dPhaInv :
	${RUN} ${DIR_BIN}/PhaInv2d_Field.x par=${DIR_CURRENT}/parfile_phainv.sh
a2dPhaInv_2dfft :
	${RUN} ${DIR_BIN}/PhaInv2d_Field_2dfft.x par=${DIR_CURRENT}/parfile_phainv.sh
a2dPhaInv_hil :
	${RUN} ${DIR_BIN}/PhaInv2d_Field_hil.x par=${DIR_CURRENT}/parfile_phainv_hil.sh
a2dPhaInv_window :
	${RUN} ${DIR_BIN}/PhaInv2d_window.x par=${DIR_CURRENT}/parfile_phainv.sh
clean:
	rm -f *.bin *.log *~
