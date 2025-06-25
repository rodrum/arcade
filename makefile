GFORT = gfortran
# locations of external necessary repos
MSISE = ./repos/NRLMSIS2.0
HWM = ./repos/HWM14
INF = ./repos/infraGA-master
# other locations
BIN = ./bin
SRC = ./src
OUT = ./output
# other files
ARCADE_OUT = run_arcade_out.txt
EIGENSEARCH_OUT = run_eigensearch_out.txt

requisites:
	mkdir -p $(BIN)
	# Data file containing DWM parameters
	cp $(HWM)/dwm07b104i.dat $(BIN)
	# Data file containing parameters for QD coordinate conversion
	cp $(HWM)/gd2qd.dat $(BIN)
	# Binary Data file containing quiet-time HWM parameters
	cp $(HWM)/hwm123114.bin $(BIN)
	# NRMLSIS2.0 Binary data file containing model parameters
	cp $(MSISE)/msis20.parm $(BIN)
	
infraga-sph: requisites
	# Range-independent case
	cd $(INF) && make infraga-sph
	cp $(INF)/bin/infraga-sph $(BIN)

infraga-sph-rngdep: requisites
	# Range-dependent case
	cd $(INF) && make infraga-sph-rngdep
	cp $(INF)/bin/infraga-sph-rngdep $(BIN)

msis-mods: requisites
	# Necessary mod files for NRMLSIS2.0
	$(GFORT) -O3 -cpp -c $(MSISE)/msis_constants.F90 -o $(MSISE)/msis_constants.mod
	$(GFORT) -O3 -cpp -c $(MSISE)/msis_init.F90 -o $(MSISE)/msis_init.mod
	$(GFORT) -O3 -cpp -c $(MSISE)/msis_gfn.F90 -o $(MSISE)/msis_gfn.mod
	$(GFORT) -O3 -cpp -c $(MSISE)/msis_tfn.F90 -o $(MSISE)/msis_tfn.mod
	$(GFORT) -O3 -cpp -c $(MSISE)/msis_dfn.F90 -o $(MSISE)/msis_dfn.mod
	$(GFORT) -O3 -cpp -c $(MSISE)/msis_calc.F90 -o $(MSISE)/msis_calc.mod

calculate-sph-nodes: requisites msis-mods
	$(GFORT) -O3 -cpp -o $(BIN)/calculate_sph_nodes $(SRC)/calculate_sph_nodes.f90 \
	$(MSISE)/alt2gph.F90 $(MSISE)/msis_constants.mod $(MSISE)/msis_init.mod \
	$(MSISE)/msis_gfn.mod $(MSISE)/msis_tfn.mod $(MSISE)/msis_dfn.mod \
	$(MSISE)/msis_calc.mod $(MSISE)/msis_gtd8d.F90 \
	$(HWM)/hwm14.f90

complete-ecmwf: requisites msis-mods
	$(GFORT) -O3 -cpp -o $(BIN)/complete_ecmwf $(SRC)/complete_ecmwf.f90 \
	$(MSISE)/alt2gph.F90 $(MSISE)/msis_constants.mod $(MSISE)/msis_init.mod \
	$(MSISE)/msis_gfn.mod $(MSISE)/msis_tfn.mod $(MSISE)/msis_dfn.mod \
	$(MSISE)/msis_calc.mod $(MSISE)/msis_gtd8d.F90 \
	$(HWM)/hwm14.f90

calculate-rngdep-nodes: requisites msis-mods
	$(GFORT) -O3 -cpp -o $(BIN)/calculate_rngdep_nodes $(SRC)/calculate_rngdep_nodes.f90 \
	$(MSISE)/alt2gph.F90 $(MSISE)/msis_constants.mod $(MSISE)/msis_init.mod \
	$(MSISE)/msis_gfn.mod $(MSISE)/msis_tfn.mod $(MSISE)/msis_dfn.mod \
	$(MSISE)/msis_calc.mod $(MSISE)/msis_gtd8d.F90 \
	$(HWM)/hwm14.f90

calculate-rngdep-nodes-ecmwf: requisites msis-mods
	$(GFORT) -O3 -cpp -o $(BIN)/calculate_rngdep_nodes_ecmwf $(SRC)/calculate_rngdep_nodes_ecmwf.f90 \
	$(MSISE)/alt2gph.F90 $(MSISE)/msis_constants.mod $(MSISE)/msis_init.mod \
	$(MSISE)/msis_gfn.mod $(MSISE)/msis_tfn.mod $(MSISE)/msis_dfn.mod \
	$(MSISE)/msis_calc.mod $(MSISE)/msis_gtd8d.F90 \
	$(HWM)/hwm14.f90

all-prep: infraga-sph infraga-sph-rngdep calculate-sph-nodes complete-ecmwf calculate-rngdep-nodes calculate-rngdep-nodes-ecmwf
	rm *mod

#===============================================================================

run-arcade: 
	mkdir -p $(OUT)/proc
	@echo "discretize.py" > $(OUT)/proc/$(ARCADE_OUT)
	@echo "-------------" >> $(OUT)/proc/$(ARCADE_OUT)
	@{ time python $(SRC)/discretize.py ; } 2>> $(OUT)/proc/$(ARCADE_OUT)
	@echo "" >> $(OUT)/proc/$(ARCADE_OUT)

	@echo "perturb_winds.py" >> $(OUT)/proc/$(ARCADE_OUT)
	@echo "----------------" >> $(OUT)/proc/$(ARCADE_OUT)
	@{ time python $(SRC)/perturb_winds.py ; } 2>> $(OUT)/proc/$(ARCADE_OUT)
	@echo "" >> $(OUT)/proc/$(ARCADE_OUT)

	@echo "plot_profiles.py" >> $(OUT)/proc/$(ARCADE_OUT)
	@echo "----------------" >> $(OUT)/proc/$(ARCADE_OUT)
	@{ time python $(SRC)/plot_profiles.py ; } 2>> $(OUT)/proc/$(ARCADE_OUT)
	@echo "" >> $(OUT)/proc/$(ARCADE_OUT)

	@echo "ARCADE_main.py" >> $(OUT)/proc/$(ARCADE_OUT)
	@echo "--------------" >> $(OUT)/proc/$(ARCADE_OUT)
	@{ time python $(SRC)/ARCADE_main.py ; } 2>> $(OUT)/proc/$(ARCADE_OUT)


eigen_search:
	mkdir -p $(OUT)/proc
	@echo "==========================================" >> $(OUT)/proc/$(EIGENSEARCH_OUT)
	@date >> $(OUT)/proc/$(EIGENSEARCH_OUT)
	@{ time python $(SRC)/eigen_search_infraGA.py ; } 2>> $(OUT)/proc/$(EIGENSEARCH_OUT)
	@echo "" >> $(OUT)/proc/$(EIGENSEARCH_OUT)

# -----
# Tests
# -----
test_discretize:
	python $(SRC)/discretize.py
test_clean:
	rm -r $(OUT)
	rm input/secs_doys_sources_stations.txt
