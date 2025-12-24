# Makefile for PMM (Parametric Model Manager) emulator programs

FC = gfortran
FFLAGS = -fPIC -fopenmp -O3 -ffree-line-length-0

# Include paths
INCS = -I.. -I../../general_modules -I../../mesh_modules -I../../pw_modules -I../../pot_modules

# Library paths
LIBS = -L.. -lcookie \
       -L../../general_modules -lgeneralmodules \
       -L../../pw_modules -lpwmodules \
       -L../../mesh_modules -lmeshmodules \
       -L../../pot_modules -lpotmodules \
       -L/opt/homebrew/opt/openblas/lib -lopenblas \
       -L/opt/homebrew/lib -lumfpack -lamd -lcholmod -lcolamd -lccolamd -lcamd -lsuitesparseconfig \
       -lblas -llapack

# Module objects
MODS = emulator_io.o param_mapping.o

# Targets
all: test_emulator_dNi emulator_train emulator_predict

# Module compilation
emulator_io.o: emulator_io.F90
	$(FC) $(FFLAGS) $(INCS) -c $<

param_mapping.o: param_mapping.F90
	$(FC) $(FFLAGS) $(INCS) -c $<

# Test program
test_emulator_dNi: test_emulator_dNi.F90
	$(FC) $(FFLAGS) $(INCS) -o $@ $< $(LIBS)

# Training program
emulator_train: emulator_train.F90 $(MODS)
	$(FC) $(FFLAGS) $(INCS) -o $@ $< $(MODS) $(LIBS)

# Prediction program
emulator_predict: emulator_predict.F90 $(MODS)
	$(FC) $(FFLAGS) $(INCS) -o $@ $< $(MODS) $(LIBS)

clean:
	rm -f test_emulator_dNi emulator_train emulator_predict *.o *.mod

.PHONY: all clean
