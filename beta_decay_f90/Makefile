FC      = gfortran
FFLAGS  = -O2 -Wall -Wextra
OBJDIR  = obj

MODS    = mod_constants.f90 mod_wave_function.f90 mod_matrix_elements.f90 \
          mod_quenching.f90 mod_shape_factor.f90 mod_phase_space.f90 mod_half_life.f90

OBJS    = $(MODS:%.f90=$(OBJDIR)/%.o)

.PHONY: all test clean

all: test

$(OBJDIR):
	mkdir -p $(OBJDIR)

$(OBJDIR)/%.o: %.f90 | $(OBJDIR)
	$(FC) $(FFLAGS) -c $< -o $@ -J$(OBJDIR)

$(OBJDIR)/test_runner.o: tests/test_runner.f90 $(OBJS) | $(OBJDIR)
	$(FC) $(FFLAGS) -c $< -o $@ -I$(OBJDIR) -J$(OBJDIR)

test: $(OBJS) $(OBJDIR)/test_runner.o
	$(FC) $(FFLAGS) -o test_runner $(OBJS) $(OBJDIR)/test_runner.o
	./test_runner

clean:
	rm -rf $(OBJDIR) test_runner
