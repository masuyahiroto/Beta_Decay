FC      = gfortran
FFLAGS  = -O2 -Wall -Wextra
SRCDIR  = src
OBJDIR  = obj

MODS    = $(SRCDIR)/mod_constants.f90 $(SRCDIR)/mod_wave_function.f90 \
          $(SRCDIR)/mod_matrix_elements.f90 $(SRCDIR)/mod_quenching.f90 \
          $(SRCDIR)/mod_shape_factor.f90 $(SRCDIR)/mod_phase_space.f90 \
          $(SRCDIR)/mod_half_life.f90

OBJS    = $(MODS:$(SRCDIR)/%.f90=$(OBJDIR)/%.o)

.PHONY: all test calc clean

all: test calc

$(OBJDIR):
	mkdir -p $(OBJDIR)

$(OBJDIR)/%.o: $(SRCDIR)/%.f90 | $(OBJDIR)
	$(FC) $(FFLAGS) -c $< -o $@ -J$(OBJDIR)

$(OBJDIR)/test_runner.o: tests/test_runner.f90 $(OBJS) | $(OBJDIR)
	$(FC) $(FFLAGS) -c $< -o $@ -I$(OBJDIR) -J$(OBJDIR)

$(OBJDIR)/beta_decay_main.o: $(SRCDIR)/beta_decay_main.f90 $(OBJS) | $(OBJDIR)
	$(FC) $(FFLAGS) -c $< -o $@ -I$(OBJDIR) -J$(OBJDIR)

test: $(OBJS) $(OBJDIR)/test_runner.o
	$(FC) $(FFLAGS) -o test_runner $(OBJS) $(OBJDIR)/test_runner.o
	./test_runner

calc: $(OBJS) $(OBJDIR)/beta_decay_main.o
	$(FC) $(FFLAGS) -o beta_decay_calc $(OBJS) $(OBJDIR)/beta_decay_main.o

clean:
	rm -rf $(OBJDIR) test_runner beta_decay_calc
