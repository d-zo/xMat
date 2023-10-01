TEMP_DIR = temp_work_dir
OUTPUT_DIR = examples
CUR_DIR = ${PWD}
TOOL_DIR = tools
TEST_DIR = tests/unit-tests
SOURCE_FILES := $(shell find ./src -name "*.f")

GNU_EXTRA_FLAGS = -O2 -flto -g -fbacktrace -fcheck=all -Wall -Wextra -Wpedantic -Wtrampolines -Wshadow -Wconversion -Wvector-operation-performance -Wstack-protector

REPRESENTATION = $(shell cat src/Math_Operations/_module.f | awk '{print $$3;}' | tr [:lower:] [:upper:])



.PHONY: default prepare all clean test examples xmat_plaxischeck xmat doc pdf doxygen

default: xmat examples

all: xmat test $(OUTPUT_DIR)/xmat_gfortran $(OUTPUT_DIR)/xmat_ifort xmat_plaxischeck doc

clean:
	-rm $(OUTPUT_DIR)/*;
	-rm xmat.f;

test: $(TEST_DIR)/xmat_fruit

examples: $(OUTPUT_DIR)/xmat_gfortran $(OUTPUT_DIR)/xmat_ifort

xmat_plaxischeck: $(OUTPUT_DIR)/plaxis_checks


xmat: xmat.f

xmat.f: $(SOURCE_FILES)
	./$(TOOL_DIR)/create_symlinks.sh;
	python $(TOOL_DIR)/unify_xmat.py;
	python $(TOOL_DIR)/adjust_xmat.py "xmat.f";

prepare:
	python $(TOOL_DIR)/component_selection.py;

doc: pdf doxygen

doxygen: html/index.html

pdf: xmat.pdf

xmat.pdf: xmat.f
	-mkdir $(TEMP_DIR);
	cp $^ $(TEMP_DIR)/;
	cp $(TOOL_DIR)/code_printer.tex $(TEMP_DIR)/$(^:.f=.tex);
	cd $(TEMP_DIR); \
	pdflatex $(^:.f=.tex); \
	pdflatex $(^:.f=.tex);
	mv $(TEMP_DIR)/$@ $(CUR_DIR)/;
	-rm -r $(TEMP_DIR);

# Setting -diag-disable 5268 prevents messages of overlong source file lines due to some LaTeX formulae in comments
$(OUTPUT_DIR)/xmat_ifort: xmat.f
	-mkdir $(TEMP_DIR);
	cp $^ $(OUTPUT_DIR)/xmat_tester.f $(TEMP_DIR)/;
	cd $(TEMP_DIR); LANG=C; ifort -free -fpp -stand f03 -diag-disable 5268 -o $(notdir $@) $^ xmat_tester.f;
	mv $(TEMP_DIR)/$(notdir $@) $@;
	-rm -r $(TEMP_DIR);

$(OUTPUT_DIR)/xmat_gfortran: xmat.f
	-mkdir $(TEMP_DIR);
	cp $^ $(OUTPUT_DIR)/xmat_tester.f $(TEMP_DIR)/;
	cd $(TEMP_DIR); LANG=C; gfortran -ffree-form -cpp -std=f2003 $(GNU_EXTRA_FLAGS) -o $(notdir $@) $^ xmat_tester.f;
	mv $(TEMP_DIR)/$(notdir $@) $@;
	-rm -r $(TEMP_DIR);

$(OUTPUT_DIR)/plaxis_checks: xmat.f
	-mkdir $(TEMP_DIR);
	cp $^ $(OUTPUT_DIR)/plaxis_checks.f $(TEMP_DIR)/;
	cd $(TEMP_DIR); LANG=C; gfortran -ffree-form -cpp -DPLAXIS_DLL -std=f2003 $(GNU_EXTRA_FLAGS) -o $(notdir $@) $^ plaxis_checks.f;
	mv $(TEMP_DIR)/$(notdir $@) $@;
	-rm -r $(TEMP_DIR);

$(TEST_DIR)/xmat_fruit: $(TEST_DIR)/xmat_mod.f $(TEST_DIR)/fruit.f90 $(TEST_DIR)/general_settings_test.f $(TEST_DIR)/debug_test.f $(TEST_DIR)/math_operations_test.f $(TEST_DIR)/xmat_fruit.f
	cd $(TEST_DIR); \
	cat $(notdir $^) > allcode.f;
	python $(TOOL_DIR)/adjust_xmat.py $(TEST_DIR)/allcode.f;
	cd $(TEST_DIR); \
	gfortran -ffree-form -cpp -std=f2003 -DREPR_$(REPRESENTATION) -o $(notdir $@) allcode.f;
	-rm $(TEST_DIR)/*.mod $(TEST_DIR)/allcode.f;

# For more thorough testing, make all functions and subroutines public in the local xmat copy
# i.e. remove default "private" restriction
$(TEST_DIR)/xmat_mod.f: xmat.f
	cp $^ $@;
	sed -i 's/^   private/!   private/g' $@;

html/index.html: tools/doxy_config xmat.f
	-mkdir $(TEMP_DIR);
	cp $^ $(TEMP_DIR)/;
	cd $(TEMP_DIR); touch fintrf.h; gfortran -ffree-form -cpp -DPLAXIS_DLL -DMATLAB_CALLING -DABQ_STD_CALLING -DABQ_EXP_CALLING -DCINTER -E xmat.f > xmat_mod.f; doxygen doxy_config;
	-rm -r html;
	cp -r $(TEMP_DIR)/html .;
	-rm -r $(TEMP_DIR);
