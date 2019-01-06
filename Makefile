OUTPUT_CK=orx_ola_fk_replacement.bc
PYTHON_SCRIPT=make_ola_high_low_ck.py
INPUT_FILES=\
naif0012.tls \
ORX_SCLKSCET.00039.tsc \
orx_v11.tf \
orx_v12.tf \
orx_v13_draft.tf \
#

test: $(OUTPUT_CK)
	python make_ola_high_low_ck.py --test

$(OUTPUT_CK): $(PYTHON_SCRIPT) $(INPUT_FILES) Makefile
	$(RM) $(OUTPUT_CK)
	python make_ola_high_low_ck.py --create

clean:
	$(RM) $(OUTPUT_CK)
