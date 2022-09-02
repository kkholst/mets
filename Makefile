
install:
	@R -q -e "devtools::install()"

check:
	@R -q -e "devtools::check()"

doc:
	R -q -e "devtools::document()"

vignette:
	@_R_FULL_VIGNETTE_=1 R -q -e "devtools::build_vignettes()"

roxy: doc

.PHONY: roxy doc vignette install
