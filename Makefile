DOCKER=nerdctl

install:
	@R -q -e "devtools::install()"

c:
	@R -q -e "devtools::check(vignettes=FALSE)"

check:
	@R -q -e "devtools::check(run_dont_test=TRUE)"

doc:
	R -q -e "devtools::document()"

test:
	R -q -e "devtools::test()"

vignette:
	@_R_FULL_VIGNETTE_=1 R -q -e "devtools::build_vignettes(clean=FALSE, quiet=FALSE)"

v:
	@R -q -e "devtools::build_vignettes(clean=FALSE, quiet=FALSE)"

roxy: doc

IMG=rocker/r-devel-ubsan-clang
dpull:
	@$(DOCKER) pull $(IMG)
## run docker, install packages, and docker commit image (e.g. -> mets)
d:
	$(DOCKER) run --cap-add SYS_PTRACE -ti --rm -v $(PWD)/../test:/data $(IMG) bash

.PHONY: c check roxy doc v vignette install dpull d
