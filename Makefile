default: run

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

.PHONY: c check roxy doc v vignette install

.PHONY: export
export:
	@mkdir -p tmp/mets
	@git archive HEAD | (cd tmp/mets; tar x)

DOCKER?=docker
IMG?=mets

.PHONY: dbuild
build:
	$(DOCKER) build . -f Dockerfile --tag $(IMG)

.PHONY: drun
drun:
	$(DOCKER) run --cap-add SYS_PTRACE --priviliged -ti --rm -v $(PWD)/tmp:/data $(IMG) bash
