R ?= /usr/bin/env R

default: run

install:
	@$(R) -q -e "devtools::install()"

c:
	@$(R) -q -e "devtools::check(vignettes=FALSE)"

check:
	@$(R) -q -e "devtools::check(run_dont_test=TRUE)"

doc:
	@$(R) -q -e "devtools::document()"

examples:
	@$(R) -q -e "devtools::run_examples()"

test:
	@$(R) -q -e 'library("mets"); tinytest::run_test_dir("inst/tinytest")'

test-loadall:
	@$(R) -q -e 'devtools::load_all("."); tinytest::test_all(".")'

slowtest:
	@$(R) -q -e 'library("mets"); tinytest::run_test_dir("inst/slowtest")'

testall: test slowtest

vignette:
	@_R_FULL_VIGNETTE_=1 $(R) -q -e "devtools::build_vignettes(clean=FALSE, quiet=FALSE)"

init: doc
	@$(R) -q -e "Rcpp::compileAttributes()"

v:
	@$(R) -q -e "devtools::build_vignettes(clean=FALSE, quiet=FALSE)"

roxy: doc

.PHONY: check doc test

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
