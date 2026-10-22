SHELL := /bin/bash
RUN := bash scripts/runtime.sh
CORES ?= 4

.PHONY: ingest build report validate view place setup
setup:
	uv venv .venv --python 3.12
	uv pip install --python .venv/bin/python -r requirements.lock

ingest:
	$(RUN) python -m shav.ingest $(if $(REFRESH),--refresh,)

build:
	$(RUN) snakemake --snakefile phylogenetic/Snakefile --cores $(CORES)

report:
	$(RUN) python -m shav.report

validate:
	$(RUN) python -m pytest
	$(RUN) mypy --cache-dir ~/.cache/shamonda-virus-phylogenomics/mypy shav
	$(RUN) python -m shav.validate

view:
	PORT=4000 HOST=localhost $(RUN) auspice view --datasetDir auspice

place:
	@test -n "$(SEQUENCES)" -a -n "$(METADATA)" || (echo 'Provide SEQUENCES=case.fasta METADATA=case.tsv' >&2; exit 2)
	$(RUN) python -m shav.place --sequences "$(SEQUENCES)" --metadata "$(METADATA)"
