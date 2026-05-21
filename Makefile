.PHONY: help build validate clean distclean test inspect adaptive-test

help:
	@echo "ngsRelate-fast — make targets"
	@echo ""
	@echo "  make build      Clone upstream, apply patch, compile both binaries"
	@echo "  make validate   Run validation pipeline (needs BEAGLE/FREQS/N/SAMPLES)"
	@echo "                  Example:"
	@echo "                    make validate BEAGLE=test.beagle.gz FREQS=test.freqs N=226 SAMPLES=samples.txt"
	@echo "  make inspect    Quick contract dump for a finished run"
	@echo "                  Example:"
	@echo "                    make inspect RES=/path/to/relatedness.res"
	@echo "  make adaptive-test  Run the adaptive scheduler pytest suite"
	@echo "  make clean      Remove build artifacts (keep bin/)"
	@echo "  make distclean  Remove everything regenerable (bin/ + build/)"

build:
	bash build.sh

validate:
	@if [ -z "$(BEAGLE)" ] || [ -z "$(FREQS)" ] || [ -z "$(N)" ] || [ -z "$(SAMPLES)" ]; then \
		echo "ERROR: need BEAGLE=, FREQS=, N=, SAMPLES="; \
		exit 1; \
	fi
	bash validate/run_validation.sh "$(BEAGLE)" "$(FREQS)" "$(N)" "$(SAMPLES)"

inspect:
	@if [ -z "$(RES)" ]; then \
		echo "ERROR: need RES=/path/to/relatedness.res"; \
		exit 1; \
	fi
	python3 scripts/contract_io.py "$(RES)"

clean:
	rm -rf build/NgsRelate-upstream/*.o build/NgsRelate-fast/*.o

distclean:
	rm -rf build/ bin/ validate/output/

adaptive-test:
	python3 -m pytest adaptive/tests/ -q
