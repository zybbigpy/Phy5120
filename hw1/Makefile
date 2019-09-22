# no assert
PYFLAG=#-O
PY_MAIN=pendulum.py

.PHONY: clean run doc

clean-pyc:
	@find . -name '*.pyc' -delete
	@find . -name '__pycache__' -type d | xargs rm -fr

clean-figs:
	@rm -rf ./figs

format:
	@yapf -i -r . --style ../yapf_style.cfg

run: format clean-figs
	python $(PYFLAG) src/$(PY_MAIN)

clean: clean-pyc clean-figs
	@echo "Clean all data."

doc: report/hw1_report.md
	cd report && make mkpdf