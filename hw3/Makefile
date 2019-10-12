# no assert
PYFLAG=-O
PY_MAIN=pendulum.py
STEP=0.01
TIME=300
ANGLE=6

.PHONY: run doc clean

clean-pyc:
	@find . -name '*.pyc' -delete
	@find . -name '__pycache__' -type d | xargs rm -fr

clean-figs:
	@rm -rf ./figs

format:
	@yapf -i -r . --style ../yapf_style.cfg

run: format
	python $(PYFLAG) src/$(PY_MAIN) -s $(STEP) -t $(TIME) -a $(ANGLE)

clean: clean-pyc clean-figs
	@echo "Clean all data."

doc: report/hw1_report.md
	@cd report && make mkpdf