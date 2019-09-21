.PHONY: clean run doc

clean-pyc:
	@find . -name '*.pyc' -delete
	@find . -name '__pycache__' -type d | xargs rm -fr

format:
	@yapf -i -r . --style ../yapf_style.cfg
run: format
	python pendulum.py 

clean: clean-pyc
	@echo "Clean all data."

doc: report/hw1_report.md
	pandoc -s -o report/hw1_report.pdf report/hw1_report.md --pdf-engine=xelatex