main.pdf: main.tex references.bib figures/*.pdf
	pdflatex -interaction=nonstopmode main.tex
	bibtex main
	pdflatex -interaction=nonstopmode main.tex
	pdflatex -interaction=nonstopmode main.tex

figures:
	python ../scripts/make_paper_figures.py --out figures

clean:
	rm -f *.aux *.bbl *.blg *.log *.out

.PHONY: figures clean
