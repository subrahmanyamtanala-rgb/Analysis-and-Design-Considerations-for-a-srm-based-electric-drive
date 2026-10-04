all: main.pdf response_to_reviewers.pdf response_to_reviewers_round2.pdf

main.pdf: main.tex references.bib figures/*.pdf
	pdflatex -interaction=nonstopmode main.tex
	bibtex main
	pdflatex -interaction=nonstopmode main.tex
	pdflatex -interaction=nonstopmode main.tex

response_to_reviewers.pdf: response_to_reviewers.tex
	pdflatex -interaction=nonstopmode response_to_reviewers.tex
	pdflatex -interaction=nonstopmode response_to_reviewers.tex

# Regenerate figures. Needs results/fea/table.npz from scripts/run_fea.py.
# The revision step must run last because it overwrites the drive figures with FE-based ones.
figures:
	python ../scripts/make_paper_figures.py --out figures
	python ../scripts/run_revision.py

clean:
	rm -f *.aux *.bbl *.blg *.log *.out

.PHONY: all figures clean

response_to_reviewers_round2.pdf: response_to_reviewers_round2.tex
	pdflatex -interaction=nonstopmode response_to_reviewers_round2.tex
	pdflatex -interaction=nonstopmode response_to_reviewers_round2.tex
