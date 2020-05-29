all: clean main-en.pdf

main-en.pdf: main-en.tex *.tex *.bib images/*
	pdflatex main-en
	bibtex main-en
	pdflatex main-en
	pdflatex main-en

clean:
	rm -f main-en.bbl main-en.blg main-en.out main-en.log main-en.toc main-en.aux
