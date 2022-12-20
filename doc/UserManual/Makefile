file=MuPhiSimManual

all:
	pdflatex $(file)
	pdflatex $(file)
	bibtex $(file)
	pdflatex $(file)
	pdflatex $(file)
	pdflatex $(file)
	@rm -rf *.aux *.bbl *.blg *.cb *.cb2 *.dvi *.idx *.ilg *.ind *.loa *.lof *.log *.lot *.maf *.mtc* *.out *.ps *.thm *.toc *.tpt *.backup *~

short: 
	pdflatex $(file)

clean:
	@rm -f *.aux *.bbl *.blg *.cb *.cb2 *.dvi *.idx *.ilg *.ind *.loa *.lof *.log *.lot *.maf *.mtc* *.out *.pdf *.ps *.thm *.toc *.tpt *.backup *~
