Complex.html: Complex.asciidoc
	# echo 'library(ascii); Asciidoc("Complex.asciidoc", width=9, height=4, prefix.string="Figs/fig", driver=cacheSweaveAsciidoc(), pdf=FALSE)' | R --vanilla
	R --vanilla < asciidoc.R
	asciidoc -a data-uri -a toc -a toclevels=4 -a num Complex.asciidoc.txt
	mv Complex.asciidoc.html Complex.html

clean: 
	\rm -r *~ *.bak *.dvi *.log *.aux *.tex *.html *.asciidoc.txt
