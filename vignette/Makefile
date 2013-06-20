RVsharing.html: RVsharing.asciidoc
	# echo 'library(ascii); Asciidoc("RVsharing.asciidoc", width=9, height=4, prefix.string="Figs/fig", driver=cacheSweaveAsciidoc(), pdf=FALSE)' | R --vanilla
	R --vanilla < asciidoc.R
	asciidoc -a data-uri -a toc -a toclevels=4 -a num RVsharing.asciidoc.txt
	mv RVsharing.asciidoc.html RVsharing.html

clean: 
	\rm -r *~ *.bak *.dvi *.log *.aux *.tex *.html *.asciidoc.txt
