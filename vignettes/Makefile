GeneDropping.html: GeneDropping.asciidoc
	R --vanilla < asciidoc.R
	asciidoc -a data-uri -a toc -a toclevels=4 -a num GeneDropping.asciidoc.txt
	mv GeneDropping.asciidoc.html GeneDropping.html

clean: 
	\rm -r *~ *.bak *.dvi *.log *.aux *.tex *.html *.asciidoc.txt

