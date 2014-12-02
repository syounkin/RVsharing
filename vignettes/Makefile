ped2trio.html: ped2trio.asciidoc
	R --vanilla < asciidoc.R
	asciidoc -a data-uri -a toc -a toclevels=4 -a num ped2trio.asciidoc.txt
	mv ped2trio.asciidoc.html ped2trio.html

clean: 
	\rm -r *~ *.bak *.dvi *.log *.aux *.tex *.html *.asciidoc.txt

