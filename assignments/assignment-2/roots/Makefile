# $MYHG/assignments/assignment-2/roots/Makefile

.PHONY: test plot clean clobber

test: testzero.py
	python testzero.py

plot: 
	@echo "Assuming .txt files are correct!"
	python makeplot.py

clean:
	rm -f *~ *.pyc

clobber: clean
	rm -f fvals.txt newtonpts.txt plot.png
