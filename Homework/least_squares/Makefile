CODE = $(filter %.cs,$^)
LIBS = $(addprefix -reference:,$(filter %.dll,$^))
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CODE)
MKLIB = mcs -target:library -out:$@ $(LIBS) $(CODE)
FUNCTION := $(shell egrep 'f1' Out.txt)

Out.txt: main.exe
	mono $< > $@

fit.gnuplot.svg: Out.txt Makefile
	echo '\
	set terminal svg background "white" ;\
	set out "$@" ;\
	set xlabel "Time [days]" ;\
	set ylabel "ln(y)" ;\
	set tics out ;\
	set grid ;\
	set key top right;\
	$(FUNCTION) ;\
	plot \
		"$<" index 1 using 1:2:3 with yerrorbars title "Measured data", f1(x) with lines title "Least-squares fit" ;\
		 ' | tee log.gpi  | gnuplot

main.exe : main.cs matrix.dll QR.dll
	$(MKEXE)

QR.dll : QR.cs matrix.dll
	$(MKLIB)

matrix.dll : matrix.cs vector.cs
	$(MKLIB)

clean:	
	$(RM) *.exe *.dll [Oo]ut* [Ll]og* data*
