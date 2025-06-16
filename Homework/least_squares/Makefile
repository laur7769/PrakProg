CODE = $(filter %.cs,$^)
LIBS = $(addprefix -reference:,$(filter %.dll,$^))
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CODE)
MKLIB = mcs -target:library -out:$@ $(LIBS) $(CODE)
FUNCTION := $(shell egrep 'f1' Out.txt)
FUNCTION2 := $(shell egrep 'f2' Out.txt)
FUNCTION3 := $(shell egrep 'f3' Out.txt)
FUNCTION4 := $(shell egrep 'f4' Out.txt)
FUNCTION5 := $(shell egrep 'f5' Out.txt)

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
	$(FUNCTION2) ;\
	$(FUNCTION3) ;\
	$(FUNCTION4) ;\
	$(FUNCTION5) ;\
	plot \
		"$<" index 1 using 1:2:3 with yerrorbars title "Measured data", f1(x) with lines title "Least-squares fit", \
		f2(x) with lines title "F=∑(ĉk-δck)fk(x)", f3(x) with lines title "F=(ĉ1+δc1)f1(x)+(ĉ2-δc2)f2(x)", \
		f4(x) with lines title "F=∑(ĉk+δck)fk(x)", f5(x) with lines title "F=(ĉ1-δc1)f1(x)+(ĉ2+δc2)f2(x)" ;\
		 ' | tee log.gpi  | gnuplot

main.exe : main.cs matrix.dll QR.dll
	$(MKEXE)

QR.dll : QR.cs matrix.dll
	$(MKLIB)

matrix.dll : matrix.cs vector.cs
	$(MKLIB)

clean:	
	$(RM) *.exe *.dll [Oo]ut* [Ll]og* data*
