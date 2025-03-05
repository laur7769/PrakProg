CODE = $(filter %.cs,$^)
LIBS = $(addprefix -reference:,$(filter %.dll,$^))
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CODE)
MKLIB = mcs -target:library -out:$@ $(LIBS) $(CODE)

Out.txt: main.exe
	mono $< > $@

interp.gnuplot.svg: Out.txt data.txt Makefile
	echo '\
	set terminal svg background "white" ;\
	set out "$@" ;\
	set xlabel "x" ;\
	set ylabel "y" ;\
	set tics out ;\
	set grid ;\
	set key top right;\
	f(x)=cos(x);\
	plot[0:9][-1.5:1.5] \
		"data.txt" index 0 using 1:2 with points title "Given data points", "data.txt" index 1 with lines title "Linear interpolation" \
		, "data.txt" index 2 using 1:2 with lines title "Quadratic splines", f(x) with lines title "f(x)=Cos(x)";\
		 ' | tee log.gpi  | gnuplot

int.gnuplot.svg: Out.txt data.txt Makefile
	echo '\
	set terminal svg background "white" ;\
	set out "$@" ;\
	set xlabel "x" ;\
	set ylabel "Integral from x[0]" ;\
	set tics out ;\
	set grid ;\
	set key top right;\
	f(x)=sin(x);\
	plot[0:9] \
		"data.txt" index 1 using 1:3 with lines title "Linear interpolation" \
		, "data.txt" index 2 using 1:3 with lines title "Quadratic splines", f(x) with lines title "Sin(x) (The definite integral from x[0] to x of cos(x))"  ;\
		 ' | tee log.gpi  | gnuplot

main.exe : main.cs matrix.dll QR.dll
	$(MKEXE)

QR.dll : QR.cs matrix.dll
	$(MKLIB)

matrix.dll : matrix.cs vector.cs
	$(MKLIB)

clean:	
	$(RM) *.exe *.dll [Oo]ut* [Ll]og* data*
