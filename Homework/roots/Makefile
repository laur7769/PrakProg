CODE = $(filter %.cs,$^)
LIBS = $(addprefix -reference:,$(filter %.dll,$^))
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CODE)
MKLIB = mcs -target:library -out:$@ $(LIBS) $(CODE)

Out.txt: main.exe
	mono $< > $@

B.gnuplot.svg: Out.txt  data_B.txt Makefile
	echo '\
	set terminal svg background "white" ;\
	set out "$@" ;\
	set xlabel "r" ;\
	set ylabel "f_0(r)" ;\
	set tics out ;\
	set grid ;\
	set key top right;\
	t(x)=x*exp(-x);\
	plot \
		"data_B.txt" index 0 using 1:2 with points title "Computed f_0(r)" \
		, t(x) with lines title "Exact f_0(r)" ;\
		 ' | tee log.gpi  | gnuplot

rmax.gnuplot.svg: Out.txt  data_B.txt Makefile
	echo '\
	set terminal svg background "white" ;\
	set out "$@" ;\
	set xlabel "rmax" ;\
	set ylabel "E_0" ;\
	set tics out ;\
	set grid ;\
	set key top right;\
	plot[1:11][-1:2.5] \
		"data_B.txt" index 1 using 1:2 with lines title "E_0 as a function of rmax, rmin=0.01" ;\
		 ' | tee log.gpi  | gnuplot

rmin.gnuplot.svg: Out.txt  data_B.txt Makefile
	echo '\
	set terminal svg background "white" ;\
	set out "$@" ;\
	set xlabel "rmin" ;\
	set ylabel "E_0" ;\
	set tics out ;\
	set grid ;\
	set key top right;\
	set logscale x ;\
	plot \
		"data_B.txt" index 2 using 1:2 with points title "E_0 as a function of rmin, rmax=8" ;\
		 ' | tee log.gpi  | gnuplot

acc.gnuplot.svg: Out.txt  data_B.txt Makefile
	echo '\
	set terminal svg background "white" ;\
	set out "$@" ;\
	set xlabel "accuracy" ;\
	set ylabel "E_0" ;\
	set tics out ;\
	set grid ;\
	set key top right;\
	set logscale x ;\
	plot \
		"data_B.txt" index 3 using 1:2 with points title "E_0 as a function of accuracy, rmin=0.01, rmax=8, eps=0.01" ;\
		 ' | tee log.gpi  | gnuplot

eps.gnuplot.svg: Out.txt  data_B.txt Makefile
	echo '\
	set terminal svg background "white" ;\
	set out "$@" ;\
	set xlabel "epsilon" ;\
	set ylabel "E_0" ;\
	set tics out ;\
	set grid ;\
	set key top right;\
	set logscale x ;\
	plot \
		"data_B.txt" index 4 using 1:2 with points title "E_0 as a function of epsilon, rmin=0.01, rmax=8, acc=0.01" ;\
		 ' | tee log.gpi  | gnuplot

main.exe : main.cs matrix.dll
	$(MKEXE)

matrix.dll : matrix.cs vector.cs QR.cs
	$(MKLIB)

clean:	
	$(RM) *.exe *.dll [Oo]ut* [Ll]og* data*
