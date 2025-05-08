CODE = $(filter %.cs,$^)
LIBS = $(addprefix -reference:,$(filter %.dll,$^))
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CODE)
MKLIB = mcs -target:library -out:$@ $(LIBS) $(CODE)

TIME = gtime --portability --append --output $@

out.times.data : main.exe
	>$@
	for N in $$(seq 0 10 1000); do \
		gtime --format "$$N %e" --output $@ --append \
		mono $< -size:$$N 1>out 2>err ;\
	done 
times.gnuplot.svg: out.times.data  Makefile
	echo '\
	set terminal svg background "white" ;\
	set out "$@" ;\
	set xlabel "N" ;\
	set ylabel "time [s]" ;\
	set tics out ;\
	set grid ;\
	set key top right;\
	f(x) = a*x*x*x ;\
	fit f(x) "out.times.data" using 1:2 via a ;\
	plot \
		"out.times.data" index 0 with points title "QR-decomposition time as a function of N", f(x) with lines title "fitted a*N^3" ;\
		 ' | tee log.gpi  | gnuplot
Out.txt : main.exe
	mono $< > $@

main.exe : main.cs matrix.dll 
	$(MKEXE)

matrix.dll : matrix.cs vector.cs
	$(MKLIB)



clean:	
	$(RM) *.exe *.dll [Oo]ut* [Ll]og*
