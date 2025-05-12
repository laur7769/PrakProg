CODE = $(filter %.cs,$^)
LIBS = $(addprefix -reference:,$(filter %.dll,$^))
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CODE)
MKLIB = mcs -target:library -out:$@ $(LIBS) $(CODE)

TIME = gtime --portability --append --output $@

Ns := $(shell seq 0 10 750)

out.times-2.txt: main.exe Makefile
	for N in $(Ns); do \
        gtime -ao log.$$N -f "$$N %e" mono $< -task:time -size:$$N  & \
	done; wait
	>$@
	for N in $(Ns); do cat log.$$N >> $@; done
	$(RM) log*

Out.txt : main.exe
	mono $< > $@

data_fr.txt: main.exe
	mono $< -rmax:10 -dr:0.025 -task:fr >> $@

data_dr.txt: main.exe
	> $@
	for number in 0.025 0.05 0.1 0.2 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.8 0.9 1 ; do \
		mono $< -rmax:10 -dr:$$number -task:dr >> $@ ;\
	done

data_rmax.txt: main.exe
	for number in 3 4 5 6 7 8 9 10 11 12 13 14 15; do \
		mono main.exe -rmax:$$number -dr:0.3 -task:rmax >> $@; \
	done

dr.gnuplot.svg : data_dr.txt Makefile
	echo 'Plot of ε[0] as a function of ∆r with rmax=10 can be seen in dr.gnuplot.svg ' >> Out.txt
	echo '\
		set terminal svg background "white" ;\
		set out "$@" ;\
		set xlabel "∆r" ;\
		set ylabel "ε[0]" ;\
		set tics out ;\
		set grid ;\
		plot \
			 "$<" using 2:1 with linespoints title "ε[0] as a function of ∆r with rmax=10" ;\
			 ' | tee log.gpi  | gnuplot

rmax.gnuplot.svg : data_rmax.txt Makefile
		echo 'Plot of ε[0] as a function of r_max with ∆r=0.3 can be seen in rmax.gnuplot.svg ' >> Out.txt
		echo '\
		set terminal svg background "white" ;\
		set out "$@" ;\
		set xlabel "rmax" ;\
		set ylabel "ε[0]" ;\
		set tics out ;\
		set grid ;\
		plot \
			 "$<" using 3:1 with linespoints title "ε[0] as a function of rmax with ∆r=0.3" ;\
			 ' | tee log.gpi  | gnuplot


fr.gnuplot.svg : data_fr.txt Makefile
	echo 'Plot with the computed first 4 eigenfunctions and analytical results can be seen in fr.gnuplot.svg' >> Out.txt
	echo 'For the eigenfunctions r_max=10 and ∆r=0.025. This was determined based on the plots in dr.gnuplot.svg and rmax.gnuplot.svg' >> Out.txt
	echo '\
	set terminal svg background "white" ;\
	set out "$@" ;\
	set xlabel "r" ;\
	set ylabel "f_k(r)" ;\
	set tics out ;\
	set grid ;\
	set key top right;\
	f1(x) = 2.0*x*exp(-x) ;\
	f2(x) = x/sqrt(2)*(1.0-x/2)*exp(-x/2.0) ;\
	f3(x) =-(2.0*x)/(3.0*sqrt(3))*(1.0-(2.0*x/3.0)+(2.0/27.0)*x**2.0)*exp(-x/3.0) ;\
	f4(x) =-(x/4.0)*(1.0-(3.0*x)/4.0+(1.0/8.0)*x**2.0-(1.0/192.0)*x**3.0)*exp(-x/4.0) ;\
	plot[0:10][-0.6:1.4] \
		"$<" index 0 with points linetype 6 lc rgb "blue" title "k=0", f1(x) with lines lc rgb "blue" title "r*R_{10}" \
		,"$<" index 1 with points linetype 6 lc rgb "red" title "k=1", f2(x) with lines lc rgb "red" title "r*R_{20}" \
		,"$<" index 2 with points linetype 6 lc rgb "purple" title "k=2", f3(x) with lines lc rgb "purple" title "-r*R_{30}" \
		,"$<" index 3 with points linetype 6 lc rgb "green" title "k=3", f4(x) with lines lc rgb "green" title "-r*R_{40}" ;\
		 ' | tee log.gpi  | gnuplot

times.gnuplot.svg:  out.times-2.txt  Makefile
	echo '\
	set terminal svg background "white" ;\
	set out "$@" ;\
	set xlabel "N" ;\
	set ylabel "time [s]" ;\
	set tics out ;\
	set grid ;\
	set key top right;\
	f(x) = a*x*x ;\
	g(x) = b*x*x*x;\
	fit f(x) "out.times-2.txt" using 1:2 via a ;\
	fit g(x) "out.times-2.txt" using 1:2 via b ;\
	plot \
		"out.times-2.txt" index 0 with points title "Measured time of diagonalization as a function og matrix size" \
		, f(x) with lines title "a*N^3 fitted to data", g(x) with lines title "b*N^2 fitted to data";\
		 ' | tee log.gpi  | gnuplot

main.exe : main.cs matrix.dll 
	$(MKEXE)

matrix.dll : matrix.cs vector.cs
	$(MKLIB)



clean:	
	$(RM) *.exe *.dll [Oo]ut* [Ll]og* data*

