CODE = $(filter %.cs,$^)
LIBS = $(addprefix -reference:,$(filter %.dll,$^))
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CODE)
MKLIB = mcs -target:library -out:$@ $(LIBS) $(CODE)

TIME = gtime --portability --append --output $@
Out.txt : main.exe
	mono $< > $@

data_fr.txt: main.exe
	mono $< -rmax:10 -dr:0.025 >> $@

data_dr.txt: main.exe
	> $@
	for number in 0.025 0.05 0.1 0.2 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.8 0.9 1 ; do \
		mono $< -rmax:10 -dr:$$number >> $@ ;\
	done

data_rmax.txt: main.exe
	for number in 3 4 5 6 7 8 9 10 11 12 13 14 15; do \
		mono main.exe -rmax:$$number -dr:0.3 >> $@; \
	done

dr.gnuplot.svg : data_dr.txt Makefile
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
	f3(x) =(2.0*x)/(3.0*sqrt(3))*(1.0-(2.0*x/3.0)+(2.0/27.0)*x**2.0)*exp(-x/3.0) ;\
	f4(x) =(x/4.0)*(1.0-(3.0*x)/4.0+(1.0/8.0)*x**2.0-(1.0/192.0)*x**3.0)*exp(-x/4.0) ;\
	plot \
		"$<" index 0 with points lc rgb "blue" title "k=0", f1(x) with lines lc rgb "blue" title "r*R_{10}" \
		,"$<" index 1 with points lc rgb "red" title "k=1", f2(x) with lines lc rgb "red" title "r*R_{20}" \
		,"$<" index 2 with points lc rgb "purple" title "k=2", f3(x) with lines lc rgb "purple" title "r*R_{30}" \
		,"$<" index 3 with points lc rgb "green" title "k=3", f4(x) with lines lc rgb "green" title "r*R_{40}" ;\
		 ' | tee log.gpi  | gnuplot

main.exe : main.cs matrix.dll 
	$(MKEXE)

matrix.dll : matrix.cs vector.cs
	$(MKLIB)



clean:	
	$(RM) *.exe *.dll [Oo]ut* [Ll]og* data*

