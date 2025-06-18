CODE = $(filter %.cs,$^)
LIBS = $(addprefix -reference:,$(filter %.dll,$^))
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CODE)
MKLIB = mcs -target:library -out:$@ $(LIBS) $(CODE)

Out.txt: main.exe
	mono $< > $@

first_ODE.gnuplot.svg: Out.txt Makefile
	echo '\
	set terminal svg background "white" ;\
	set out "$@" ;\
	set xlabel "x" ;\
	set ylabel "y" ;\
	set tics out ;\
	set grid ;\
	set key top right;\
	f(x) = 5*exp(-x) ;\
	plot \
		"Out.txt" index 1 using 1:2 with points title "Computed values" \
		, f(x)  with lines title "Analytical solution f(x)=5exp(-x)";\
		 ' | tee log.gpi  | gnuplot

osc_fric.gnuplot.svg: Out.txt Makefile
	echo '\
	set terminal svg background "white" ;\
	set out "$@" ;\
	set xlabel "x" ;\
	set ylabel "y" ;\
	set tics out ;\
	set grid ;\
	set key top right;\
	plot \
		"Out.txt" index 2 using 1:2 with points title "theta(t)", "Out.txt" index 2 using 1:3 with points title "omega(t)" ;\
		 ' | tee log.gpi  | gnuplot

planet.gnuplot.svg: Out.txt data_B.txt Makefile
	echo '\
	set terminal svg background "white" ;\
	set out "$@" ;\
	set xlabel "x" ;\
	set ylabel "y" ;\
	set tics out ;\
	set grid ;\
	set key top right;\
	plot[-2:2][-1:3] \
		"data_B.txt" index 0 using (1/$$2)*cos($$1):(1/$$2)*sin($$1) with points pt 6 ps 1.5 title "ε=0, u(0)=1.0, du/dφ(0)=0.0", "data_B.txt" index 1 using (1/$$2)*cos($$1):(1/$$2)*sin($$1) with lines linewidth 2 title "ε=0, u(0)=1.0, du/dφ(0)=-0.5" \
		,"data_B.txt" index 2 using (1/$$2)*cos($$1):(1/$$2)*sin($$1) with lines title "ε=0.01, u(0)=1.0, du/dφ(0)=-0.5";\
		 ' | tee log.gpi  | gnuplot

23.gnuplot.svg: Out.txt data_C.txt Makefile
	echo '\
	set terminal svg background "white" ;\
	set out "$@" ;\
	set xlabel "x" ;\
	set ylabel "y" ;\
	set tics out ;\
	set grid ;\
	set key top right;\
	plot \
		"data_C.txt" index 0 using 1:2 with points title "Numerical results using a 23 stepper", (1.0/3.0)*x**3 with lines title "Analytical result: y=1/3x^3" ;\
		 ' | tee log.gpi  | gnuplot

three_body.gnuplot.svg: data_C.txt Makefile
	echo '\
	set terminal svg background "white" ;\
	set out "$@" ;\
	set xlabel "x" ;\
	set ylabel "y" ;\
	set tics out ;\
	set grid ;\
	set key top right;\
	plot[-1.5:1.5][-1:1] \
		"data_C.txt" index 1 using 2:3 with lines linewidth 3 title "Body 1", "data_C.txt" index 1 using 4:5 with lines linewidth 3 title "Body 2", "data_C.txt" index 1 using 6:7 with lines linewidth 3 title "Body 3", \
		"data_C.txt" index 2 with points pt 6 ps 2 title "Initial positions" ;\
		 ' | tee log.gpi  | gnuplot
main.exe : main.cs matrix.dll
	$(MKEXE)

matrix.dll : matrix.cs vector.cs splines.cs
	$(MKLIB)

clean:	
	$(RM) *.exe *.dll [Oo]ut* [Ll]og* data*
