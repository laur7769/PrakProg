CODE = $(filter %.cs,$^)
LIBS = $(addprefix -reference:,$(filter %.dll,$^))
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CODE)
MKLIB = mcs -target:library -out:$@ $(LIBS) $(CODE)

Out.txt: main.exe
	mono main.exe > out.txt

A.gnuplot.svg: data_A.txt  Makefile
	echo '\
	set terminal svg background "white" ;\
	set out "$@" ;\
	set xlabel "x" ;\
	set ylabel "y" ;\
	set tics out ;\
	set grid ;\
	set key top right;\
	plot \
		"data_A.txt" index 0 with lines title "Cos(5*x-1)*Exp(-x*x)" \
		, "data_A.txt" index 1 with points title "Neural network";\
		 ' | tee log.gpi  | gnuplot

B.gnuplot.svg: data_B.txt  Makefile
	echo '\
	set terminal svg background "white" ;\
	set out "$@" ;\
	set xlabel "x" ;\
	set ylabel "y" ;\
	set tics out ;\
	set grid ;\
	set key top right;\
	df(x) = -5*sin(5*x-1)*exp(-x*x)-2*x*exp(-x*x)*cos(5*x-1) ;\
	ddf(x) = exp(-x*x)*((4*x*x-27)*cos(5*x-1)+20*x*sin(5*x-1)) ;\
	F(x) = -0.5*exp(-x*x) ;\
	plot[-1:1][-30:50] \
		"data_B.txt" index 0 with points title "Neural derivative" \
		, "data_B.txt" index 1 with points title "Neural second derivative"\
		,"data_B.txt" index 2 with points title "Neural anti-derivative"\
		,df(x) with lines title "dg(x)" \
		, ddf(x) with lines title "ddg(x)" \
		, F(x) with lines title "G(x)" ;\
		 ' | tee log.gpi  | gnuplot

C.gnuplot.svg:  data_C.txt  Makefile
	echo '\
	set terminal svg background "white" ;\
	set out "$@" ;\
	set xlabel "x" ;\
	set ylabel "y" ;\
	set tics out ;\
	set grid ;\
	set key top right;\
	plot \
		"data_C.txt" index 0 with points title "Numerical solution" \
		, "data_C.txt" index 1 with points title "Neural network response";\
		 ' | tee log.gpi  | gnuplot
main.exe : main.cs matrix.dll
	$(MKEXE)

matrix.dll : matrix.cs vector.cs QR.cs
	$(MKLIB)

clean:	
	$(RM) *.exe *.dll [Oo]ut* [Ll]og* data* fit*
