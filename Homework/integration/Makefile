CODE = $(filter %.cs,$^)
LIBS = $(addprefix -reference:,$(filter %.dll,$^))
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CODE)
MKLIB = mcs -target:library -out:$@ $(LIBS) $(CODE)

Out.txt: main.exe
	mono $< > $@

erf.gnuplot.svg: Out.txt  data_A.txt Makefile
	echo '\
	set terminal svg background "white" ;\
	set out "$@" ;\
	set xlabel "x" ;\
	set ylabel "erf(x)" ;\
	set tics out ;\
	set grid ;\
	set key top right;\
	plot[-5:5][-1.4:1.4] \
		"data_A.txt" index 1 using 1:2 with lines title "Computed values for the error function" \
		, "data_A.txt" index 0  with points title "Tabulated values";\
		 ' | tee log.gpi  | gnuplot

acc.gnuplot.svg: Out.txt  data_A.txt Makefile
	echo '\
	set terminal svg background "white" ;\
	set out "$@" ;\
	set xlabel "acc" ;\
	set ylabel "Absolute difference" ;\
	set logscale y;\
	set logscale x;\
	set format x "10^{%T}";\
	set tics out ;\
	set grid ;\
	set key top right;\
	plot[0.0000001:1][0.0000000001:0.01] \
		"data_A.txt" index 2 using 1:2 with points title "Difference between computed and tabulated values for z=1 with eps=0";\
		 ' | tee log.gpi  | gnuplot

main.exe : main.cs matrix.dll
	$(MKEXE)

matrix.dll : matrix.cs vector.cs
	$(MKLIB)

clean:	
	$(RM) *.exe *.dll [Oo]ut* [Ll]og* data*
