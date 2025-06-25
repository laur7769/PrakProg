CODE = $(filter %.cs,$^)
LIBS = $(addprefix -reference:,$(filter %.dll,$^))
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CODE)
MKLIB = mcs -target:library -out:$@ $(LIBS) $(CODE)

Out.txt: main.exe
	mono main.exe < higgs.txt 1> out.txt

Higgs.gnuplot.svg: Out.txt  Higgs.txt fit.txt  Makefile
	echo '\
	set terminal svg background "white" ;\
	set out "$@" ;\
	set xlabel "E [GeV]" ;\
	set ylabel "σ(E) [Certain units]" ;\
	set tics out ;\
	set grid ;\
	set key top right;\
	plot[100:160][-5:15] \
		"Higgs.txt" index 0 with yerrorbars title "Experimental data" \
		, "fit.txt" index 0 with lines title "Fit", "fit.txt" index 1 with lines title "Fit central difference approximation" ;\
		 ' | tee log.gpi  | gnuplot



main.exe : main.cs matrix.dll
	$(MKEXE)

matrix.dll : matrix.cs vector.cs QR.cs
	$(MKLIB)

clean:	
	$(RM) *.exe *.dll [Oo]ut* [Ll]og* data* fit*
