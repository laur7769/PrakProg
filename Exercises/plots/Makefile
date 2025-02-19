CODE = $(filter %.cs,$^)
LIBS = $(addprefix -reference:,$(filter %.dll,$^))
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CODE)
MKLIB = mcs -target:library -out:$@ $(LIBS) $(CODE)

Out.erf.gnuplot.svg : Out.txt Makefile
	echo '\
		set terminal svg background "white" ;\
		set out "$@" ;\
		set key bottom right ;\
		set xlabel "x" ;\
		set ylabel "erf(x)" ;\
		set tics out ;\
		set grid ;\
		set title "Error function" ;\
		plot [0:2.5][0:1.2]\
			 "$<" index 0 with lines title "Error function approximation" \
			 ,"erf.txt" with points pointtype 6 title "Table values";\
			 ' | tee log.gpi  | gnuplot 

Out.sgamma.gnuplot.svg : Out.txt Makefile
	echo '\
		set terminal svg background "white" ;\
		set out "$@" ;\
		set key bottom right ;\
		set xlabel "x" ;\
		set ylabel "Γ(x)" ;\
		set tics out ;\
		set grid ;\
		set title "Gamma function" ;\
		plot [-5:5][-5:5]\
			 "$<" index 1 with lines title "Stirling approximation" \
			 ,"Gamma.txt" with points pointtype 6 title "Table values";\
			 ' | tee log2.gpi  | gnuplot

Out.lngamma.gnuplot.svg : Out.txt Makefile
	echo '\
		set terminal svg background "white" ;\
		set out "$@" ;\
		set key bottom right ;\
		set xlabel "x" ;\
		set ylabel "ln(Γ(x))" ;\
		set tics out ;\
		set grid ;\
		set title "Gamma function" ;\
		plot [0:5]\
			 "$<" index 2 with lines title " lnGamma" ;\
			 ' | tee log2.gpi  | gnuplot

Out.txt: main.exe
	mono $< > $@

sfuns.dll : sfuns.cs ; $(MKLIB)
main.exe : main.cs sfuns.dll; $(MKEXE)

clean:	
	$(RM) *.exe *.dll [Oo]ut* [Ll]og*
