CODE = $(filter %.cs,$^)
LIBS = $(addprefix -reference:,$(filter %.dll,$^))
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CODE)
MKLIB = mcs -target:library -out:$@ $(LIBS) $(CODE)

Out.txt: main.exe
	mono $< > $@

MCB.gnuplot.svg: data_B.txt data_A.txt Makefile
	echo '\
	set terminal svg background "white" ;\
	set out "$@" ;\
	set xlabel "N" ;\
	set ylabel "error" ;\
	set tics out ;\
	set grid ;\
	set key top right;\
	t(x)=a/sqrt(x);\
	fit t(x) "data_A.txt" index 0 using 1:3 via a;\
	plot[0:100000][0:0.03] \
		"data_B.txt" index 0 using 1:3 with lines lw 0.5 title "Actual error, Halton MC" \
		,"data_B.txt" index 0 using 1:2 with lines lw 0.5 title "Estimated error, Halton MC", "data_A.txt" index 0 using 1:2 with lines lw 0.5 title "Estimated error, plain MC" ;\
		 ' | tee log.gpi  | gnuplot

MCA.gnuplot.svg: data_A.txt Makefile
	echo '\
	set terminal svg background "white" ;\
	set out "$@" ;\
	set xlabel "N" ;\
	set ylabel "error" ;\
	set tics out ;\
	set grid ;\
	set key top right;\
	t(x)=a/sqrt(x);\
	fit t(x) "data_A.txt" index 0 using 1:3 via a;\
	plot[0:100000][0:0.03] \
		"data_A.txt" index 0 using 1:3 with lines lw 0.5 title "Actual error, plain MC" \
		,"data_A.txt" index 0 using 1:2 with lines lw 0.5 title "Estimated error, plain MC", ;\
		 ' | tee log.gpi  | gnuplot

MCC.gnuplot.svg: data_C.txt data_A.txt Makefile
	echo '\
	set terminal svg background "white" ;\
	set out "$@" ;\
	set xlabel "N" ;\
	set ylabel "error" ;\
	set tics out ;\
	set grid ;\
	set key top right;\
	t(x)=a/sqrt(x);\
	fit t(x) "data_C.txt" index 0 using 1:3 via a;\
	plot[0:100000][0:0.02] \
		"data_C.txt" index 0 using 1:3 with lines lw 0.5 title "Actual error, stratified sampling" \
		,"data_C.txt" index 0 using 1:2 with lines lw 0.5 title "Estimated error, stratified sampling", "data_A.txt" index 0 using 1:2 with lines lw 0.5 title "Estimated error, plain MC" ;\
		 ' | tee log.gpi  | gnuplot

main.exe : main.cs matrix.dll
	$(MKEXE)

matrix.dll : matrix.cs vector.cs
	$(MKLIB)

clean:	
	$(RM) *.exe *.dll [Oo]ut* [Ll]og* data* *svg
