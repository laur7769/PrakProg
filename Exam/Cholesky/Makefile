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
	g(x) = 5.26945e-09*x*x*x;\
	plot \
		"out.times.data" index 0 with points title "QR-decomposition time as a function of N", g(x) with lines title "Fitted function decribing time of QR decomposition"  ;\
		 ' | tee log.gpi  | gnuplot

Out.txt: main.exe
	mono $< > $@



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
