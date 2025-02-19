CODE = $(filter %.cs,$^)
LIBS = $(addprefix -reference:,$(filter %.dll,$^))
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CODE)
MKLIB = mcs -target:library -out:$@ $(LIBS) $(CODE)

all: plot.test.svg

plot.test.svg: test.gpi
	gnuplot $<		

Out.txt: main.exe
	mono $< > $@

sfuns.dll : sfuns.cs ; $(MKLIB)
main.exe : main.cs sfuns.dll; $(MKEXE)

clean:	
	$(RM) *.exe *.dll [Oo]ut* [Ll]og*
