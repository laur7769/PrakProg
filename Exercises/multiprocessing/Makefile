CODE = $(filter %.cs,$^)
LIBS = $(addprefix -reference:,$(filter %.dll,$^))
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CODE)
MKLIB = mcs -target:library -out:$@ $(LIBS) $(CODE)

N=1e8
TIME = gtime --portability --append --output $@
Out.txt : main.exe main2.exe main3.exe Makefile
	>$@
	echo '1. Manual multithreading' >>$@
	$(TIME) mono $< -terms:$N -threads:1 >>$@
	$(TIME) mono $< -terms:$N -threads:2 >>$@
	$(TIME) mono $< -terms:$N -threads:3 >>$@
	$(TIME) mono $< -terms:$N -threads:4 >>$@
	$(TIME) mono main2.exe -terms:$N >>$@ 
	$(TIME) mono main3.exe -terms:$N >>$@

main.exe : main.cs; $(MKEXE)

main2.exe: main2.cs; $(MKEXE)

main3.exe: main3.cs; $(MKEXE)

clean:	
	$(RM) *.exe *.dll [Oo]ut* [Ll]og*
