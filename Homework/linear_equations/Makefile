CODE = $(filter %.cs,$^)
LIBS = $(addprefix -reference:,$(filter %.dll,$^))
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CODE)
MKLIB = mcs -target:library -out:$@ $(LIBS) $(CODE)

TIME = gtime --portability --append --output $@
Out.txt : main.exe
	mono $< > $@

main.exe : main.cs matrix.dll 
	$(MKEXE)

matrix.dll : matrix.cs vector.cs
	$(MKLIB)



clean:	
	$(RM) *.exe *.dll [Oo]ut* [Ll]og*
