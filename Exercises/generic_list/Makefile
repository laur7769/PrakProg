CODE = $(filter %.cs,$^)
DLLS = $(addprefix -reference:, $(filter %.dll,$^))
MKLIB = mcs -out:$@ -target:library $(DLLS) $(CODE)
MKEXE = mcs -out:$@ -target:exe $(DLLS) $(CODE)
Out.txt: main.exe
	cat inputfile.txt | mono main.exe > $@ 2> Error.txt
main.exe: main.cs
	$(MKEXE)
