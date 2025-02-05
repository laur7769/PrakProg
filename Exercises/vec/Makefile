CODE = $(filter %.cs,$^)
DLLS = $(addprefix -reference:, $(filter %.dll,$^))
MKLIB = mcs -out:$@ -target:library $(DLLS) $(CODE)
MKEXE = mcs -out:$@ -target:exe $(DLLS) $(CODE)

Out.txt: main.exe; mono $< > $@
main.exe: main.cs vec.dll; $(MKEXE)
vec.dll: vec.cs; $(MKLIB)
