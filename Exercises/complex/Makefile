CODE = $(filter %.cs,$^)
DLLS = $(addprefix -reference:, $(filter %.dll,$^))
MKLIB = mcs -out:$@ -target:library $(DLLS) $(CODE)
MKEXE = mcs -out:$@ -target:exe $(DLLS) $(CODE)
Out.txt: main.exe
	mono main.exe > $@ 
main.exe : main.cs cmath.dll
	mcs -reference:cmath.dll -target:exe -out:main.exe main.cs

cmath.dll : cmath.cs complex.cs
	mcs -target:library -out:./cmath.dll $^
