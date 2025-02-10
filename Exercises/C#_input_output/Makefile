CODE = $(filter %.cs,$^)
DLLS = $(addprefix -reference:, $(filter %.dll,$^))
MKLIB = mcs -out:$@ -target:library $(DLLS) $(CODE)
MKEXE = mcs -out:$@ -target:exe $(DLLS) $(CODE)
Out.txt: main.exe
	echo 1 2 3 4 5 | mono main.exe -numbers:1,2,3,4,5 -input:inputfile.txt -output:outputfile.txt > $@ 2> Error.txt

main.exe: main.cs
	$(MKEXE)
