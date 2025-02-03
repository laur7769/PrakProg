Out.txt : main.exe
	mono $< > $@
	#Mono main.exe > Out.txt

main.exe : main.cs sfuns.dll
	mcs -target:exe -reference:sfuns.dll -out:main.exe main.cs

sfuns.dll : sfuns.cs
	mcs -target:library -out:sfuns.dll sfuns.cs
	# Vi specificerer at det skal være et library, ellers ledes der efter en Main funktion, da defualt er en executable

clean:
	$(RM) sfuns.dll main.exe Out.txt

test:
	echo $(RM)

