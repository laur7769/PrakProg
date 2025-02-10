Out.txt : Makefile main.exe input.txt
	echo "hello" > $@
	echo 1,2,3,4,5,6 | sed -e "s/,/\./g" >> $@
	echo "hello world from stdin" | mono main.exe 23,23 44,44 55, 55 1>> $@ 2> log
	cat input.txt | mono main.exe  1>> $@ 2>> log
	mono main.exe <input.txt 1>> $@ 2>> log

input.txt : Makefile
	echo "these go $@" > $@

main.exe : main.cs
	mcs -target:exe -out:$@ $<

test:
	echo $@
	echo '$@'
	echo "$@"


