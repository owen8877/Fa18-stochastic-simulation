cc = g++

main: main.cpp
	$(cc) $< -lboost_program_options -o $@ -Ofast

clean:
	rm ./main

run: main
	./main

.PHONY: run clean