cc = g++

main: main.cpp
	$(cc) $< -lboost_program_options -o $@ -O3

main_bcc: main_bcc.cpp
	$(cc) $< -lboost_program_options -o $@ -O3

clean:
	rm ./main ./main_bcc

.PHONY: run clean