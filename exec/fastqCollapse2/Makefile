x3: fastqCollapse
	./fastqCollapse -i test_1.fq test_2.fq -sigLength 50 -o test -summary summary.txt

clean:
	rm fastqCollapse

fastqCollapse: fastqCollapse.cpp
	g++ fastqCollapse.cpp -static -o fastqCollapse -lz
