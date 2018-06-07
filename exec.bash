for (( i = 0; i < 20; i++ )); do
	{ time mpirun --hostfile hostfile ./ray 1 4 1 ; } 2>> saida-1-4-16.txt
done