for (( i = 0; i < 20; i++ )); do
	{ time mpirun --hostfile hostfile ./ray 2 4 1 ; } 2>> saida-2-4-4.txt
done