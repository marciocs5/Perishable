for heur in {0,1,2,3}
do
	for m in {2,5,7} 
	do
		for gamma in {1,3,5} 
		do
			for file in DATA/* 
			do
				timeout 30m ./main RESULTS_5.txt $file $gamma $m 1 $heur
			done
		done
	done
done
for heur in {1,2,3}
do
	for m in {2,5,7} 
	do
		for gamma in {1,3,5} 
		do
			for file in DATA/* 
			do
				timeout 30m ./main RESULTS_5.txt $file $gamma $m 0 $heur
			done
		done
	done
done
