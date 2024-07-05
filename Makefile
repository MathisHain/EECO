all:
	matlab -nodisplay -nosplash -nodesktop -r "try, run('EnsembleRun.m'); catch, end, exit"