#PBS -S /bin/bash
#PBS -N bam-index
#PBS -q medium
#PBS -l nodes=1
#PBS -l cput=72:00:00

module add samtools
mkdir err/
mkdir logs/

cd $PBS_O_WORKDIR

for bamfile in *.bam; do
	echo $bamfile
	# samtools sort "${bamfile}" "${bamfile}.sorted"
	# samtools index ${bamfile}.sorted > err/$PBS_JOBID > logs/$PBS_JOBID
	samtools index $bamfile > err/$PBS_JOBID > logs/$PBS_JOBID
done
exit 0
