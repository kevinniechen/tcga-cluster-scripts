#PBS -S /bin/bash
#PBS -N separate-dif-len-bams
#PBS -q medium
#PBS -l nodes=1
#PBS -l cput=24:00:00

module add samtools
mkdir err/
mkdir logs/

cd $PBS_O_WORKDIR

touch bam-information.txt

for f in *.bam; do
	if [$(samtools view $f | head -n 10000000 | cut -f 10 | perl -ne 'chomp;print length($_) . "\n"' | sort | uniq -c | wc -l) -gt 1]; then
		rm $f
		rm "${f}.bai"
	else
		readlength=$(samtools view $f | head -n 10000000 | cut -f 10 | perl -ne 'chomp;print length($_) . "\n"' | sort | uniq -c | awk '{ print $2 }')
		readcount=$(samtools view $f | head -n 10000000 | cut -f 10 | perl -ne 'chomp;print length($_) . "\n"' | sort | uniq -c | awk '{ print $1 }')
		printf "${f} ${readlength} ${readcount}\n" >> bam-information.txt

		mkdir -p "bams-${readlength}"
		mv $f "bams-${readlength}"
		mv "${f}.bai" "bams-${readlength}"
	fi
done
exit 0
