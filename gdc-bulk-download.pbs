#PBS -S /bin/bash
#PBS -N gdc-bulk-download
#PBS -q large
#PBS -l nodes=1
#PBS -l cput=72:00:00

cd $PBS_O_WORKDIR

/ifs/projects/GRCBL-NGS/slowtemp/local/gdc-client download -m manifest.xml -t token.txt
