### Introduction

### Sample workflow
#### Login to the Moab HPC cluster
Open your terminal and log in with your NIH username and password:
```ssh USERNAME@moab.ncifcrf.gov
```

Then, open up an interactive job (so that the commands/scripts that you run
won't slow down the head node of the server):
```qsub -I
```
Enter in `exit` if you ever want
to leave it.

Change to the lab directory:
```cd /ifs/projects/GRCBL-NGS/slowtemp
```

Create and go inside a project directory:
```mkdir project1
   cd project1
```
