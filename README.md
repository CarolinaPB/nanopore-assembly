# snakemake-template
Template directory for creating a snakemake pipeline

## To use this repository as a template:
Click the "use this template" button

## Other instructions
Install `conda` if you don't have it

### Create conda environment

Recommended - give the profile a name related to your pipeline (ex: polish-assembly)

```
conda create --name <env> --file requirements.txt
```

(by creating an environment from requirements.txt you'll be creating and environment that already has snakemake)
### Activate environment
```
conda activate <env>
```

### To deactivate the environment (if you want to leave the conda environment)
```
conda deactivate
```

### Create hpc config file ([good example](https://www.sichong.site/2020/02/25/snakemake-and-slurm-how-to-manage-workflow-with-resource-constraint-on-hpc/))

Necessary for snakemake to prepare and send jobs.   
Recommended - give the profile a name related to this pipeline (ex: polish-assembly)

#### Start with creating the directory
```
mkdir -p ~/.config/snakemake/<profile name>
```

#### Add config.yaml to that directory and add the specifications:
```
jobs: 10
cluster: "sbatch -t 1:0:0 --mem=16000 -c 16 --job-name={rule} --exclude=fat001,fat002,fat101,fat100 --output=logs_slurm/{rule}.out --error=logs_slurm/{rule}.err"

use-conda: true
```
(change the options between square brackets)

## How to run

First it's good to always make a dry run: shows if there are any problems with the rules and we can use it to look at the commands and verify that all the fields are in the correct place

Dry run (prints execution plan and commands that will be run)
```
snakemake -np 
```
Run in the HPC 
```
snakemake --profile <profile name>
```

Other flags:
- --forceall : run all the steps, even if it's not needed
- --rerun-incomplete : rerun incomplete steps
- -R [rulename] : run this specific rule
- --max-jobs-per-second \<N> : sometimes there are some problems with the job timings/ many jobs being submitted at once so it's good to choose a low number

--------
If running the rules using slurm it's important that the logs_slurm directory has been created beforehand
