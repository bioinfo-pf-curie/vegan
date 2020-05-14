# NF_VEGAN : steps
# Preprocessing – main.nf (based on GATK best practices)
# Germline variant calling – germlineVC.nf
# Somatic variant calling – somaticVC.nf (optional)
# Annotation – annotate.nf (optional)
# Reporting (multiqc)


# exemple de lancement (a modifier avec vos path) du script qui effectue les etapes : de creation du repertoire d'exec
# qui contient l'ensemble des composants necessaires au lancement du pipeline,
# creation des variables d'environnements, des repertoires resultats, configs, script, 
# creation des scripts bash de lancement et d'enchainements des etapes du pipeline.
#
# Update with your nextflow or specifics tools installation
 export PATH=/mnt/beegfs/common/apps/graphviz/graphviz-2.44.0/bin:$PATH
###  ABACUS
# (peut facilement être utilisé sur un autre cluster : calcsub, ...)
# modes de lancement : automatique 
# singularity sur le cluster (mode nominal) 
export SUBMIT="srun";export RACINE_PIPELINES_DIR=$(pwd);bash ${RACINE_PIPELINES_DIR}/scripts/run_pipeline.sh -p NF_VEGAN-1.0 -r DREAM_insilico_1 -g hg19_base -b /mnt/beegfs/data/annotations/pipelines/nf-vegan/hg19_base -s "--maxMemory '120.GB' --input /mnt/beegfs/data/dataset/EUCANCan/dream/tsv/dream-fastqs-set1_spit.tsv --singularityImagePath /mnt/beegfs/data/containers/singularity/nf-vegan/images --annotation_cache true --snpeff_cache /mnt/beegfs/data/annotations/pipelines/nf-vegan/hg19_base/databases/snpEff_v4_3/hg19 --tools manta,mutect2,haplotypecaller,ascat,snpeff --skipQC versions --params.annotateTools manta,mutect2" -u singularity,cluster -m /mnt/beegfs/data/metrics/nf-vegan/results -q dev -o /mnt/beegfs/data/tmp/NGS_RUN_TEMP 2>&1 > /mnt/beegfs/data/tmp/logs/NF_VEGAN-1.0_$$.log &

