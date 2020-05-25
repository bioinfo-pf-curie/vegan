#!/bin/bash
#set -eu
#
#  Copyright (c) 2019-2020, U900, Institut Curie
#  Copyright (c) 2019-2020, Philippe La Rosa (Equipe HPC)
#
# run_pipeline.sh : This file is part of NGS=>Nextflow Pipeline
#
# Objectifs :
# Generateur d'architecture et d'environnement (espace d'exec) et lancement du pipeline via le WorkFlow
# Manager NextFlow. 
#
# Effectue la creation des repertoires et fichiers qui sont necessaires au lancement du pipeline via le
# script nextflow sur le cluster de calcul ou en local, pour les differentes étapes du pipeline. 
#
# Contexte d'usages :
# Dans le cadre d'un deploiement dans l'espace cible d'installation canonique selon implementation.
# Ce script sera lancé à partir d'un des noued de soummission par n'importe quel user ayant d'acces aux
# data (fastq.gz, references, annotations, ...)
#
# Parametres Entree :
# -p : le nom du projet (exemple : NF_VEGAN)
# -r : le nom du run (exemple : challengeSet1 )
# -s : params du pipemine (exemple : )         
# -i : chemin complet du fichier target.bed
# -u : profiles (exemple : singularity,cluster)         
# -g : genomes (exemple : hg19_base)         
# -b : chemin complet du repertoire genomesbase          
# -o : racine chemin du repertoire ou seront sauvegarde l'ensemble des
#      fichiers pour toutes les etapes : d'integration initiale et d'analyses 
# -q : nom de la file (exemple : dev) 
# -m : chemin complet du repertoire resultats des metriques 
# [-h] : usage
#
# Sorties :
# Soumission au cluster via qsub, ou en local  :
# affichage console (stdout) exemple pour le run TEST du projet NF_VEGAN :
#### scripts/run_sarek.sh for projet NF_VEGAN-2.3 env: and run:TEST by plarosa
#1) create_env
#2) create_nxf_work_dir and configurations files
#3) install_nxf_script
#4) nxf_pipeline_nfcore-sarek
#WORK_DIR = /data/tmp/NGS_RUN_TEMP/NF_VEGAN-2.3_TEST_1571140311631
#RACINE_PIPELINES_DIR = /bioinfo/users/plarosa/projects/GIT/hackathon_intel_genci/script/pipelines/Sarek/ScLifeLab
#CONFIG_NXF_PATH = /data/tmp/NGS_RUN_TEMP/NF_VEGAN-2.3_TEST_1571140311631/conf
#LOCAL_SCRIPTS (lancemment et enchainement automatique sur le cluster): 
#SCRIPT_STEP_MAPPING = /data/tmp/NGS_RUN_TEMP/NF_VEGAN-2.3_TEST_1571140311631/run_mapping.sh
#SCRIPT_STEP_ANNOTATE = /data/tmp/NGS_RUN_TEMP/NF_VEGAN-2.3_TEST_1571140311631/run_annotate.sh

NF_VEGAN_DEBUG=${NF_VEGAN_DEBUG:=0}; [[ "$NF_VEGAN_DEBUG" == 'x' ]] && set -x

if [[ $TERM && $TERM != 'dumb' ]]
  then
    if command -v tput &>/dev/null
      then
        GREEN=$(tput setaf 2; tput bold)
        YELLOW=$(tput setaf 3)
        RED=$(tput setaf 1)
        NORMAL=$(tput sgr0)
    fi
fi

function echo_red() {
    >&2 echo -e "$RED$*$NORMAL"
}

function echo_green() {
    echo -e "$GREEN$*$NORMAL"
}

function echo_yellow() {
    >&2 echo -e "$YELLOW$*$NORMAL"
}

function die() {
  echo_red "$*"
  exit 1
}

### usage ###
function usage (){
    echo -e "\nUsage: $0"
    echo -e "\n [Options]"    
    echo -e "\t-p : nom du projet (exemple : NF_VEGAN)"     
    echo -e "\t-r : nom du run (exemple : 10X)"  
    echo -e "\t-s : chemin complet sampleplan"         
    echo -e "\t-u : profiles (exemple : singularityPath,cluster)"         
    echo -e "\t-g : genomes (exemple : hg19)"         
    echo -e "\t-b : chemin complet du repertoire genomesbase "         
    echo -e "\t-o : racine chemin du repertoire ou seront sauvegarde l'ensemble des fichiers pour toutes les etapes : d'integration initiale et d'analyses" 
    echo -e "\t-q : nom de la file (exemple : dev)" 
    echo -e "\t-m : chemin complet du repertoire resultats des metriques" 
    echo -e "\t-h : usage"          
    echo -e "\n\n [Example]: \n\t# export SUBMIT=qsub -N NF_VEGAN-2.5.1;export EXECUTOR=pbs;export RACINE_PIPELINES_DIR=$(pwd);bash scripts/run_sarek.sh -p NF_VEGAN-2.5.1 -r BENCH_TOOLS -g GRCh38 -b /data/kdi_prod/.kdi/project_workspace_0/1430/acl/01.00/public_data/annotations/GRCh38 -s --input Sarek-data/HPC-bench/tsv/HPC-bench-test.tsv -u singularityPath,cluster -i /data/kdi_prod/.kdi/project_workspace_0/1430/acl/01.00/public_data/annotations/illumina/TargetRegions.bed  -m /data/kdi_prod/.kdi/project_workspace_0/1430/acl/01.00/results/test/sarek -q dev -o /data/tmp/NGS"
    exit 1
}

((!$#)) && echo "Il n'y a pas d'arguments!!" && usage

if [[ ($# < 17) || ($# > 18) ]]
then
    echo "Nombre d'arguments incorrect ($#) !!"
    usage
fi

while getopts "p:r:s:b:g:o:u:q:m:h" optionName; do
case "$optionName" in

p) project_name="$OPTARG";;
r) run="$OPTARG";;
s) sampledata="$OPTARG";;
u) profile="$OPTARG";;
g) genome="$OPTARG";;
b) genomebase="$OPTARG";;
o) racine_work_dir="$OPTARG";;
q) queue="$OPTARG";;
m) results_dir="$OPTARG";;
h) usage;;
*) usage;;
esac
done

function pipeline_create_env() {
   echo " 1) create_env"
   # path par defaut du repertoire racine d'exec du pipeline
   RACINE_WORK_DIR_DEF="/data/tmp/NGS_RUN_TEMP"
   [[ ! $racine_work_dir ]] && racine_work_dir=${RACINE_WORK_DIR_DEF}
   WORK_DIR=${WORK_DIR:="${racine_work_dir}/${project_name}_${run}_${name_end}"}
   LOGNAME=${LOGNAME:="inconnu"}
   [[ ! $RACINE_PIPELINES_DIR ]] && RACINE_PIPELINES_DIR="/bioinfo/pipelines/nf-vegan/dev"
   NXF_DIR=./
   NXF_NAME=${NXF_NAME:="${script_nf}"}
   NXF_BIN_DIR=${NXF_BIN_DIR:="/bioinfo/local/build/Centos/nextflow/nextflow-19.04.0.5069"}
   CONFIGS_DIR=${CONFIGS_DIR:="conf"}
   CONFIGS_NXF_DIR=${CONFIGS_NXF_DIR:="nextflow"}
   SCRIPT_STEP_MAPPING=${SCRIPT_STEP_MAPPING:="${WORK_DIR}/run_vegan_mapping.sh"}
   SCRIPT_STEP_ANNOTATE=${SCRIPT_STEP_ANNOTATE:="${WORK_DIR}/run_vegan_annotate.sh"}
   SCRIPT_STEP_COPIE_RESULTS=${SCRIPT_STEP_COPIE_RESULTS:="${WORK_DIR}/run_copie_results.sh"}
   BIN_DIR=${BIN_DIR:="bin"}
   SCRIPTS_DIR=${SCRIPTS_DIR:="scripts"}
   ASSETS_DIR=${ASSETS_DIR:="assets"}
   LOCAL_LOG_DIR=${LOCAL_LOG_DIR:="LOG"}
   PIPELINES_PATH=${PIPELINES_PATH:="/bioinfo/pipelines"}
   COMM="bash"
   SUBMIT=${SUBMIT:="qsub"}
}


function pipeline_create_nxf_work_dir() {
   echo " 2) create_nxf_work_dir and configurations files"
   # creation architecture
   mkdir -p ${WORK_DIR}
   # copie des elements dans l'espace de travail
   cp -r ${RACINE_PIPELINES_DIR}/${CONFIGS_DIR} ${WORK_DIR}
   cp -r ${RACINE_PIPELINES_DIR}/${ASSETS_DIR} ${WORK_DIR}
   cp ${RACINE_PIPELINES_DIR}/nextflow.config ${WORK_DIR}
   # creation scripts de lancement de nextflow
   cat <<NF_VEGAN_MAPPING > ${SCRIPT_STEP_MAPPING}
# NF_VEGAN Mapping, Recalibrate, VariantCalling, Annotate.
cd  ${WORK_DIR}
# lancement du pipeline sur le cluster 
nextflow run main.nf -resume -profile ${profile} ${sampledata} --step mapping --genome ${genome} --genomeAnnotationPath  ${genomebase} -c nextflow.config && ${COMM} ${SCRIPT_STEP_ANNOTATE}
NF_VEGAN_MAPPING

# creation scripts de lancement de nextflow pour etape suivante
   cat <<NF_VEGAN_ANNOTATE > ${SCRIPT_STEP_ANNOTATE}
# NF_VEGAN Annotate.
cd  ${WORK_DIR}
# lancement du pipeline sur le cluster pour etape suivante
nextflow run main.nf -resume -profile ${profile} ${sampledata} --step annotate --genome ${genome} --genomeAnnotationPath  ${genomebase} -c nextflow.config && ${COMM} ${SCRIPT_STEP_COPIE_RESULTS}
NF_VEGAN_ANNOTATE

# creation script de copie des metriques  
   cat <<COPIE > ${SCRIPT_STEP_COPIE_RESULTS}
cd  ${WORK_DIR}/results/
cp -r Reports/MultiQC ${results_dir}/${name_end} 
cp -r pipeline_info ${results_dir}/${name_end} 
COPIE

}

function pipeline_install_nxf_script() {
   echo " 3) install_nxf_script"
   [[ -e ${RACINE_PIPELINES_DIR}/${NXF_NAME} ]] || die "${RACINE_PIPELINES_DIR}/${NXF_NAME} n'est pas accessible"
   cp  ${RACINE_PIPELINES_DIR}/*.nf ${WORK_DIR}
   cp -r ${RACINE_PIPELINES_DIR}/${BIN_DIR} ${WORK_DIR}
   cp -r ${RACINE_PIPELINES_DIR}/${SCRIPTS_DIR} ${WORK_DIR}
}

function nxf_pipeline_launch() {
   echo " 4) nxf_pipeline_nf-vegan"
   
   echo_green "WORK_DIR = ${WORK_DIR}"
   echo_green "RACINE_PIPELINES_DIR = ${RACINE_PIPELINES_DIR}"
   echo_green "CONFIG_NXF_PATH = ${WORK_DIR}/${CONFIGS_DIR}"
   echo_green "LOCAL_SCRIPTS (lancemment et enchainement automatique sur le cluster): "
   echo_green "SCRIPT_STEP_MAPPING = ${SCRIPT_STEP_MAPPING}"
   echo_green "SI OK => SCRIPT_STEP_ANNOTATE = ${SCRIPT_STEP_ANNOTATE}"
   echo_green "SI OK => SCRIPT_STEP_COPIE_RESULTS = ${SCRIPT_STEP_COPIE_RESULTS}"
   d=$(date)
   cd ${WORK_DIR}
   chmod +x ${SCRIPT_STEP_MAPPING}
   ${COMM} ${SCRIPT_STEP_MAPPING}
}


function nxf_date() {
    local ts=$(date +%s%3N); [[ $ts == *3N ]] && date +%s000 || echo $ts
}

function get_date() {
    local ts=$(date +%d-%m-%y) && echo $ts
}
name_end=$(nxf_date)
date_run=$(get_date)

##
echo_yellow "#### $0 for projet ${project_name} and run:${run} by ${LOGNAME}"
#
# creation des variables d'environnement nécessaires au lancement du pipeline 
pipeline_create_env

#
# creation du repertoire d'exec : doit contenir l'ensemble des composants necessaire au lancement du pipeline
# creation repertoires resultats, configs, script, creation script bash de lancement 
# de nextflow pointant sur les parametres.
pipeline_create_nxf_work_dir

# ajout lien vers le script nxf du pipeline dans WORK_DIR
pipeline_install_nxf_script

#
# lancement de nexflow pour le pipeline pour la première étape qui enchaine la suivante sur succes
nxf_pipeline_launch


