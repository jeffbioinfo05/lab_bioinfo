#!/bin/bash
# Executa o pipeline do projeto no diretório atual.
# Lê pipeline_type de config.yaml e chama o Snakefile correto.
# Uso: run_pipeline.sh [flags extras do snakemake]
# Exemplos:
#   run_pipeline.sh
#   run_pipeline.sh --dry-run
#   run_pipeline.sh --cores 16
#   run_pipeline.sh --rerun-incomplete

set -e

if [ ! -f "config.yaml" ]; then
    echo "ERRO: config.yaml não encontrado."
    echo "Execute este script a partir do diretório do projeto."
    exit 1
fi

PIPELINE_TYPE=$(python3 -c "import yaml; c=yaml.safe_load(open('config.yaml')); print(c['pipeline_type'])")
PROJECT_NAME=$(python3  -c "import yaml; c=yaml.safe_load(open('config.yaml')); print(c.get('project_name',''))")
SNAKEFILE="$HOME/lab/pipelines/$PIPELINE_TYPE/Snakefile"

if [ ! -f "$SNAKEFILE" ]; then
    echo "ERRO: Pipeline '$PIPELINE_TYPE' não encontrado."
    echo "Esperado em: $SNAKEFILE"
    echo "Pipelines disponíveis:"
    ls "$HOME/lab/pipelines/"
    exit 1
fi

echo "========================================"
echo "Pipeline : $PIPELINE_TYPE"
echo "Projeto  : $PROJECT_NAME"
echo "Snakefile: $SNAKEFILE"
echo "========================================"
echo ""

snakemake \
    --snakefile "$SNAKEFILE" \
    --configfile config.yaml \
    --profile "$HOME/lab/profiles/default" \
    "$@"
