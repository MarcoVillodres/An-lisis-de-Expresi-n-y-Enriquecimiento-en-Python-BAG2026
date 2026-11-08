##Detener el script  en caso de error
set -e
echo "1-. Activación del entorno en carpeta"
##Al abrir el archivo del main se inicia en una terminal sin conda por lo que es necesario especificar conda con source e iniciarlo en la termianl nueva creada con main
source $(conda info --base)/etc/profile.d/conda.sh
##Una vez en conda, se activa el entorno
conda activate ./.venv

echo "2-. Ejecutando script"

python src/01_run_deseq2.py

echo "3-. Ejecución realizada con éxito"