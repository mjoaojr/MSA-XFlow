import sys
import json
import re
from pathlib import Path
from datetime import datetime

# Verifica se o diretório foi passado como argumento
if len(sys.argv) < 2:
    print("Uso correto: python tempos.py <caminho_do_diretorio>")
    sys.exit(1)

# Captura o caminho do diretório
diretorio = Path(sys.argv[1])

# Verifica se o diretório realmente existe
if not diretorio.is_dir():
        print(f"Erro: O diretório '{diretorio}' não existe.")
        sys.exit(1)

# Busca apenas arquivos que terminam com .date
arquivos = list(diretorio.glob('*.date'))

mintime = datetime(2980, 1, 1)
tempos = []
vet = []
# Percorre todos os arquivos do diretório
for arquivo in arquivos:
        if arquivo.is_file():  # Garante que é um arquivo, não uma pasta
#            print(f"\nLendo o arquivo: {arquivo.name}")
            try:
                # Abre e lê o arquivo (com tratamento para evitar erros de codificação)
                with open(arquivo, 'r', encoding='utf-8') as f:
                    job = f.readline().strip()
                    starttime_s = f.readline().strip()
                    if job == "lcs":
                        print (f"start = {starttime_s}")

                    starttime_s = re.sub(r'(\d),(\d)', r'\1.\2', starttime_s)
                    starttime_s = re.sub(r'-03:00', '', starttime_s)
                    starttime_s = starttime_s if len(starttime_s) <= 25 else starttime_s[:25]
                    endtime_s = f.readline().strip()
                    endtime_s = re.sub(r'(\d),(\d)', r'\1.\2', endtime_s)
                    endtime_s = re.sub(r'-03:00', '', endtime_s)
                    endtime_s = endtime_s if len(endtime_s) <= 25 else endtime_s[:25]

                    starttime = datetime.strptime(starttime_s,"%Y-%m-%dT%H:%M:%S.%f")
                    mintime = starttime if starttime < mintime else mintime
                    endtime = datetime.strptime(endtime_s,"%Y-%m-%dT%H:%M:%S.%f")

                    vet.append({"job":job, "starttime":starttime_s, "endtime":endtime_s})
            except Exception as e:
                print(f"Não foi possível ler o arquivo. Erro: {e}")


for it in vet:
#    print (datetime.strptime(it['starttime'], "%Y-%m-%dT%H:%M:%S.%f") - mintime)
    tempos.append({"job":it['job'], "starttime":(datetime.strptime(it['starttime'], "%Y-%m-%dT%H:%M:%S.%f") - mintime).total_seconds(), "endtime":(datetime.strptime(it['endtime'], "%Y-%m-%dT%H:%M:%S.%f") - mintime).total_seconds()})

#print (mintime)
#ordenado = sorted(tempos, key=lambda x: x['starttime'])
ordenado = tempos
#outputs = json.dumps(ordenado, indent=4)
outputs = json.dumps(ordenado)
print(outputs)
