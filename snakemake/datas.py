import json
import sys

# Abre e carrega o arquivo JSON
with open(sys.argv[1], 'r', encoding='utf-8') as arquivo:
    dados = json.load(arquivo)

# Acessa a chave 'starttime'
start_time = dados.get('starttime')
end_time = dados.get('endtime')
rule = dados.get('rule')
arq = dados.get('input')[0]
if rule == "align":
    mode = arq.split('.')[1]
    gt = arq.split('.')[2]
else:
    if rule == "aligncons" or rule == "alignconsupgma" or rule == "alignconsslmin":
        arq = dados.get('input')[1]
        mode = arq.split('.')[1]
        gt = arq.split('.')[2]
    else:
        if rule == "nj" or rule == "upgma" or rule == "slmin" or rule == "slmax":
            arq = dados.get('input')[1]
            mode = arq.split('.')[1]
            gt = "-"
        else:
            mode = "-"
            gt = "-"

#print("rule;mode;gt;starttime;endtime")
print(f"{rule};{mode};{gt};{start_time};{end_time}")
#if len(rule) > 8:
#    print(f"{rule}\t{mode}\t{gt}\t{start_time}\t{end_time}")
#else:
#    print(f"{rule}\t\t{mode}\t{gt}\t{start_time}\t{end_time}")
