
### Código ###
from Bio import Entrez
from Bio import SeqIO

Entrez.email="mmayen@lcg.unam.mx"

# Buscamos A-beta lactamase en Gene y obtenemos la lista de ids de los genes encontrados #

# Se realiza la búsqueda para saber en qué bases de datos se encuentran resultados
termino = "B beta-lactamase"
organismo = "Acinetobacter baumannii"


### Búsqueda del término ingresado en la base de datos 'gene' y recuperación de la lista de IDs ###
handle = Entrez.egquery(term=termino)
record = Entrez.read(handle)
handle.close()

# Iteramos sobre las bases de datos regresadas por la primera búsqueda
for resultados in record['eGQueryResult']:
  # Si existen resultados en la base de datos 'gene' guardamos los IDs
  if (resultados['Count'] != '0') and (resultados['Status'] == 'Ok') and (resultados['DbName']=='gene'): 
    handle2 = Entrez.esearch(db='gene', term=termino)
    record2 = Entrez.read(handle2)
    IDlist = record2['IdList']
    handle2.close()
print(IDlist)
