'''
La función crear término recibe un diccionario con los términos a buscar (keys)
y los organismos en los que se deseabuscar dicho término (values), y regresa una
lista de términos para su búsqueda con  
'''

def crear_termino(info):
  # Craemos una lista vacia para guardar los términos creados
  terminos = []
  # Iteramos sobre cada termino y sus organismos a buscar a buscar
  for (term, organismos) in info.items():
    organismo = str(organismos).split(',')
    # Iteramos sobre los organismos a buscar en un mismo término
    for o in organismo:
      # Escribimos el termino buscado y el organismo en el que se desea buscar 
      cadena = '(' + term + ' AND ' + o + '[Orgn])'
      # Agregamos cada termino completo a la lista de términos
      terminos.append(cadena)

  return(terminos)


informacion = {'A beta-lactamase':'Acinetobacter baumannii,Pseudomonas aeruginosa',
        'B beta-lactamase': 'Acinetobacter baumannii,Pseudomonas aeruginosa'}

crear_termino(informacion)


'''
La función id_list toma un término a buscar en distintas bases de datos y un mail
de contacto para realizar la búsqueda con E-utilities, y regresa una lista de IDs
corresponientes a la base de datos 'gene' asociados al término buscado.
'''

def id_list(termino, mail):
  # Importamos librerías 
  from Bio import Entrez
  from Bio import SeqIO
  # Agregamos un e-mail de contacto
  Entrez.email = mail

  # Realizamos la busqueda con egquerry del término ingresado
  handle = Entrez.egquery(term=termino)
  record = Entrez.read(handle)
  handle.close()

  # Iteramos sobre las bases de datos regresadas por la primera búsqueda
  for resultados in record['eGQueryResult']:
    # Si existen resultados en la base de datos 'gene' guardamos los IDs
    if (resultados['Count'] != '0') and (resultados['Status'] == 'Ok') and (resultados['DbName']=='gene'): 
      # Realizamos la búsqueda de nuestro término en la base de datos 'gene' y obtenemos los IDs
      handle2 = Entrez.esearch(db='gene', term=termino)
      record2 = Entrez.read(handle2)
      IDlist = record2['IdList']
      handle2.close()

  return(IDlist)

mail = 'mmayen@lcg.unam.mx'
term = '(A beta-lactamase AND Acinetobacter baumannii[Orgn])'

id_list(term, mail)


### Búsqueda de los genes correspondientes a la lista de IDs y creación del archivo multifasta #


#Abrimos un archivo en modo escritura 
path = "C:\\Users\\Melissa\\Downloads\\archivo_fastas"
archivo = open(path, "w")

# Creamos 2 listas para guardar los IDs y las secuencias de los genes encontrados
Gene_names = []
Gene_sequences = []

# Busqueda de los nombres y secuencias de la lista de IDs
for id in IDlist:
  print(id)
  handle = Entrez.efetch(db='gene', id=id, rettype='gb', retmode='text')
  # leemos archivo genebank
  record = SeqIO.read(handle,'genbank')
  handle.close()

  # Se guardan los nombres de los genes
  Gene_names.append(record.name)
  # Se guardan las secuencias de los genes
  Gene_sequences.append(record.seq)

# Se escriben los ids y las secuencias de la búsqueda de la lista de IDs el archivo multifasta 
for i in range(0, len(Gene_names)):
  archivo_fastas.write("> " + Gene_names[i] + "\n")
  archivo_fastas.write(Gene_sequences[i] + "\n\n")

archivo_calidad.close()
