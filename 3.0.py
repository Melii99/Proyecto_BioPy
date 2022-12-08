'''
DOCUMENTACIÓN
'''

### Funciónes ###

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

'''
La función id_list toma un término a buscar en distintas bases de datos y un mail
de contacto para realizar la búsqueda con E-utilities; y regresa una lista de IDs
corresponientes a la base de datos 'gene' asociados al término buscado.
'''
def crear_idlist(termino, mail):
  # Importamos librerías 
  from Bio import Entrez
  from Bio import SeqIO
  # Agregamos un e-mail de contacto
  Entrez.email = mail

  # Craemos una lista vacía 
  lista = []

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
      handle2.close()
      lista = record2['IdList'] 
      print('holi')

  return(lista)
    
'''
La funcion archivo_fastas toma una lista de IDs correspondientes a la base de 
datos 'gene', un mail, y una ruta para la descarga del archivo fasta; y
realiza su búsqueda y obtiene los nombres (IDs) de los genes y su
secuencia, regresando un archivo multifasta correspondiente a los genes buscados
'''
def archivo_fastas(IDlist, mail, ruta):

  # Importamos librerías 
  from Bio import Entrez
  from Bio import SeqIO
  # Agregamos un e-mail de contacto
  Entrez.email = mail

  #Abrimos un archivo en modo escritura con la ruta proporcionada por el usuario
  archivo = open(ruta, "w")

  # Creamos 2 listas para guardar los IDs y las secuencias de los genes encontrados
  Gene_names = []
  Gene_sequences = []

  # Iteramos sobre cada ID de la lista de IDs
  for id in IDlist:
    # Realizamos la busqueda de cada ID en la base de datos 'gene'
    handle = Entrez.efetch(db='gene', id=id, rettype='gb', retmode='text')
    record = SeqIO.read(handle,'genbank')
    handle.close()
    # Se guardan los nombres de los genes
    Gene_names.append(record.name)
    # Se guardan las secuencias de los genes
    Gene_sequences.append(record.seq)

  # Se escriben los IDs y las secuencias recuperadas en el archivo 
  for i in range(0, len(Gene_names)):
    archivo_fastas.write("> " + Gene_names[i] + "\n")
    archivo_fastas.write(Gene_sequences[i] + "\n\n")
  archivo_fastas.close()

'''
La función crear path toma un término y una ruta para la descarga de archivos y 
crea una nueva ruta para la descarga de un archivo específico correspondiente al 
término ingresado
'''
def crear_path(termino, path):
  cadena = termino.replace(' AND ', '_')
  cadena = cadena.replace('[Orgn]', '')
  cadena = cadena.replace('(', '')
  cadena = cadena.replace(')', '')
  cadena = path + '\\' + cadena
  return(cadena)

'''
La función crear_termkey toma un término (usado para la búsqueda) y regresa un 
nuevo término que se usará como key del diccionario term_count
'''
def crear_termkey(termino):
  cadena = termino.replace(' AND ', '_')
  cadena = cadena.replace('[Orgn]', '')
  cadena = cadena.replace(' beta-lactamase', '')
  cadena = cadena.replace('(', '')
  cadena = cadena.replace(')', '')
  return(cadena)


### Búsqueda de términos ###

# Ingresamos nuestro diccionario con los términos y organismos que deseamos buscar
informacion = {'A beta-lactamase':'Escherichia coli,Klebsiella pneumoniae,Proteus vulgaris,Acinetobacter baumannii,Pseudomonas aeruginosa',
        'B beta-lactamase':'Escherichia coli,Klebsiella pneumoniae,Proteus vulgaris,Acinetobacter baumannii,Pseudomonas aeruginosa',
        'C beta-lactamase':'Escherichia coli,Klebsiella pneumoniae,Proteus vulgaris,Acinetobacter baumannii,Pseudomonas aeruginosa',
        'D beta-lactamase':'Escherichia coli,Klebsiella pneumoniae,Proteus vulgaris,Acinetobacter baumannii,Pseudomonas aeruginosa'}

# Agregamos el e-mail que se usará para las búsquedas con E-utilities
mail = 'mmayen@lcg.unam.mx'
# Agregamos el path donde se agregarán los archivos fasta de nuestras búsquedas 
# (los nombres de cada archivo serán diferentes acorde al término buscado)
path = 'C:\\Users\\Melissa\\Downloads'

# Creamos los términos que vamos a necesitar para nuestra búsqueda
terminos = crear_termino(informacion)

# Creamos un diccionario vacío para agregar la cuenta de genes encontrados en cada término
term_count = {}

# Iteramos sobre todos los términos a buscar
for term in terminos:
  
  # Cambiamos la ruta de descarga para cada archivo (nombre del archivo)
  path_archivo = crear_path(term, path)

  # Hacemos la búsqueda de cada término 
  ID_list = crear_idlist(term, mail)

  # Creamos el término necesario para usar como key en el diccionario
  termkey = crear_termkey(term)

  # Agregamos items al diccionario term_count
  term_count[termkey] = len(ID_list)
 
  # Hacemos la busqueda de la lista de IDs en la base de datos 'gene' y
  # Cremos el archivo multifasta correspondiente a nuestro término buscado
  archivo_fastas(ID_list, mail, path_archivo)

### Gráficas ###

# Importación de librerías
import matplotlib.pyplot as plt

# Creamos 2 listas vacías para agregar las especies y el tipo de beta-lactamasa que contienen
# y el numero que se encontró en cada una
species = []
number = []

for key,value in dicc.items():
  species.append(key)
  number.append(value)

# Graficamos el numero de genes asociados a un cierto tipo de beta-lactamasa en cada especie
plt.plot(number,species,'^k:')


