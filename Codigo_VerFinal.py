"""
NAME
       Codigo_VerFinal.py
VERSION
        [4.0]
AUTHOR
        Melissa Mayén Quiroz
DESCRIPTION
        Dentro del programa se encuentran varias funciones: 
        crear_termino, crear_idslist, archivo_fastas, crear_path y crear_termkey. Cada 
        una de ellas es necesaria para completar un paso del pipline general del código.
        Primero: se ingresó un diccionario que contiene los términos de nuestra búsqueda
        (en este caso, se buscaron genes asociados a 4 tipos de beta-lactamasas en 5 especies
        de basterias; el diccionario es modificable para realizar una búsqueda diferente).
        Además agregamos el e-mail que se usará para las búsquedas con E-utilities (en este
        caso se usó el de la autora: mmayen@lcg.unam.mx) y una ruta para las descargas de los
        archivos fasta creados (en este caso se usó la ruta del equipo de la autora:
        'C:\\Users\\Melissa\\Downloads')
        Segundo: con la función crear_terminos creamos los términos para nuestra búsqueda así
        como un diccionario vacío para guardar el conteo de los genes correspondientes.
        Tercero: iteramos sobre cáda término y realizamos la búsqueda de genes asociados al
        mismo con egquery y e search obteniendo una lista de IDs (correspondientes a 'gene').
        Cuarto: dentro del mismo bucle realizamos el conteo del numero de genes asociados al
        término correspondiente y se agregan al diccionario term_count (término:numero de genes).
        Quinto: dentro del mismo bucle tomamos la lista de IDs de genes, buscamos su nombre y
        secuencia y creamos un archivo multifasta de todos los genes encontrados asociados al 
        término buscado (para cada término se crea un archivo diferente y con la función). *En
        este paso se producen algunos errores de conexión en algunos casos por lo que dentro de
        éste código se encuentra en la línea 228 como comentario.
        Sexto: creamos una gráfica general correspondiente al número de genes asociados a cada
        término.
        Séptino: clasificamos el numero de genes encontrados por tipo de beta-lactamasa y creamos
        4 gráficas que indican el número de genes por cada especie. Cada gráfica corresponde a 
        un tipo de beta-lactamasa (A, B, C, D).
        Octavo: clasificamos el numero de genes encontrados por especie y creamos 5 gráficas que 
        indican el número de genes por tipo de beta-lactamasa. Cada gráfica corresponde a una
        especie bacteriana diferente.
        Noveno: Creamos una gráfica del numero total de genes encontrados por tipo de b-lactamasa.
        Décimo: Creamos una gráfica del numero total de genes encontrados por especie.
        Los archivos generados (FASTA y gráficas) se encuentran dentro del repositorio del proyecto:
        https://github.com/Melii99/Python_class_2022/tree/main/Data/archivo_articulos 
  
        *Para el correcto funcionamiento del programa es necesario tener instalado Biopython
        *El e-mail así como el diccionario proporcionado es modificable para la parte de 'código'
        (no para la parte de 'gráficas')
        *Para la generación de archivos la ruta (path) escrita en este programa no funcionará
        en otros equipos por lo que deberá modificarse.
        *Dentro de las funciones crear_idlist y archivos_fasta se realiza importación de librerías:
        Entrez y SeqIO
        *Dentro de la parte de 'gráficas' se realiza importación de librerías: matplotlib.pyplot
INPUT
OUTPUT
        - Archivos FASTA
        - Gráficas
EXAMPLES
        Input
        Output
        >NC_000913.3:c663963-662752 Escherichia coli str. K-12 substr. MG1655, complete genome
        ATGAATACCATTTTTTCCGCTCGTATCATGAAGCGCCTGGCGCTCACCACGGCTCTTTGCACAGCCTTTA
        TCTCTGCTGCACATGCCGATGACCTGAATATCAAAACTATGATCCCGGGTGTACCGCAGATCGATGCGGA
        GTCCTACATCCTGATTGACTATAACTCCGGCAAAGTGCTCGCCGAACAGAACGCAGATGTCCGCCGCGAT
        CCTGCCAGCCTGACCAAAATGATGACCAGTTACGTTATCGGCCAGGCAATGAAAGCCGGTAAATTTAAAG
        AAACTGATTTAGTCACTATCGGCAACGACGCATGGGCCACCGGTAACCCGGTGTTTAAAGGTTCTTCGCT
        GATGTTCCTCAAACCGGGCATGCAGGTTCCGGTTTCTCAGCTGATCCGCGGTATTAACCTGCAATCGGGT
        AACGATGCTTGTGTCGCCATGGCCGATTTTGCCGCTGGTAGCCAGGACGCTTTTGTTGGCTTGATGAACA
        ...
        (Gráficas)
GITHUB 
        https://github.com/Melii99/Proyecto_BioPy/blob/main/Codigo_VerFinal.py
"""

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

  # Creamos una lista vacía 
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
      lista = record2['IdList']
      handle2.close()


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

  # Hacemos la búsqueda de cada término 
  lista = crear_idlist(term, mail)
  ID_list = lista

  # Creamos el término necesario para usar como key en el diccionario
  termkey = crear_termkey(term)

  # Agregamos items al diccionario term_count
  term_count[termkey] = len(ID_list)

  # Cambiamos la ruta de descarga para cada archivo (nombre del archivo)
  path_archivo = crear_path(term, path)
 
  # Hacemos la busqueda de la lista de IDs en la base de datos 'gene' y
  # Cremos el archivo multifasta correspondiente a nuestro término buscado
  #archivo_fastas(ID_list, mail, path_archivo)

### Gráficas ###

# Importación de librerías
import matplotlib.pyplot as plt

# Creamos 2 listas vacías para agregar las especies y el tipo de beta-lactamasa que contienen
# y el numero que se encontró en cada una
species = []
number = []

for key,value in term_count.items():
  species.append(key)
  number.append(value)

# Grafica de numero de genes asociados a un cierto tipo de beta-lactamasa en cada especie #
plt.plot(number,species,'^k:')
plt.title("Numero de genes asociados a cada término")
plt.show()

# Clasificamos el numero de genes encontrados por tipo de beta-lactamasa #
A_num = []
B_num = []
C_num = []
D_num = []
A_sp = []
B_sp = []
C_sp = []
D_sp = []

for i in range(0,len(species)):
  splitstr = species[i].split('_')

  if splitstr[0] == 'A':
    A_num.append(number[i])
    A_sp.append(splitstr[1])
  if splitstr[0] == 'B':
    B_num.append(number[i])
    B_sp.append(splitstr[1])
  if splitstr[0] == 'C':
    C_num.append(number[i])
    C_sp.append(splitstr[1])
  if splitstr[0] == 'D':
    D_num.append(number[i])
    D_sp.append(splitstr[1]) 

# Graficamos por tipo de beta-lactamasa #

# Gráfica de numero de genes asociados a A beta-lactamasa por especie
plt.plot(A_num, A_sp,'or')
plt.title("A beta-lactamase")
plt.show()

# Gráfica de numero de genes asociados a B beta-lactamasa por especie
plt.plot(B_num, B_sp,'or')
plt.title("B beta-lactamase")
plt.show()

# Gráfica de numero de genes asociados a B beta-lactamasa por especie
plt.plot(C_num, C_sp,'or')
plt.title("C beta-lactamase")
plt.show()

# Gráfica de numero de genes asociados a B beta-lactamasa por especie
plt.plot(D_num, D_sp,'or')
plt.title("D beta-lactamase")
plt.show()

# Clasificamos el número total de genes por tipo de beta_lactamasa 
tipo_lac = ['A','B','C','D']
genes_tipolac = [sum(A_num), sum(B_num), sum(C_num), sum(D_num)]

# Clasificamos el numero de genes encontrados por especie #
ec_num = []
kp_num = []
pv_num = []
ab_num = []
pa_num = []
ec_bl = []
kp_bl = []
pv_bl = []
ab_bl = []
pa_bl = []

for i in range(0,len(species)):
  splitstr = species[i].split('_')

  if splitstr[1] == 'Escherichia coli':
    ec_num.append(number[i])
    ec_bl.append(splitstr[0])
  if splitstr[1] == 'Klebsiella pneumoniae':
    kp_num.append(number[i])
    kp_bl.append(splitstr[0])
  if splitstr[1] == 'Proteus vulgaris':
    pv_num.append(number[i])
    pv_bl.append(splitstr[0])
  if splitstr[1] == 'Acinetobacter baumannii':
    ab_num.append(number[i])
    ab_bl.append(splitstr[0]) 
  if splitstr[1] == 'Pseudomonas aeruginosa':
    pa_num.append(number[i])
    pa_bl.append(splitstr[0]) 

# Graficamos por especie #

# Gráfica de numero de genes de E. coli 
plt.plot(ec_bl, ec_num,'ob')
plt.title('Escherichia coli')
plt.show()

# Gráfica de numero de genes de K. pneumoniae 
plt.plot(kp_bl, kp_num,'ob')
plt.title('Klebsiella pneumoniae')
plt.show()
# Gráfica de numero de genes de P. vulgaris 
plt.plot(pv_bl, pv_num,'ob')
plt.title('Proteus vulgaris')
plt.show()

# Gráfica de numero de genes de A. baumannii 
plt.plot(ab_bl, ab_num,'ob')
plt.title('Acinetobacter baumannii')
plt.show()

# Gráfica de numero de genes de P. aeuruginosa 
plt.plot(pa_bl, pa_num,'ob')
plt.title('Pseudomonas aeruginosa')
plt.show()

# Clasificamos el número total de genes por tipo de beta_lactamasa #
tipo_lac = ['A','B','C','D']
genes_tipolac = [sum(A_num), sum(B_num), sum(C_num), sum(D_num)]

# Gráfica de numero de genes asociados a cada tipo de beta-lactamasa
plt.bar(tipo_lac, genes_tipolac, color=['green', 'green', 'green', 'green'])
plt.title("Número total de genes por tipo de beta-lactamasa")
plt.show()

# Clasificamos el número total de genes por especie #
especie = ['E.coli','K.pneumoniae','P.vulgaris','A.baumannii', 'P.aeruginosa']
genes_especie = [sum(ec_num), sum(kp_num), sum(pv_num), sum(ab_num), sum(pa_num)]

# Gráfica de numero de genes asociados a cada tipo de beta-lactamasa
plt.bar(especie,genes_especie, color=['cyan', 'cyan', 'cyan', 'cyan','cyan'])
plt.title("Número total de genes por especie")
plt.show()