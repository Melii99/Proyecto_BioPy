# README

Nuestro proyecto tiene como objetivo Genes asociados a la resistencia contra antibióticos: Beta-lactámicos tipo A, B, C y D



Dentro del repositorio se encuentran todas las versiónes del código, incluyendo la versión final

        Dentro de la versión final se encuentran varias funciones: crear_termino, crear_idslist, 
        archivo_fastas, crear_path y crear_termkey. Cada una de ellas es necesaria para 
        completar un paso del pipline general del código.
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

Además se encuentra una carpeta para las gráficas generadas y los archivos FASTA generados.

