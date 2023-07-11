
<img align="right" src="https://visitor-badge.laobi.icu/badge?page_id=zumrudu-anka.zumrudu-anka">

<h1 align="center">
  <a href="https://git.io/typing-svg">
    <img src="https://readme-typing-svg.herokuapp.com/?lines=Hello,+There!+;This+is+Angel+Loayza...;Nice+to+meet+you!&center=true&size=30">
  </a>
</h1>

> Antes de, el repositorio a GitHub se encuentra en el siguiente link: https://github.com/angel452/BallStarTree . Ahora si, continuemos :)

# 💫 Instrucciones:
- Una vez clonado el repositorio, o descargado la carpeta, ejecutar el comando
```bash
cmake CMakeLists.txt
```
> Esto les creara algunos archivos necesarios para la compilacion del proyecto y tambien creara una carpeta llamada "src/", aqui estaran los archivos .cpp y .h que contienen el codigo fuente del proyecto.
- Luego, compilar el proyecto con el comando
```bash
make
```
> Esto crara un ejecutable llamado "BallStarTree"
- Finalmente, ejecutar el programa con el comando
```bash
./BallStarTree
```
# 📊 Sobre el proyecto:
- Lectura de archivos .csv
> Si se desea colocar un nuevo dataSet, se debe editar la linea 594 del archivo main.cpp.

> Necesariamente, debemos mantener el mismo formato de los archivos .csv para evitar problemas. Es decir: estar separado por comas y tener en las 2 ultimas columnas un id y una nombre. El numero de filas no importa.

>  Por ultimo, si cambiamos el dataSet, debemos verificar la linea 575, pues aqui controlamos la candidad de dimensiones que tendran los puntos, nodos, etc. Para este dataSet, se tienen 15 datos por columna, por lo que en el codigo fuente se debe colocar 14 (linea 575). La razon es que la ultima columna (el nombre) la guardamos en otro vector.

- Creacion del arbol
> Se puede cambiar el numero maximo de puntos/datos por nodo, en la linea 574 del archivo main.cpp. Esto controlara hasta cuanto el arbol se dividira hacia abajo.

- Knn Vecinos mas cercanos
> Tambien podemos controlar cuantos vecinos cercanos queremos ver. Esto se puede editar en la linea 645 del archivo main.cpp.

> **IMPORTANTE:** Por la misma logica usada para implementar esta funcion, el programa esta limitado a solo mostrar un numero menor de vecinos respecto a M (linea 574). Es decir; si en nuestro M es 500, solo podremos mostrar 499 vecinos cercanos o menos (< M), si le tratamos de pasar un numero mas alto, el programa abortara su ejecucion.

- Muestra de resultados
> Tras la ejecucion del programa, en la misma terminal se mostraran 2 cuadros con el tiempo en nanosegundos que le tomo indexar (insertar y crear la estructura), y obtener los knn Vecinos mas cercanos.
