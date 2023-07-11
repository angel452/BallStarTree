// ######################## MAIN LIBRARIES #######################
#include <iostream>
#include <vector>

#include <sstream>
#include <string>
#include <fstream>

#include <chrono>
#include <algorithm>
#include <cmath>

// ######################## PCA LIBRARIES #######################
#include "pca.h"

using namespace std;

vector<string> allNamesSongs;
vector<vector<float>> allKnnDistances;

template <int ndim>
struct Point{
    float puntos[ndim];

    Point(){
        for(int i = 0; i < ndim; i++) {
            puntos[i] = 0;
        }
    }

    Point( float _puntos[] ){
        for(int i = 0; i < ndim; i++){
            puntos[i] = _puntos[i];
        }
    }

    void printPointOfNode(){
        for(int i = 0; i < ndim - 1; i++){
            cout << puntos[i] << "  ";
        }
        cout << " ID: " << puntos[ndim - 1] ;
        //cout << endl;
    }

    void printPoint(){
        for(int i = 0; i < ndim; i++){
            cout << puntos[i] << "  ";
        }
        cout << endl;
    }

    float &operator[](int posicion){
        return puntos[posicion];
    }
};

bool compareFirstColumn(const std::vector<float>& a, const std::vector<float>& b) {
    return a[0] < b[0];
}

template <typename T, int ndim>
struct Node{
    int nPuntos; // numero de puntos en el nodo
    int M; // Maximo de modos por nodo
    bool isLeaf; // si es hoja o no

    vector< Point<ndim> > grupoRecords;
    vector< Node<T, ndim> > grupoHijos;

    Point<ndim-1> puntoCentral;
    float radio;

    Node( int _maxRecords = 0 ){
        nPuntos = 0;
        M = _maxRecords;
        isLeaf = true;

        puntoCentral;
        radio = 0;
    }

    // ################### PRINT FUNCTIONS ###################################
    void printPointsNode( Node<T, ndim> *ptr ){
        for(int i = 0; i < ptr->nPuntos; i++){
            //cout << i+1 << ". " << '\t';
            cout << "( " ;
            ptr->grupoRecords[i].printPointOfNode();
            cout << ")";
        }
        cout << endl;
    }

    // ################### OPERATIONS ###################################
    // Distancia entre un record completo y un punto central
    float getEuclideanDistance(Point<ndim> punto1, Point<ndim-1> punto2){
        float res = 0;
        for(int i = 0; i < ndim - 1; i++){
            res = res + pow( punto2[i]-punto1[i] , 2);
        }
        return sqrt(res);
    }

    // Distancia entre un record y un punto central
    float getEuclideanDistance2(Point<ndim-1> punto1, Point<ndim-1> punto2){
        float res = 0;
        for(int i = 0; i < ndim - 1; i++){
            res = res + pow( punto2[i]-punto1[i] , 2);
        }
        return sqrt(res);
    }

    // Distancia entre un record y un record entero
    float getEuclideanDistance3(Point<ndim-1> punto1, Point<ndim> punto2){
        float res = 0;
        for(int i = 0; i < ndim - 1; i++){
            res = res + pow( punto2[i]-punto1[i] , 2);
        }
        return sqrt(res);
    }

    // ####################### ALGORITMOS DE ORDENAMIENTO ###########################################
    /*
     * - Estos 2 olgoritmos 3 algoritmos a continuacion, fueron una de las primeras maneras de ordenar una matrix proyectada
     * - En la actualidad se usa otro metodo mas eficiente, pero queda como referencia de como se hizo en un principio
     * - Ademas son algoritmos que funcionan a la perfeccion, pero son muy lentos para ordenar una matriz de 35144 x n datos
    */
     int particion(Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> &reprojectPoints, int inicio, int final, Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> &OriginalPoints, Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> &OriginalIds){
        float pivote = reprojectPoints(inicio, 0); // El pivote sera el primer elemento
        int aux = inicio + 1;
        // Recorremos desde el pivote hasta el final
        for(int i = aux; i <= final; i++){
            if(reprojectPoints(i, 0) < pivote){
                // Intercambiamos los valores de la matriz de puntos proyectados
                reprojectPoints.row(aux).swap(reprojectPoints.row(i));
                // Interacambiamos los valores de la matriz de puntos originales
                OriginalPoints.row(aux).swap(OriginalPoints.row(i));
                // Intercambiamos los valores de la matriz de ids originales
                OriginalIds.row(aux).swap(OriginalIds.row(i));
                aux++;
            }
        }
        reprojectPoints.row(inicio).swap(reprojectPoints.row(aux-1));
        OriginalPoints.row(inicio).swap(OriginalPoints.row(aux-1));
        OriginalIds.row(inicio).swap(OriginalIds.row(aux-1));
        return aux-1;
    }

    void quickSort(Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> &reprojectPoints, int inicio, int final, Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> &OriginalPoints, Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> &OriginalIds ){
        if(inicio < final){
            int pivote = particion(reprojectPoints, inicio, final, OriginalPoints, OriginalIds);
            quickSort(reprojectPoints, inicio, pivote-1, OriginalPoints, OriginalIds);
            quickSort(reprojectPoints, pivote+1, final, OriginalPoints, OriginalIds);
        }
    }

    Eigen::MatrixXf sortMatrix(const Eigen::MatrixXf &original_matrix) {

        Eigen::MatrixXf sorted_matrix(original_matrix.rows(), original_matrix.cols());
        Eigen::VectorXf sorted_rows = original_matrix.col(0);
        Eigen::VectorXf::Index min_row;

        for ( int i = 0; i < original_matrix.rows(); i++) {
            sorted_rows.minCoeff(&min_row);
            sorted_matrix.row(i) = original_matrix.row(min_row);
            sorted_rows(min_row) = std::numeric_limits<double>::max();
        }

        return sorted_matrix;
    }

    // ################### MAIN FUNCTIONS ###################################
    void getCircleInfo(Node<T, ndim> *ptr){
        // --> Obtener el punto central
        /*
         * LOGICA:
         * - Primero obtener el maximo y minimo valor en la dimension i.
         * - Luego ubicar el punto central en la mitad entre el minimo y el maximo
         * - Repetir para cada dimension
         */

        for(int i = 0; i < ndim - 1; i++){
            float pntoCentralAux;
            float minPoint = ptr->grupoRecords[0][i]; // Almacenamos el minimo valor del punto en la dimension i
            float maxPoint = ptr->grupoRecords[0][i]; // Almacenamos el maximo valor del punto en la dimension i
            for(int j = 1; j < ptr->nPuntos; j++){
                float minPointAux = ptr->grupoRecords[j][i];
                float maxPointAux = ptr->grupoRecords[j][i];

                if(minPointAux < minPoint){
                    minPoint = minPointAux;
                }
                if(maxPointAux > maxPoint){
                    maxPoint = maxPointAux;
                }
            }
            pntoCentralAux = (maxPoint + minPoint)/2; // Ubicamos la mitad entre el minimo y el maximo
            ptr->puntoCentral[i] = pntoCentralAux;
        }

        // --> Obtener el radio
        /*
         * LOGICA:
         * - Sacar la distancia el record[i] y el punto central
         * - Guardamos el mayor de estos valores
         * - Este valor es el radio
         */

        float maxDistance = getEuclideanDistance( ptr->grupoRecords[0], ptr->puntoCentral);
        for(int i = 1; i < ptr->nPuntos; i++){
            float maxDistanceAux = getEuclideanDistance( ptr->grupoRecords[i],ptr->puntoCentral);
            if(maxDistanceAux > maxDistance){ // Guardamos maxima la distancia entre el record y el punto central
                maxDistance = maxDistanceAux;
            }
        }
        ptr->radio = maxDistance;
    }

    void insertToNode( Point<ndim> record, Node<T, ndim> *ptr ){
        // LOGICA:
        /*
         * Insertamos el record a donde apunta el puntero
         * Aumentamos el numero de puntos
         * Marcamos que ya no es hoja
         */
        ptr->grupoRecords.push_back(record);
        ptr->nPuntos++;
        ptr->isLeaf = false;
    }

    pair<Node, Node> splitPCA(Node<T, ndim> *ptr){
        // LOGICA:
        /*
         * - Guardamos en la matriz pca_data_matrix los puntos SIN el id
         * - Guardamos en la matriz pca_ids_matrix los ids de los puntos
         * - Calculamos PCA, esto me devolvera una matriz con los puntos proyectados
         * - Ordenamos estos puntos de menor a mayor
         * - Insertamos en el hijo 1 los datos desde la posicion 0 hasta n/2
         * - Insertamos en el hijo 2 los datos desde la posicion n/2 hasta n
         */

        pair<Node, Node> resultSplit;

        // --> Variables:
        const int nPuntosAux = ptr->nPuntos;
        Eigen::MatrixXf pca_data_matrix(nPuntosAux, ndim-1); // Guardamos los puntos sin el id
        Eigen::MatrixXf pca_ids_matrix(nPuntosAux, 1); // Guardamos los ids de los puntos

        // --> Insertamos los puntos en una matriz Eigen:
        for(int i = 0; i < nPuntosAux; i++){
            for(int j = 0; j < ndim - 1; j++){
                pca_data_matrix(i, j) = ptr->grupoRecords[i][j];
            }
            // * Insertamos los ids de los puntos en una matriz Eigen:
            pca_ids_matrix(i, 0) = ptr->grupoRecords[i][ndim - 1];
        }

        // --> PCA:
        pca_t<float> pca;
        pca.set_input(pca_data_matrix);
        pca.compute(); // Calculamos PCA

        Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> reprojection = pca.reprojection2();

        // TEST 1. (Usando quicksort)
        // quickSort(reprojection, 0, (ptr->nPuntos)-1, pca_data_matrix, pca_ids_matrix);

        // --> Juntamos las matrices reprojection, pca_data_matrix y pca_ids_matrix en una sola matriz:
        Eigen::MatrixXf combinate_Matrix(ptr->nPuntos, reprojection.cols() + pca_data_matrix.cols() + pca_ids_matrix.cols());
        combinate_Matrix << reprojection, pca_data_matrix, pca_ids_matrix;

        // --> Lo pasamos a un vector de vectores:
        vector<vector<float>> combinate_Matrix_Vector;
        for(int i = 0; i < combinate_Matrix.rows(); i++){
            vector<float> auxVector;
            for(int j = 0; j < combinate_Matrix.row(i).size(); j++){
                auxVector.push_back(combinate_Matrix.row(i)[j]);
            }
            combinate_Matrix_Vector.push_back(auxVector);
        }

        // --> Ordenamos en base a la primera columna:
        sort(combinate_Matrix_Vector.begin(), combinate_Matrix_Vector.end(), compareFirstColumn);

        // --> Separar los puntos en los hijos
        Node<T, ndim> leftChild;
        Node<T, ndim> *ptrLeftChild;
        Node<T, ndim> auxNode1(M);
        leftChild = auxNode1;
        ptrLeftChild = &leftChild;

        Node<T, ndim> rightChild;
        Node<T, ndim> *ptrRightChild;
        Node<T, ndim> auxNode2(M);
        rightChild = auxNode2;
        ptrRightChild = &rightChild;

        // 1. Insertar puntos en el hijo izquierdo
        for(int i = 0; i < (ptr->nPuntos)/2; i++){
            Point<ndim> auxPoint;
            int contadorAux = 0;
            for(int j = ndim-1; j < (ndim-1)+(ndim-1); j++){ // Desde la posicion 13 hasta 25
                auxPoint[contadorAux] = combinate_Matrix_Vector[i][j];
                contadorAux++;
            }
            // * Insertar el id del punto en la ultima posicion
            auxPoint[ndim - 1] = combinate_Matrix_Vector[i][(ndim-1)*2]; // Ultima posicion: (i,26)

            ptrLeftChild->grupoRecords.push_back(auxPoint);
            ptrLeftChild->nPuntos++;
        }

        // 2. Insertar puntos en el hijo derecho
        for(int i = (ptr->nPuntos)/2; i < ptr->nPuntos; i++){
            Point<ndim> auxPoint;
            int contadorAux = 0;
            for(int j = ndim-1; j < (ndim-1)+(ndim-1); j++){ // Desde la posicion 13 hasta 25
                auxPoint[contadorAux] = combinate_Matrix_Vector[i][j]; // Ultima posicion: (i,26)
                contadorAux++;
            }
            // * Insertar el id del punto
            auxPoint[ndim - 1] = combinate_Matrix_Vector[i][(ndim-1)*2];

            ptrRightChild->grupoRecords.push_back(auxPoint);
            ptrRightChild->nPuntos++;
        }

        resultSplit.first = leftChild;
        resultSplit.second = rightChild;

        return resultSplit;
    }

    void createEstructure(Node *ptr){
        // --> Print number of points in node
        //cout << "Points in node... " << ptr->nPuntos << endl << endl;
        //printPointsNode(ptr);

        // 1. Sacamos el radio y el punto central del nodo actual (ptr)
        getCircleInfo(ptr);
        //cout << endl;
        //cout << "El punto central es: "; ptr->puntoCentral.printPoint();
        //cout << "El radio es: " << ptr->radio << endl;
        //cout << "El numero de puntos dentro es: " << ptr->nPuntos << endl;
        //cout << endl;

        // 2. Split
        if((ptr->nPuntos)/2 >= M){
            pair<Node, Node> resultSplit = splitPCA(ptr);

            // --> Insertamos los hijos en el vector de hijos del al nodo actual (ptr)
            ptr->grupoHijos.push_back(resultSplit.first);
            ptr->grupoHijos.push_back(resultSplit.second);

            // --> Creamos punteros a los hijos
            Node<T, ndim> *ptrLeftChild = &ptr->grupoHijos[0];
            Node<T, ndim> *ptrRightChild = &ptr->grupoHijos[1];

            ptrLeftChild->isLeaf = false;
            createEstructure(ptrLeftChild);

            ptrRightChild->isLeaf = false;
            createEstructure(ptrRightChild);
        }
        else{
            // Ya que no se puede dividir mas, consideramos a este nodo un nodo hoja
            ptr->isLeaf = true;
        }
    }

    // ############################## KNN SEARCH ############################
    void getDistancesOfGroup(Node<T, ndim> *ptr, Point <ndim-1> recordBase){
        for(int i = 0; i < ptr->grupoRecords.size(); i++){
            float distance = getEuclideanDistance3(recordBase, ptr->grupoRecords[i]);

            //cout << "La distancia entre el record y el punto " << i << " es: " << distance << endl;

            // --> Guardamos la distancia y su id en "allKnnDistances"
            vector<float> auxVector;
            auxVector.push_back(distance);
            auxVector.push_back(ptr->grupoRecords[i][ndim-1]);
            allKnnDistances.push_back(auxVector);
        }
    }

    void knnSearch( Point<ndim-1> recordBase , Node<T, ndim> *ptr, int Nvecinos){
        // LOGICA:
        /*
         * OJO: Hacemos una busqueda recursiva hasta llegar a un nodo hoja
         * - Empezando en el root, vemos en que hijo se encuentra el recordBase
         * - Si esta en el hijo 1, llamamos a la funcion recursivamente pero con un puntero a su hijo 1
         * - Si esta en el hijo 2, llamamos a la funcion recursivamente pero con un puntero a su hijo 2
         * - Si entramos en un nodo hoja, calculamos distancias con los puntos que contiene dicho nodo
         */

        if(ptr->isLeaf){
            getDistancesOfGroup(ptr, recordBase);

            // --> Ordenamos el vector de distancias
            sort(allKnnDistances.begin(), allKnnDistances.end());

            return;
        }
        else{
            // --> Obtenemos la distancia entre el record y el punto central de un hijo
            float distance1 = getEuclideanDistance2(recordBase, ptr->grupoHijos[0].puntoCentral);
            if(distance1 <= ptr->grupoHijos[0].radio){
                // * El record se encuentra en el hijo 1
                Node<T, ndim> *ptrChild;
                ptrChild = &ptr->grupoHijos[0];
                knnSearch(recordBase, ptrChild, Nvecinos);
            }
            else{
                // * El record se encuentra en el hijo 2
                Node<T, ndim> *ptrChild;
                ptrChild = &ptr->grupoHijos[1];
                knnSearch(recordBase, ptrChild, Nvecinos);
            }
        }
    }

    // ############################# KNN SEARCH V2 ##########################
    float logicDistMinOfNode(Point<ndim-1> recordBase, Node<T, ndim> *ptr){
        if(ptr->isLeaf){
            float distMin1 = getEuclideanDistance2(recordBase, ptr->puntoCentral) - ptr->radio;
            float distMin2 = 0.0;
            return max(distMin1, distMin2);
        }
        else{
            float distMin1 = getEuclideanDistance2(recordBase, ptr->puntoCentral) - ptr->radio;
            float distmin2 = logicDistMinOfNode(recordBase, ptr->grupoHijos[0]);
            return max(distMin1, distmin2);
        }
    }

    vector< pair< Point<ndim-1>,float > > knnSearchV2( Node<T, ndim> *ptr, Point<ndim-1> recordBase, int Nvecinos, vector< pair< Point<ndim-1>,float > > p_in){

        vector< pair< Point<ndim-1>,float > > p_out;

        // --> Obtenemos la distancia maxima del nodo actual
        float distanciaMaximaActual;
        if(p_in.size() < Nvecinos){
            distanciaMaximaActual = numeric_limits<float>::infinity();
        }
        else{
            float dist_Max = 0;
            for(int i = 0; i < p_in.size(); i++){
                float dist_Max_aux = getEuclideanDistance2(p_in[i].first, recordBase);
                if(dist_Max_aux > dist_Max){
                    dist_Max = dist_Max_aux;
                }
            }
            distanciaMaximaActual = p_in[p_in.size()-1].second;
        }

        // --> Obtenemos la distancia minia del nodo
        float distanciaMinima = logicDistMinOfNode(recordBase, ptr);

        // --> Recorrido del arbol
        if(distanciaMinima > distanciaMaximaActual){
            return p_in;
        }
        else if(ptr->isLeaf == true){
            p_out = p_in;
            for(int i = 0; i < ptr->grupoRecords.size(); i++){
                float dist_aux = getEuclideanDistance(ptr->grupoRecords[i], recordBase);
                if(dist_aux < distanciaMaximaActual){
                    pair< Point<ndim-1>,float > aux_pair;
                    aux_pair.first = ptr->grupoRecords[i];
                    aux_pair.second = dist_aux;
                    p_out.push_back(aux_pair);
                }
            }
        }
        else if(ptr->isLeaf == false){
            float dist_aux = getEuclideanDistance2(ptr->grupoHijos[0].puntoCentral, recordBase);
            if(dist_aux < distanciaMaximaActual){
                p_out = knnSearchV2(&ptr->grupoHijos[0], recordBase, Nvecinos, p_in);
            }
            else{
                p_out = p_in;
            }

            dist_aux = getEuclideanDistance2(ptr->grupoHijos[1].puntoCentral, recordBase);
            if(dist_aux < distanciaMaximaActual){
                p_out = knnSearchV2(&ptr->grupoHijos[1], recordBase, Nvecinos, p_out);
            }
        }
        return p_out;
    }
};

template <typename T, int ndim>
class BallStarTree {
private:
    Node<T, ndim> root;
    Node<T, ndim> *ptr;

public:
    BallStarTree(int _M){
        Node<T, ndim> aux(_M);
        root = aux;

        ptr = &root;
    }

    void insertToRoot( Point<ndim> newRecord ){
        root.insertToNode(newRecord, ptr);
    }

    void _createEstructure(){
        root.createEstructure(ptr);
    }

    int getIdOfSong( string song_name ){
        for(int i = 1; i <= allNamesSongs.size(); i++){
            if(allNamesSongs[i-1] == song_name){
                return i;
            }
        }
        cout << "El nombre de la cancion no se encontro en el DataSet" << endl;
        return -1;
    }

    Point<ndim - 1> getPoints(int idSong ){
        Point<ndim - 1> auxPoint;
        for(int i = 0; i < ptr->grupoRecords.size(); i++){
            if(ptr->grupoRecords[i][ndim-1] == idSong){ // Si se encontro...
                for(int j = 0; j < ndim - 1; j++){
                    auxPoint[j] = ptr->grupoRecords[i][j];
                }
                return auxPoint;
            }
        }
        cout << "El id de la cancion no se encontro en el DataSet" << endl;
        return auxPoint;
    }

    void _knnSearch( string song_name, int Nvecinos ){
        // LOGICA
        /*
         * - Buscar que indice tiene el nombre de la cancion
         * - Obtener la informacion del id buscado
         * - Buscar los Nvecinos mas cercanos
         * - Mostrar los Nvecinos mas cercanos
         */

        // --> Buscar que indice tiene el nombre de la cancion
        int idSong = getIdOfSong(song_name);
        // * En caso no se encuentre la cancion
        if(idSong == -1){
            cout << "Abortando..." << endl;
            return;
        }

        cout << "La cancion: <<" << song_name << ">>, tiene el id " << idSong << endl;

        // --> Obtener la informacion del id buscado
        Point<ndim-1> recordBase;
        recordBase = getPoints(idSong);

        // --> Buscar los Nvecinos mas cercanos
        // * Try 1 (Simple)
        root.knnSearch(recordBase, ptr, Nvecinos);

        // * Try 2 (Mas seguro)
        //vector< pair< Point<ndim-1>,float > > p_in;
        //vector< pair< Point<ndim-1>,float > > knnList = root.knnSearchV2(ptr, recordBase, Nvecinos, p_in);
    }
};

int main()
{
    // ------- MAIN INFO -------
    const int M = 17590; // 17590
    const int ndim = 14;
    // ------- TEST INFO -------
    //const int M = 2;
    //const int ndim = 4;
    // -------------------------

    cout << endl << "############################################### " << endl;
    cout << "               BALL * TREE " << endl;
    cout << "############################################## " << endl;

    // ---------------- CONFIGURACION PARA OBTENER EL TIEMPO -------------------
    auto start_Index = std::chrono::steady_clock::now();
    // -------------------------------------------------------------------------

    BallStarTree<Point<ndim>, ndim> test1(M);

    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INDEXADO Y CREACION DE LA ESTRUCTURA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    // ########################## LECTURA .SVG ##################################
    string filename = "finalDataSet.csv";
    //string filename = "ej1_main.csv"; // TEST FILE

    ifstream in(filename);
    string line;

    while(getline(in, line)){
        stringstream ss(line);
        string token;
        float recordAux[ndim];

        // * Guardamos unicamente los datos que tienen informacion numerica y el ID
        for(int i = 0; i < ndim; i++){
            getline(ss, token, ',');
            recordAux[i] = stof(token);
        }

        // * Guardamos el nombre de la cancion en un vector aparte
        getline(ss, token, ',');
        allNamesSongs.push_back(token);

        // * Guardamos el record en el arbol
        Point<ndim> record(recordAux);
        test1.insertToRoot(record);
    }
    // ######################### ESTRUCTURA #####################################
    // * Creamos la estructura del arbol
    test1._createEstructure();
    // ##########################################################################

    // ---------------- CONFIGURACION PARA OBTENER EL TIEMPO -------------------
    auto end_Index = std::chrono::steady_clock::now();
    // -------------------------------------------------------------------------

    // ######################### MOSTRAMOS ESTADISTICAS ############################
    cout << endl << endl;
    cout << "################## TIEMPO DE INDEXACION ############################ " << endl;
    // --> Tiempo
    cout << "##" << '\t' << "Execution Time: "
                        << chrono::duration_cast<chrono::nanoseconds>(end_Index - start_Index).count()
                        << "ns" << '\t' << '\t' << '\t' << '\t' << "##" << endl;
    cout << "#################################################################### " << endl;

    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% KNN VECINOS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    cout << endl << endl << "Iniciando busqueda del vecino mas cercano ..." << endl << endl;

    // ---------------- CONFIGURACION PARA OBTENER EL TIEMPO -------------------
    auto start_knn = std::chrono::steady_clock::now();
    // -------------------------------------------------------------------------

    int Nvecinos = 50; // Tiene que ser menor que M
    //int Nvecinos = 50; // TEST NVECINOS

    if(Nvecinos > M){
        cout << "El numero de vecinos no puede ser mayor que M (Linea 573) " << endl;
        return 0;
    }

    string song_name = "CHOPSTICK"; // CHOPSTICK
    test1._knnSearch(song_name, Nvecinos);

    // ---------------- CONFIGURACION PARA OBTENER EL TIEMPO -------------------
    auto end_knn = std::chrono::steady_clock::now();
    // -------------------------------------------------------------------------

    // --> Mostramos los primeros N vecinos mas cercanos
    cout << endl << "Los " << Nvecinos << " vecinos mas cercanos son: " << endl;
    for(int i = 0; i < Nvecinos; i++){
        // --> Buscamos el nombre de la cancion
        float idSong2 = allKnnDistances[i][1];
        string songName = allNamesSongs[idSong2-1];
        cout << i+1 << ". " << songName << '\t' << '\t' << '\t' << "Distance: " << allKnnDistances[i][0] << endl;
    }

    // ######################### MOSTRAMOS ESTADISTICAS ############################
    cout << endl << endl;
    cout << "################## TIEMPO DE CONSULTA KNN ######################### " << endl;
    // --> Tiempo
    cout << "##" << '\t' << "Execution Time: "
                        << chrono::duration_cast<chrono::nanoseconds>(end_knn - start_knn).count()
                        << "ns" << '\t' << '\t' << '\t' << '\t' << "##" << endl;
    cout << "################################################################## " << endl << endl;

    return 0;
}