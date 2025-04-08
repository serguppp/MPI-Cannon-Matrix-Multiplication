#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

#define N 4 // Rozmiar macierzy
#define PP 2   // Pierwiastek z liczby procesów
#define P 4   // Liczba procesów

float rawA[N][N], rawB[N][N]; // Macierze A i B z pliku
float distrA[N][N], distrB[N][N]; //Macierze A i B po poczatkowej dystrybucji

float a[N / PP][N / PP], b[N / PP][N / PP], c[N / PP][N / PP]; 
float aa[N / PP][N / PP], bb[N / PP][N / PP]; 

float(*psa)[N / PP], (*psb)[N / PP], (*pra)[N / PP], (*prb)[N / PP];

float CSek[N][N], Cglob[N][N];  

int rank, np, koniec;

double startwtime1, startwtime2, endwtime;

//Wczytywanie plików
void loadFile(FILE **file, const char *path) {
    printf("Proces %d sprawdza obecność pliku tekstowego z danymi: %s\n", rank, path);
    *file = fopen(path, "r"); 
    if (*file == NULL) {
        printf("Proces %d wskazuje błąd otwarcia pliku \"%s\"\n", rank, path);
        koniec = 1;
        MPI_Bcast(&koniec, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Finalize();
        exit(0);
    } else {
        printf("Proces %d poprawnie otworzył plik \"%s\"\n", rank, path);
        koniec = 0;
        MPI_Bcast(&koniec, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }
}

void loadMatrix(FILE **file, const char *path, const char id){
    printf("Proces %d odczytuje tablice z pliku \"%s\"\n", rank, path);

    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++) {
            if (id=='A'){
                fscanf(*file, "%f", &rawA[i][j]);
                //printf("Wczytano a[%d][%d] %f\n", i, j, A[i][j]);
            }
            else if (id=='B'){
                fscanf(*file, "%f", &rawB[i][j]);
                //printf("Wczytano b[%d][%d] %f\n", i, j, B[i][j]);
            }

        }
}
// Funkcja TESTUJĄCAA do wypisywania tablicy 
void printMatrix(float matrix[N][N], int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            printf("%f ", matrix[i][j]); 
        }
        printf("\n");
    }
}

// Funkcja zajmująca się wstępną inicjalizacją macierzy
void initMatrix(){
    int new_j_A;
    int new_i_B;
    for (int i = 0; i<N; i++){
        for (int j=0; j<N; j++){
            new_j_A = (j - i + N) % N;
            new_i_B = (i - j + N) % N;            
            distrA[i][new_j_A] = rawA[i][j];
            distrB[new_i_B][j] = rawB[i][j];
        }
    }
}

void initC(){
    for (int i = 0; i < N / PP; i++)
        for (int j = 0; j < N / PP; j++)
             c[i][j] = 0;
    //printf("Utworzono tablice C\n");
}

// Funkcja wypisująca blok procesów dla macierzy a, b i c
void printBlock(float matrix[N / PP][N / PP], const char* name) {
    printf("Proces %d, Macierz %s:\n", rank, name);
    for (int i = 0; i < N / PP; i++) {
        for (int j = 0; j < N / PP; j++) {
            printf("%f ", matrix[i][j]);
        }
        printf("\n");
    }
}

int main(int argc, char** argv) {
    FILE *fileA, *fileB;
    const char *pathA = "matrix1.txt";
    const char *pathB = "matrix2.txt";

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    int row = rank / PP;
    int col = rank % PP;
    int mod = 0;
    int tag = 101;

    MPI_Status statSend[4], statRecv[4];
    MPI_Request reqSend[4], reqRecv[4];
    MPI_Comm cartComm;


    //Sprawdzamy, czy liczba uruchomionych procesów MPI (np) jest zgodna z oczekiwaną liczbą procesów (P)
    // i wykonuje odpowiednią akcję, jeśli liczba procesów jest nieprawidłowa.
    if (np != P) {
        if (rank == 0)
            printf("Ilość procesorów MPI niezgodna z oczekiwaną liczbą procesów. Uruchom mpiexec -n %d program\n", P);
        MPI_Finalize();
        exit(0);
    }

    if (rank == 0)
        printf("Obliczenia metodą Cannona dla macierzy %d x %d\n", N, N);

    if (rank == 0) startwtime1 = MPI_Wtime();

    //wczytanie danych przez proces rank=0
	if (rank == 0){
        loadFile(&fileA, pathA);        
        loadFile(&fileB, pathB);    
    }
	else{
		MPI_Bcast(&koniec, 1, MPI_INT, 0, MPI_COMM_WORLD);
		if (koniec) { 
            MPI_Finalize(); exit(0); 
        }
	}

    if (rank == 0) {
            // Odczyt danych wejściowych z plików
            loadMatrix(&fileA, pathA, 'A');
            loadMatrix(&fileB, pathB, 'B');


            // Wstępna dystrybucja macierzy
            initMatrix();

            // Rozesłanie tablic do pozostałych procesów
            // Iterujemy przez wszystkie procesy
            for (int p = 0; p < P; p++){
                int i = p % PP; // numer kolumny
                int j = p / PP; // numer wiersza

                float subA[N / PP][N / PP];
                float subB[N / PP][N / PP];

                // Dla każdego bloku a i b kopiujemy odpowiedni fragment z A i B:
                 // cos tu jest nie tak
                for (int x = 0; x < N/PP; x++){
                    for (int y = 0; y < N/PP; y++){
                        subA[x][y] = distrA[j * (N / PP) + x][i * (N / PP) + y];
                        subB[x][y] = distrB[j * (N / PP) + x][i * (N / PP) + y];
                    }
                }

                // Wysłanie bloków do procesów
                if (p != 0){
                    MPI_Isend(subA, N * N / PP / PP, MPI_FLOAT, p, tag, MPI_COMM_WORLD, &reqSend[0]);
                    MPI_Isend(subB, N * N / PP / PP, MPI_FLOAT, p, tag, MPI_COMM_WORLD, &reqSend[1]);
                }
                // proces 0 zapisuje swój fragment do lokalnych a i b, zamiast go wysyłać.
                else{
                    memcpy(a, subA, sizeof(subA));
                    memcpy(b, subB, sizeof(subB));
                }
            }
    } else {
        // Procesy różne od rank=0 odbierają tablicę a i b
        MPI_Irecv(a, N * N / PP / PP, MPI_FLOAT, 0, tag, MPI_COMM_WORLD, &reqRecv[0]);
        MPI_Irecv(b, N * N / PP / PP, MPI_FLOAT, 0, tag, MPI_COMM_WORLD, &reqRecv[1]);
        
        MPI_Wait(&reqRecv[0], MPI_STATUS_IGNORE);
        MPI_Wait(&reqRecv[1], MPI_STATUS_IGNORE);
    }

    //Zerowanie tablic aa i bb    
    memset(aa, 0, sizeof(aa));
    memset(bb, 0, sizeof(bb));

    //Inicjalziacja tablicy wynikowej
    initC();

    if (rank == 0) startwtime2 = MPI_Wtime();
    //Algorytm Cannona
    pra = aa; prb = bb; psa = a; psb = b;
    
    for (int kk = 0; kk < PP; kk++) {   
        for (int i = 0; i < N / PP; i++)
            for (int j = 0; j < N / PP; j++)
                for (int k = 0; k < N / PP; k++)
                    c[i][j] += psa[i][k] * psb[k][j];

        int left = ((col - 1 + PP) % PP) + row * PP;
        int right = ((col + 1) % PP) + row * PP;
        int up = ((row - 1 + PP) % PP) * PP + col;
        int down = ((row + 1) % PP) * PP + col;

        MPI_Irecv(pra, N * N / PP / PP, MPI_FLOAT, right, tag, MPI_COMM_WORLD, &reqRecv[2]); 
        MPI_Irecv(prb, N * N / PP / PP, MPI_FLOAT, down, tag, MPI_COMM_WORLD, &reqRecv[3]);
    
        MPI_Isend(psa, N * N / PP / PP, MPI_FLOAT, left, tag, MPI_COMM_WORLD, &reqSend[2]);
        MPI_Isend(psb, N * N / PP / PP, MPI_FLOAT, up, tag, MPI_COMM_WORLD, &reqSend[3]);

        MPI_Wait(&reqSend[2], MPI_STATUS_IGNORE);
        MPI_Wait(&reqSend[3], MPI_STATUS_IGNORE);
        MPI_Wait(&reqRecv[2], MPI_STATUS_IGNORE);
        MPI_Wait(&reqRecv[3], MPI_STATUS_IGNORE);

        MPI_Barrier(MPI_COMM_WORLD);

        mod = (mod + 1) % 2; // Zmiana wartości mod przed warunkiem

        if (mod == 1) {  
            pra = a; prb = b; psa = aa; psb = bb;
        } else {
            pra = aa; prb = bb; psa = a; psb = b;
        }

    }

    MPI_Gather(c, N * N / P, MPI_FLOAT, Cglob, N * N / P, MPI_FLOAT, 0, MPI_COMM_WORLD);

    //printf("Ukonczono mnozenie macierzy dla rank %d\n", rank);
    //Koniec zliczania czasow
    if (rank == 0) {
        endwtime = MPI_Wtime();
        printf("Czas przetwarzania: %f sekund\n", endwtime - startwtime1);
        printf("Czas obliczeń: %f sekund\n", endwtime - startwtime2);
    }

    //Test poprawności wyników:
    if (rank == 0) {
        printf("Macierz A do sek: \n");
        printMatrix(rawA, N, N);
        printf("Macierz B do sek: \n");
        printMatrix(rawB, N, N);
        // Obliczenia sekwencyjne mnożenia tablic CSek=A*B
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++) {
                CSek[i][j] = 0;
                for (int k = 0; k < N; k++) {
                    CSek[i][j] += rawA[i][k] * rawB[k][j];
                }
            }
    
        // Odbiór wyników obliczeń równoległych do globalnej tablicy wynikowej Cglob
        // Rozesłanie wyników do procesów (w razie potrzeby)
        // Porównanie poprawności obliczeń (Csek, Cglob) przy uwzględnieniu progu poprawności
        int errors = 0;
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                if (fabs(CSek[i][j] - Cglob[i][j]) > 1e-6) {
                    //printf("Błąd na pozycji (%d, %d): Csek = %f, Cglob = %f\n", i, j, CSek[i][j], Cglob[i][j]);
                    errors++;
                }
            }
        }
        if (errors == 0) {
            printf("Wyniki są poprawne!\n");
        } else {
            printf("Znaleziono %d błędów w wynikach.\n", errors);
        }

        printf("Macierz wynikowa Cglob:\n");
        printMatrix(Cglob, N, N);
        printf("Macierz wynikowa Csek:\n");
        printMatrix(CSek, N, N);
    }

    if (rank == 0) {
        printf("Program MPI zakończony poprawnie!\n");
    }

    MPI_Finalize();
    return 0;
}
