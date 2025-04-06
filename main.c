#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

#define N 3 // Rozmiar macierzy
#define PP 3    // Pierwiastek z liczby procesów
#define P 9     // Liczba procesów

float A[N][N], B[N][N];
float a[N / PP][N / PP], b[N / PP][N / PP], c[N / PP][N / PP]; 
float aa[N / PP][N / PP], bb[N / PP][N / PP]; // cos tu jest nie tak - do czego one sluza, czy maja byc puste?
float(*psa)[N / PP], (*psb)[N / PP], (*pra)[N / PP], (*prb)[N / PP];

float CSek[N][N], Cglob[N][N];   // cos tu jest nie tak a konkretnie ich obliczenia sekwencyjne mnozenia tablic
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
                fscanf(*file, "%f", &A[i][j]);
                printf("Wczytano a[%d][%d] %f\n", i, j, A[i][j]);
            }
            else if (id=='B'){
                fscanf(*file, "%f", &B[i][j]);
                printf("Wczytano b[%d][%d] %f\n", i, j, B[i][j]);
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
    float tmpA[N][N];
    float tmpB[N][N];
    int new_j_A;
    int new_i_B;

    memcpy(tmpA, A, sizeof(A));  // Kopiowanie całej tablicy A do tmpA
    memcpy(tmpB, B, sizeof(B));  // Kopiowanie całej tablicy A do tmpA

    for (int i = 0; i<N; i++){
        for (int j=0; j<N; j++){
            new_j_A = (j - i + PP) % PP;
            new_i_B = (i - j + PP) % PP;
            A[i][new_j_A] = tmpA[i][j];
            B[new_i_B][j] = tmpB[i][j];
        }
    }

}

void initC(){
    int i = rank % N; // kolumna w siatce procesów
    int j = rank / N; // wiersz w siatce procesów

    for (int i = 0; i < N / PP; i++)
        for (int j = 0; j < N / PP; j++)
             c[i][j] = 0;
    printf("Utworzono tablice C\n");
}

int main(int argc, char** argv) {
    FILE *fileA, *fileB;
    const char *pathA = "matrix1.txt";
    const char *pathB = "matrix2.txt";

    int row = rank / PP;
    int col = rank % PP;
    int mod = 0;
    int tag = 101;
    
    MPI_Status statSend[4], statRecv[4];
    MPI_Request reqSend[4], reqRecv[4];

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    if (rank == 0)
        printf("Obliczenia metodą Cannona dla macierzy %d x %d\n", N, N);

    if (rank == 0) startwtime1 = MPI_Wtime();

    //wczytanie danych przez proces rank=0
	if (rank == 0){
        loadFile(&fileA, pathA);        
        loadFile(&fileB, pathA);    
    }
	else{
		MPI_Bcast(&koniec, 1, MPI_INT, 0, MPI_COMM_WORLD);
		if (koniec) { 
            MPI_Finalize(); exit(0); 
        }
	}

    //Sprawdzamy, czy liczba uruchomionych procesów MPI (np) jest zgodna z oczekiwaną liczbą procesów (P)
    // i wykonuje odpowiednią akcję, jeśli liczba procesów jest nieprawidłowa.
    if (np != P) {
        if (rank == 0)
            printf("Ilość procesorów MPI niezgodna z oczekiwaną liczbą procesów. Uruchom mpiexec -n %d program\n", P);
        MPI_Finalize();
        exit(0);
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
                        subA[x][y] = A[i * (N / PP) + x][j * (N / PP) + y]; // i * (N / PP) + x – określa rzeczywisty wiersz w A lub B
                        printf("subA[%d][%d] = %f \n", x, y, subA[x][y]);
                        subB[x][y] = B[i * (N / PP) + x][j * (N / PP) + y]; // j * (N / PP) + y – określa rzeczywistą kolumnę w A lub B
                        printf("subB[%d][%d] = %f \n", x, y, subB[x][y]);
                    }
                }

                // Wysłanie bloków do procesów
                if (p != 0){
                    //MPI_Send(&subA, (N / PP) * (N / PP), MPI_FLOAT, p, 0, MPI_COMM_WORLD);
                    //MPI_Send(&subB, (N / PP) * (N / PP), MPI_FLOAT, p, 1, MPI_COMM_WORLD);
                    MPI_Isend(&subA, N * N / PP / PP, MPI_FLOAT, p, tag, MPI_COMM_WORLD, &reqSend[0]);
                    MPI_Isend(&subB, N * N / PP / PP, MPI_FLOAT, p, tag, MPI_COMM_WORLD, &reqSend[1]);
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
        
        //MPI_Wait(&reqRecv[0], MPI_STATUS_IGNORE);
        //MPI_Wait(&reqRecv[1], MPI_STATUS_IGNORE);

        MPI_Barrier(MPI_COMM_WORLD);
        
    }

    // Przygotowanie lokalnej tablicy wynikowej      
    initC();

    if (rank == 0) startwtime2 = MPI_Wtime();

    //Algorytm Cannona
    pra = aa; prb = bb; psa = a; psb = b;
    for (int kk = 0; kk < PP; kk++) {
        for (int i = 0; i < N / PP; i++)
            for (int k = 0; k < N / PP; k++)
                for (int j = 0; j < N / PP; j++)
                    c[i][j] += psa[i][k] * psb[k][j];

        int right = (col+1)%N + row; //Sąsiad po prawo
        int up =  (row - 1 + N) % N * N + col; // Sąsiad góraa
        int left = (col - 1 + N) % N + row; // Sąsiad lewo
        int down = (row + 1) % N * N + col; // Sąsiad niżej

        MPI_Irecv(pra, N * N / PP / PP, MPI_FLOAT, right, tag, MPI_COMM_WORLD, &reqRecv[2]); 
        printf("%d PROBUJE ODEBRAC Z PRAWEJ OD %d!!!\n", kk, right);
        MPI_Irecv(prb, N * N / PP / PP, MPI_FLOAT, down, tag, MPI_COMM_WORLD, &reqRecv[3]);
        printf("%d PROBUJE ODEBRAC Z DOLU OD %d!!!\n", kk, down);
    
        MPI_Isend(psa, N * N / PP / PP, MPI_FLOAT, left, tag, MPI_COMM_WORLD, &reqSend[2]);
        printf("%d PROBUJE WYSLAC DO LEFT DO %d!!!\n", kk, left);
        MPI_Isend(psb, N * N / PP / PP, MPI_FLOAT, up, tag, MPI_COMM_WORLD, &reqSend[3]);
        printf("%d PROBUJE WYSLAC DO GORY DO %d!!!\n", kk, up);
            //MPI_Wait(&reqRecv[2], &statRecv[2]);
        //MPI_Wait(&reqRecv[3], &statRecv[3]);
    
        MPI_Barrier(MPI_COMM_WORLD);

        if (mod = ((mod + 1) % 2)) {
            pra = a; prb = b; psa = aa; psb = bb;
        } else {
            pra = aa; prb = bb; psa = a; psb = b;
        }

    }

    printf("Ukonczono mnozenie macierzy dla rank %d\n", rank);
    //Koniec zliczania czasow
    if (rank == 0) {
        endwtime = MPI_Wtime();
        printf("Czas przetwarzania: %f sekund\n", endwtime - startwtime1);
        printf("Czas obliczeń: %f sekund\n", endwtime - startwtime2);
    }

    //Test poprawności wyników:

    if (rank == 0) {
        // Obliczenia sekwencyjne mnożenia tablic CSek=A*B
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++) {
                CSek[i][j] = 0;
                for (int k = 0; k < N; k++) {
                    CSek[i][j] += a[i][k] * b[k][j];
                }
            }
    
        // Odbiór wyników obliczeń równoległych do globalnej tablicy wynikowej Cglob
        // Rozesłanie wyników do procesów (w razie potrzeby)
        // Porównanie poprawności obliczeń (Csek, Cglob) przy uwzględnieniu progu poprawności
        int errors = 0;
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                if (fabs(CSek[i][j] - Cglob[i][j]) > 1e-6) {
                    errors++;
                }
            }
        }
        if (errors == 0) {
            printf("Wyniki są poprawne!\n");
        } else {
            printf("Znaleziono %d błędów w wynikach.\n", errors);
        }

        //printMatrix(CSek, 3, 3);
    }

    MPI_Finalize();
    return 0;
}
