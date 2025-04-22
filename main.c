#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

#define N 2000 // Rozmiar macierzy
#define PP 5   // Pierwiastek z liczby procesów
#define P 25   // Liczba procesów

// --- Globalne wskaźniki do macierzy (dynamicznie alokowane) ---

// Macierze A i B wczytane przez proces 0
double **rawA = NULL, **rawB = NULL;
// Macierze A i B po poczatkowej dystrybucji (przesunięciu) obliczone przez proces 0
double **distrA = NULL, **distrB = NULL;

// Lokalne bloki macierzy dla każdego procesu
double **a = NULL, **b = NULL, **c = NULL;
// Bufory pomocnicze do odbierania danych w pętli Cannona
double **aa = NULL, **bb = NULL;

// Wskaźniki do buforów wysyłania/odbierania w pętli Cannona
double **psa = NULL, **psb = NULL, **pra = NULL, **prb = NULL;

// Macierz wynikowa obliczona sekwencyjnie (dla weryfikacji)
double **cSek = NULL;
// Globalna macierz wynikowa zbierana w procesie 0
double **cGlob = NULL; 

// Zmienne pomocnicze: numer procesu, liczba procesów, flaga końca programu
int rank, np, finalize;

// Zmienne służące do pomiarów czasów cGlob i cSek
double start_cSek, end_cSek;
double start_cGlob, end_cGlob;

/// --- Dynamiczna alokacja/dealokacja tablic ---

// Alokacja tablicy
double **allocateMatrix(int rows, int cols) {
    double *data = (double *)malloc(rows * cols * sizeof(double));
    if (data == NULL) {
        perror("Nie udało się zaalokować bloku danych macierzy");
        return NULL;
    }
    double **matrix = (double **)malloc(rows * sizeof(double *));
    if (matrix == NULL) {
        perror("Nie udało się zaalokować wskaźników wierszy macierzy");
        free(data);  // Zwolnij blok danych, jeśli alokacja wskaźników się nie powiedzie
        return NULL;
    }
    for (int i = 0; i < rows; i++) {
        matrix[i] = &(data[i * cols]);
    }
    return matrix;
}

// Zwolnienie pamięci macierzy
void freeMatrix(double **matrix) {
    if (matrix != NULL) {
        if (matrix[0] != NULL) {
            free(matrix[0]); // Zwolnij ciągły blok danych
        }
        free(matrix); // Zwolnij wskaźniki wierszy
    }
}

// Funkcja wczytywania plików 
void loadFile(FILE **file, const char *path) {
    if (rank == 0) { // Tylko proces 0 otwiera pliki
        *file = fopen(path, "r");
        if (*file == NULL) {
            perror("Błąd otwarcia pliku");
            finalize = 1; // Ustaw flagę błędu
        } else {
            finalize = 0; // Brak błędu
        }
    }
    // Rozgłoś status otwarcia pliku do wszystkich procesów
    MPI_Bcast(&finalize, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (finalize == 1) { // Jeśli był błąd, zakończ
        if (rank == 0) printf("Zakończenie programu z powodu błędu pliku.\n");
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
}

// Funkcja wczytująca dane z pliku do macierzy
void loadMatrix(FILE **file, const char *path, const char id, double **matrix, int rows, int cols) {
    if (rank == 0) { // Tylko proces 0 wczytuje dane
        // printf("Wczytywanie macierzy %c z pliku %s...\n", id, path);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                if (fscanf(*file, "%lf", &matrix[i][j]) != 1) {
                     if (feof(*file)) {
                        fprintf(stderr, "Błąd: Niespodziewany koniec pliku %s przy (%d, %d). Za mało danych.\n", path, i, j);
                    } else if (ferror(*file)) {
                        perror("Błąd odczytu danych z pliku");
                    } else {
                        fprintf(stderr, "Błąd: Nieprawidłowy format danych w pliku %s przy (%d, %d).\n", path, i, j);
                    }
                    fclose(*file);
                    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE); // Zakończ wszystkie procesy
                }
            }
        }
        fclose(*file); // Zamknij plik po wczytaniu
        // printf("Wczytywanie macierzy %c zakończone.\n", id);
    }
}

// Funkcja testowa wypisująca macierz
void printMatrix(double **matrix, int rows, int cols, const char* name) {
     printf("Macierz %s (%dx%d):\n", name, rows, cols);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            printf("%6.1f ", matrix[i][j]);
        }
        printf("\n");
    }
}

// Funkcja inicjalizująca macierze
void initMatrix(int blockSize) {
    // Bufory pomocnicze 
    double **tempA = NULL, **tempB = NULL;
    if (rank == 0) {
        tempA = allocateMatrix(N, N);
        tempB = allocateMatrix(N, N);
        if (!tempA || !tempB) {
             fprintf(stderr, "Rank 0: Błąd alokacji tempA/tempB w initMatrix.\n");
             MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }

        // Kopiujemy początkowe dane do buforów 
        memcpy(tempA[0], rawA[0], N * N * sizeof(double));
        memcpy(tempB[0], rawB[0], N * N * sizeof(double));

        // printf("Inicjalizacja (przesuwanie) macierzy A i B...\n");
        // Przesuwanie bloków zgodnie z metodą Cannona
        for (int blockRow = 0; blockRow < PP; blockRow++) {
            for (int blockCol = 0; blockCol < PP; blockCol++) {
                // Obliczamy nowe indeksy bloków po przesunięciu
                int newColA = (blockCol - blockRow + PP) % PP; // A przesuwamy w lewo o blockRow
                int newRowB = (blockRow - blockCol + PP) % PP; // B przesuwamy w górę o blockCol

                // Skopiuj całe bloki
                for (int i = 0; i < blockSize; i++) {
                    double *destA_row_start = &distrA[blockRow * blockSize + i][newColA * blockSize];
                    double *srcA_row_start = &tempA[blockRow * blockSize + i][blockCol * blockSize];
                    double *destB_row_start = &distrB[newRowB * blockSize + i][blockCol * blockSize];
                    double *srcB_row_start = &tempB[blockRow * blockSize + i][blockCol * blockSize];

                    memcpy(destA_row_start, srcA_row_start, blockSize * sizeof(double));
                    memcpy(destB_row_start, srcB_row_start, blockSize * sizeof(double));
                }
            }
        }
        // printf("Inicjalizacja macierzy A i B zakończona.\n");

        // Zwolnij tymczasowe bufory
        freeMatrix(tempA);
        freeMatrix(tempB);
    }
}


// Funkcja wypisująca lokalny blok macierzy
void printBlock(double **matrix, int rows, int cols, const char* name) {
    printf("Proces %d, Macierz %s (%dx%d):\n", rank, name, rows, cols);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            printf("%6.1f ", matrix[i][j]);
        }
        printf("\n");
    }
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
    MPI_Comm_size(MPI_COMM_WORLD, &np); 

    FILE *fileA = NULL, *fileB = NULL;
    const char *pathA = "matrix1.txt";
    const char *pathB = "matrix2.txt";

    // --- Sprawdzanie podstawowych parametrów ---
    if (N % PP != 0) {
        if (rank == 0) {
            fprintf(stderr, "Błąd: Rozmiar macierzy N (%d) musi być podzielny przez PP (%d).\n", N, PP);
        }
        MPI_Finalize();
        return EXIT_FAILURE;
    }
     if (np != P) {
        if (rank == 0)
            printf("Ilość procesorów MPI (%d) niezgodna z oczekiwaną liczbą procesów P (%d). Uruchom mpiexec -n %d program\n", np, P, P);
        MPI_Finalize();
        return EXIT_FAILURE; 
    }
    if (PP * PP != P) {
         if (rank == 0)
            fprintf(stderr, "Błąd: P (%d) musi być kwadratem PP (%d).\n", P, PP);
        MPI_Finalize();
        return EXIT_FAILURE;
    }

    int row = rank / PP; // Wiersz procesu w siatce PPxPP
    int col = rank % PP; // Kolumna procesu w siatce PPxPP
    int tag = 101; // Tag dla wiadomości MPI
    int blockSize = N / PP; // Rozmiar boku bloku lokalnego`
    int localSize = blockSize * blockSize; // Liczba elementów w bloku lokalnym

    // --- Alokacja lokalnych macierzy dla wszyzstkich procesów ---
    a = allocateMatrix(blockSize, blockSize);
    b = allocateMatrix(blockSize, blockSize);
    c = allocateMatrix(blockSize, blockSize);
    aa = allocateMatrix(blockSize, blockSize); 
    bb = allocateMatrix(blockSize, blockSize); 
    if (!a || !b || !c || !aa || !bb) {
        fprintf(stderr, "Proces %d: Błąd alokacji macierzy lokalnych (a,b,c,aa,bb).\n", rank);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    // --- Rank 0: Wczytywanie danych, lokalna dystrybucja tablocy, rozsyłanie bloków między procesami ---
    if (rank == 0)
    {
        //printf("Proces 0 rozpoczyna pracę...\n");
        //printf("Alokacja pamięci dla macierzy globalnych (rank 0)...\n");
        rawA = allocateMatrix(N, N);
        rawB = allocateMatrix(N, N);
        distrA = allocateMatrix(N, N);
        distrB = allocateMatrix(N, N);
        cGlob = allocateMatrix(N, N); 
        cSek = allocateMatrix(N, N); 

        if (!rawA || !rawB || !distrA || !distrB || !cGlob || !cSek) {
             fprintf(stderr, "Rank 0: Błąd alokacji macierzy globalnych.\n");
             MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }
        //printf("Alokacja pamięci zakończona.\n");

        // Wczytanie plików
        loadFile(&fileA, pathA);
        loadFile(&fileB, pathB);

        // Wczytanie danych do rawA i rawB
        loadMatrix(&fileA, pathA, 'A', rawA, N, N);
        loadMatrix(&fileB, pathB, 'B', rawB, N, N);

        // Dystrybucja poczatkowa blokow 
        initMatrix(blockSize);

        // --- Rozesłanie odpowiednich bloków z distrA i distrB ---
        // printf("Rozsyłanie bloków początkowych do procesów...\n");
        for (int p = 0; p < P; p++){
            int target_proc_row = p / PP; // wiersz procesu docelowego
            int target_proc_col = p % PP; // kolumna procesu docelowego

            // Alokacja buforów pomocniczów do rozsyłania  bloków
            double **subA = allocateMatrix(blockSize, blockSize);
            double **subB = allocateMatrix(blockSize, blockSize);
             if (!subA || !subB) {
                 fprintf(stderr, "Rank 0: Błąd alokacji subA/subB w pętli dystrybucji.\n");
                 MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
            }

            // Kopiowanie odpowiedniego bloku z przesuniętych macierzy distrA/distrB
            for (int x = 0; x < blockSize; x++){
                for (int y = 0; y < blockSize; y++){
                    subA[x][y] = distrA[target_proc_row * blockSize + x][target_proc_col * blockSize + y];
                    subB[x][y] = distrB[target_proc_row * blockSize + x][target_proc_col * blockSize + y];
                }
            }

            // Wyślij bloki do procesów 'p' lub skopiuj lokalnie dla p=0
            if (p != 0){
                MPI_Send(subA[0], localSize, MPI_DOUBLE, p, tag, MPI_COMM_WORLD);
                MPI_Send(subB[0], localSize, MPI_DOUBLE, p, tag + 1, MPI_COMM_WORLD); // Użyj innego tagu dla B
            }
            else { 
                 memcpy(a[0], subA[0], localSize * sizeof(double));
                 memcpy(b[0], subB[0], localSize * sizeof(double));
            }

            // Zwalnianie bufoów
            freeMatrix(subA);
            freeMatrix(subB);
        }
        // printf("Rozsyłanie bloków zakończone.\n");

        // Zwalnianie tablic, gdyż nie będą już używane

        freeMatrix(distrA); distrA = NULL;
        freeMatrix(distrB); distrB = NULL;


    } else { // --- Rank !=0 - odbiór danych ---
        //  Odbierz do ciągłego bloku danych (a[0] i b[0])
        MPI_Recv(a[0], localSize, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(b[0], localSize, MPI_DOUBLE, 0, tag + 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    // --- Synchronizacja przed rozpoczęciem obliczeń---
    MPI_Barrier(MPI_COMM_WORLD);
    
    // Zerowanie macierzy wynikowej c  
    memset(c[0], 0, localSize * sizeof(double)); 

    if (rank == 0) {
        //printf("Rozpoczynanie algorytmu Cannona...\n");
        start_cGlob = MPI_Wtime(); // Rozpoczęcie mierzenia czasu
    }

    // --- Algorytm Cannona ---

    pra = aa; prb = bb; // Bufory odbiorcze
    psa = a; psb = b;   // Bufory do wysłania 

    for (int kk = 0; kk < PP; kk++) {
        // Krok 1: Mnożenie lokalnych bloków C += A * B 
        for (int i = 0; i < blockSize; i++) {
            for (int k = 0; k < blockSize; k++) {
                for (int j = 0; j < blockSize; j++) {
                    c[i][j] += psa[i][k] * psb[k][j];
                }
            }
        }

        // Krok 2: Przesunięcie bloków 
        int left = (col - 1 + PP) % PP + row * PP;
        int right = (col + 1) % PP + row * PP;
        int up = col + ((row - 1 + PP) % PP) * PP;
        int down = col + ((row + 1) % PP) * PP;

        // Wyślij obecny blok A (psa) w lewo, odbierz nowy z prawej do bufora odbiorczego A (pra)
        MPI_Sendrecv(psa[0], localSize, MPI_DOUBLE, left, tag,
                     pra[0], localSize, MPI_DOUBLE, right, tag,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // Wyślij obecny blok B (psb) w górę, odbierz nowy z dołu do bufora odbiorczego B (prb)
        MPI_Sendrecv(psb[0], localSize, MPI_DOUBLE, up, tag + 1, 
                     prb[0], localSize, MPI_DOUBLE, down, tag + 1, 
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // Krok 3: Zamiana wskaźników buforów (Pointer Swap)
        // Odebrane bloki (w pra, prb) stają się blokami do mnożenia/wysłania (psa, psb) w następnej iteracji.
        // Bloki, które właśnie wysłaliśmy (były w psa, psb) stają się buforami odbiorczymi (pra, prb).
        double **temp_ptr;
        temp_ptr = psa; psa = pra; pra = temp_ptr; 
        temp_ptr = psb; psb = prb; prb = temp_ptr; 
    }

    // --- Zbieranie wyników do procesu 0 ---
    //  cGlob jest [N][N], ale MPI_Gather oczekuje liniowego bufora, dlatego wymaga będzie rekonstrukcja tej macierzy.
    // Zbieramy lokalne wyniki (c[0]) do globalnej macierzy wynikowej (cGlob[0]) dla ranku 0

    MPI_Gather(c[0], localSize, MPI_DOUBLE,           // Dane wysyłane: lokalny blok c
               (rank == 0) ? cGlob[0] : NULL,         // Dane odbierane: globalna macierz cGlob (traktowana liniowo)
               localSize, MPI_DOUBLE,                 // Liczba i typ danych odbieranych
               0, MPI_COMM_WORLD);                    // Proces root - my używamy rank 0




    // --- Rank 0: Rekonstrukcja macierzy i weryfikacja ---
    if (rank == 0) {
        end_cGlob = MPI_Wtime(); // Zakończenie pomiaru czasy
        // printf("Algorytm Cannona zakończony.\n");
        // printf("Rekonstrukcja macierzy wynikowej cGlob...\n")
        freeMatrix(rawA); 
        freeMatrix(rawB);

    }
    // --- Czyszczenie danych ---
    // Zwalniamy macierze dla wszystkich procesów
    freeMatrix(a);
    freeMatrix(b);
    freeMatrix(c);
    freeMatrix(aa);
    freeMatrix(bb);
    if (rank == 0){
        // Alokacja pomocniczej tablicy do rekonstrukcji cGlob;
        double **Ctemp = allocateMatrix(N, N);
         if (!Ctemp) {
             fprintf(stderr, "Rank 0: Błąd alokacji Ctemp.\n");
             MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
         }

        // --- Rekonstrukcja cGlob ---
        // MPI_Gather umieszcza bloki kolejno w pamięci według rangi (0, 1, 2, ...)
        double *gathered_data = cGlob[0];  // Wskaźnik na początek zebranych danych
        int current_pos = 0; // Pozycja w liniowym buforze 

        for (int p = 0; p < P; p++) {
            // Oblicz lewy górny róg bloku (proc_row, proc_col) dla procesu 'p'
            int proc_row = p / PP;
            int proc_col = p % PP;

            // Skopiuj dane z liniowego gathered_data do Ctemp na odpowiednią pozycję bloku
            for (int i = 0; i < blockSize; i++) { // i - wiersz w bloku
                for (int j = 0; j < blockSize; j++) { // j - kolumna w bloku
                    int global_row = proc_row * blockSize + i;
                    int global_col = proc_col * blockSize + j;
                    Ctemp[global_row][global_col] = gathered_data[current_pos++];
                }
            }
        }
        // Skopiuj poprawnie zrekonstruowaną macierz z Ctemp do cGlob
        memcpy(cGlob[0], Ctemp[0], N * N * sizeof(double));
        freeMatrix(Ctemp); 

        // printf("Rekonstrukcja zakończona.\n");

        // --- Weryfikacja ---
        printf("Czas obliczeń: %f sekund\n", end_cGlob - start_cGlob);

        // printf("Rozpoczynanie obliczeń sekwencyjnych dla weryfikacji...\n");
        start_cSek = MPI_Wtime();
        for (int i = 0; i < N; i++){
            for (int k = 0; k < N; k++) {
                for (int j = 0; j < N; j++) {
                    cSek[i][j] += rawA[i][k] * rawB[k][j];
                }
            }
        }
        end_cSek = MPI_Wtime();
        printf("Czas obliczeń sekwencyjnych (cSek): %f sekund\n", end_cSek - start_cSek);

        // Porównanie wyników cGlob (z MPI) i cSek (sekwencyjnie)
        // printf("Porównywanie wyników cGlob i cSek...\n");
        int errors = 0;
        double max_diff = 0.0;
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                double diff = fabs(cSek[i][j] - cGlob[i][j]);
                if (diff > 1e-9) { 
                    errors++;
                    if (diff > max_diff) max_diff = diff;
                }
            }
        }
        if (errors == 0) {
            printf("Weryfikacja zakończona: Wyniki są POPRAWNE!\n");
        } else {
            printf("Weryfikacja zakończona: Znaleziono %d błędów w wynikach. Maksymalna różnica: %e\n", errors, max_diff);
        }

        freeMatrix(cGlob);
        freeMatrix(cSek);
    } 

    MPI_Finalize();
    return 0;
}