#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

#define N 2000 // Rozmiar macierzy
#define PP 2// Pierwiastek z liczby procesów
#define P 4// Liczba procesów

// Macierze A i B wczytane przez proces 0
double rawA[N][N], rawB[N][N];
// Macierze A i B po poczatkowej dystrybucji (przesunięciu) obliczone przez proces 0
double distrA[N][N], distrB[N][N];

// Lokalne bloki macierzy dla każdego procesu
double a[N / PP][N / PP], b[N / PP][N / PP], c[N / PP][N / PP];
// Bufory pomocnicze do odbierania danych w pętli Cannona
double aa[N / PP][N / PP], bb[N / PP][N / PP];

// Wskaźniki do buforów wysyłania/odbierania w pętli Cannona
double(*psa)[N / PP], (*psb)[N / PP], (*pra)[N / PP], (*prb)[N / PP];

// Macierz wynikowa obliczona sekwencyjnie (dla weryfikacji)
double CSek[N][N];
// Globalna macierz wynikowa zbierana w procesie 0
double Cglob[N][N];

int rank, np, finalize; // Zmienna 'finalize' do obsługi błędów plików

double startwtime1, startwtime2, endwtime;

// Funkcja wczytywania plików 
void loadFile(FILE **file, const char *path) {
    if (rank == 0) { // Tylko proces 0 otwiera pliki
        *file = fopen(path, "r");
        if (*file == NULL) {
            perror("Błąd otwarcia pliku");
            finalize = 1; // Ustaw flagę błędu
        } else {
            printf("Proces 0 poprawnie otworzył plik \"%s\"\n", path);
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

// Funkcja wczytywania danych z pliku do macierzy 
void loadMatrix(FILE **file, const char *path, const char id) {
    if (rank == 0) { // Tylko proces 0 wczytuje dane
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                if (id == 'A') {
                    if (fscanf(*file, "%lf", &rawA[i][j]) != 1) {
                         fprintf(stderr, "Błąd odczytu danych z pliku %s\n", path);
                         fclose(*file);
                         MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE); // Zakończ wszystkie procesy
                    }
                } else if (id == 'B') {
                     if (fscanf(*file, "%lf", &rawB[i][j]) != 1) {
                         fprintf(stderr, "Błąd odczytu danych z pliku %s\n", path);
                         fclose(*file);
                         MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE); // Zakończ wszystkie procesy
                    }
                }
            }
        }
        fclose(*file); // Zamknij plik po wczytaniu
    }
}

// Funkcja TESTUJĄCA do wypisywania macierzy
void printMatrix(double matrix[N][N], int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            printf("%6.1f ", matrix[i][j]);
        }
        printf("\n");
    }
}

// Funkcja zajmująca się wstępną inicjalizacją macierzy

void initMatrix() {
    int blockSize = N / PP; // Rozmiar jednego bloku

    double tempA[N][N], tempB[N][N]; // Bufory pomocnicze do przesunięć bloków

    // Kopiujemy początkowe dane do buforów
    memcpy(tempA, rawA, sizeof(tempA));
    memcpy(tempB, rawB, sizeof(tempB));

    // Przesuwanie bloków zgodnie z metodą Cannona
    for (int blockRow = 0; blockRow < PP; blockRow++) { 
        for (int blockCol = 0; blockCol < PP; blockCol++) { 
            // Obliczamy nowe indeksy bloków po przesunięciu
            int newColA = (blockCol - blockRow + PP) % PP; // A przesuwamy w lewo o blockRow
            int newRowB = (blockRow - blockCol + PP) % PP; // B przesuwamy w górę o blockCol

            // Kopiujemy **całe bloki**
            for (int i = 0; i < blockSize; i++) {
                memcpy(
                    &distrA[blockRow * blockSize + i][newColA * blockSize], // Miejsce docelowe w distrA
                    &tempA[blockRow * blockSize + i][blockCol * blockSize], // Źródło w tempA
                    blockSize * sizeof(double) // Kopiuj cały wiersz bloku
                );

                memcpy(
                    &distrB[newRowB * blockSize + i][blockCol * blockSize], 
                    &tempB[blockRow * blockSize + i][blockCol * blockSize], 
                    blockSize * sizeof(double)
                );
            }
        }
    }
    
}

// Funkcja wypisująca blok procesów dla macierzy a, b i c
void printBlock(double matrix[N / PP][N / PP], const char* name) {
    printf("Proces %d, Macierz %s:\n", rank, name);
    for (int i = 0; i < N / PP; i++) {
        for (int j = 0; j < N / PP; j++) {
            printf("%f ", matrix[i][j]);
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
    int row = rank / PP; // Wiersz procesu w siatce PPxPP
    int col = rank % PP; // Kolumna procesu w siatce PPxPP
    int mod = 0;
    int tag = 101; // Tag dla wiadomości MPI
    int blockSize = N / PP; // Rozmiar boku bloku lokalnego
    int localSize = blockSize * blockSize; // Liczba elementów w bloku lokalnym

    MPI_Status statSend[6], statRecv[6];
    MPI_Request reqSend[6], reqRecv[6];

    //Sprawdzamy, czy liczba uruchomionych procesów MPI (np) jest zgodna z oczekiwaną liczbą procesów (P)
    // i wykonuje odpowiednią akcję, jeśli liczba procesów jest nieprawidłowa.
    if (np != P) {
        if (rank == 0)
            printf("Ilość procesorów MPI niezgodna z oczekiwaną liczbą procesów. Uruchom mpiexec -n %d program\n", P);
        MPI_Finalize();
        exit(0);
    }

    if (rank == 0)
    {
        printf("Obliczenia metodą Cannona dla macierzy %d x %d\n", N, N);
        startwtime1 = MPI_Wtime();
        // wczytanie danych 
        loadFile(&fileA, pathA);        
        loadFile(&fileB, pathB);    

        // Odczyt danych wejściowych z plików
        loadMatrix(&fileA, pathA, 'A');
        loadMatrix(&fileB, pathB, 'B');

        //Dystrybucja poczatkowa blokow
        initMatrix();

       /* printf("Dystrybuowana macierz a\n");
        printMatrix(distrA,N,N);
        printf("Dystrybuowana macierz b\n");
        printMatrix(distrB,N,N);
        */

        // Rozesłanie odpowiednich bloków z distrA i distrB do pozostałych procesów
        for (int p = 0; p < P; p++){
            int row = p % PP; // numer kolumny bloku docelowego
            int col = p / PP; // numer wiersza bloku docelowego

           // Tymczasowe bufory na bloki do wysłania/skopiowania
           double subA[blockSize][blockSize];
           double subB[blockSize][blockSize];

            // Kopiowanie odpowiedniego bloku z przesuniętych macierzy distrA/distrB
            for (int x = 0; x < blockSize; x++){
                for (int y = 0; y < blockSize; y++){
                    subA[x][y] = distrA[col * blockSize + x][row * blockSize + y];
                    subB[x][y] = distrB[col * blockSize + x][row * blockSize + y];
                }
            }
            // Wyślłanie bloków do procesów 'p' lub skopiowanie lokalnie dla p=0
            if (p != 0){
                MPI_Send(subA, localSize, MPI_DOUBLE, p, tag, MPI_COMM_WORLD);
                MPI_Send(subB, localSize, MPI_DOUBLE, p, tag + 1, MPI_COMM_WORLD); // Użyj innego tagu dla B
            }
            // Proces 0 kopiuje swoje bloki bezpośrednio
            else{
                memcpy(a, subA, sizeof(subA));
                memcpy(b, subB, sizeof(subB));
            }
            //printf("Proces %d dostanie\n", p);
            //printBlock(subA,"subA");
            //printBlock(subB,"subB");
        }
    } else {
        // Procesy różne od rank=0 odbierają tablicę a i b
        MPI_Recv(a, localSize, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(b, localSize, MPI_DOUBLE, 0, tag + 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }


    // Zerowanie macierzy pomocniczych i wynikowej 
    memset(aa, 0, sizeof(aa));
    memset(bb, 0, sizeof(bb));
    memset(c, 0, sizeof(c));

    if (rank == 0) startwtime2 = MPI_Wtime();

    //Algorytm Cannona

    pra = aa; prb = bb; // Bufory odbiorcze
    psa = a; psb = b;   // Bufory do wysłania (początkowo)
    
    for (int kk = 0; kk < PP; kk++) {   
        // Krok 1: Mnożenie lokalnych bloków
        for (int i = 0; i < blockSize; i++)
            for (int j = 0; j < blockSize; j++)
                for (int k = 0; k < blockSize; k++)
                    c[i][j] += psa[i][k] * psb[k][j];

        // Krok 2: Przesunięcie bloków (komunikacja)
        int left = (col - 1 + PP) % PP + row * PP; 
        int right = (col + 1) % PP + row * PP;
        int up = col + ((row - 1 + PP) % PP) * PP;
        int down = col + ((row + 1) % PP) * PP;
        
        // Wyślij obecny blok A w lewo, odbierz nowy z prawej do bufora odbiorczego A (pra)
        MPI_Sendrecv(psa, localSize, MPI_DOUBLE, left, tag,
            pra, localSize, MPI_DOUBLE, right, tag,
            MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // Wyślij obecny blok B w górę, odbierz nowy z dołu do bufora odbiorczego B (prb)
        MPI_Sendrecv(psb, localSize, MPI_DOUBLE, up, tag + 1, // Inny tag dla B
                    prb, localSize, MPI_DOUBLE, down, tag + 1, // Inny tag dla B
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // Krok 3: Zamiana wskaźników buforów 
        // Teraz odebrane bloki (w pra, prb) stają się blokami do wysłania (psa, psb) w następnej iteracji
        // a bloki, które właśnie wysłaliśmy (w psa, psb) stają się buforami odbiorczymi (pra, prb)
        double(*temp_ptr)[N / PP];
        temp_ptr = psa; psa = pra; pra = temp_ptr;
        temp_ptr = psb; psb = prb; prb = temp_ptr;

    }

    // Zbieranie wyników do procesu 0 
    // Uwaga: Cglob jest [N][N], ale MPI_Gather oczekuje liniowego bufora.
    // Można zbierać do tymczasowego bufora liniowego lub bezpośrednio do Cglob traktując go jak liniowy.
    MPI_Gather(c, localSize, MPI_DOUBLE,       // Dane wysyłane: lokalny blok c
        Cglob, localSize, MPI_DOUBLE,          // Dane odbierane: globalna macierz Cglob (traktowana liniowo)
        0, MPI_COMM_WORLD);                   // Odbiera tylko proces 0

    // Proces 0 rekonstruuje macierz Cglob z zebranych bloków
    if (rank == 0) {
        double Ctemp[N][N]; // Bufor tymczasowy na poprawnie ułożoną macierz
        int current_pos = 0; // Pozycja w liniowym buforze Cglob

        for (int p = 0; p < P; p++) {
            // Oblicz współrzędne bloku (proc_row, proc_col) dla procesu 'p'
            int proc_row = p / PP;
            int proc_col = p % PP;

            // Skopiuj dane z liniowego Cglob do Ctemp w odpowiednie miejsce
            for (int i = 0; i < blockSize; i++) { // i - wiersz wewnątrz bloku
                for (int j = 0; j < blockSize; j++) { // j - kolumna wewnątrz bloku
                    int global_row = proc_row * blockSize + i;
                    int global_col = proc_col * blockSize + j;
                    // Używamy wskaźnika do Cglob traktowanego jako double*
                    Ctemp[global_row][global_col] = ((double*)Cglob)[current_pos++];
                }
            }
        }
        // Skopiuj poprawnie ułożoną macierz z Ctemp do Cglob
        memcpy(Cglob, Ctemp, sizeof(Ctemp));
    }

    //printf("Ukonczono mnozenie macierzy dla rank %d\n", rank);
    //finalize zliczania czasow
    if (rank == 0) {
        endwtime = MPI_Wtime();
        printf("Czas przetwarzania: %f sekund\n", endwtime - startwtime1);
        printf("Czas obliczeń: %f sekund\n", endwtime - startwtime2);
    
        // Weryfikacja wyników przez porównanie z obliczeniami sekwencyjnymi
        // Obliczenia sekwencyjne mnożenia tablic CSek=A*B
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++) {
                CSek[i][j] = 0;
                for (int k = 0; k < N; k++) {
                    CSek[i][j] += rawA[i][k] * rawB[k][j];
                }
            }
    
        // Porównanie wyników Cglob (z MPI) i CSek (sekwencyjnie)
        int errors = 0;
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                double diff = fabs((double)CSek[i][j] - (double)Cglob[i][j]);
                if (fabs(CSek[i][j] - Cglob[i][j])>1e-5) {
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

        //printf("Macierz wynikowa Cglob:\n");
        // printMatrix(Cglob, N, N);
        // printf("Macierz wynikowa Csek:\n");
        // printMatrix(CSek, N, N);
        printf("Program MPI zakończony poprawnie!\n");
    }

    MPI_Finalize();
    return 0;
}