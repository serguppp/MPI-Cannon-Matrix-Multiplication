#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    FILE *plik, *plik2;

    // Określenie ścieżki w zależności od systemu operacyjnego
    #ifdef _WIN32  // Windows
        const char *path1 = "matrix1.txt";
        const char *path2 = "matrix2.txt";
    #else  // Linux
        const char *path1 = "./matrix1.txt";
        const char *path2 = "./matrix2.txt";
    #endif

    // Otwieranie plików
    plik = fopen(path1, "w");
    if (plik == NULL) {
        perror("Błąd otwarcia pliku");
        exit(EXIT_FAILURE);
    }

    plik2 = fopen(path2, "w");
    if (plik2 == NULL) {
        perror("Błąd otwarcia pliku");
        exit(EXIT_FAILURE);
    }

    printf("Pliki otwarte poprawnie, zapisuję dane\n");

    // Zapis danych do plików
    for (int i = 0;  i < 4008; i++) {
        for (int j = 0; j < 4008; j++) {
            fprintf(plik, "%10.2lf ", (double)i + j);
            fprintf(plik2, "%10.2lf ", (double)i + j);
        }
        fprintf(plik, "\n");
        fprintf(plik2, "\n");
    }

    // Zamknięcie plików
    fclose(plik);
    fclose(plik2);
    printf("Pliki zostały zapisane i zamknięte poprawnie\n");

    MPI_Finalize();

    return 0;
}
