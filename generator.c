#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv) {
    FILE *plik, *plik2;
    int i, j;

    plik = fopen("C:/Users/sergi/Desktop/mpi_projekt/MPI-Cannon-Matrix-Multiplication/matrix1.txt", "w");
    if (plik == NULL) {
        perror("Błąd otwarcia pliku");
        exit(EXIT_FAILURE);
    }

    plik2 = fopen("C:/Users/sergi/Desktop/mpi_projekt/MPI-Cannon-Matrix-Multiplication/matrix2.txt", "w");
    if (plik2 == NULL) {
        perror("Błąd otwarcia pliku");
        exit(EXIT_FAILURE);
    }
    printf("Pliki otwarte poprawnie, zapisuje dane\n");

    for (i = 0; i < 192; i++) {
        for (j = 0; j < 192; j++){
            fprintf(plik, "%6.1f ", (float)i + j);
            fprintf(plik2, "%6.1f ", (float)i + j);
        }
        fprintf(plik, "\n");
        fprintf(plik2, "\n");
    }

    fclose(plik);
    fclose(plik2);
    printf("Pliki został zapisane i zamknięte poprawnie\n");

    return 0;
}
