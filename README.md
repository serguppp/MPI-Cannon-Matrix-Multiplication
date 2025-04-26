# MPI-Cannon-Matrix-Multiplication

## Opis
Mnożenie macierzy metodą Cannona z użyciem wielu procesorów.

## Użycie

1. Sklonuj repozytorium:
   ```sh
   git clone https://github.com/serguppp/MPI-Cannon-Matrix-Multiplication
   ```

2. Skonfiguruj generator macierzy:
   ```sh
   mpicc generator.c -o generator
   ```

3. Uruchom generator macierzy:
   ```sh
   mpirun -np 1 ./generator
   ```

4. Skompiluj główny program:
   ```sh
   mpicc main.c -O3 -o main
   ```

5. Uruchom mnożenie macierzy, określając liczbę procesorów `P`:
   ```sh
   mpirun -np P main
   ```
