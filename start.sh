#!/bin/bash

for i in $(seq 1 10)
do
  echo "Uruchomienie numer: $i"
  mpirun  --oversubscribe -np 25 ./main
  echo "-------------------------------------"
done

echo "Zakończono 10 uruchomień."