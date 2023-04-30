#include <iostream>
#include "Macierz.h"

//#define N 905
#define N 905
#define A1 5
#define A2 -1
#define A3 -1

using std::cout;

int main()
{
    Macierz * A = new Macierz(N, N);
    
    // utworzono macierz zer "A" - o rozmiarah N x N
    A->UstawDiagonale(A1, A2, A3, N);
    //A->Drukuj();

    Macierz* B = new Macierz(N); 
    B->V[B->n - 1] = sin(N * 2);
    //B->Drukuj();

    Macierz* X = new Macierz(N);
    //X->Drukuj();   
    
    Macierz* Wynik = A->Jacobi(B);
    Wynik->Drukuj();

    return 0;
}
