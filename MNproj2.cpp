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

    Macierz* B = new Macierz(N); 
    for(int i = 0; i<N; i++) B->V[i] = sin(i * 2);
    

    cout << "\n\n\n=============== Jacobi =================\n";
    Macierz* WynikJ = A->Jacobi(B);
    //WynikJ->Drukuj();


    cout << "\n\n\n=============== GaussSeidl =================\n";
    Macierz* WynikGS = A->GaussSeidl(B);
    //WynikGS->Drukuj();


    cout << "\n\n\n=============== FaktoryzacjaLU =================\n";
    Macierz* WynikLU = A->faktoryzacjaLU(B);
    //WynikGS->Drukuj();

    return 0;
}
