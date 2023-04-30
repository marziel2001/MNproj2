#include "Macierz.h"
#include <iostream>

using std::cout;

double Macierz::docelowaNorma = 1e-9;

Macierz::Macierz(int rows, int cols)
{
    this->n = rows;

    typ = macierz;

    M = new double* [rows];
    for (int i = 0; i < rows; i++)
    {
        this->M[i] = new double[cols];
    }

    this->InicjalizujMacierz(rows, cols, 0);
}
Macierz::Macierz(int rows)
{
    typ = wektor;
    this->n = rows; 

    V = new double[rows];

    this->InicjalizujWektor(rows, 0);
}
Macierz::~Macierz()
{
    if (typ == macierz)
    {
      
        for (int i = 0; i < this->n; i++)
        {
            delete[] this->M[i];
        }

        delete[] this->M;
    }
    else if (typ == wektor)
    {
        delete[] this->V;
    }
}

void Macierz::InicjalizujMacierz(int rows, int cols, double wartosc)
{
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            this->M[i][j] = wartosc;
        }
    }
}

void Macierz::InicjalizujWektor(int rows, double wartosc)
{
    for (int i = 0; i < rows; i++)
    {
        this->V[i] = wartosc;
    }
}

//void Macierz::StworzMacierzeLUD(Macierz*& L, Macierz*& U, Macierz*& D)
//{
//    for (int i = 0; i < this->n; i++)
//    {
//        for (int j = 0; j < this->n; j++)
//        {
//            if (i > j)
//            {
//                L->M[i][j] = this->M[i][j];
//            }
//            else if (i < j)
//            {
//                U->M[i][j] = this->M[i][j];
//            }
//            else
//            {
//                D->M[i][j] = this->M[i][j];
//            }
//        }
//    }
//}

Macierz*& Macierz::Jacobi(Macierz*& B)
{
    int iteracje = 0;
    double norm = 0.0;
    Macierz * Xprev = new Macierz(this->n);
    Macierz * X = new Macierz(this->n);
    X->InicjalizujWektor(this->n, 1);

    Macierz * Residuum = new Macierz(this->n);
    Macierz* A = this;

    Xprev->InicjalizujWektor(this->n, 1);

    do
    {
        iteracje++;
        for (int i = 0; i < A->n; i++)
        {
            double suma = 0;

            for (int j = 0; j < A->n; j++)
            {
                if (j != i)
                {
                    suma -= A->M[i][j] * Xprev->V[j];
                }
            }      
            X->V[i] = (B->V[i] + suma) / A->M[i][i];

        }
        // liczymy residuum i norm
        Residuum = Macierz::MnozeniePrzezWektor(A, X);

        for (int i = 0; i < A->n; i++)
        {
            Residuum->V[i] -= B->V[i];
        }

        for (int i = 0; i < this->n; i++)
        {
            Xprev->V[i] = X->V[i];
        }

        double suma = 0.0;
        for (int i = 0; i < this->n; i++)
        {
            suma += pow(Residuum->V[i], 2);
        }
        norm = sqrt(suma);

        if (iteracje > 10000) break;
    } while (norm > Macierz::docelowaNorma);

    cout << "iteracje: " << iteracje << "\n";

    return X;
    ///exp1 = -D\(L + U)
    ///exp2 = D\b;
    /// 
    /// norma to bêdzie pierwiastek z sumy kwadratów wektora residuum
}



void Macierz::GaussSeidl()
{
    cout << "Hello GS";

}

Macierz*& Macierz::Dodawanie(Macierz*& A, Macierz*& B)
{
    Macierz* C = new Macierz(A->n, A->n);

    for (int i = 0; i < A->n; i++)
    {
        for (int j = 0; j < A->n; j++)
        {
            C->M[i][j] = B->M[i][j] + A->M[i][j];
        }
    }

    return C;
}
Macierz*& Macierz::MnozenieMacierzy(Macierz*& A, Macierz*& B)
{
    Macierz* C = new Macierz(A->n, A->n);


    for (int j = 0; j < B->n; j++)
    {  
        for (int i = 0; i < A->n; i++)
        {
            double s = 0;

            for (int k = 0; k < B->n; k++)
{
                s += A->M[i][k] * B->M[k][j];
            }

            C->M[i][j] = s;
        }
    }

    return C;
}
Macierz*& Macierz::MnozeniePrzezWektor(Macierz*& A, Macierz*& V)
{
    Macierz* C = new Macierz(A->n);


        for (int i = 0; i < A->n; i++)
{
            double s = 0;

            for (int k = 0; k < V->n; k++)
            {
                s += A->M[i][k] * V->V[k];
            }

            C->V[i] = s;
        }

    return C;
}
Macierz*& Macierz::RazyMinusJeden()
{

    Macierz* C;

    if (this->typ == macierz)
    {   
        C = new Macierz(this->n, this->n);

        for (int i = 0;i<this->n; i++)
        {
            for (int j = 0; j < this->n; j++)
            {
                C->M[i][j] = -this->M[i][j];
            }
        }
    }
    else
    {
        C = new Macierz(this->n);

        for (int j = 0; j < this->n; j++)
        {
            C->V[j] = -this->V[j];
        }
    }

    return C;
}

void Macierz::faktoryzacjaLU()
{

}


void Macierz::Drukuj()
{
    if (typ == macierz)
    {
        cout << "\n[\n ";
        for (int i = 0; i < this->n; i++)
        {
            cout << "\t[";
            for (int j = 0; j < this->n; j++)
            {
                cout << M[i][j] << ", ";
            }
            cout << "];\n";
        }
        cout << "]\n";
    }
    else if (typ == wektor)
    {
        cout << "\n[ ";
        for (int i = 0; i < this->n; i++)
        {
            cout << V[i] << ", ";
            
        }
        cout << " ]\n";
    }
}

void Macierz::UstawDiagonale(double a1, double a2, double a3, int rozmiar)
{
    for (int i = 0; i < rozmiar; i++)
    {
        for (int j = 0; j < rozmiar; j++)
        {
            if (i == j)
            {
                M[i][j] = a1;
            }
            if (i == j - 1 || i == j + 1)
            {
                M[i][j] = a2;
            }
            if (i == j - 2 || i == j + 2)
            {
                M[i][j] = a3;
            }
        }
    }
}
