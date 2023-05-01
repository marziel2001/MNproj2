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
Macierz::Macierz(const Macierz& const N)
{
    this->n = N.n;
    this->M = new double* [this->n];

    for (int i = 0; i < this->n; i++)
    {
        this->M[i] = new double[this->n];
    }
    
    for (int i = 0; i < this->n; i++)
    {
        for (int j = 0; j < this->n; j++)
        {
            int a = N.M[i][j];
            this->M[i][j] = a;
        }
    }

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

Macierz*& Macierz::Jacobi(Macierz *B)
{
    int iteracje = 0;
    double norm = 0.0;
    Macierz * Xprev = new Macierz(this->n);
    Macierz * X = new Macierz(this->n);
    X->InicjalizujWektor(this->n, 1);
    Xprev->InicjalizujWektor(this->n, 1);

    Macierz* residuum;
    Macierz* A = this;

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
        for (int i = 0; i < this->n; i++) Xprev->V[i] = X->V[i];

        residuum = Macierz::Residuum(A, X, B);
        norm = Macierz::LiczNormeRes(residuum);
        
    } while (norm > Macierz::docelowaNorma);

    cout << "iteracje: " << iteracje << "\n";
    cout << "norma residuum dla Jacobiego: " << norm << "\n";

    return X;
}

Macierz*& Macierz::GaussSeidl(Macierz*& B)
{
    int iteracje = 0;
    double norm = 0.0;
    Macierz* Xminus1 = new Macierz(this->n);
    Macierz* X = new Macierz(this->n);
    X->InicjalizujWektor(this->n, 1);

    Macierz* residuum;
    Macierz* A = this;

    Xminus1->InicjalizujWektor(this->n, 1);

    do
    {
        iteracje++;
        for (int i = 0; i < A->n; i++)
        {
            double suma = 0;

            for (int j = 0; j <i; j++)
            {
                suma -= A->M[i][j] * X->V[j];
            }
            for (int j = i+1; j < A->n; j++)
            {
                suma -= A->M[i][j] * Xminus1->V[j];
            }

            X->V[i] = (B->V[i] + suma) / A->M[i][i];
        }
        for (int i = 0; i < this->n; i++) Xminus1->V[i] = X->V[i];

        residuum = Macierz::Residuum(A, X, B);
        norm = Macierz::LiczNormeRes(residuum);

    } while (norm > Macierz::docelowaNorma);

    cout << "iteracje Gauss Seidel: " << iteracje << "\n";
    cout << "norma residuum dla GaussSeidel: " << norm << "\n";

    return X;
}

Macierz*& Macierz::faktoryzacjaLU(Macierz*& B)
{
    
    Macierz* U = new Macierz(*this);
    Macierz* L = new Macierz(this->n, this->n);


    for (int i = 0; i < this->n; i++) L->M[i][i] = 1;

   // U->Drukuj();
    //L->Drukuj();

    //int iteracje = 0;
    double norm = 0.0;

    Macierz* X = new Macierz(this->n);

    X->InicjalizujWektor(this->n, 0);

    Macierz* residuum;
    


    // tworzenie macierzy L i U
    for (int k = 0; k < this->n - 1; k++)
    {



        for (int j = k + 1; j < this->n; j++)
        {



            L->M[j][k] = U->M[j][k] / U->M[k][k];

            for (int i = k; i < this->n; i++)
            {
                U->M[j][i] = U->M[j][i] - L->M[j][k] * U->M[k][i];
            }
        }
    }



    Macierz* Y = new Macierz(this->n);


    for (int i = 0; i < this->n; i++)
    {
        double suma = 0;

        for (int j = 0; j < i; j++)
        {
            suma -= L->M[i][j] * Y->V[j];
        }

        Y->V[i] = (B->V[i] + suma) / L->M[i][i];
    }

    

    for (int i = this->n - 1; i >=0; i--)
    {
        double suma = 0;

        for (int j = i+1; j < U->n; j++)
        {
            suma -= U->M[i][j] * X->V[j];
        }

        X->V[i] = (Y->V[i] + suma) / U->M[i][i];
    }

    
    Macierz* A = this;

    residuum = Macierz::Residuum(A, X, B);
    norm = Macierz::LiczNormeRes(residuum);

    cout << "norma residuum dla faktoryzacji LU: " << norm << "\n";
    return X;
}

Macierz*& Macierz::Residuum(Macierz*& A, Macierz*& X, Macierz*& B)
{
    Macierz * residuum = Macierz::MnozeniePrzezWektor(A, X);
    for (int i = 0; i < residuum->n; i++) residuum->V[i] -= B->V[i];
    return residuum;
}
double Macierz::LiczNormeRes(Macierz*& Residuum)
{
    double suma = 0.0;
    for (int i = 0; i < this->n; i++)
    {
        suma += pow(Residuum->V[i], 2);
    }
    return sqrt(suma);
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
