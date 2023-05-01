#pragma once

class Macierz
{
	enum typ
	{
		macierz, 
		wektor
	};

public:

	double** M = nullptr;
	double* V = nullptr;

	int n = 0;
	static double docelowaNorma;

	Macierz(int rows, int cols);
	Macierz(int rows);
	Macierz(const Macierz& const N);
	~Macierz();

	void Drukuj();
	void UstawDiagonale(double a1, double a2, double a3, int rozmiar);
	//void StworzMacierzeLUD(Macierz*& L, Macierz*& U, Macierz*& D);
	void InicjalizujMacierz(int rows, int cols, double wartosc);
	void InicjalizujWektor(int rows, double wartosc);

	// Dzialania
	static Macierz*& Dodawanie(Macierz*& A, Macierz*& B);
	static Macierz*& MnozeniePrzezWektor(Macierz*& A, Macierz*& V);
	static Macierz*& MnozenieMacierzy(Macierz*& A, Macierz*& B);
	Macierz*& RazyMinusJeden();

	Macierz*& Jacobi( Macierz* B);
	Macierz*& GaussSeidl(Macierz*& B);
	Macierz*& faktoryzacjaLU(Macierz*& B);
	static Macierz*& Residuum(Macierz*& A, Macierz*& X, Macierz*& B);

	double LiczNormeRes(Macierz*& R);

private:
	typ typ = macierz;



};
