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

	Macierz(int rows, int cols);
	Macierz(int rows);
	~Macierz();

	void Drukuj();
	void UstawDiagonale(double a1, double a2, double a3, int rozmiar);
	void StworzMacierzeLUD(Macierz*& L, Macierz*& U, Macierz*& D);
	void InicjalizujMacierz(int rows, int cols, double wartosc);
	void InicjalizujWektor(int rows, double wartosc);

	// Dzialania
	static void Dodawanie(Macierz*& A, Macierz*& B);

	void Jacobi(Macierz*& B, Macierz*& X, Macierz*& R);
	void GaussSeidl();
	void faktoryzacjaLU();

private:
	typ typ = macierz;



};
