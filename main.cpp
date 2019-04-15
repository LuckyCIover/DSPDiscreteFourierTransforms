#define _USE_MATH_DEFINES

#include <algorithm> 
#include <iostream>
#include <fstream>
#include <complex>
#include <string>
#include <vector>
#include <math.h>
#include <limits>
#include <string>
#include <time.h> 

using namespace std;

using z_type = complex<double>;

vector <z_type> parseArray(char * fileName)
{
	ifstream inputFile(fileName);

	vector <z_type> inputDataVector;

	double Re, Im;

	if (inputFile.is_open())
	{
		while (inputFile >> Re >> Im)
			inputDataVector.push_back(z_type(Re,Im));
	}
	else
	{
		cout << "Error opening input file(s)" << endl;
		exit(1);
	}

	inputFile.close();

	return inputDataVector;
}

void writeOutputToFile(char * fileName, vector <z_type> inputVectorX)
{
	size_t N = inputVectorX.size();

	ofstream outputFile(fileName);

	if (outputFile.is_open())
	{
		outputFile.precision(numeric_limits<double>::max_digits10 + 1);
		for (size_t n = 0; n < N; n++)
			outputFile << inputVectorX[n].real() << " " << inputVectorX[n].imag() << endl;
	}
	else
	{
		cout << "Error opening output file" << endl;
		exit(1);
	}

	outputFile.close();
}

vector <z_type> DFT(vector <z_type> inputVectorX)
{
	size_t N = inputVectorX.size();

	vector <z_type> outputVectorY(N, 0);

	clock_t start = clock();

	for (size_t k = 0; k < N; k++)
	{
		z_type tempSum = 0;
		for (size_t j = 0; j < N; j++)
		{
			tempSum += inputVectorX[j] * exp(-(2 * M_PI * 1i * (double)k * (double)j)/(double)N);
		}
		outputVectorY[k] = (1 / sqrt((double)N)) * tempSum;
	}

	clock_t end = clock();

	double seconds = (double)(end - start) / CLOCKS_PER_SEC;

	cout.precision(numeric_limits<double>::max_digits10 + 1);
	cout << "Computing time: " << seconds << " seconds" << endl;

	return outputVectorY;
}

vector <z_type> IDFT(vector <z_type> inputVectorY)
{
	size_t N = inputVectorY.size();

	vector <z_type> outputVectorX(N, 0);

	for (size_t k = 0; k < N; k++)
	{
		z_type tempSum = 0;
		for (size_t j = 0; j < N; j++)
		{
			tempSum += inputVectorY[j] * exp((2 * M_PI * 1i * (double)k * (double)j) / (double)N);
		}
		outputVectorX[k] = (1 / sqrt((double)N)) * tempSum;
	}

	return outputVectorX;
}

vector <z_type> FFT_wop(vector <z_type> inputVectorX)
{
	int n = (int)log2((double)inputVectorX.size());

	size_t N = inputVectorX.size();
	
	vector <z_type> outputVectorY(N, 0);

	z_type omega;

	clock_t start = clock();

	for (int k = 1; k < n + 1; k++)
	{
		int jLimit = (int)pow(2, k - 1);
		int lLimit = (int)pow(2, n - k);
		for (int j = 0; j < jLimit; j++)
		{
			omega = exp(-(2 * j * M_PI * 1i ) / pow(2, k));
			for (int l = 0; l < lLimit; l++)
			{
				outputVectorY[j * (int)pow(2, n - k) + l] = (inputVectorX[j * (int)pow(2, n - k + 1) + l] +
					omega * inputVectorX[j * (int)pow(2, n - k + 1) + l + (int)pow(2, n - k)]) / sqrt(2);

				outputVectorY[(int)pow(2, n - 1) + j * (int)pow(2, n - k) + l] = (inputVectorX[j * (int)pow(2, n - k + 1) + l] -
					omega * inputVectorX[j * (int)pow(2, n - k + 1) + l + (int)pow(2, n - k)]) / sqrt(2);
			}
		}

		inputVectorX = outputVectorY;
	}

	clock_t end = clock();

	double seconds = (double)(end - start) / CLOCKS_PER_SEC;

	cout.precision(numeric_limits<double>::max_digits10 + 1);
	cout << "Computing time: " << seconds << " seconds" << endl;

	return outputVectorY;
}

vector <z_type> IFFT_wop(vector <z_type> inputVectorY)
{
	int n = (int)log2(inputVectorY.size());

	size_t N = inputVectorY.size();

	vector <z_type> outputVectorX(N, 0);

	z_type omega;

	for (int k = 1; k < n + 1; k++)
	{
		int jLimit = (int)pow(2, k - 1);
		int lLimit = (int)pow(2, n - k);
		for (int j = 0; j < jLimit; j++)
		{
			omega = exp((2 * j * M_PI * 1i) / pow(2, k));
			for (int l = 0; l < lLimit; l++)
			{
				outputVectorX[j * (int)pow(2, n - k) + l] = (inputVectorY[j * (int)pow(2, n - k + 1) + l] +
					omega * inputVectorY[j * (int)pow(2, n - k + 1) + l + (int)pow(2, n - k)]) / sqrt(2);

				outputVectorX[(int)pow(2, n - 1) + j * (int)pow(2, n - k) + l] = (inputVectorY[j * (int)pow(2, n - k + 1) + l] -
					omega * inputVectorY[j * (int)pow(2, n - k + 1) + l + (int)pow(2, n - k)]) / sqrt(2);
			}
		}

		inputVectorY = outputVectorX;
	}

	return outputVectorX;
}

vector <z_type> FFT_wp(vector <z_type> inputVectorX)
{
	int n = (int)log2(inputVectorX.size());

	size_t N = inputVectorX.size();

	vector <z_type> outputVectorY(N, 0);

	z_type omega;

	for (int k = 1; k < n + 1; k++)
	{
		int jLimit = (int)pow(2, k - 1);
		int lLimit = (int)pow(2, n - k);
		for (int j = 0; j < jLimit; j++)
		{
			for (int l = 0; l < lLimit; l++)
			{
				omega = exp(-(2 * l * M_PI * 1i) / pow(2, k));

				outputVectorY[j * (int)pow(2, k) + l] = (inputVectorX[j * (int)pow(2, k) + l] +
					omega * inputVectorX[j * (int)pow(2, k) + l + (int)pow(2, k - 1)]) / sqrt(2);

				outputVectorY[(int)pow(2, k - 1) + j * (int)pow(2, k) + l] = (inputVectorX[j * (int)pow(2, k) + l] -
					omega * inputVectorX[j * (int)pow(2, k) + l + (int)pow(2, k - 1)]) / sqrt(2);
			}
		}

		inputVectorX = outputVectorY;
	}

	return outputVectorY;
}

vector <z_type> IFFT_wp(vector <z_type> inputVectorY)
{
	int n = (int)log2(inputVectorY.size());

	size_t N = inputVectorY.size();

	vector <z_type> outputVectorX(N, 0);

	z_type omega;

	for (int k = 1; k < n + 1; k++)
	{
		int jLimit = (int)pow(2, k - 1);
		int lLimit = (int)pow(2, n - k);
		for (int j = 0; j < jLimit; j++)
		{
			for (int l = 0; l < lLimit; l++)
			{
				omega = exp((2 * l * M_PI * 1i) / pow(2, k));

				outputVectorX[j * (int)pow(2, k) + l] = (inputVectorY[j * (int)pow(2, k) + l] +
					omega * inputVectorY[j * (int)pow(2, k) + l + (int)pow(2, k - 1)]) / sqrt(2);

				outputVectorX[(int)pow(2, k - 1) + j * (int)pow(2, k) + l] = (inputVectorY[j * (int)pow(2, k) + l] -
					omega * inputVectorY[j * (int)pow(2, k) + l + (int)pow(2, k - 1)]) / sqrt(2);
			}
		}

		inputVectorY = outputVectorX;
	}

	return outputVectorX;
}

vector<z_type> directConvolution(vector<z_type> X, vector<z_type> Y)
{
	size_t N = X.size();
	size_t M = Y.size();

	vector<z_type> U(N + M - 1, 0);

	for (size_t n = 0; n < N + M - 1; n++)
	{
		for (size_t m = 0; m < N + M - 1; m++)
		{
			if (n - m < 0 || n - m >= M || m < 0 || m >= N)
				U[n] += 0;
			else
				U[n] += X[m] * Y[n - m];
		}
	}
		
	return U;
}

vector<z_type> FFTConvolution(vector<z_type> X, vector<z_type> Y)
{
	size_t L = X.size();
	size_t M = Y.size();

	int N = (int)pow(2, ceil(log2(max(L, M))));

	X.resize(2 * N, 0);
	Y.resize(2 * N, 0);

	vector<z_type> SX = FFT_wop(X);
	vector<z_type> SY = FFT_wop(Y);
	vector<z_type> SU(2 * N, 0);
	vector<z_type> U(2 * N, 0);

	for (int k = 0; k < 2 * N; k++)
		SU[k] = sqrt(2 * N) * SX[k] * SY[k];

	U = IFFT_wop(SU);
	U.resize(L + M - 1);

	return U;
}

int main(int argc, char *argv[])
{
	if (argc < 2)
	{
		cout << "Program usage template: ./dsp.exe operation inputfile1 inputfile2(for CONV/FCONV) outputfile" << endl;
		cout << "Operations:" << endl;
		cout << "DFT IDFT FFT IFFT CONV FCONV" << endl;
		return 1;
	}
	else
	{
		vector <z_type> X, Y;
		vector <z_type> A, B, C;

		if (string(argv[1]) == "DFT")
		{
			if (argc != 4)
			{
				cout << "Program usage template: ./dsp.exe DFT inputfile outputfile" << endl;
				exit(1);
			}
			else
			{
				cout << "DFT choosed" << endl;
				X = parseArray(argv[2]);
				Y = DFT(X);
				writeOutputToFile(argv[3], Y);
			}
		}

		if (string(argv[1]) == "IDFT")
		{
			if (argc != 4)
			{
				cout << "Program usage template: ./dsp.exe IDFT inputfile outputfile" << endl;
				exit(1);
			}
			else
			{
				cout << "IDFT choosed" << endl;
				Y = parseArray(argv[2]);
				X = IDFT(Y);
				writeOutputToFile(argv[3], X);
			}
		}

		if (string(argv[1]) == "FFT")
		{
			if (argc != 4)
			{
				cout << "Program usage template: ./dsp.exe FFT inputfile outputfile" << endl;
				exit(1);
			}
			else
			{
				cout << "FFT choosed" << endl;
				X = parseArray(argv[2]);
				Y = FFT_wop(X);
				writeOutputToFile(argv[3], Y);
			}
		}

		if (string(argv[1]) == "IFFT")
		{
			if (argc != 4)
			{
				cout << "Program usage template: ./dsp.exe IFFT inputfile outputfile" << endl;
				exit(1);
			}
			else
			{
				cout << "IFFT choosed" << endl;
				Y = parseArray(argv[2]);
				X = IFFT_wop(Y);
				writeOutputToFile(argv[3], X);
			}
		}

		if (string(argv[1]) == "CONV")
		{
			if (argc != 5)
			{

				cout << "Program usage template: ./dsp.exe CONV inputfile1 inputfile2 outputfile" << endl;
				exit(1);
			}
			else
			{
				cout << "CONV choosed" << endl;
				A = parseArray(argv[2]);
				B = parseArray(argv[3]);
				C = directConvolution(A, B);
				writeOutputToFile(argv[4], C);
			}
		}
		
		if (string(argv[1]) == "FCONV")
		{
			if (argc != 5)
			{
				cout << "Program usage template: ./dsp.exe FCONV inputfile1 inputfile2 outputfile" << endl;
				exit(1);
			}
			else
			{
				cout << "FCONV choosed" << endl;
				A = parseArray(argv[2]);
				B = parseArray(argv[3]);
				C = FFTConvolution(A, B);
				writeOutputToFile(argv[4], C);
			}
		}

		if (string(argv[1]) == "DFT_FFT_TIME_TEST")
		{
			if (argc != 4)
			{
				cout << "Program usage template: ./dsp.exe DFT_FFT_TIME_TEST inputfile outputfile" << endl;
				exit(1);
			}
			else
			{
				// Maksim, code here:

				//
			}
		}

		if (string(argv[1]) == "CONV_FCONV_TIME_TEST")
		{
			if (argc != 5)
			{
				cout << "Program usage template: ./dsp.exe CONV_FCONV_TIME_TEST inputfile1 inputfile2 outputfile" << endl;
				exit(1);
			}
			else
			{
				// Maksim, code here:

				//
			}
		}

		cout << "Check your result in the output file" << endl;
	}

	return 0;
}

