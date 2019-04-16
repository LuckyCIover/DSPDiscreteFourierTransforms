#define _USE_MATH_DEFINES

#include <algorithm> 
#include <iostream>
#include <fstream>
#include <iomanip>
#include <complex>
#include <string>
#include <vector>
#include <math.h>
#include <limits>
#include <string>
#include <time.h> 
#include <windows.h>

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
		cout << "Program usage template: .\\dsp.exe operation inputfile1 inputfile2(for CONV/FCONV) outputfile" << endl;
		cout << "Operations:" << endl;
		cout << "DFT IDFT FFT IFFT CONV FCONV" << endl;
		return 1;
	}
	else
	{
		vector <z_type> X, Y;
		vector <z_type> A, B, C, C1, C2;

		if (string(argv[1]) == "DFT")
		{
			if (argc != 4)
			{
				cout << "Program usage template: .\\dsp.exe DFT inputfile outputfile" << endl;
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
				cout << "Program usage template: .\\dsp.exe IDFT inputfile outputfile" << endl;
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
				cout << "Program usage template: .\\dsp.exe FFT inputfile outputfile" << endl;
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
				cout << "Program usage template: .\\dsp.exe IFFT inputfile outputfile" << endl;
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

				cout << "Program usage template: .\\dsp.exe CONV inputfile1 inputfile2 outputfile" << endl;
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
				cout << "Program usage template: .\\dsp.exe FCONV inputfile1 inputfile2 outputfile" << endl;
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
				cout << "Program usage template: .\\dsp.exe DFT_FFT_TIME_TEST inputfile outputfile" << endl;
				exit(1);
			}
			else
			{
				vector<z_type> Y_DFT, Y_FFT;
				vector<double> DFT_time, FFT_time;

				cout << "The results of the temporary comparison of DFT and FFT will appear below" << endl;
				cout << "Please, wait..." << endl;

				X = parseArray(argv[2]);

				size_t N = X.size();
				size_t n_max = (size_t)log2(N);
				size_t N_current = N;

				LARGE_INTEGER begin_DFT, end_DFT, begin_FFT, end_FFT, freq1, freq2;

				cout << "---------------------------------------" << endl;

				for (size_t n = 0; n < n_max; n++)
				{
					QueryPerformanceCounter(&begin_DFT);
					Y_DFT = DFT(X);
					QueryPerformanceCounter(&end_DFT);
					QueryPerformanceFrequency(&freq1);

					QueryPerformanceCounter(&begin_FFT);
					Y_FFT = FFT_wop(X);
					QueryPerformanceCounter(&end_FFT);
					QueryPerformanceFrequency(&freq2);

					double norm = 0.0;

					for (size_t k = 0; k < C1.size(); k++)
					{
						if (norm < abs(Y_DFT[k] - Y_FFT[k]))
							norm = abs(Y_DFT[k] - Y_FFT[k]);
					}

					double seconds_DFT = (double)(end_DFT.QuadPart - begin_DFT.QuadPart) / freq1.QuadPart;
					double seconds_FFT = (double)(end_FFT.QuadPart - begin_FFT.QuadPart) / freq2.QuadPart;

					cout.precision(numeric_limits<double>::max_digits10 + 1);
					cout << "N: " << N_current << endl;
					cout << "Err_norm: " << norm << endl;
					cout << "DFT_time: " << seconds_DFT << " seconds" << endl;
					cout << "FFT_time: " << seconds_FFT << " seconds" << endl;
					cout << "---------------------------------------" << endl;

					DFT_time.insert(DFT_time.begin(), seconds_DFT);
					FFT_time.insert(FFT_time.begin(), seconds_FFT);

					N_current = N_current / 2;

					X.resize(N_current);
				}

				ofstream outputFile(argv[3]);

				if (outputFile.is_open())
				{
					outputFile.precision(5);
					outputFile.setf(ios::fixed);
					for (size_t n = 0; n < n_max; n++)
						outputFile << DFT_time[n] << " " << FFT_time[n] << endl;
				}
				else
				{
					cout << "Error opening output file" << endl;
					exit(1);
				}

				outputFile.close();
			}
		}

		if (string(argv[1]) == "CONV_FCONV_TIME_TEST")
		{
			if (argc != 6)
			{
				cout << "Program usage template: .\\dsp.exe CONV_FCONV_TIME_TEST 0/1(both change size) inputfile1 inputfile2 outputfile" << endl;
				exit(1);
			}
			else
			{
				vector<double> CONV_time, FCONV_time;

				cout << "The results of the temporary comparison of DFT and FFT will appear below" << endl;
				cout << "Please, wait..." << endl;

				A = parseArray(argv[3]);
				B = parseArray(argv[4]);

				size_t L = A.size();
				size_t M = B.size();
				size_t L_current = L;
				size_t M_current;

				LARGE_INTEGER begin_СONV, end_CONV, begin_FCONV, end_FCONV, freq1, freq2;

				if (atoi(argv[2]) == 1)
					cout << "Both 1st and 2nd are changing size: " << endl;
				else
					cout << "1st array is fixed, 2nd is changing size: " << endl;

				cout << "---------------------------------------" << endl;

				for (M_current = M; M_current > M / 10; M_current -= M / 10)
				{
					QueryPerformanceCounter(&begin_СONV);
					C1 = directConvolution(A, B);
					QueryPerformanceCounter(&end_CONV);
					QueryPerformanceFrequency(&freq1);

					QueryPerformanceCounter(&begin_FCONV);
					C2 = FFTConvolution(A, B);
					QueryPerformanceCounter(&end_FCONV);
					QueryPerformanceFrequency(&freq2);

					double norm = 0.0;

					for (size_t k = 0; k < C1.size(); k++)
					{
						if (norm < abs(C2[k] - C1[k]))
							norm = abs(C2[k] - C1[k]);
					}

					double seconds_CONV = (double)(end_CONV.QuadPart - begin_СONV.QuadPart) / freq1.QuadPart;
					double seconds_FCONV = (double)(end_FCONV.QuadPart - begin_FCONV.QuadPart) / freq2.QuadPart;

					cout.precision(numeric_limits<double>::max_digits10 + 1);
					cout << "L: " << L_current << endl;
					cout << "M: " << M_current << endl;
					cout << "Err_norm: " << norm << endl;
					cout << "CONV_time: "  << seconds_CONV << " seconds" << endl;
					cout << "FCONV_time: " << seconds_FCONV << " seconds" << endl;
					cout << "---------------------------------------" << endl;

					CONV_time.insert(CONV_time.begin(), seconds_CONV);
					FCONV_time.insert(FCONV_time.begin(), seconds_FCONV);

					if (atoi(argv[2]) == 1)
					{
						L_current = L_current - L / 10;
						A.resize(L_current);
					}

					B.resize(M_current);
				}

				ofstream outputFile(argv[5]);

				if (outputFile.is_open())
				{
					outputFile.precision(5);
					outputFile.setf(ios::fixed);
				
                    for (size_t m = 0; m < CONV_time.size(); m++)
                    {
                        if (atoi(argv[2]) == 1)
                            outputFile << L_current + (m + 1) * (L / 10) << " " << M_current + (m + 1) * (M / 10) << " " << CONV_time[m] << " " << FCONV_time[m] << endl;
                        else
                            outputFile << L_current << " " << M_current + (m + 1) * (M / 10) << " " << CONV_time[m] << " " << FCONV_time[m] << endl;
                    }
				}
				else
				{
					cout << "Error opening output file" << endl;
					exit(1);
				}

				outputFile.close();
			}
		}

		cout << "Check your result in the output file" << endl;
	}

	system("pause");
	return 0;
}
