#ifndef nfsglft_h
#define nfsglft_h

void ndsglft (int B, int M, double** x, double* f_real, double* f_imag, double*** f_hat_real, double*** f_hat_imag);
void nfsglft (int B, int M, double** x, double* f_real, double* f_imag, double*** f_hat_real, double*** f_hat_imag, int nfft_sigma, int nfft_q);

void ndsglft_adjoint (int B, int M, double** x, double* f_real, double* f_imag, double*** f_hat_real, double*** f_hat_imag);
void nfsglft_adjoint (int B, int M, double** x, double* f_real, double* f_imag, double*** f_hat_real, double*** f_hat_imag, int nfft_sigma, int nfft_q);

void infsglft_golub_van_loan (int B, int M, double** x, double* f_real, double* f_imag, double*** f_hat_real, double*** f_hat_imag, double*** f_hat_real_ground_truth, double*** f_hat_imag_ground_truth, int n_iter, int nfft_sigma, int nfft_q, char* preferences);
void infsglft_bjoerck (int B, int M, double** x, double* f_real, double* f_imag, double*** f_hat_real, double*** f_hat_imag, double*** f_hat_real_ground_truth, double*** f_hat_imag_ground_truth, int n_iter, int nfft_sigma, int nfft_q, char* preferences);

long double R (int n, int l, long double r);

double M_P (int l, int m, double xi);
void Y (int l, int m, double theta, double phi, double* result);

#endif