#include <stdio.h>
void updateEq(double* uNext, double* u, double* uPrev, double* coeffs, int Nx, int Ny)
{
 for (int l = 2; l < Nx - 2; ++l) { for (int m = 2; m < Ny - 2; ++m) { uNext[l + m * Nx] = coeffs[0] * u[l + m * Nx] + coeffs[1] * (u[l + (m-2) * Nx] + u[l + (m+2) * Nx] + u[l+2 + m * Nx] + u[l-2 + m * Nx]) + coeffs[2] * (u[l-1 + (m-1) * Nx] + u[l+1 + (m+1) * Nx] + u[l+1 + (m-1) * Nx] + u[l-1 + (m+1) * Nx]) + coeffs[3] * (u[l + (m-1) * Nx] + u[l + (m+1) * Nx] + u[l+1 + m * Nx] + u[l-1 + m * Nx]) + coeffs[4] * uPrev[l + m * Nx] + coeffs[5] * (uPrev[l + (m-1) * Nx] + uPrev[l + (m+1) * Nx] + uPrev[l+1 + m * Nx] + uPrev[l-1 + m * Nx]);} }
}