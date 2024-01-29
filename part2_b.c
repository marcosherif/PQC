#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MODULUS 97
#define ROOT_W 22
#define INV_ROOT_W 33


int mod_pow(int base, int exp, int mod) {
    int result = 1;
    base = base % mod;
    while (exp > 0) {
        if (exp % 2 == 1) {
            result = (result * base) % mod;
        }
        exp = exp >> 1;
        base = (base * base) % mod;
    }
    return result;
}


void ntt(int *a, int n, int root) {
    int m, i, j, k, w, omega;
    for (m = n; m > 1; m /= 2) {
        omega = mod_pow(root, (MODULUS - 1) / m, MODULUS);
        w = 1;
        for (k = 0; k < n; k += m) {
            for (j = k; j < k + m / 2; j++) {
                i = j + m / 2;
                int u = a[j];
                int v = (long long)w * a[i] % MODULUS;
                a[j] = (u + v) % MODULUS;
                a[i] = (u - v + MODULUS) % MODULUS;
            }
            w = (long long)w * omega % MODULUS;
        }
    }
}


void intt(int *a, int n, int inv_root) {
    int m, i, j, k, w, omega_inv;
    for (m = 2; m <= n; m *= 2) {
        omega_inv = mod_pow(inv_root, (MODULUS - 1) / m, MODULUS);
        w = 1;
        for (k = 0; k < n; k += m) {
            for (j = k; j < k + m / 2; j++) {
                i = j + m / 2;
                int u = a[j];
                int v = (long long)w * a[i] % MODULUS;
                a[j] = (u + v) % MODULUS;
                a[i] = (u - v + MODULUS) % MODULUS;
                a[j] = (long long)a[j] * inv_root % MODULUS;
                a[i] = (long long)a[i] * inv_root % MODULUS;
            }
            w = (long long)w * omega_inv % MODULUS;
        }
    }
}


void polynomial_multiply(int *a, int *b, int *result, int n) {
    int i;

    
    ntt(a, n, ROOT_W);
    ntt(b, n, ROOT_W);

    
    for (i = 0; i < n; i++) {
        result[i] = (long long)a[i] * b[i] % MODULUS;
    }

    
    intt(result, n, INV_ROOT_W);
    
    
    int inv_n = mod_pow(n, MODULUS - 2, MODULUS);
    for (i = 0; i < n; i++) {
        result[i] = (long long)result[i] * inv_n % MODULUS;
    }
}


int main() {
    int a[] = {1, 2, 3, 4};
    int b[] = {5, 6, 7, 8};
    int n = 4;
    int result[4];

    polynomial_multiply(a, b, result, n);

    printf("Result of polynomial multiplication: ");
    for (int i = 0; i < n; i++) {
        printf("%d ", result[i]);
    }
    printf("\n");

    return 0;
}
