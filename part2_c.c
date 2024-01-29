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


void permute_bitreverse(int *src, int *dst, int n) {
    int i, j;
    for (i = 0; i < n; i++) {
        dst[i] = 0;
        for (j = 0; j < n; j++) {
            if ((i & (1 << j)) != 0) {
                dst[i] |= (1 << (n - 1 - j));
            }
        }
    }
}


void ntt_fast(int *data, int n, int forward) {
    int *x = (int *)malloc(n * sizeof(int));
    permute_bitreverse(data, x, n);

    int *factors;
    int w;
    if (forward) {
        factors = (int *)malloc(n * sizeof(int));
        w = ROOT_W;
        for (int i = 0; i < n; i++) {
            factors[i] = mod_pow(w, i, MODULUS);
        }
    } else {
        factors = (int *)malloc(n * sizeof(int));
        w = INV_ROOT_W;
        for (int i = 0; i < n; i++) {
            factors[i] = mod_pow(w, i * (n - 1), MODULUS);
        }
    }

    
    for (int s = 1; s <= log2(n); s++) {
        int m = 1 << s;
        for (int b = 0; b < n; b += m) {
            int factor = 1;
            for (int op = 0; op < m / 2; op++) {
                int a = x[b + op];
                int b_val = (long long)x[b + op + m / 2] * factor % MODULUS;
                x[b + op] = (a + b_val) % MODULUS;
                x[b + op + m / 2] = (a - b_val + MODULUS) % MODULUS;
                factor = (long long)factor * factors[op] % MODULUS;
            }
        }
    }

    
    if (forward) {
        for (int i = 0; i < n; i++) {
            x[i] = (long long)x[i] * factors[i] % MODULUS;
        }
    }

    
    for (int i = 0; i < n; i++) {
        data[i] = x[i];
    }

    free(x);
    free(factors);
}


int main() {
    int data[] = {1, 2, 3, 4};
    int n = 4;

    printf("Original data: ");
    for (int i = 0; i < n; i++) {
        printf("%d ", data[i]);
    }
    printf("\n");

    ntt_fast(data, n, 1);  

    printf("Transformed data: ");
    for (int i = 0; i < n; i++) {
        printf("%d ", data[i]);
    }
    printf("\n");


    return 0;
}
