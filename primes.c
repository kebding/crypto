/* primes.c
This program checks if the input number is a prime number.
If the number is prime, it also checks whether it is a "safe prime."
If the number is a safe prime, it finds the Sophie Germain prime of the input.
Finally, if the number is a safe prime, it find which numbers [1, 20] are
full-period generators mod the input.

The default value is because this script was written for a school assignment
that asked to find the above prime properties of that value.

Since this program is designed to handle very large numbers, it uses the GNU MP
(GMP) library. To compile, make sure to have GMP installed and to include the
flag -lgmp when running gcc.
*/

#include <stdio.h>
#include <gmp.h>

#define MILLER_RABIN_REPS 32 // GMP recommends 15-50

int main(int argc, char** argv) {

    if (argc < 2) {
        printf("USAGE: %s <number>\n", argv[0]);
        return -1;
    }

    // set up the large numbers. mpz_t is GMP's integer (z) type
    mpz_t n, q, temp1, temp2;
    mpz_init(n); // input number
    mpz_init(q); // input's Sophie Germain prime (if input is a safe prime)
    mpz_init(temp1);
    mpz_init(temp2);

    // get the number to calculate from argv[1]
    if (mpz_set_str(n, argv[1], 0) < 0) {
        printf("error reading input\n");
        return -1;
    }

    // determine if the input is prime
    if (!mpz_probab_prime_p(n, MILLER_RABIN_REPS)) {
        printf("input is not prime\n");
        mpz_clears(n, q, temp1, temp2, NULL);
        return 0;
    }
    else {
        printf("input is prime\n");
    }

    // make sure the input isn't 2. if so, it is an unsafe prime so we're done.
    // if it isn't 2, we know the input is odd, which helps below.
    if (mpz_cmp_ui(n, 2) == 0) {
        printf("input is not a safe prime\n");
        mpz_clears(n, q, temp1, temp2, NULL);
        return 0;
    }

    // determine if the input is a safe prime by checking if q=(n-1)/2 is prime
    mpz_sub_ui(temp1, n, 1);
    mpz_divexact_ui(q, temp1, 2);
    if (!mpz_probab_prime_p(q, MILLER_RABIN_REPS)) {
        printf("input is not a safe prime\n");
        mpz_clears(n, q, temp1, temp2, NULL);
        return 0;
    }
    else {
        printf("input is a safe prime. its Sophie Germain prime is:\n");
        mpz_out_str(stdout, 10, q);
        printf("\n");
    }

    // determine which elements of [1, 20] are full-period generators mod n
    // Lagrange's Theorem tells us the order of any element must be {1, 2, q, n-1}
    int i;
    for (i = 1; i <= 20; i++) {
        mpz_set_ui(temp1, i);

        // temp2 = (temp1 ^ X) mod n
        mpz_powm_ui(temp2, temp1, 1, n);
        if (mpz_cmp_ui(temp2, 1) == 0)
            continue;
        mpz_powm_ui(temp2, temp1, 2, n);
        if (mpz_cmp_ui(temp2, 1) == 0)
            continue;
        mpz_powm(temp2, temp1, q, n);
        if (mpz_cmp_ui(temp2, 1) == 0)
            continue;

        // if we made it here, none of {1, 2, q} are the order of i,
        // so the order of i must be n-1, meaning i is a full-period generator
        printf("%d is a full-period generator\n", i);
    }

    mpz_clears(n, q, temp1, temp2, NULL);
    return 0;
}
