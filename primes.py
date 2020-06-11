#!/usr/bin/env python3
''' primes.py

This program checks if the input number is a prime number.
If the number is prime, it also checks whether it is a "safe prime."
If the number is a safe prime, it finds the Sophie Germain prime of the input.
Finally, if the number is a safe prime, it find which numbers {1...20} are
full-period generators mod the input.

The default value is because this script was written for a school assignment
that asked to find the above prime properties of that value.

This script depends on the miller_rabin package, which can be installed with
pip3 install miller_rabin
'''

import math
from miller_rabin import miller_rabin

def get_prime_factors(n):
    # first, verify that n is a valid input
    if type(n) is not int:
        raise TypeError("n must be an integer")
    if n < 2:
        raise ValueError("n must be greater than 1")

    # get all the prime factors of n
    factors = []
    # check 2 as a factor first so we can skip all even numbers after
    while (n % 2 == 0):
        factors.append(2)
        n //= 2
    # now check all odd numbers up to sqrt(n)
    i = 3
    # use a while loop to save memory compared to a for loop with range()
    # (which consequently saves a lot of time)
    while i < int(math.sqrt(n)):
        while (n % i == 0):
            factors.append(i)
            n //= i
        i += 2
    factors.append(n)
    return factors

def is_prime(n):
    # first, verify that n is a valid input
    if type(n) is not int:
        raise TypeError("n must be an integer")
    if n < 2:
        raise ValueError("n must be greater than 1")

    return miller_rabin(n)

def is_safe_prime(n):
    # first, verify that n is a valid input
    if type(n) is not int:
        raise TypeError("n must be an integer")
    if n < 2:
        raise ValueError("n must be greater than 1")

    # handle the cases of n = 2, 3 to ensure safe input to is_prime
    if n == 2 or n == 3:
        return False

    q = int((n - 1) / 2)    # int cast will not floor value since n is odd
    return miller_rabin(q)


def main():
    import argparse as ap

    parser = ap.ArgumentParser(description=
            "Find some prime number properties of the input")
    parser.add_argument("input", type=int, metavar="<input>", 
            help="number to find prime properties of")
    parser.add_argument("-f", "--factor", action="store_true", dest="factor",
            default=False, help="get the prime factors of the input")

    args = parser.parse_args()
    n = args.input

    if args.factor:
        factors = get_prime_factors(n)
        if len(factors) > 1:
            print("the prime factors of {} are:".format(n))
            print(factors)
            return  # no need to compute if it's a safe prime
        else:
            print("{} is prime".format(n))

    else:
        # check if n is prime
        primality = is_prime(n)
        if primality is True:
            print("{} is prime".format(n))
        elif primality is False:
            print("{} is not prime".format(n))
            return

    # if n is a prime, check if it is a safe prime
    primality = is_safe_prime(n)
    if primality is True:
        print("{} is a safe prime".format(n))
        print("The Sophie Germain prime is {}".format(int((n - 1) / 2)))
    elif primality is False:
        print("{} is not a safe prime".format(n))
        return
    else:
        print("invalid input: input must be an integer greater than 1")
        return


if __name__ == "__main__":
    main()
