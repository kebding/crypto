#!/usr/bin/env python3
''' rsa.py

This script takes an input file and encrypts or decrypts it using the RSA
encryption algorithm using Euler's totient.
The input file, key, output file, and mode of operation are supplied as
command-line arguments. The default mode of operation is encryption.
'''

import math
import random

# import primality test function. if not found, import miller-rabin and
# build the primality test. if miller-rabin cannot be imported, abort.
try:
    from primes import is_prime
except ImportError:
    try:
        from miller_rabin import miller_rabin
    except ImportError:
        import sys
        sys.exit("error: could not import primes or miller-rabin")
    def is_prime(n):
        if type(n) is not int:
            raise TypeError("n must be an integer")
        if n < 2:
            raise ValueError("n must be greater than 1")
        return miller_rabin(n)

# This function returns the greatest common divisor of its inputs.
def gcd(a, b):
    while b != 0:
        (a, b) = (b, a % b)
    return a


# This function returns the greatest common divisor of its inputs using the
# extended Euclidean algorithm. It returns the gcd and the Bezout coefficients
# (i.e. the x and y in the equation ax + by == gcd(a, b) ) for use in 
# multiplicative inverse calculations.
# returns: (gcd, x, y)
def extended_gcd(a, b):
    if type(a) is not int or type(b) is not int:
        raise TypeError("inputs must be two integers")
    r = b
    r_prev = a
    x = 0
    x_prev = 1
    y = 1
    y_prev = 0

    while r != 0:
        quotient = r_prev // r
        (r_prev, r) = (r, r_prev - (quotient * r))
        (x_prev, x) = (x, x_prev - (quotient * x))
        (y_prev, y) = (y, y_prev - (quotient * y))

    return (r_prev, x_prev, y_prev)


# This function returns the multiplicative inverse of a mod b.
# Since ax + by == gcd(a, b) == 1, ax - 1 == (-y)b, which means that
# ax is congruent to 1 with regards to mod b. therefore, the 
# multiplicative inverse of a mod b is x.
def multiplicative_inverse(a, b):
    (gcd, x, y) = extended_gcd(a, b)
    if gcd != 1:  # error
        return None
    return x


# This function takes an input bytes object and key and returns
# a bytes object for the corresponding ciphertext. Endianness is assumed
# to be little.
def rsa_encrypt(plaintext, key):
    # since the key contains a key and a modulus, unpack it
    exponent, modulus = key

    # since the message must be shorter than the modulus, break the message
    # into modulus-sized blocks, encrypt/decrypt the blocks, and string the
    # output blocks together as the output text.

    # read 1 less byte than the modulus' length in case the plaintext's int
    # value is greater than the modulus
    block_len = math.ceil((modulus.bit_length()) / 8) - 1
    block_num = 0
    ciphertext = bytes()
    while(block_len * block_num < len(plaintext)):
        block = plaintext[(block_len * block_num):(block_len * (block_num + 1))]
        block_num += 1
        plaintext_int = int.from_bytes(block, 'little')
        ciphertext_int = pow(plaintext_int, exponent, modulus)
        # write 1 more byte than block_len since the result can be up to the
        # modulus' value, which cannot be represented in block_len bytes (see
        # block_len calculation above)
        ciphertext_bytes = ciphertext_int.to_bytes(block_len + 1, 'little')
        ciphertext += ciphertext_bytes

    return ciphertext


# This function takes an input bytes object and key and returns
# a bytes object for the corresponding plaintext. Endianness is assumed
# to be little. 
def rsa_decrypt(ciphertext, key):
    # since the key contains a key and a modulus, unpack it
    exponent, modulus = key

    # since the message must be shorter than the modulus, break the message
    # into modulus-sized blocks, encrypt/decrypt the blocks, and string the
    # output blocks together as the output text.

    # unlike in rsa_encrypt, we know each block's value is less than the value
    # of the modulus and that the ciphertext bytes were written in
    # modulus-length blocks, so we can read a full modulus-length bytes for the
    # block
    block_len = math.ceil((modulus.bit_length()) / 8)
    block_num = 0
    plaintext = bytes()
    while(block_len * (block_num + 1) < len(ciphertext)):
        block = ciphertext[(block_len * block_num):(block_len * (block_num + 1))]
        block_num += 1
        ciphertext_int = int.from_bytes(block, 'little')
        plaintext_int = pow(ciphertext_int, exponent, modulus)
        plaintext_bytes = plaintext_int.to_bytes(block_len - 1, 'little')
        plaintext += plaintext_bytes
    # handle final block separately to strip zero-padding from encryption
    block = ciphertext[(block_len * block_num):(block_len * (block_num + 1))]
    ciphertext_int = int.from_bytes(block, 'little')
    plaintext_int = pow(ciphertext_int, exponent, modulus)
    plaintext_len = math.ceil(plaintext_int.bit_length() / 8)
    plaintext_bytes = plaintext_int.to_bytes(plaintext_len, 'little')
    plaintext += plaintext_bytes

    return plaintext


# This function generates a pair of keys from two input prime numbers. Each key
# is a tuple of (exponent, modulus). The keys are returned as a tuple of the
# form (public, private).
def generate_keys(p, q):
    if not is_prime(p) or not is_prime(q) or p == q:
        raise ValueError("p and q must be distinct prime numbers")
    n = p * q
    phi_n = (p - 1) * (q - 1)
    # use public exponent of 2^(16) + 1 since it is the most common value. if
    # it is not coprime to phi_n, find a new number
    e = pow(2, 16) + 1
    g = gcd(e, phi_n)
    while g != 1:
        # try a new odd number 3 <= n < phi_n
        e = random.randrange(3, phi_n, 2)
        g = gcd(e, phi_n)

    # get the multiplicative inverse of e mod phi
    d = multiplicative_inverse(e, phi_n)

    return ((e, n), (d, n))


if __name__ == "__main__":
    import argparse as ap

    parser = ap.ArgumentParser(description=
            "RSA encryption/decryption of files")
    parser.add_argument("input", type=str,
            help="file to encrypt/decrypt")
    parser.add_argument("output", type=str,
            help="file to write the output to")
    parser.add_argument("p", type=int,
            help="first prime number to generate a key with")
    parser.add_argument("q", type=int,
            help="second prime number to generate a key with")
    parser.add_argument("-d", "--decrypt", dest="decrypt", action="store_true",
            default=False, help="decrypt the input instead of encrypting")

    args = parser.parse_args()

    public_key, private_key = generate_keys(args.p, args.q)

    infile = open(args.input, 'rb')
    input_data = infile.read()
    infile.close()

    if args.decrypt:
        output_data = rsa_decrypt(input_data, private_key)
    else:
        output_data = rsa_encrypt(input_data, public_key)

    outfile = open(args.output, 'wb')
    outfile.write(output_data)
    outfile.close()
