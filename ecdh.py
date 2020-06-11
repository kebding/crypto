''' ecdh.py
This file contains a library for elliptic-curve Diffie-Hellman key exchange.
'''

import elliptic_curve as ec
from math import ceil, sqrt

class ECDH():
    def __init__(self, curve, g):
        ''' ECDH.__init__
        initializes the elliptic-curve Diffie-Hellman instance
        inputs:
            elliptic_curve.EC curve
            elliptic_curve.Point g
        returns:
            none
        '''
        if not isinstance(curve, ec.EC):
            raise TypeError("curve must be an elliptic curve object")
        if not isinstance(g, ec.Point):
            raise TypeError("g must be a Point object")
        self.curve = curve
        self.g = g


    def generate_public_key(self, private_key):
        ''' ECDH.generate_public_key
        this function generates an ECDH public key from an input private key
        inputs:
            int private_key
        returns:
            elliptic_curve.Point public_key
        '''
        if not (type(private_key) is int):
            raise TypeError("private key must be an integer")
        if not (private_key > 0):
            raise ValueError("invalid private key value")

        public_key = self.curve.mult(self.g, private_key)
        return public_key


    def calculate_shared_secret(self, private_key, public_key):
        ''' ECDH.calculate_shared_secret
        this function calculates the shared secret for an ECDH key exchange
        inputs:
            int private_key
            elliptic_curve.Point public_key
        returns:
            elliptic_curve.Point shared_secret
        '''
        if not (type(private_key) is int):
            raise TypeError("private key must be an integer")
        if not (private_key > 0):
            raise ValueError("invalid private key value")
        if not isinstance(public_key, ec.Point):
            raise TypeError("public key must be a Point object")
        if not self.curve.has_point(public_key):
            raise ValueError("public key not on curve")

        shared_secret = self.curve.mult(public_key, private_key)
        return shared_secret
