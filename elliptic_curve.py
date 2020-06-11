''' elliptic_curve.py

This program defines a library for elliptic curve operations.
'''

from math import sqrt
import bigfloat


def extended_gcd(a, b):
    ''' extended_gcd
    This function returns the greatest common divisor of its inputs using the
    extended Euclidean algorithm. It returns the gcd and the Bezout coefficients
    (i.e. the x and y in the equation ax + by == gcd(a, b) ) for use in
    multiplicative inverse calculations.

    inputs:
        int a
        int b
    returns:
        (gcd, x, y)
    '''
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


def modinv(a, b):
    ''' modinv
    This function returns the multiplicative inverse of a mod b.
    Since ax + by == gcd(a, b) == 1, ax - 1 == (-y)b, which means that
    ax is congruent to 1 with regards to mod b. therefore, the
    multiplicative inverse of a mod b is x.

    inputs:
        int a
        int b
    returns:
        x = multiplicative inverse of a mod b
    '''
    (gcd, x, y) = extended_gcd(a, b)
    if gcd != 1:
        return None
    return x % b  # mod is to handle negative results



# A class for an integer point on a Cartesian plane.
class Point():
    def __init__(self, x = 0, y = 0):
        if not (type(x) is int and type(y) is int):
            raise TypeError("x-value must be integers")
        self.x = x
        self.y = y

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            if self.x == other.x and self.y == other.y:
                return True
            else:
                return False
        return False

    def __str__(self):
        return "Point({},{})".format(self.x, self.y)




# A class for an elliptic curve.
class EC():
    def __init__(self, a, b, mod):
        ''' EC.__init__
        initializes the elliptic curve, requiring that a and b are less than the
        modulus and that a and b are positive numbers
        inputs:
            int a
            int b
            int mod
        returns:
            none
        '''
        a > 0 and b > 0 and a < mod and b < mod and mod > 2
        assert ((4 * (a**3)) + (27 * (b**2)) != 0)
        self.a = a
        self.b = b
        self.mod = mod
        self.identity_element = Point(0, 0)

    def __str__(self):
        return "Curve(a={}, b={}, mod={}".format(self.a, self.b, self.mod)

    def add(self, p1, p2):
        ''' EC.add
        this function adds two points on an elliptic curve
        conceptually, it creates a line from p1 to y and returns the negation of the
        point hit by continuing the line beyond y.
        if p1 and p2 are the same point, a tangent of the curve at p1 is used as the
        line for the computation.

        inputs:
            Point p1
            Point p2
        returns:
            Point (p1 + p2)
        '''
        # sum of two identity elements is zero
        if p1 == p2 and p1.y == 0:
            return self.identity_element
        # sum of two points on a vertical line is zero
        if p1.x == p2.x and p1.y != p2.y:
            return self.identity_element
        # sum of a point and the identity is the point
        if p1 == self.identity_element:
            return p2
        if p2 == self.identity_element:
            return p1

        # else, we have two non-zero points to add
        if p1 != p2:
            s = ((p2.y - p1.y) * modinv(p2.x - p1.x, self.mod)) % self.mod
        else:
            s = ((3 * (p1.x**2) + self.a) * modinv(2 * p1.y, self.mod)) \
                    % self.mod

        ret_x = (s**2 - p1.x - p2.x) % self.mod
        ret_y = (s * (p1.x - ret_x) - p1.y) % self.mod

        return Point(ret_x, ret_y)


    def mult(self, p, n):
        ''' EC.mult
        this function computes p * n for a point on an elliptic curve. instead of
        the naive iterative addition method, this function uses the more efficient
        double-and-add method.

        inputs:
            Point p
            int n
        returns:
            Point (p * n)
        '''
        if n == 0:
            return self.identity_element
        # get the bitwise representation of n as an iterable string
        n_bits = bin(n)[2:]
        #print("n={}, n_bits={}".format(n, n_bits))

        pn = Point(p.x, p.y)
        #print("starting point = {}".format(pn))
        for i in n_bits[1:]:
            pn = self.add(pn, pn)  # double
            #print("doubled: new point = {}".format(pn))
            if i == '1':
                pn = self.add(pn, p)  # add
                #print("added: new point = {}".format(pn))
            else:
                #print("no add")
                pass

        return pn


    def points_at(self, x):
        ''' EC.points_at
        this function returns the points on the curve at the given x-value.

        inputs:
            int x
        returns:
            (Point(x, y(x)), Point(x, -y(x)))
        '''
        if x >= self.mod:
            raise ValueError("x-value must be less than curve's modulus")

        # y^2 = x^3 + ax + b mod n
        y_squared = ((x**3) + (self.a * x) + self.b) % self.mod
        y = sqrt(y_squared)
        return (Point(x, y), Point(x, -y))


    def has_point(self, p):
        ''' EC.has_point
        this function returns True if the input point is on the curve and False
        otherwise by checking if plugging the point's x-value into the curve
        equation results in the same y value as the point (mod n).

        inputs:
            Point p
        returns:
            True if p is on the curve, False if not
        '''
        if not isinstance(p, Point):
            raise TypeError("input must be a Point object")

        if p.x >= self.mod:
            return False

        # compute y^2, the left-hand-side of the curve equation, for the
        # input point's x-value
        curve_y_sq = ((p.x**3) + (self.a * p.x) + self.b) % self.mod
        # compute the point's y^2 mod n
        p_y_sq = (p.y**2) % self.mod

        return curve_y_sq == p_y_sq


    def hasses_bound(self):
        ''' EC.hasses_bound
        this function returns a tuple indicating the interval of Hasse's Theorem
        for the curve with a modulus of p:
            p + 1 - 2*sqrt(p) <= num points on curve <= p + 1 + 2*sqrt(p)
        in order to support finding the bounds for large curves, this function
        utilizes the bigfloat libary and returns a bigfloat object.
        inputs:
            none
        returns:
            tuple (lowerbound, upperbound)
        '''
        mod_sqrt = bigfloat.sqrt(self.mod)
        lower = self.mod + 1 - (2 * mod_sqrt)
        upper = self.mod + 1 + (2 * mod_sqrt)
        return (lower, upper)


    def negative(p):
        ''' EC.negative
        returns the negative of the input point
        inputs:
            Point p
        returns:
            Point (p.x, -p.y)
        '''
        if not self.has_point(p):
            raise ValueError("point must be on the curve")
        return Point(p.x, ((-p.y) % self.mod))
