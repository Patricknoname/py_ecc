from __future__ import absolute_import

from py_ecc.fields import (
    bls12_377_FQ as FQ,
    bls12_377_FQ2 as FQ2,
    bls12_377_FQ12 as FQ12,
    bls12_377_FQP as FQP,
)
from py_ecc.fields.field_properties import (
    field_properties,
)

from py_ecc.typing import (
    Field,
    GeneralPoint,
    Point2D,
)


field_modulus = field_properties["bls12_377"]["field_modulus"]
curve_order = 8444461749428370424248824938781546531375899335154063827935233455917409239041

# Curve order should be prime
assert pow(2, curve_order, curve_order) == 2
# Curve order should be a factor of field_modulus**12 - 1
assert (field_modulus ** 12 - 1) % curve_order == 0

# Curve is y**2 = x**3 + 1
b = FQ(1)
# Twisted curve over FQ**2
b2 = FQ2((0, 0x010222f6db0fd6f343bd03737460c589dc7b4f91cd5fd889129207b63c6bf8000dd39e5c1ccccccd1c9ed9999999999a))
# Extension curve over FQ**12; same b value as over FQ
b12 = FQ12((4,) + (0,) * 11)

# Generator for curve over FQ
G1 = (
    FQ(0x008848defe740a67c8fc6225bf87ff5485951e2caa9d41bb188282c8bd37cb5cd5481512ffcd394eeab9b16eb21be9ef),  # noqa: E501
    FQ(0x01914a69c5102eff1f674f5d30afeec4bd7fb348ca3e52d96d182ad44fb82305c2fe3d3634a9591afd82de55559c8ea6),  # noqa: E501
)
# Generator for twisted curve over FQ2
G2 = (
    FQ2([
        0x024aa2b2f08f0a91260805272dc51051c6e47ad4fa403b02b4510b647ae3d1770bac0326a805bbefd48056c8c121bdb8,  # noqa: E501
        0x0ce5d527727d6e118cc9cdc6da2e351aadfd9baa8cbdd3a76d429a695160d12c923ac9cc3baca289e193548608b82801,  # noqa: E501
    ]),
    FQ2([
        0x13e02b6052719f607dacd3a088274f65596bd0d09920b61ab5da61bbdc7f5049334cf11213945d57e5ac7d055d042b7e,  # noqa: E501
        0x0606c4a02ea734cc32acd2b02bc28b99cb3e287e85a763af267492ab572e99ab3f370d275cec1da1aaa9075ff05f79be,  # noqa: E501
    ]),
)
# Point at infinity over FQ
Z1 = None
# Point at infinity for twisted curve over FQ2
Z2 = None


# Check if a point is the point at infinity
def is_inf(pt: GeneralPoint[Field]) -> bool:
    return pt is None


# Check that a point is on the curve defined by y**2 == x**3 + b
def is_on_curve(pt: Point2D[Field], b: Field) -> bool:
    if is_inf(pt):
        return True
    x, y = pt
    return y**2 - x**3 == b


assert is_on_curve(G1, b)
#assert is_on_curve(G2, b2)


# Elliptic curve doubling
def double(pt: Point2D[Field]) -> Point2D[Field]:
    if is_inf(pt):
        return pt
    x, y = pt
    m = 3 * x**2 / (2 * y)
    newx = m**2 - 2 * x
    newy = -m * newx + m * x - y
    return (newx, newy)


# Elliptic curve addition
def add(p1: Point2D[Field], p2: Point2D[Field]) -> Point2D[Field]:
    if p1 is None or p2 is None:
        return p1 if p2 is None else p2
    x1, y1 = p1
    x2, y2 = p2
    if x2 == x1 and y2 == y1:
        return double(p1)
    elif x2 == x1:
        return None
    else:
        m = (y2 - y1) / (x2 - x1)
    newx = m**2 - x1 - x2
    newy = -m * newx + m * x1 - y1
    assert newy == (-m * newx + m * x2 - y2)
    return (newx, newy)


# Elliptic curve point multiplication
def multiply(pt: Point2D[Field], n: int) -> Point2D[Field]:
    if n == 0:
        return None
    elif n == 1:
        return pt
    elif not n % 2:
        return multiply(double(pt), n // 2)
    else:
        return add(multiply(double(pt), int(n // 2)), pt)


def eq(p1: GeneralPoint[Field], p2: GeneralPoint[Field]) -> bool:
    return p1 == p2


# "Twist" a point in E(FQ2) into a point in E(FQ12)
w = FQ12([0, 1] + [0] * 10)


# Convert P => -P
def neg(pt: Point2D[Field]) -> Point2D[Field]:
    if pt is None:
        return None
    x, y = pt
    return (x, -y)


def twist(pt: Point2D[FQP]) -> Point2D[FQ12]:
    if pt is None:
        return None
    _x, _y = pt
    # Field isomorphism from Z[p] / x**2 to Z[p] / x**2 - 2*x + 2
    xcoeffs = [_x.coeffs[0] - _x.coeffs[1], _x.coeffs[1]]
    ycoeffs = [_y.coeffs[0] - _y.coeffs[1], _y.coeffs[1]]
    # Isomorphism into subfield of Z[p] / w**12 - 2 * w**6 + 2,
    # where w**6 = x
    nx = FQ12([xcoeffs[0]] + [0] * 5 + [xcoeffs[1]] + [0] * 5)
    ny = FQ12([ycoeffs[0]] + [0] * 5 + [ycoeffs[1]] + [0] * 5)
    # Divide x coord by w**2 and y coord by w**3
    return (nx / w ** 2, ny / w**3)


G12 = twist(G2)
# Check that the twist creates a point that is on the curve
#assert is_on_curve(G12, b12)
p=double(G1)
assert(is_on_curve(p,b))
p=add(G1,p)
assert(is_on_curve(p,b))
x=multiply(G1,3)
assert(eq(x,p))

xx=multiply(x, curve_order+1)
assert(eq(xx,x))