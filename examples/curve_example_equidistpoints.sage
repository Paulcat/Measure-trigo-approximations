# This Sage code computes equispaced points for the upper quadrant on our example curve
#
#     0 = 4*cos(x1)*cos(x2) + 4*cos(x1) + 4*cos(x2) - 1
#
# The starting point of the interval is excluded, the end point is included.
#
# Copyright (c) 2021 Markus Wageringel

prec=80
c1 = 5/8; c2 = 0

from sage.all import ComplexBallField, fast_callable, sqrt
i = SR.I()
t1 = SR.var('t1')
assume(t1 > 0)
T2 = acos(-1 + ZZ(2)*c1 / (cos(t1) + 1 + c2))
T1 = acos(-1 - c2 + ZZ(2)*c1 / (cos(t1) + 1))  # t1, t2 coordinates swapped

KK = ComplexBallField(prec)
gam = vector([exp(i*t1), exp(i * T2)])
gam1 = vector([gk.derivative(t1) for gk in gam])
gam1norm = gam1.norm()

fc = fast_callable(gam1norm, vars=[t1], domain=KK)

t1_a = KK(-acos(-ZZ(1)/2*c2 - 1 + ZZ(1)/2*sqrt(8*c1+c2**2)))
t1_b = KK(acos(-ZZ(1)/2*c2 - 1 + ZZ(1)/2*sqrt(8*c1+c2**2)))

f = lambda t,_: fc(t)

# recursive binary search
# - offset is the area to be accounted for left of point a
def equidist_points(f, a, b, num_points, area=None, offset=KK.zero()):
    assert num_points > 0
    c = (a + b) / 2
    if a.overlaps(b):
        if num_points > 1:
            raise ValueError("cannot distribute %s points between %s and %s" % (num_points, a, b))
        yield c
        return

    f_ac = KK.integral(f, a, c)
    f_cb = KK.integral(f, c, b)
    f_oac = offset + f_ac
    if area is None:
        area = (f_oac + f_cb) / num_points

    num_left = f_oac / area
    if num_left.contains_integer():
        num_left += KK(0.5)
    num_left = floor(num_left)
    num_right = num_points - num_left

    if num_left > 0:
        yield from equidist_points(f, a, c, num_left, area=area, offset=offset)
    if num_right > 0:
        offset_new = f_oac - num_left * area
        yield from equidist_points(f, c, b, num_right, area=area, offset=offset_new)


num_points = 3000 // 4
%time points = list(RDF(p) for p in equidist_points(f, t1_a, t1_b, num_points))
with open(f"points_curve_M{num_points*4}_quarter.txt", "w") as fname:
    num = fname.write(','.join(str(a) for a in points) + ';\n')

num_points = 5000 // 4
%time points = list(RDF(p) for p in equidist_points(f, t1_a, t1_b, num_points))
with open(f"points_curve_M{num_points*4}_quarter.txt", "w") as fname:
    num = fname.write(','.join(str(a) for a in points) + ';\n')
