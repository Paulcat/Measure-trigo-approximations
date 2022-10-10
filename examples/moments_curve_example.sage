# This Sage code computes moments up to high precision for our example curve
#
#     0 = 4*cos(x1)*cos(x2) + 4*cos(x1) + 4*cos(x2) - 1
#
# Dependencies: https://github.com/mwageringel/momentproblems
#
# Copyright (c) 2021 Markus Wageringel

from momentproblems.moment_functionals import mk_moment_trigonometric_curve_example
mom = mk_moment_trigonometric_curve_example(5/8, 0, prec=80)

d = 100
alphas = [(a,b) for a in range(d+1) for b in range(d+1)]

for k in range(d+1):
    # chunks for progress
    alpha_chunk = alphas[(d+1)*k:(d+1)*(k+1)]
    print("precomputing chunk k=%s : %s..%s" % (k, alpha_chunk[0], alpha_chunk[-1]))
    mom.precompute([(α,) for α in alpha_chunk], num_processes=20)

err = max(mom(α).rad() for α in alphas)
print('maximum error radius: %s' % err)

# note that all moments are real
M = matrix(RDF, d+1, d+1, {α:mom(α).real() for α in alphas}, sparse=False)

print("saving sobj file")
save((err, M), f"moments_curve_n{d}.sobj")

print("saving txt file")
with open(f"moments_curve_n{d}.txt", "w") as f:
    for row in M.rows():
        num = f.write(','.join(str(a) for a in row) + ';\n')
