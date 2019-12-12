"""
Test if the ellipcorr lib has been installed correctly?

See https://github.com/GeoscienceAustralia/ellip-corr/

print(ellipcorr.ellipticity_corr('sP', 65, 124, 45, 39), "=? -0.38976147770881653")
print(ellipcorr.ellipticity_corr('P', 65, 124, 45, 39), "=?  -0.3774765431880951")
print(ellipcorr.ellipticity_corr('pPKiKP', 65, 124, 45, 39),"=? -0.7800000309944153" )
"""

import ellipcorr
x = ellipcorr.ellipticity_corr('sP', 65, 124, 45, 39)
assert (x== -0.38976147770881653)

x = ellipcorr.ellipticity_corr('P', 65, 124, 45, 39)
assert (x== -0.3774765431880951)


print("ellipcorr testing completed")
