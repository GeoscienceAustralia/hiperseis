"""
Test if the ellipcorr lib has been installed correctly?

1. Install ellip-corr into a python environment. See guide https://github.com/GeoscienceAustralia/ellip-corr/

2. export  ELLIPCORR=~/Githubz/hiperseis/ellip-corr/

3. run this script in the envionment where ellipcorr has been installed: 
    source ~/Venvs/hiperseis/bin/activate

"""

import ellipcorr

x = ellipcorr.ellipticity_corr('sP', 65, 124, 45, 39)
assert (x== -0.38976147770881653)

x = ellipcorr.ellipticity_corr('P', 65, 124, 45, 39)
assert (x== -0.3774765431880951)

print(ellipcorr.ellipticity_corr('sP', 65, 124, 45, 39), "=? -0.38976147770881653")
print(ellipcorr.ellipticity_corr('P', 65, 124, 45, 39), "=?  -0.3774765431880951")
print(ellipcorr.ellipticity_corr('pPKiKP', 65, 124, 45, 39),"=? -0.7800000309944153" )

print("ellipcorr testing completed")
