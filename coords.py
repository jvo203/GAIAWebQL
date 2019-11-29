import numpy as np
import astropy.coordinates as coord
import astropy.units as u

c1 = coord.ICRS(ra=89.014303*u.degree, dec=13.924912*u.degree,
                distance=(37.59*u.mas).to(u.pc, u.parallax()),
                pm_ra_cosdec=372.72*u.mas/u.yr,
                pm_dec=-483.69*u.mas/u.yr,
                radial_velocity=0.37*u.km/u.s)
print(c1)

gc1 = c1.transform_to(coord.Galactocentric)
print(gc1)

gc1.representation_type = 'cylindrical'
print(gc1)

print("V_R:",gc1.d_rho.to(u.km/u.s))
print("V_Phi:",gc1.d_phi)
#print(gc1.d_phi.to(u.km/u.s))
print("V_Z:",gc1.d_z.to(u.km/u.s))
