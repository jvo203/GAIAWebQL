/*#include <pybind11/pybind11.h>
namespace py = pybind11;*/

#include <stdio.h>
#include <math.h>

#include "PJMCoords.h"

int main()
{
    double alpha = 1.553592668717866;
    double delta = 0.24303556245046865;

    double ra = 89.014303;
    double dec = 13.924912;
    double d = 26.60281989890928;

    double pmra = 372.72;          //[mas/yr]
    double pmdec = -483.69;        //[mas/yr]
    double radial_velocity = 0.37; //[km/s]

    double rICRS[3] = {d * cos(alpha) * cos(delta), d * sin(alpha) * cos(delta), d * sin(delta)};

    //py::scoped_interpreter guard{};
    double R[3][3] = {{-0.0548755604162154, -0.8734370902348850, -0.4838350155487132}, {+0.4941094278755837, -0.4448296299600112, +0.7469822444972189}, {-0.8676661490190047, -0.1980763734312015, +0.4559837761750669}};
    double theta = 0.003253017785496385;
    double H[3][3] = {{cos(theta), 0.0, sin(theta)}, {0.0, 1.0, 0.0}, {-sin(theta), 0.0, cos(theta)}};
    double dGC = 8300.0;

    double rGC[3] = {0.0, 0.0, 0.0};
    double r[3] = {-dGC, 0.0, 0.0};
    //r = R * r_icrs - dGC * xGC

    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            r[i] += R[i][j] * rICRS[j];

    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            rGC[i] += H[i][j] * r[j];

    for (int i = 0; i < 3; i++)
        printf("rGC[%d] = %f\t", i, rGC[i]);
    printf("\n");

    //the final solution
    OmniCoords OC;
    OC.change_sol_pos(8.3, 0.027);

    /*
        position = (r,a,d), velocity = (vr,ma*,md)
        units: as defined in Units.h (kpc,radian,kpc/Myr,radian/Myr)
        NOTE:  ma* = cos(d) * da/dt = cos(d) * ma
        The equatorial coordinates (distance, right ascension, declination)
    */

    //"\tHEQ\t\tHeliocentric EQuatorial polar (r,a,d,vr,mu_a*,mu_d)\n\n"
    /*coords[1] *=Units::degree_i;
  coords[2] *=Units::degree_i;
  coords[3] *=Units::kms_i;
  coords[4] *=Units::masyr_i;
coords[5] *=Units::masyr_i;*/
    /*coords[1] *=Units::degree_i;
  coords[2] *=Units::degree_i;
  coords[3] *=Units::kms_i;
  coords[4] *=Units::masyr_i;
coords[5] *=Units::masyr_i;*/

    vec6 sHEQ, sGCA, sGCY;

    sHEQ[0] = d / 1000.0;      //[kpc]
    sHEQ[1] = ra;              //[deg]
    sHEQ[2] = dec;             //[deg]
    sHEQ[3] = radial_velocity; //[km/s]
    sHEQ[4] = pmra;
    sHEQ[5] = pmdec;
    //vec6 GCAfromHEQ(vec6 sHEQ) {take_HEQ(sHEQ); return give_GCA();}

    OC.take_HEQ_units(sHEQ);
    sGCA = OC.give_GCA_units();

    double X = sGCA[0];
    double Y = sGCA[1];
    double RR = sqrt(X * X + Y * Y);
    std::cout << "X: " << X << "\tY: " << Y << "\tR: " << RR << std::endl;

    std::cout << "to GCA:" << std::endl;
    for (int i = 0; i < 6; i++)
        std::cout << sGCA[i] << std::endl;

    OC.take_HEQ_units(sHEQ);
    sGCY = OC.give_GCY_units();

    std::cout << "to GCY:" << std::endl;
    for (int i = 0; i < 6; i++)
        std::cout << sGCY[i] << std::endl;
}