
def SSNonResXSec( kl, kt, c2, cg, c2g ):
  AAX = [2.09078, 10.1517, 0.282307, 0.101205, 1.33191, -8.51168, -1.37309, 2.82636, 1.45767, -4.91761, -0.675197, 1.86189, 0.321422, -0.836276, -0.568156]
  SM_XSec_BR = 33.45*0.0026

  tot = 0
  tot += AAX[0]*kt*kt*kt*kt
  tot += AAX[1]*c2*c2
  tot += ( AAX[2]*kt*kt + AAX[3]*cg*cg )*kl*kl
  tot += AAX[4]*c2g*c2g
  tot += ( AAX[5]*c2 + AAX[6]*kt*kl )*kt*kt
  tot += ( AAX[7]*kt*kl + AAX[8]*cg*kl )*c2
  tot += AAX[9]*c2*c2g
  tot += ( AAX[10]*cg*kl + AAX[11]*c2g )*kt*kt
  tot += ( AAX[12]*kl*cg + AAX[13]*c2g )*kt*kl
  tot += AAX[14]*cg*c2g*kl

  return tot*SM_XSec_BR

