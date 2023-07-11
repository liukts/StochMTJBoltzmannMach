# ===========================================================
#   Various resistance values for 
#   memristive devices from the literature (references included)
#   which have been used to implement cross bar arrays.
# ===========================================================
from dataclasses import dataclass

@dataclass(slots=True)
class RRAM_dev:
    LRS   : float
    HRS   : float
    # amplification for re-input to crossbar array
    scale : float
    ref   : str
    def __str__(self):
        return f"{self.ref}"

#   scale is gleaned from simulation data
#   non 1 scale values are working for MaxCut3
HfO2      = RRAM_dev(10e1    ,7*10e6 ,{"Max Sat": 1e11,
                                       "Max Cut": 1e12},"Lee,2008")
ZnO       = RRAM_dev(2*10e-1 ,2*10e6 ,1    ,"Kim 2009")
MnO2      = RRAM_dev(2*10e-1 ,10e4   ,1    ,"Yang 2009")
MnO_Ta2O5 = RRAM_dev(5*10e3  ,10e9   ,1    ,"Hu Q 2018")
test      = RRAM_dev(1000    ,3000   ,1e13 ,"N/A")
test2     = RRAM_dev(1000    ,7000   ,3.5e7 ,"N/A")
