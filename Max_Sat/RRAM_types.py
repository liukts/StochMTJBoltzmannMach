# ===========================================================
#   Various resistance values for 
#   memristive devices from the literature (references included)
#   which have been used to implement cross bar arrays.
# ===========================================================
from dataclasses import dataclass

#   add slots=True to dataclass if python 3.11
@dataclass()
class RRAM_dev:
    LRS   : float
    HRS   : float
    # amplification for re-input to crossbar array, may vary for prob
    Max_Sat_amp : float
    Max_Cut_amp : float
    ref   : str
    def __str__(self):
        return f"{self.ref}"

#   scale is acquired through simulation data... very tediously
HfO2      = RRAM_dev(10e1    ,7*10e6 ,1e11      ,1e12   ,"Lee,2008")
#working with no sigmoid HfHfO2    = RRAM_dev(10e5    ,10e6   ,2.5e15    ,0    ,"intrinsic switching variability")
#FIXME: WRONG RESITANCE VALUES, should be 1e5, 1e6
HfHfO2    = RRAM_dev(1e5    ,1e6   ,3.8e13    ,0    ,"intrinsic switching variability")
MTJ_INC   = RRAM_dev(1000    ,3000   ,1.25e12   ,0   ,"N/A")

'''
ZnO       = RRAM_dev(2*10e-1 ,2*10e6 ,1    ,"Kim 2009")
MnO2      = RRAM_dev(2*10e-1 ,10e4   ,1    ,"Yang 2009")
MnO_Ta2O5 = RRAM_dev(5*10e3  ,10e9   ,1    ,"Hu Q 2018")
test2     = RRAM_dev(1000    ,7000   ,3.5e7 ,"N/A")
'''
