# -*- coding: utf-8 -*-
"""
Created on Thu Mar 25 10:18:24 2021

@author: Han Lab
"""

import os, sys

def makedir(dr):
    if not os.path.exists(dr): os.mkdir(dr)
    return 

daynum = str(input("What day of recording is it? (e.g. d4) "))
#animals we are recording from
ans = ["e156", "e157", "e158"]
#internal HD src
src = r"Z:\\"
#iterate thru the animal directories
for an in ans:
    andir = os.path.join(src, an)
    #make day of recording dir
    recordir = os.path.join(andir, daynum); makedir(recordir)
    #make the rest of the folders
    #behavior
    makedir(os.path.join(recordir, "behavior"))
    makedir(os.path.join(recordir, "behavior", "vr"))
    makedir(os.path.join(recordir, "behavior", "clampex"))
    #videos
    makedir(os.path.join(recordir, "eye"))
    makedir(os.path.join(recordir, "tail"))

print("\n*************made folders!*************\n")