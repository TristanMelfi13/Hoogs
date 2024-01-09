import math
import numpy as np
import time
import pandas as pd
start_time = time.time()
CurrentFrame = pd.read_csv("Structures/Duplex0.gro")


#  For GC Base Pairs; need Dists for Resis: G O6 <-->  C N4(H41/2)
#                                           G N1(H1) <--> C N1
#                                           G N2(H21/2) <--> C O2
#=============================================================#
#For AT Base Pairs; Need Dists for Resis: T O4 <--> A N6(H61/2)
#                                         T N1(H3) <--> A N1
class Residue:
    Name = ""
    Atom_Locations = {}

    def __init__(self, name):
        self.Name = name

    def AddAtomLoc(self, AtomName, Xcor, Ycor, Zcor):
        Location = [float(Xcor), float(Ycor), float(Zcor)]
        self.Atom_Locations[AtomName] = Location

Name_List = []
Atom_Data = CurrentFrame.values
#  print(Atom_Data)
def distance(x1, y1, z1, x2, y2, z2):
    return math.sqrt(math.pow(x2 - x1, 2) + math.pow(y2 - y1, 2) + math.pow(z2 - z1, 2))
def mag(x, y, z):
    return math.sqrt(math.pow(x, 2) + math.pow(y, 2) + math.pow(z, 2))
def CreateObjects(Atom_Data, Seq_Length, Strands):
    Iter = 1
    Resi_List = []
    Seen_End = 0
    Current_List_Pos = 0
    while len(Resi_List) != 24:
        # Iter Accesses the ith atom information, [0] gets the data as a string,
        # .split() removes whitespace, [0] accesses the full resi name, [2:3] Gets just the C/A/G/T
        # Note that the first resi should be n - 2, & the last n + 1
        # When looking at just the name, run Name[0:2]
        Resi_Name = Atom_Data[Iter][0].split()[0]
        Resi_List.append(Residue(Resi_Name))
        for i in range(Iter, 1000):
            #  Atom_Data contains all the information from our .gro file...
            Resi_For_Check = Atom_Data[i][0].split()  # Get the resi info for this part of the loop
            if Resi_For_Check[0] != Resi_Name:
                Iter = i
                Current_List_Pos += 1
                break
            else:
                if len(Resi_List) == 0:
                    Resi_List[0].AddAtomLoc("{}_{}".format(Resi_Name, Resi_For_Check[1]), Resi_For_Check[3], Resi_For_Check[4], Resi_For_Check[5])
                else:
                    Resi_List[len(Resi_List) - 1].AddAtomLoc("{}_{}".format(Resi_Name, Resi_For_Check[1]), Resi_For_Check[3], Resi_For_Check[4], Resi_For_Check[5])
    return Resi_List
def FindHBondsCtoG(Cytosine, Guanine):
    # Cytosine Residues

    H41 = Cytosine.Atom_Locations.get("{}_H41".format(Cytosine.Name))
    H42 = Cytosine.Atom_Locations.get("{}_H42".format(Cytosine.Name))
    N3 = Cytosine.Atom_Locations.get("{}_N3".format(Cytosine.Name))
    O2 = Cytosine.Atom_Locations.get("{}_O2".format(Cytosine.Name))

    # Guanine Residues
    H1 = Guanine.Atom_Locations.get("{}_H1".format(Guanine.Name))
    O6 = Guanine.Atom_Locations.get("{}_O6".format(Guanine.Name))
    H21 = Guanine.Atom_Locations.get("{}_H21".format(Guanine.Name))
    H22 = Guanine.Atom_Locations.get("{}_H22".format(Guanine.Name))
    # Distances in Angstroms

    Dist_O6_H41 = distance(H41[0], H41[1], H42[2], O6[0], O6[1], O6[2]) * 10
    Dist_O6_H42 = distance(H42[0], H42[1], H42[2], O6[0], O6[1], O6[2]) * 10
    Dist_H1_N3 = distance(H1[0], H1[1], H1[2], N3[0], N3[1], N3[2]) * 10
    Dist_H21_O2 = distance(H21[0], H21[1], H21[2], O2[0], O2[1], O2[2]) * 10
    Dist_H22_O2 = distance(H22[0], H22[1], H22[2], O2[0], O2[1], O2[2])


    AllGood = 0

    if Dist_O6_H41 < 3 or Dist_O6_H42 < 3:
        AllGood += 1
    if Dist_H1_N3 < 3:
        AllGood += 1
    if Dist_H21_O2 < 3 or Dist_H22_O2 < 3:
        AllGood += 1
    #  print("Dist: {} -- {}: {} || {}".format(Cytosine.Name, Guanine.Name, Dist_H1_N3, AllGood == 3))
    return AllGood == 3
def FindPlaneAngleCtoG(Cytosine, Guanine):
    # Jus for sanity
    # Note Cytosine is our cytosine
    # Resi 2 is our guanine
    C_N1 = (Cytosine.Atom_Locations.get("{}_N1".format(Cytosine.Name)))
    C_N4 = Cytosine.Atom_Locations.get("{}_N4".format(Cytosine.Name))
    C_O2 = Cytosine.Atom_Locations.get("{}_O2".format(Cytosine.Name))

    C_V1 = [C_N4[0] - C_N1[0], C_N4[1] - C_N1[1], C_N4[2] - C_N1[2]]
    C_V2 = [C_O2[0] - C_N1[0], C_O2[1] - C_N1[1], C_O2[2] - C_N1[2]]
    C_Normal = np.cross(C_V1, C_V2)



    G_C8 = Guanine.Atom_Locations.get("{}_C8".format(Guanine.Name))
    G_O6 = Guanine.Atom_Locations.get("{}_C4".format(Guanine.Name))
    G_N2 = Guanine.Atom_Locations.get("{}_C5".format(Guanine.Name))

    G_V1 = [G_C8[0] - G_O6[0], G_C8[1] - G_O6[1], G_C8[2] - G_O6[2]]
    G_V2 = [G_N2[0] - G_O6[0], G_N2[1] - G_O6[1], G_N2[2] - G_O6[2]]
    G_Normal = np.cross(G_V1, G_V2)

    Angle_Rads = (math.acos(np.dot(C_Normal, G_Normal) / (mag(C_Normal[0], C_Normal[1], C_Normal[2]) * mag(G_Normal[0], G_Normal[1], G_Normal[2]))))
    return Angle_Rads * (180 / np.pi)
def FindHBondsAtoT(Adenine, Thymine):
    #  Adenine
    H61 = Adenine.Atom_Locations.get("{}_H61".format(Adenine.Name))
    H62 = Adenine.Atom_Locations.get("{}_H62".format(Adenine.Name))
    N1 = Adenine.Atom_Locations.get("{}_N1".format(Adenine.Name))

    #  Thymine
    H3 = Thymine.Atom_Locations.get("{}_H3".format(Thymine.Name))
    O4 = Thymine.Atom_Locations.get("{}_O4".format(Thymine.Name))

    Dist_H61_O4 = distance(H61[0], H61[1], H61[2], O4[0], O4[1], O4[2]) * 10
    Dist_H62_O4 = distance(H62[0], H62[1], H62[2], O4[0], O4[1], O4[2]) * 10
    Dist_N1_H3 = distance(N1[0], N1[1], N1[2], H3[0], H3[1], H3[2]) * 10

    AllGood = 0

    if Dist_H61_O4 < 3 or Dist_H62_O4 < 3:
        AllGood += 1
    if Dist_N1_H3 < 3:
        AllGood += 1

    return AllGood == 2
def FindPlaneAngleAtoT(Adenine, Thymine):
    # Adenine Plane
    C8 = Adenine.Atom_Locations.get("{}_C8".format(Adenine.Name))
    N6 = Adenine.Atom_Locations.get("{}_N6".format(Adenine.Name))
    C2 = Adenine.Atom_Locations.get("{}_N3".format(Adenine.Name))
    A_V1 = [C8[0] - N6[0], C8[1] - N6[1], C8[2] - N6[2]]
    A_V2 = [C2[0] - N6[0], C2[1] - N6[1], C2[2] - N6[2]]
    A_normal = np.cross(A_V1, A_V2)


    # Thymine Plane
    N1 = Thymine.Atom_Locations.get("{}_N1".format(Thymine.Name))
    O4 = Thymine.Atom_Locations.get("{}_O4".format(Thymine.Name))
    O2 = Thymine.Atom_Locations.get("{}_O2".format(Thymine.Name))
    T_V1 = [N1[0] - O4[0], N1[1] - O4[1], N1[2] - O4[2]]
    T_V2 = [O2[0] - O4[0], O2[1] - O4[1], O2[2] - O4[2]]
    T_normal = np.cross(T_V1, T_V2)

    Angle_Rads = math.acos((np.dot(A_normal, T_normal) / (mag(A_normal[0], A_normal[1], A_normal[2]) * mag(T_normal[0], T_normal[1], T_normal[2]))))
    return Angle_Rads * (180 / np.pi)
def FindWatsAtoT(Adenine, ResiList):
    for i in range(len(ResiList)):
        if "T" in ResiList[i].Name:
            Thymine = ResiList[i]
            HBonds_Good = FindHBondsAtoT(Adenine, Thymine)
            Plane_Angle = FindPlaneAngleAtoT(Adenine, Thymine)
            if (Plane_Angle < 20 or Plane_Angle > 160) and HBonds_Good:
                print("Effective Watson Crick Pair Between: {} <===> {}".format(Adenine.Name, Thymine.Name))
def FindWatsCtoG(Cytosine, ResiList):
    #  We pass in a C, look for G's
    #  For GC Base Pairs; need Dists for Resis: G O6 <-->  C N4(H41/2)
    #                                           G N1(H1) <--> C N1
    #                                           G N2(H21/2) <--> C O2
    for i in range(len(ResiList)):
        if "G" in ResiList[i].Name:
            Guanine = ResiList[i]
            HBonds_Good = FindHBondsCtoG(Cytosine, Guanine)
            Plane_Angle = FindPlaneAngleCtoG(Cytosine, Guanine)
            if HBonds_Good and (Plane_Angle > 160 or Plane_Angle < 20):
                print("Effective Watson-Crick Pair Between: {} -- {}".format(Cytosine.Name, Guanine.Name))
def FindWats(ResiList):
    for i in range(len(ResiList)):
        Current_Resi = ResiList[i]
        if "C" in Current_Resi.Name:
            FindWatsCtoG(Current_Resi, ResiList)
        if "A" in Current_Resi.Name:
            FindWatsAtoT(Current_Resi, ResiList)

#

Current_Residue_Locations = CreateObjects(Atom_Data, Seq_Length=12, Strands=2)


FindWats(Current_Residue_Locations)
#  print(Current_Residue_Locations[0].Atom_Locations.get("1DC_O5'"))

end_time = time.time()
print("Execution time: {}".format(round(end_time - start_time, 2)))
