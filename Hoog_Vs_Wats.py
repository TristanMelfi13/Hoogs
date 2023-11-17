import pandas as pd
CurrentFrame = pd.read_csv("Structures/Duplex0.gro")


class Residue:
    Name = ""
    Atom_Locations = {}

    def __init__(self, name):
        self.Name = name

    def AddAtomLoc(self, AtomName, Xcor, Ycor, Zcor):
        Location = [Xcor, Ycor, Zcor]
        self.Atom_Locations[AtomName] = Location


Name_List = []

Atom_Data = CurrentFrame.values
print(Atom_Data)

def CreateObjects(Atom_Data, Seq_Length, Strands):
    Iter = 1
    Resi_List = []
    Seen_End = 0
    while Seen_End != Strands:
        # Iter Accesses the ith atom information, [0] gets the data as a string,
        # .split() removes whitespace, [0] accesses the full resi name, [2:3] Gets just the C/A/G/T
        # Note that the first resi should be n - 2, & the last n + 1
        Nucleo_Name = Atom_Data[Iter][0].split()[0]
        Resi_List.append(Residue(Nucleo_Name))
        for i in range(Iter, 10000):
            Current_Atom = Atom_Data[i][0].split()
            if Current_Atom[0] != Nucleo_Name:
                Iter = i
                break
            else:
                Resi_List[len(Resi_List) - 1].AddAtomLoc(Current_Atom[1], Current_Atom[3], Current_Atom[4], Current_Atom[5])
        print(Atom_Data[370][0].split()[0][0:2])
        break
        if Nucleo_Name[0:2] == str(Strands):
            Seen_End += 1
    print(len(Resi_List[0].Atom_Locations))
    return Resi_List


Current_Residue_Locations = CreateObjects(Atom_Data, Seq_Length=12,Strands=2)
