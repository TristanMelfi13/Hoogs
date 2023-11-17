import pandas as pd

Seq_Length = 12
Num_Strands = 2
CurrentFrame = pd.read_csv("Structures/Duplex0.gro")


class Nucleotide:
    Name = ""
    Atom_Locations = {}

    def __init__(self, name):
        self.Name = name


Nucleotide_List_For_CurrentFrame = []
Name_List = []

Atom_Data = CurrentFrame.values
print(Atom_Data)

Seen_End = 0
Iter = 32

while Seen_End != Num_Strands:
    # Iter Accesses the ith atom information, [0] gets the data as a string,
    # .split() removes whitespace, [0] accesses the full resi name, [2:3] Gets just the C/A/G/T
    # Note that the first resi should be n - 2, & the last n + 1
    Nucleo_Dict = {"C": 29, "A": 31, "G": 32, "T": 31}
    Nucleo_Name = Atom_Data[Iter][0].split()[0]
    Alter = 0
    if "1D" == Nucleo_Name[0:2]:
        Alter = -2
    elif str(Seq_Length) == Nucleo_Name[0:2]:
        Alter = 1
    Atom_Count = Nucleo_Dict.get(Nucleo_Name[-1]) + Alter
    break

