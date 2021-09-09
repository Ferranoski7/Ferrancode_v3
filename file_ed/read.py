from molecule import Molecule


def read_xyz(sysl):  # returns atm_l with ['Symbol',x,y,z]
    # ------------------  Open .bands file -------------------------
    bfile=sysl+'.xyz'
    print(' ')
    print('reading the {} file'.format(bfile))
    print(' ')


    f=open(bfile, 'r')

    # ------------------  read Header ------------------------------

    nat=int(f.readline())
    print("Number of atoms  {:4d}".format(nat))
    title =f.readline().strip()
    print(title)

    # ------------------  read Geometry ------------------------------
    sym_l=[]
    coord_a=[]
    for atom in range (nat):
        line=f.readline().strip()
        line_l=line.split()
        sym_l.append(line_l[0])
        coord_a.append([float(line_l[1]),float(line_l[2]),float(line_l[3])])
    atm_l=Molecule(sym_l, coord_a, sysl)
    f.close()
    return atm_l

def read_txt(sysl): #returns atm_l with ['Symbol',x,y,z]
    # ------------------  Open .bands file -------------------------

    bfile=sysl+'.txt'
    print(' ')
    print('reading the {} file'.format(bfile))
    print(' ')


    f=open(bfile, 'r')

    # ------------------  read Header ------------------------------

    for i in range(3):
        header =f.readline().strip()

    # ------------------  read Geometry ------------------------------



    line=f.readline().strip()
    line_l=line.split()
    sym_l=[]
    coord_a=[]
    while len(line_l)!=0:
        sym_l.append(line_l[0])
        coord_a.append([float(line_l[1]), float(line_l[2]), float(line_l[3])])
        line=f.readline().strip()
        line_l=line.split()
    atm_l = Molecule(sym_l, coord_a, sysl)
    f.close()
    return atm_l
