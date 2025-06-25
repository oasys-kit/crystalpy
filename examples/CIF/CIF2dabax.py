#
# from a CIF file:
# i) CIF_GetCrystal(cif_filename) generatesa dictionary (xraylib format)
# ii) get_dabax_text(CIF_GetCrystal(cif_filename)) generates text for a new entry for DABAX Crystals.dat
#
# It uses the library dans_diffraction https://github.com/DanPorter/Dans_Diffraction
#

import Dans_Diffraction as dif
import re

import numpy
import xraylib
from dabax.common_tools import bragg_metrictensor

def no_digits(text):
    return ''.join([c for c in text if not c.isdigit()])

def to_sentence_case(text):
    # Split into sentences
    sentences = re.split('([.!?] *)', text)
    # Capitalize each sentence part
    sentences = [s.capitalize() for s in sentences]
    return ''.join(sentences)

def CIF_GetCrystal(entry_name='YB66'):
    """
    parse a complex crystal structure from CIF file into a dictionary (like xraylib.Crystal_GetCrystal('Si'))

    it has an additional fields for each atom: the charge

    return a dictionary containing crystal information
    """

    try:
        xtl = dif.Crystal(entry_name)
        flag_found = True
    except:
        flag_found = False

    if not flag_found:
        raise (Exception("CIF filename not loaded: %s" % entry_name))

    # print(xtl)

    uvw, type1, label, occupancy, uiso, mxmymz = xtl.Structure.get()

    Zs = []

    for i in range(len(type1)):
        try:
            Z = xraylib.SymbolToAtomicNumber(to_sentence_case(no_digits(label[i])))
        except:
            Z = -1
        Zs.append(Z)

    cryst = {'name': entry_name}  # returned dictionary like that one created by xraylib.Crystal_GetCrystal(descriptor)

    cryst['a']     = xtl.Cell.a
    cryst['b']     = xtl.Cell.b
    cryst['c']     = xtl.Cell.c
    cryst['alpha'] = xtl.Cell.alpha
    cryst['beta']  = xtl.Cell.beta
    cryst['gamma'] = xtl.Cell.gamma

    volume = bragg_metrictensor(xtl.Cell.a    ,
                                xtl.Cell.b    ,
                                xtl.Cell.c    ,
                                xtl.Cell.alpha,
                                xtl.Cell.beta ,
                                xtl.Cell.gamma,
                                RETURN_VOLUME=1)

    cryst['volume'] = volume
    cryst['n_atom'] = len(type1)
    atom = []

    for i in range(len(type1)):
        if 1: # cell_data.shape[0] == 5:  # standard 5 columns
            # not here, this info is not in the dabax file
            # s = symbol_to_from_atomic_number(int(cell_data[0,i]))
            atom.append({
                'Zatom': int(Zs[i]),
                'fraction': occupancy[i],
                'x': uvw[i][0],
                'y': uvw[i][1],
                'z': uvw[i][2],
                'charge': 0.0, })
        else:  # 6 columns (charge)
            # 'AtomicName' required to compatible my current code
            # s = symbol_to_from_atomic_number(int(cell_data[0,i]))
            # if cell_data[5, i] != 0:  #charged
            #     s = s + f'%+.6g'%cell_data[5, i]
            pass

    cryst['atom'] = atom
    cryst['cpointer'] = None



    if 0:  # found Anisotropic coefficients in the header, process it
        ANISO_KEY = "UANISO_COFF"  # prefix for a line with anisotropic coefficients

        AnisoItem = {'Name': '       ',
                     'start': 0,
                     'end': 0,
                     'beta11': 0.0,
                     'beta22': 0.0,
                     'beta33': 0.0,
                     'beta12': 0.0,
                     'beta13': 0.0,
                     'beta23': 0.0}

        a = sorted(a, key=lambda x: int(x[1][0]),
                   reverse=False)  # sort 'Start' ascendant, avoid order changed by the SpecFile
        n = 0
        Aniso = []
        for x in a:  # tuple('UANISO_COFF_B1',[1 96 0.00038 0.00044 0.00039 0 0 0])
            AnisoItem['Name'] = x[0][len(ANISO_KEY) + 1:]  # get site atom name starting from 13th character 'B1', etc
            AnisoItem['start'] = int(x[1][0])
            AnisoItem['end'] = int(x[1][1])
            AnisoItem['beta11'] = float(x[1][2])
            AnisoItem['beta22'] = float(x[1][3])
            AnisoItem['beta33'] = float(x[1][4])
            AnisoItem['beta12'] = float(x[1][5])
            AnisoItem['beta13'] = float(x[1][6])
            AnisoItem['beta23'] = float(x[1][7])
            Aniso.append(AnisoItem.copy())
            n = n + 1
        cryst['Aniso'] = Aniso  # if having key 'Ansio' when there is anisotropic data,otherwise no
        cryst['n_aniso'] = n
    else:  # create a dummy Aniso to compatible with xraylib
        AnisoItem = {'Name': '       ',
                     'start': 0,
                     'end': 0,
                     'beta11': 0.0,
                     'beta22': 0.0,
                     'beta33': 0.0,
                     'beta12': 0.0,
                     'beta13': 0.0,
                     'beta23': 0.0}
        cryst['Aniso'] = [AnisoItem.copy()]
        cryst['n_aniso'] = 1

    return cryst

def get_dabax_text(xx):

    txt = "\n"
    txt += "#S %d %s\n" % (200, cif_filename)
    txt += "#UCELL %f %f %f %f %f %f\n" % (
                                            xx["a"],
                                            xx["b"],
                                            xx["c"],
                                            xx["alpha"],
                                            xx["beta"],
                                            xx["gamma"],
                                            )
    txt += "#UTEMP 298\n"
    txt += "#UCOMMENT %s\n" % "http://www.crystallography.net/cod/2102909.html"
    txt += "#N  %d\n" % (xx["n_atom"])
    txt += "#L  AtomicNumber  Fraction  X  Y  Z\n"

    for at in xx["atom"]:
        txt += "%d %f %f %f %f\n" % (at["Zatom"], at["fraction"], at["x"], at["y"], at["z"])
    txt += "\n"
    return txt


def Dans_Diffraction_example1(cif_filename=""):
    """
    Dans_Diffraction Examples
    Build your very own crystal structure,
    using lattice parameters and atomic coordinates
    """
    import matplotlib.pyplot as plt  # Plotting
    import sys, os

    cf = os.path.dirname(__file__)
    sys.path.insert(0, os.path.join(cf, '..'))

    if cif_filename == "":
        xtl = dif.Crystal()

        xtl.name = 'Oh! What a Lovely Crystal'
        xtl.new_cell([2.8,2.8,6.0,90,90,90])
        #xtl.new_atoms(u=[0,0.5], v=[0,0.5], w=[0,0.25], type=['Na','O'], label=['Na1','O1'], occupancy=None, uiso=None, mxmymz=None)
        xtl.Atoms.changeatom(0,u=0,   v=0,   w=0,    type='Na',label='Na1',occupancy=None, uiso=None, mxmymz=None) # there is an Fe ion added by default
        xtl.Atoms.addatom(u=0.5, v=0.5, w=0.25, type='O', label='O1', occupancy=None, uiso=None, mxmymz=None)
        #xtl.Symmetry.addsym('x,y,z+1/2')
        xtl.Symmetry.load_spacegroup(7)
        xtl.generate_structure() # apply symmetry to atomic positions.
    else:
        xtl = dif.Crystal(cif_filename)

    print(xtl.info())
    xtl.Plot.plot_crystal()
    plt.show()

    # Plot Powder pattern:
    xtl.Plot.simulate_powder(energy_kev=17)
    plt.show()

if __name__ == "__main__":

    cif_filename = 'BiFeO3.cif'

    # plot structure and powder pattern
    # Dans_Diffraction_example1(cif_filename)

    # get the dict in xraylib format
    xx = CIF_GetCrystal(cif_filename)
    for key in xx.keys():
        print(key, xx[key])

    # print text for DABAX Crystals.dat
    print("------------------------%s----------------------------------"  % get_dabax_text(xx))

