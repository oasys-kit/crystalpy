
#
#
# This example shows the diffraction by a Si 111 crystal calculated in its simplest implementation:
#
#
#    - calculate_simple_diffraction()
#      Uses a crystal setup and calculates the complex transmitivity and reflectivity
#
#
import numpy

try:
    from xoppylib.crystals.tools import bragg_calc2, run_diff_pat
    import xraylib
except:
    raise Exception("This example needs xoppylib and must be installed")
from dabax.dabax_xraylib import DabaxXraylib

from crystalpy.diffraction.GeometryType import BraggDiffraction, BraggTransmission, LaueDiffraction, LaueTransmission
from crystalpy.diffraction.DiffractionSetupXraylib import DiffractionSetupXraylib
from crystalpy.diffraction.Diffraction import Diffraction
from crystalpy.util.Photon import Photon

from srxraylib.plot.gol import plot

def plot_crystal_sketch(v0,vH,H):

    v0_h = v0.components()[1]
    v0_v = v0.components()[2]

    vH_h = vH.components()[1]
    vH_v = vH.components()[2]

    H_h = H.components()[1]
    H_v = H.components()[2]

    plot_crystal_sketch_components(v0_h,v0_v,vH_h ,vH_v ,H_h ,H_v)

def plot_crystal_sketch_components(v0_h,v0_v,vH_h ,vH_v ,H_h ,H_v):

    import matplotlib.pyplot as plt
    from matplotlib.path import Path
    import matplotlib.patches as patches

    hshift = 0.2
    vshift = 0.0

    plt.figure(1, figsize=(6,6))
    ax = plt.subplot(111)
    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('none')
    ax.xaxis.set_ticklabels([])
    ax.yaxis.set_ticklabels([])


    # draw axes

    plt.xlim([-2.2,2.2])
    plt.ylim([-2.2,2.2])

    ax.annotate("",
                xy    =(0.0,0.0), xycoords='data',
                xytext=(2.0,0.0), textcoords='data',
                arrowprops=dict(arrowstyle="<-",connectionstyle="arc3"),)
    plt.text(2, 0,"$x_2$", color='k')

    ax.annotate("",
                xy    =(0.0,0.0), xycoords='data',
                xytext=(0.0,2.0), textcoords='data',
                arrowprops=dict(arrowstyle="<-",connectionstyle="arc3"),)
    plt.text(0, 2,"$x_3$", color='k')


    # draw vectors

    ax.annotate("",
                xy    =(-v0_h,-v0_v), xycoords='data',
                xytext=(0.0,0.0), textcoords='data',
                arrowprops=dict(arrowstyle="<-",connectionstyle="arc3",color='red'),)
    plt.text(-v0_h+hshift,-v0_v+vshift, r"$\vec k_0$", color='r')


    ax.annotate("",
                xy    =(0,0), xycoords='data',
                xytext=(vH_h,vH_v), textcoords='data',
                arrowprops=dict(arrowstyle="<-",connectionstyle="arc3",color='red'),)
    plt.text(vH_h+hshift,vH_v+vshift, r"$\vec k_H$", color='r')


    ax.annotate("",
                xy    =(0,0), xycoords='data',
                xytext=(H_h,H_v), textcoords='data',
                arrowprops=dict(arrowstyle="<-",connectionstyle="arc3",color='blue'),)
    plt.text(H_h+hshift,H_v+vshift, r"$\vec H$", color='b')


    # draw Bragg plane

    ax.annotate("",
                xy    =(0,0), xycoords='data',
                xytext=( -H_v*1.5, H_h*1.5), textcoords='data',
                arrowprops=dict(arrowstyle="-",connectionstyle="arc3",color='green'),)

    ax.annotate("",
                xy    =(0,0), xycoords='data',
                xytext=(H_v*1.5,-H_h*1.5), textcoords='data',
                arrowprops=dict(arrowstyle="-",connectionstyle="arc3",color='green'),)

    # draw crystal
    #
    x1 = -0.8
    y1 = -0.1
    x2 =  0.8
    y2 =  0.0

    verts = [
        (x1,y1), # left, bottom
        (x2,y1), # left, top
        (x2,y2), # right, top
        (x1,y2), # right, bottom
        (x1,y1), # ignored
        ]


    codes = [Path.MOVETO,
             Path.LINETO,
             Path.LINETO,
             Path.LINETO,
             Path.CLOSEPOLY,
             ]

    path = Path(verts, codes)
    patch = patches.PathPatch(path, facecolor='orange', lw=2)
    ax.add_patch(patch)

    plt.show()

#
def calculate_with_crystalpy(bragg_or_laue=0,  #
                            diffracted_or_transmitted=0, #
                            crystal_name           = "Si",    # string
                            thickness              = 1e-2,    # meters
                            miller_h               = 1,       # int
                            miller_k               = 1,       # int
                            miller_l               = 1,       # int
                            asymmetry_angle        = 0.0,     # radians
                            energy                 = 8000.0,  # eV
                            angle_deviation_min    = -100e-6, # radians
                            angle_deviation_max    = 100e-6,  # radians
                            angle_deviation_points = 500,
                            ):

    if bragg_or_laue == 0:
        if diffracted_or_transmitted == 0:
            geometry_type = BraggDiffraction()
        elif diffracted_or_transmitted == 1:
            geometry_type = BraggTransmission()
        else:
            raise Exception("Bad geometry type")
    elif bragg_or_laue == 1:
        if diffracted_or_transmitted == 0:
            geometry_type = LaueDiffraction()
        elif diffracted_or_transmitted == 1:
            geometry_type = LaueTransmission()
        else:
            raise Exception("Bad geometry type")
    else:
        raise Exception("Bad geometry type")

    # Create a diffraction setup.

    print("\nCreating a diffraction setup...")
    diffraction_setup = DiffractionSetupXraylib(geometry_type          = geometry_type,
                                               crystal_name           = crystal_name,
                                               thickness              = thickness,
                                               miller_h               = miller_h,
                                               miller_k               = miller_k,
                                               miller_l               = miller_l,
                                               asymmetry_angle        = asymmetry_angle,
                                               azimuthal_angle        = 0.0)

    angle_step = (angle_deviation_max-angle_deviation_min)/angle_deviation_points

    #
    # gets Bragg angle needed to create deviation's scan
    #
    bragg_angle = diffraction_setup.angleBragg(energy)

    print("Bragg angle for E=%f eV is %f deg"%(energy,bragg_angle*180.0/numpy.pi))


    # Create a Diffraction object (the calculator)
    diffraction = Diffraction()

    # initialize arrays for storing outputs
    deviations = numpy.zeros(angle_deviation_points)
    intensityS = numpy.zeros(angle_deviation_points)
    intensityP = numpy.zeros(angle_deviation_points)

    # k_0_unitary = diffraction_setup.incomingPhotonDirection(energy, 0.0)
    # photon_0 = Photon(energy_in_ev=energy,direction_vector=k_0_unitary)
    # k_H_unitary = diffraction_setup._
    # print(">>>>>>>>>>>>>>>>>>>>>>>>k_0: ",k_0_unitary._components )

    # plot_crystal_sketch(k_0_unitary,k_0_unitary,)


    for ia in range(angle_deviation_points):
        deviation = angle_deviation_min + ia * angle_step

        # angle = deviation  + bragg_angle + asymmetry_angle
        # # calculate the components of the unitary vector of the incident photon scan
        # # Note that diffraction plane is YZ
        # yy = numpy.cos(angle)
        # zz = - numpy.abs(numpy.sin(angle))
        # photon = Photon(energy_in_ev=energy,direction_vector=Vector(0.0,yy,zz))

        k_unitary = diffraction_setup.incomingPhotonDirection(energy, deviation)

        # or equivalently
        # k_0_unitary = diffraction_setup.incomingPhotonDirection(energy, 0.0)
        # k_unitary = k_0_unitary.rotateAroundAxis( Vector(1.0,0.0,0.0), -deviation)

        photon = Photon(energy_in_ev=energy,direction_vector=k_unitary)

        # perform the calculation
        coeffs = diffraction.calculateDiffractedComplexAmplitudes(diffraction_setup,photon)

        # store results
        deviations[ia] = deviation
        intensityS[ia] = numpy.abs(coeffs['S'])**2
        intensityP[ia] = numpy.abs(coeffs['P'])**2

    return deviations,intensityS,intensityP


def calculate_with_xoppy(bragg_or_laue=0,  #
                            diffracted_or_transmitted=0, #
                            crystal_name           = "Si",    # string
                            thickness              = 1e-2,    # meters
                            miller_h               = 1,       # int
                            miller_k               = 1,       # int
                            miller_l               = 1,       # int
                            asymmetry_angle        = 0.0,     # radians
                            energy                 = 8000.0,  # eV
                            angle_deviation_min    = -100e-6, # radians
                            angle_deviation_max    = 100e-6,  # radians
                            angle_deviation_points = 500,
                            ):

    MILLER_INDEX_H = miller_h
    MILLER_INDEX_K = miller_k
    MILLER_INDEX_L = miller_l
    TEMPER = 1.0
    MOSAIC = 0


    if bragg_or_laue == 0:
        if diffracted_or_transmitted == 0:
            GEOMETRY = 0
        elif diffracted_or_transmitted == 1:
            GEOMETRY = 2
        else:
            raise Exception("Bad geometry type")
    elif bragg_or_laue == 1:
        if diffracted_or_transmitted == 0:
            GEOMETRY = 1
        elif diffracted_or_transmitted == 1:
            GEOMETRY = 3
        else:
            raise Exception("Bad geometry type")
    else:
        raise Exception("Bad geometry type")


    SCAN = 2
    UNIT = 0 # rad
    SCANFROM = angle_deviation_min
    SCANTO = angle_deviation_max
    SCANPOINTS = angle_deviation_points
    ENERGY = energy
    ASYMMETRY_ANGLE = asymmetry_angle*180.0/numpy.pi
    THICKNESS = thickness*100
    MOSAIC_FWHM = 0.0
    RSAG = 0.0
    RMER = 0.0
    ANISOTROPY = 0
    POISSON = 0.0
    CUT = 0
    FILECOMPLIANCE = ""


    for file in ["diff_pat.dat","diff_pat.gle","diff_pat.par","diff_pat.xop","xcrystal.bra"]:
        try:
            os.remove(os.path.join(locations.home_bin_run(),file))
        except:
            pass


    if (GEOMETRY == 1) or (GEOMETRY == 3):
        if ASYMMETRY_ANGLE == 0.0:
            print("xoppy_calc_xcrystal: WARNING: In xcrystal the asymmetry angle is the angle between Bragg planes and crystal surface,"+
                  "in BOTH Bragg and Laue geometries.")


    descriptor = crystal_name
    if SCAN == 3: # energy scan
        emin = SCANFROM - 1
        emax = SCANTO + 1
    else:
        emin = ENERGY - 100.0
        emax = ENERGY + 100.0

    print("Using crystal descriptor: ",descriptor)



    bragg_dictionary = bragg_calc2(
                        descriptor=descriptor,
                        hh=MILLER_INDEX_H, kk=MILLER_INDEX_K, ll=MILLER_INDEX_L,
                        temper=float(TEMPER),
                        emin=emin, emax=emax, estep=5.0, fileout="xcrystal.bra",
                        ANISO_SEL=0,
                        do_not_prototype = 0,  # 0=use site groups (recommended), 1=use all individual sites
                        verbose = False,
                        material_constants_library = xraylib,
                        )

    run_diff_pat(
        bragg_dictionary,
        preprocessor_file="xcrystal.bra",
        descriptor=descriptor,
        MOSAIC=MOSAIC,
        GEOMETRY=GEOMETRY,
        SCAN=SCAN,
        UNIT=UNIT,
        SCANFROM=SCANFROM,
        SCANTO=SCANTO,
        SCANPOINTS=SCANPOINTS,
        ENERGY=ENERGY,
        ASYMMETRY_ANGLE=ASYMMETRY_ANGLE,
        THICKNESS=THICKNESS,
        MOSAIC_FWHM=MOSAIC_FWHM,
        RSAG=RSAG,
        RMER=RMER,
        ANISOTROPY=ANISOTROPY,
        POISSON=POISSON,
        CUT=CUT,
        FILECOMPLIANCE=FILECOMPLIANCE,
    )

    #show calculated parameters in standard output
    txt_info = open("diff_pat.par").read()
    for line in txt_info:
        print(line,end="")



    try:
        calculated_data =  numpy.loadtxt("diff_pat.dat", skiprows=5)

    except Exception as e:
        raise Exception("Error loading diff_pat.dat :" + str(e))

    print(calculated_data.shape)

    deviations = calculated_data[:,0].copy()
    intensityS = calculated_data[:,6].copy()
    intensityP = calculated_data[:,5].copy()

    return deviations,intensityS,intensityP


def input_cases(case):
    input_dict = {}
    if case == "bragg_symmetric":

        # # Bragg symmetric
        input_dict["bragg_or_laue"]             = 0
        input_dict["diffracted_or_transmitted"] = 0
        input_dict["crystal_name"]              = "Si"    # string
        input_dict["thickness"]                 = 1e-2    # meters
        input_dict["miller_h"]                  = 1       # int
        input_dict["miller_k"]                  = 1       # int
        input_dict["miller_l"]                  = 1       # int
        input_dict["asymmetry_angle"]           = 0*numpy.pi/180  # radians
        input_dict["energy"]                    = 8000.0             # eV
        input_dict["angle_deviation_min"]       = -100e-6            # radians
        input_dict["angle_deviation_max"]       = 300e-6             # radians
        input_dict["angle_deviation_points"]    = 500
    elif case == "laue_symmetric":
        # Laue symmetric
        input_dict["bragg_or_laue"]             = 1
        input_dict["diffracted_or_transmitted"] = 0
        input_dict["crystal_name"]              = "Si"    # string
        input_dict["thickness"]                 = 10e-6   # meters
        input_dict["miller_h"]                  = 1       # int
        input_dict["miller_k"]                  = 1       # int
        input_dict["miller_l"]                  = 1       # int
        input_dict["asymmetry_angle"]           = 90.0 * numpy.pi/180  # radians
        input_dict["energy"]                    = 8000.0             # eV
        input_dict["angle_deviation_min"]       = -100e-6            # radians
        input_dict["angle_deviation_max"]       = 300e-6             # radians
        input_dict["angle_deviation_points"]    = 500

    if case == "bragg_asymmetric":

        # # Bragg symmetric
        input_dict["bragg_or_laue"]             = 0
        input_dict["diffracted_or_transmitted"] = 0
        input_dict["crystal_name"]              = "Si"    # string
        input_dict["thickness"]                 = 1e-2    # meters
        input_dict["miller_h"]                  = 1       # int
        input_dict["miller_k"]                  = 1       # int
        input_dict["miller_l"]                  = 1       # int
        input_dict["asymmetry_angle"]           = 10*numpy.pi/180  # radians
        input_dict["energy"]                    = 8000.0             # eV
        input_dict["angle_deviation_min"]       = -100e-6            # radians
        input_dict["angle_deviation_max"]       = 300e-6             # radians
        input_dict["angle_deviation_points"]    = 500
    elif case == "laue_asymmetric":
        # Laue symmetric
        input_dict["bragg_or_laue"]             = 1
        input_dict["diffracted_or_transmitted"] = 0
        input_dict["crystal_name"]              = "Si"    # string
        input_dict["thickness"]                 = 10e-6   # meters
        input_dict["miller_h"]                  = 1       # int
        input_dict["miller_k"]                  = 1       # int
        input_dict["miller_l"]                  = 1       # int
        input_dict["asymmetry_angle"]           = 80.0 * numpy.pi/180  # radians
        input_dict["energy"]                    = 8000.0             # eV
        input_dict["angle_deviation_min"]       = -100e-6            # radians
        input_dict["angle_deviation_max"]       = 300e-6             # radians
        input_dict["angle_deviation_points"]    = 500

    return input_dict
#
# main
#
if __name__ == "__main__":


    if True:
        for case in ["bragg_symmetric","laue_symmetric","bragg_asymmetric","laue_asymmetric"]:
            input_dict = input_cases(case)




            angle1,intS1, intP1 = calculate_with_crystalpy(
                                    bragg_or_laue             =  input_dict["bragg_or_laue"],
                                    diffracted_or_transmitted =  input_dict["diffracted_or_transmitted"],
                                    crystal_name              =  input_dict["crystal_name"],
                                    thickness                 =  input_dict["thickness"],
                                    miller_h                  =  input_dict["miller_h"],
                                    miller_k                  =  input_dict["miller_k"],
                                    miller_l                  =  input_dict["miller_l"],
                                    asymmetry_angle           =  input_dict["asymmetry_angle"],
                                    energy                    =  input_dict["energy"],
                                    angle_deviation_min       =  input_dict["angle_deviation_min"],
                                    angle_deviation_max       =  input_dict["angle_deviation_max"],
                                    angle_deviation_points    =  input_dict["angle_deviation_points"],)

            # plot(angle1,intS1,legend=['S-pol crystalpy'])

            angle2,intS2, intP2 = calculate_with_xoppy(
                                    bragg_or_laue             =  input_dict["bragg_or_laue"],
                                    diffracted_or_transmitted =  input_dict["diffracted_or_transmitted"],
                                    crystal_name              =  input_dict["crystal_name"],
                                    thickness                 =  input_dict["thickness"],
                                    miller_h                  =  input_dict["miller_h"],
                                    miller_k                  =  input_dict["miller_k"],
                                    miller_l                  =  input_dict["miller_l"],
                                    asymmetry_angle           =  input_dict["asymmetry_angle"],
                                    energy                    =  input_dict["energy"],
                                    angle_deviation_min       =  input_dict["angle_deviation_min"],
                                    angle_deviation_max       =  input_dict["angle_deviation_max"],
                                    angle_deviation_points    =  input_dict["angle_deviation_points"],)



            plot(angle1,intS1,angle2,intS2,legend=['S-pol crystalpy','S-pol xoppy'],title=case)


            # v0_h = -7.5128581375702067E-002
            # v0_v = -0.99717385458127339
            # vH_h = 0.41165131819335760
            # vH_v = -0.91134142462069123
            # H_h = 0.98480775301220813
            # H_v = 0.17364817766693044
            f = open("diff_pat.gle",'r')
            txt = f.readlines()
            f.close()
            lines = txt[5:11]
            print(lines)
            for line in lines:
                exec(line.strip())
            plot_crystal_sketch_components(v0_h,v0_v,vH_h ,vH_v ,H_h ,H_v)