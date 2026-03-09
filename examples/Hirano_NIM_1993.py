# Hirano et al. Nucl. Instr. and Meth. A 336 (1993) 343-353

import numpy

from crystalpy.diffraction.GeometryType import BraggDiffraction, BraggTransmission
from crystalpy.diffraction.GeometryType import LaueDiffraction, LaueTransmission
from crystalpy.diffraction.DiffractionSetupXraylib import DiffractionSetupXraylib
from crystalpy.diffraction.Diffraction import Diffraction

from crystalpy.util.Vector import Vector
from crystalpy.util.Photon import Photon
from crystalpy.util.ComplexAmplitudePhoton import ComplexAmplitudePhoton
from crystalpy.util.ComplexAmplitudePhotonBunch import ComplexAmplitudePhotonBunch

from srxraylib.plot.gol import plot, plot_show

def matplotlib_combine(old_fig, old_axes, grid=1, show=0):
    import matplotlib.pyplot as plt

    # 1. Create your new "Master" 2x2 grid
    fig_master, axs_master = plt.subplots(2, 2, figsize=(12, 10))
    new_axes = axs_master.flatten()  # Makes it easier to loop (0, 1, 2, 3)

    # 2. Put your existing axes in a list
    # old_axes = [ax1_1, ax1_2, ax1_3, ax1_4]

    # 3. Loop through and "clone" the content
    for i, old_ax in enumerate(old_axes):
        target_ax = new_axes[i]
        if grid: target_ax.grid(True)

        # Clone all Lines (plot, step, etc.)
        for line in old_ax.get_lines():
            target_ax.plot(line.get_xdata(), line.get_ydata(),
                           label=line.get_label(),
                           color=line.get_color(),
                           linestyle=line.get_linestyle())

        # Clone all Collections (scatter plots)
        for coll in old_ax.collections:
            target_ax.add_collection(coll)  # This is a bit more advanced

        # Copy labels and titles
        target_ax.set_title(old_ax.get_title())
        target_ax.set_xlabel(old_ax.get_xlabel())
        target_ax.set_ylabel(old_ax.get_ylabel())
        if old_ax.get_legend():
            target_ax.legend()

    # 4. Close the old individual windows to save RAM
    for i in range(4):
        plt.close(old_fig[i])

    plt.tight_layout()
    if show: plt.show()

#
def calculate_fig3(calculation_method=0, calculation_strategy_flag=0):

    # Create a diffraction setup.

    print("\nCreating a diffraction setup...")
    diffraction_setup = DiffractionSetupXraylib(geometry_type          = BraggDiffraction(),  # GeometryType object
                                               crystal_name           = "Ge",                             # string
                                               thickness              = 1e-2,                             # meters
                                               miller_h               = 2,                                # int
                                               miller_k               = 2,                                # int
                                               miller_l               = 0,                                # int
                                               asymmetry_angle        = 0,#10.0*numpy.pi/180.,                              # radians
                                               azimuthal_angle        = 0.0)                              # radians                            # int


    energy                 = 10000.0                           # eV
    angle_deviation_min    = -100e-6                          # radians
    angle_deviation_max    = 100e-6                           # radians
    angle_deviation_points = 500

    angle_step = (angle_deviation_max-angle_deviation_min)/angle_deviation_points

    #
    # gets Bragg angle needed to create deviation's scan
    #
    bragg_angle = diffraction_setup.angleBragg(energy)
    bragg_angle_corrected = diffraction_setup.angleBraggCorrected(energy)
    half_darwin_s = diffraction_setup.darwinHalfwidthS(energy)

    print("Bragg angle for E=%f eV is %f deg"%(energy, bragg_angle_corrected * 180.0 / numpy.pi))
    print("Bragg angle for E=%f eV is %f deg" % (energy, bragg_angle * 180.0 / numpy.pi))
    print("Darwin width S for E=%f eV is %f deg" % (energy, 2 * half_darwin_s))



    # initialize arrays for storing outputs
    deviations = numpy.zeros(angle_deviation_points, dtype=float)
    intensityS = numpy.zeros(angle_deviation_points, dtype=float)
    intensityP = numpy.zeros(angle_deviation_points, dtype=float)
    phaseS = numpy.zeros(angle_deviation_points, dtype=float)
    phaseP = numpy.zeros(angle_deviation_points, dtype=float)

    for ia in range(angle_deviation_points):
        deviation = angle_deviation_min + ia * angle_step
        angle = deviation  + bragg_angle_corrected

        # calculate the components of the unitary vector of the incident photon scan
        # Note that diffraction plane is YZ
        yy = numpy.cos(angle)
        zz = - numpy.abs(numpy.sin(angle))
        photon = Photon(energy_in_ev=energy, direction_vector=Vector(0.0,yy,zz))

        # perform the calculation
        coeffs = Diffraction.calculateDiffractedComplexAmplitudes(diffraction_setup, photon,
                                                                  calculation_method=calculation_method,
                                                                  calculation_strategy_flag=calculation_strategy_flag)
        R_S = complex(coeffs['S'])
        R_P = complex(coeffs['P'])
        # store results
        deviations[ia] = deviation / half_darwin_s
        intensityS[ia] = numpy.abs(R_S) ** 2
        intensityP[ia] = numpy.abs(R_P) ** 2
        phaseS[ia] = numpy.angle(R_S)
        phaseP[ia] = numpy.angle(R_P)

    # plot results

    fig1, ax1 = plot(
        deviations, intensityS ** 4,
        deviations, intensityP ** 4,
        xtitle="deviation angle [normalised]",
        ytitle="Reflectivity",
        title="BraggDiffraction Fig 3; Ge220 * 4 reflections",
        legend=["Sigma-polarization","Pi-polarization"],
        xrange=[-2, 2],
        grid=1,
        show=0,
    )

    fig2, ax2 = plot(
        deviations, numpy.degrees(phaseP - phaseS) * 4,
        deviations, numpy.degrees(phaseP - phaseS) * 4,
        xtitle="deviation angle [normalised]",
        ytitle="Phase difference [deg]",
        title="BraggDiffraction Fig 3; PhaseP - PhaseS; Ge220 * 4 reflections",
        xrange=[-2, 2], yrange=[-180, 180],
        grid=1,
        show=0,
    )

    return fig1, fig2, ax1, ax2

#
def calculate_fig4(calculation_method=0, calculation_strategy_flag=0):

    # Create a diffraction setup.

    print("\nCreating a diffraction setup...")
    diffraction_setup = DiffractionSetupXraylib(geometry_type          = BraggTransmission(),  # GeometryType object
                                               crystal_name           = "Si",                             # string
                                               thickness              = 50e-6,                             # meters
                                               miller_h               = 2,                                # int
                                               miller_k               = 2,                                # int
                                               miller_l               = 0,                                # int
                                               asymmetry_angle        = 0,#10.0*numpy.pi/180.,                              # radians
                                               azimuthal_angle        = 0.0)                              # radians                            # int


    energy                 = 10000.0                           # eV
    angle_deviation_min    = -100e-6                          # radians
    angle_deviation_max    = 100e-6                           # radians
    angle_deviation_points = 1000

    angle_step = (angle_deviation_max-angle_deviation_min)/angle_deviation_points

    #
    # gets Bragg angle needed to create deviation's scan
    #
    bragg_angle = diffraction_setup.angleBragg(energy)
    bragg_angle_corrected = diffraction_setup.angleBraggCorrected(energy)
    half_darwin_s = diffraction_setup.darwinHalfwidthS(energy)

    print("Bragg angle for E=%f eV is %f deg"%(energy, bragg_angle_corrected * 180.0 / numpy.pi))
    print("Bragg angle for E=%f eV is %f deg" % (energy, bragg_angle * 180.0 / numpy.pi))
    print("Darwin width S for E=%f eV is %f deg" % (energy, 2 * half_darwin_s))



    # initialize arrays for storing outputs
    deviations = numpy.zeros(angle_deviation_points, dtype=float)
    intensityS = numpy.zeros(angle_deviation_points, dtype=float)
    intensityP = numpy.zeros(angle_deviation_points, dtype=float)
    phaseS = numpy.zeros(angle_deviation_points, dtype=float)
    phaseP = numpy.zeros(angle_deviation_points, dtype=float)

    for ia in range(angle_deviation_points):
        deviation = angle_deviation_min + ia * angle_step
        angle = deviation  + bragg_angle_corrected

        # calculate the components of the unitary vector of the incident photon scan
        # Note that diffraction plane is YZ
        yy = numpy.cos(angle)
        zz = - numpy.abs(numpy.sin(angle))
        photon = Photon(energy_in_ev=energy, direction_vector=Vector(0.0,yy,zz))

        # perform the calculation
        coeffs = Diffraction.calculateDiffractedComplexAmplitudes(diffraction_setup, photon,
                                                                  calculation_method=calculation_method,
                                                                  calculation_strategy_flag=calculation_strategy_flag)
        R_S = complex(coeffs['S'])
        R_P = complex(coeffs['P'])
        # store results
        deviations[ia] = deviation / half_darwin_s
        intensityS[ia] = numpy.abs(R_S) ** 2
        intensityP[ia] = numpy.abs(R_P) ** 2
        phaseS[ia] = numpy.angle(R_S)
        phaseP[ia] = numpy.angle(R_P)

    # plot results

    fig1, ax1 = plot(
        deviations, intensityS,
        deviations, intensityP,
        xtitle="deviation angle [normalised]",
        ytitle="Reflectivity",
        title="BraggTransmission Fig 4; Si220",
        legend=["Sigma-polarization","Pi-polarization"],
        xrange=[-5, 5], yrange=[0,1],
        grid=1,
        show=0,
    )

    phaseSu = numpy.unwrap(phaseS)
    phasePu = numpy.unwrap(phaseP)

    fig2, ax2 = plot(
        deviations, numpy.degrees(phasePu - phaseSu),
        deviations, numpy.degrees(phasePu - phaseSu),
        xtitle="deviation angle [normalised]",
        ytitle="Phase difference [deg]",
        title="BraggTransmission Fig 4; PhaseP - PhaseS; Si220",
        xrange=[-5, 5], yrange=[-180, 180],
        grid=1,
        show=0,
    )

    return fig1, fig2, ax1, ax2

#
def calculate_fig6(calculation_method=0, calculation_strategy_flag=0):

    # Create a diffraction setup.

    print("\nCreating a diffraction setup...")
    diffraction_setup = DiffractionSetupXraylib(geometry_type=LaueDiffraction(),  # GeometryType object
                                                crystal_name="Si",  # string
                                                thickness=50e-6,  # meters
                                                miller_h=2,  # int
                                                miller_k=2,  # int
                                                miller_l=0,  # int
                                                asymmetry_angle=numpy.pi / 2,
                                                # 10.0*numpy.pi/180.,                              # radians
                                                azimuthal_angle=0.0)  # radians                            # int

    energy = 10000.0  # eV
    angle_deviation_min = -100e-6  # radians
    angle_deviation_max = 100e-6  # radians
    angle_deviation_points = 1000

    angle_step = (angle_deviation_max - angle_deviation_min) / angle_deviation_points

    #
    # gets Bragg angle needed to create deviation's scan
    #
    bragg_angle = diffraction_setup.angleBragg(energy)
    bragg_angle_corrected = diffraction_setup.angleBraggCorrected(energy)
    half_darwin_s = diffraction_setup.darwinHalfwidthS(energy)

    print("Bragg angle for E=%f eV is %f deg" % (energy, bragg_angle_corrected * 180.0 / numpy.pi))
    print("Bragg angle for E=%f eV is %f deg" % (energy, bragg_angle * 180.0 / numpy.pi))
    print("Darwin width S for E=%f eV is %f deg" % (energy, 2 * half_darwin_s))

    # initialize arrays for storing outputs
    deviations = numpy.zeros(angle_deviation_points, dtype=float)
    intensityS = numpy.zeros(angle_deviation_points, dtype=float)
    intensityP = numpy.zeros(angle_deviation_points, dtype=float)
    phaseS = numpy.zeros(angle_deviation_points, dtype=float)
    phaseP = numpy.zeros(angle_deviation_points, dtype=float)

    for ia in range(angle_deviation_points):
        deviation = angle_deviation_min + ia * angle_step
        angle = deviation + bragg_angle_corrected

        # calculate the components of the unitary vector of the incident photon scan
        # Note that diffraction plane is YZ
        yy = numpy.cos(angle)
        zz = - numpy.abs(numpy.sin(angle))
        photon = Photon(energy_in_ev=energy, direction_vector=Vector(0.0, yy, zz))
        photon = ComplexAmplitudePhoton(energy_in_ev=energy,
                                        direction_vector=diffraction_setup.vectorIncomingPhotonDirection(
                                            energy,deviation,
                                            angle_center_flag=2,# 0=Absolute angle, 1=Theta Bragg Corrected, 2=Theta Bragg
                                        ),
                                        Esigma=1,
                                        Epi=1)

        # perform the calculation
        coeffs = Diffraction.calculateDiffractedComplexAmplitudes(diffraction_setup, photon,
                                                                  calculation_method=calculation_method,
                                                                  calculation_strategy_flag=calculation_strategy_flag)
        R_S = complex(coeffs['S'])
        R_P = complex(coeffs['P'])
        # store results
        deviations[ia] = deviation / half_darwin_s
        intensityS[ia] = numpy.abs(R_S) ** 2
        intensityP[ia] = numpy.abs(R_P) ** 2
        phaseS[ia] = numpy.angle(R_S)
        phaseP[ia] = numpy.angle(R_P)

    # plot results

    fig1, ax1 = plot(
        deviations, intensityS,
        deviations, intensityP,
        xtitle="deviation angle [normalised]",
        ytitle="Reflectivity",
        title="LaueDiffraction Fig 6; Si220",
        legend=["Sigma-polarization", "Pi-polarization"],
        xrange=[-5, 5], yrange=[0, 1],
        grid=1,
        show=0,
    )

    phaseSu = numpy.unwrap(phaseS)
    phasePu = numpy.unwrap(phaseP)

    fig2, ax2 = plot(
        deviations, numpy.degrees(phasePu - phaseSu),
        deviations, numpy.degrees(phasePu - phaseSu),
        xtitle="deviation angle [normalised]",
        ytitle="Phase difference [deg]",
        title="Laue Diffraction Fig 6; PhaseP - PhaseS; Si220",
        xrange=[-5, 5], yrange=[-180, 180],
        grid=1,
        show=0,
    )

    return fig1, fig2, ax1, ax2

#
def calculate_fig9(calculation_method=0, calculation_strategy_flag=0):

    # Create a diffraction setup.

    print("\nCreating a diffraction setup...")
    diffraction_setup = DiffractionSetupXraylib(geometry_type=LaueTransmission(),  # GeometryType object
                                                crystal_name="Si",  # string
                                                thickness=150e-6,  # meters
                                                miller_h=2,  # int
                                                miller_k=2,  # int
                                                miller_l=0,  # int
                                                asymmetry_angle=numpy.pi / 2,
                                                # 10.0*numpy.pi/180.,                              # radians
                                                azimuthal_angle=0.0)  # radians                            # int

    energy = 10000.0  # eV
    angle_deviation_min = -100e-6  # radians
    angle_deviation_max = 100e-6  # radians
    angle_deviation_points = 1000

    angle_step = (angle_deviation_max - angle_deviation_min) / angle_deviation_points

    #
    # gets Bragg angle needed to create deviation's scan
    #
    bragg_angle = diffraction_setup.angleBragg(energy)
    bragg_angle_corrected = diffraction_setup.angleBraggCorrected(energy)
    half_darwin_s = diffraction_setup.darwinHalfwidthS(energy)

    print("Bragg angle for E=%f eV is %f deg" % (energy, bragg_angle_corrected * 180.0 / numpy.pi))
    print("Bragg angle for E=%f eV is %f deg" % (energy, bragg_angle * 180.0 / numpy.pi))
    print("Darwin width S for E=%f eV is %f deg" % (energy, 2 * half_darwin_s))

    # initialize arrays for storing outputs
    deviations = numpy.zeros(angle_deviation_points, dtype=float)
    intensityS = numpy.zeros(angle_deviation_points, dtype=float)
    intensityP = numpy.zeros(angle_deviation_points, dtype=float)
    phaseS = numpy.zeros(angle_deviation_points, dtype=float)
    phaseP = numpy.zeros(angle_deviation_points, dtype=float)

    for ia in range(angle_deviation_points):
        deviation = angle_deviation_min + ia * angle_step
        angle = deviation + bragg_angle_corrected

        # calculate the components of the unitary vector of the incident photon scan
        # Note that diffraction plane is YZ
        yy = numpy.cos(angle)
        zz = - numpy.abs(numpy.sin(angle))
        photon = Photon(energy_in_ev=energy, direction_vector=Vector(0.0, yy, zz))
        photon = ComplexAmplitudePhoton(energy_in_ev=energy,
                                        direction_vector=diffraction_setup.vectorIncomingPhotonDirection(
                                            energy,deviation,
                                            angle_center_flag=2,# 0=Absolute angle, 1=Theta Bragg Corrected, 2=Theta Bragg
                                        ),
                                        Esigma=1,
                                        Epi=1)

        # perform the calculation
        coeffs = Diffraction.calculateDiffractedComplexAmplitudes(diffraction_setup, photon,
                                                                  calculation_method=calculation_method,
                                                                  calculation_strategy_flag=calculation_strategy_flag)
        R_S = complex(coeffs['S'])
        R_P = complex(coeffs['P'])
        # store results
        deviations[ia] = deviation / half_darwin_s
        intensityS[ia] = numpy.abs(R_S) ** 2
        intensityP[ia] = numpy.abs(R_P) ** 2
        phaseS[ia] = numpy.angle(R_S)
        phaseP[ia] = numpy.angle(R_P)

    # plot results

    fig1, ax1 = plot(
        deviations, intensityS,
        deviations, intensityP,
        xtitle="deviation angle [normalised]",
        ytitle="Reflectivity",
        title="Laue Transmission Fig 9; Si220",
        legend=["Sigma-polarization", "Pi-polarization"],
        xrange=[-5, 5], yrange=[0, 1],
        grid=1,
        show=0,
    )

    # phaseSu = numpy.unwrap(phaseS)
    # phasePu = numpy.unwrap(phaseP)

    angle_deg = numpy.degrees(numpy.unwrap(phaseP - phaseS))
    # Force result back into [-180, 180] range
    final_angle_deg = (angle_deg + 180) % 360 - 180
    fig2, ax2 = plot(
        deviations, final_angle_deg,
        xtitle="deviation angle [normalised]",
        ytitle="Phase difference [deg]",
        title="Laue Transmission Fig 9; PhaseP - PhaseS; Si220",
        xrange=[-5, 5], yrange=[-180, 180],
        grid=1,
        show=0,
    )

    return fig1, fig2, ax1, ax2

#
# main
#
if __name__ == "__main__":

    calculation_method = 1         # 0=Zachariasen, 1=Guigay
    calculation_strategy_flag = 0  # 0=mpmath 1=numpy 2=numpy-truncated

    fig1_1, fig2_1, ax1_1, ax2_1 = calculate_fig3(calculation_method=0, calculation_strategy_flag=0)
    fig1_2, fig2_2, ax1_2, ax2_2 = calculate_fig4(calculation_method=0, calculation_strategy_flag=0)
    fig1_3, fig2_3, ax1_3, ax2_3 = calculate_fig6(calculation_method=0, calculation_strategy_flag=0)
    fig1_4, fig2_4, ax1_4, ax2_4 = calculate_fig9(calculation_method=0, calculation_strategy_flag=0)

    individual_plots = 0

    if individual_plots:
        pass
    else:
        matplotlib_combine([fig1_1, fig1_2, fig1_3, fig1_4], [ax1_1, ax1_2, ax1_3, ax1_4], show=0)
        matplotlib_combine([fig2_1, fig2_2, fig2_3, fig2_4], [ax2_1, ax2_2, ax2_3, ax2_4], show=0)
    plot_show()

