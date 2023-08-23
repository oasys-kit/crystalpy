import numpy
from crystalpy.util.calc_xcrystal import calc_xcrystal_angular_scan, calc_xcrystal_energy_scan
from crystalpy.util.calc_xcrystal import calc_xcrystal_alphazachariasen_scan, calc_xcrystal_double_scan

if __name__ == "__main__":
    from srxraylib.plot.gol import set_qt
    set_qt()

    calculation_method = 1
    is_thick = 0
    use_transfer_matrix = 0

    if True:
        calc_xcrystal_angular_scan(material_constants_library_flag=0, do_plot=True, calculation_method=calculation_method,
                                                                           is_thick=is_thick,
                                                                           use_transfer_matrix=use_transfer_matrix)

    if True:
        calc_xcrystal_angular_scan(material_constants_library_flag=0, geometry_type_index=1, thickness=10e-6,
                                   asymmetry_angle=numpy.radians(90), do_plot=True, calculation_method=calculation_method,
                                                                           is_thick=is_thick,
                                                                           use_transfer_matrix=use_transfer_matrix)

        calc_xcrystal_energy_scan(material_constants_library_flag=0, do_plot=True, calculation_method=calculation_method,
                                                                           is_thick=is_thick,
                                                                           use_transfer_matrix=use_transfer_matrix)

        calc_xcrystal_alphazachariasen_scan(do_plot=1, calculation_method=calculation_method,
                                                                           is_thick=is_thick,
                                                                           use_transfer_matrix=use_transfer_matrix)

    if True:
        calc_xcrystal_double_scan(        material_constants_library_flag=0,
            crystal_name="Si",
            thickness=1e-2,
            miller_h=1,
            miller_k=1,
            miller_l=1,
            asymmetry_angle=0.0,
            energy_min=8000,
            energy_max=8010,
            energy_points=1,
            angle_deviation_min=-100e-6,
            angle_deviation_max=100e-6,
            angle_deviation_points=200,
            angle_center_flag=2, # 0=Absolute angle, 1=Theta Bragg Corrected, 2=Theta Bragg
            calculation_method=calculation_method,
            is_thick=is_thick,
            use_transfer_matrix=use_transfer_matrix,
            geometry_type_index=0,
            do_plot=1,)

        calc_xcrystal_double_scan(        material_constants_library_flag=0,
            crystal_name="Si",
            thickness=1e-2,
            miller_h=1,
            miller_k=1,
            miller_l=1,
            asymmetry_angle=0.0,
            energy_min=7990,
            energy_max=8010,
            energy_points=100,
            angle_deviation_min=0,
            angle_deviation_max=100e-6,
            angle_deviation_points=1,
            angle_center_flag=2, # 0=Absolute angle, 1=Theta Bragg Corrected, 2=Theta Bragg
            calculation_method=calculation_method,
            is_thick=is_thick,
            use_transfer_matrix=use_transfer_matrix,
            geometry_type_index=0,
            do_plot=1,)

    if True:
        calc_xcrystal_double_scan(        material_constants_library_flag=0,
            crystal_name="Si",
            thickness=0.010,
            miller_h=1,
            miller_k=1,
            miller_l=1,
            asymmetry_angle=0.0,
            energy_min=7990,
            energy_max=8010,
            energy_points=150,
            angle_deviation_min=-100e-6,
            angle_deviation_max=100e-6,
            angle_deviation_points=150,
            angle_center_flag=2, # 0=Absolute angle, 1=Theta Bragg Corrected, 2=Theta Bragg
            calculation_method=calculation_method,
            is_thick=is_thick,
            use_transfer_matrix=use_transfer_matrix,
            geometry_type_index=0,
            do_plot=1,)