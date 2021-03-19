"""
    read polynom from csv-file (exported from Detectors.mdb) by DRGen-gw and plot graph
    parameters: <path_to_CalcResults.csv> [<point_num>]
    without point_num outputs list of energies, for which response was calculated
"""

import sys
import csv
import matplotlib.pyplot as plt
import numpy as np

MC2 = 511.008
MC2x2 = 1022.016

def poly(values, coeffs):
    """
        calculates polynom y from x (list or array) and coefficients
    """
    return np.array(
        [np.sum(np.array([coeffs[i] * x**i for i in range(len(coeffs))])) 
            for x in values
        ]
    )

def compton_backscat(energy):
    """
        calculate energy for Compton backscatering (wait energy in keV)
    """
    return MC2 / (2 + MC2 / energy)

def get_energy_bounds(energy):
    """
        energy bounds for polynom zones. energy in keV
    """
    e_back_scat = compton_backscat(energy)
    e_compton_edge = energy - e_back_scat
    if energy <= MC2x2:
        res = np.array([0, e_back_scat, e_compton_edge, energy])
    else:
        res = np.array([0, e_back_scat, e_compton_edge, 
                energy - MC2x2, energy - 4.0/3 * MC2, energy - MC2, energy])
    res.sort()
    return res
    
def read_det_energies(csv_file_name):
    """
        read energy list for what response was calculated
    """
    header = {}
    energies = []
    with open(csv_file_name) as f:
        reader = csv.reader(f)
        is_first = True
        for row in reader:
            # read header
            if is_first:
                for i in range(len(row)):
                    header[row[i]] = i
                is_first = False
                if "energy" not in header or "point" not in header:
                    print("There is no point or energy column in file", csv_file_name)
                    sys.exit()
                continue
            # read energies
            energy = float(row[header["energy"]])
            energies.append(energy)
    return energies
    
def read_a_polys(csv_file_name):
    """
        read all polynom coefficients for all energies
    """
    header = {}
    a_poly = []
    with open(csv_file_name) as f:
        reader = csv.reader(f)
        is_first = True
        for row in reader:
            # read header
            if is_first:
                for i in range(len(row)):
                    header[row[i]] = i
                is_first = False
                continue
            # E -> a-poly, E -> b-poly
            energy = float(row[header["energy"]])
            i0 = header["a00"]
            a = np.zeros(6*9)
            for i in range(6*9):
                a[i] = float(row[i0+i])
            a = a.reshape((6,9))
            a_poly.append((energy, a))
    return a_poly
    
def test_get_energy_bounds():
    # small energy
    e = 100
    res = get_energy_bounds(e)
    e_bs = MC2 / (2 + MC2 / e)
    expected = np.array([0, e - e_bs, e_bs, e])
    assert (res == expected).all()
    # big energy
    e = 2000
    res = get_energy_bounds(e)
    e_bs = MC2 / (2 + MC2 / e)
    expected = np.array([0, e_bs, e -MC2x2, e - 4/3*MC2, e - MC2, e - e_bs, e])
    assert (res == expected).all()
    

if __name__ == "__main__":
    test_get_energy_bounds()
    if len(sys.argv) <= 1:
        print("Plot detector response function\nplot_calc_result_poly: <path_to_CalcResults.csv> [<point_num>]")
        sys.exit()

    csv_file_name = sys.argv[1]
    point_num = int(sys.argv[2]) if len(sys.argv) > 2 else -1
    
    if point_num < 1:
        energies = read_det_energies(csv_file_name)
        print("Found energies:")
        for i,e in enumerate(energies):
            print(i+1, '\t', e)
        sys.exit()
    
    # read csv, header -- column names, poly -- list of (energy, np.array(poly_coeffs))
    a_poly = read_a_polys(csv_file_name)
    point_num -= 1
    if point_num >= len(a_poly):
        print("Point index is out of range. Max index can be", len(a_poly)+1)
        sys.exit()
        
    energy, polys = a_poly[point_num]
    NUM_POINTS = 20*6
    if energy < 1.022016:
        polys = polys[3:]
        NUM_POINTS = 20*3
    bounds = get_energy_bounds(energy*1000)
    assert len(polys) + 1 == len(bounds)
    
    res = np.array([])
    egrid = np.array([])
    for i, poly_zone in enumerate(polys):
        num_points = round((bounds[i+1] - bounds[i])/(energy*1000) * NUM_POINTS)
        assert num_points <= NUM_POINTS
        grid = np.linspace(0,1,num_points)
        egrid = np.concatenate((egrid, (bounds[i+1] - bounds[i]) * grid + bounds[i]))
        p = poly(grid, poly_zone)
        res = np.concatenate((res, p))
    
    # plot poly
    plt.plot(e_grid, res, 'r')
    plt.title('Compton for energy {} keV'.format(energy*1000))
    plt.xlabel('energy, MeV') 
    plt.ylabel('a')
    plt.show()
            
    
    
