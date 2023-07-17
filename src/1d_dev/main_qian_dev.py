# usr/bin/env python = 3

import numpy as np

class ADSolver():
    """
    Time is in hours.
    """
    def __init__(self) -> None:
        pass
    
    def cal_tran(self, ki, kr, vel, conc, dt, dx):
        """
        The function calculates the advection of a concentration field using given parameters.
        
        :param ki: ki is a list of values representing the reaction rate constants at different
        positions in a system
        :param kr: The parameter "kr" represents the reaction rate coefficient
        :param vel: The parameter "vel" represents the velocity of the fluid or medium in which the
        transport process is occurring
        :param conc: The parameter "conc" represents the concentration of a substance at different
        points in space. It is a list or array containing the concentration values at each point
        :param dt: dt is the time step size, which represents the amount of time between each iteration
        of the calculation. It determines how much time is advanced in each iteration
        :param dx: The parameter "dx" represents the spatial step size or the distance between two
        adjacent grid points in the spatial domain. It is used to calculate the spatial derivative of
        the concentration
        :return: the value of `adv_c`.
        """
        ad1 = ki[2:] * vel[2:] * conc[2:] - ki[1:-1] * vel[1:-1] * conc[1:-1]
        ad2 = kr[1:-1] * vel[1:-1] * conc[1:-1] - kr[0:-2] * vel[0:-2] * conc[0:-2]
        adv_c = (-ad1 + ad2) * (dt / dx)
        return adv_c

    def cal_diff(self, diff_coff, conc, dt, dx):
        diff_c = (diff_coff[2:] * conc[2:] - 2 * diff_coff[1:-1] * conc[1:-1] + diff_coff[0:-2] * conc[0:-2]) * (dt / dx ** 2)
        return diff_c

    def cal_reaction():
        pass

    def wqpsl():
        """
        Point source load.
        """
        pass

    def wqwet():
        """
        Loadings during wet days
        """
        #TODO:maybe this is not implemented.
        pass

    def bdc():
        """
        Boundary conditions
        """
        pass

def read_sheet(workbook, name):
    """
    The function `read_sheet` reads data from a specified sheet in an Excel workbook and returns it as a
    numpy array.
    
    :param workbook: The workbook parameter is the Excel workbook object that contains the sheet you
    want to read. This object represents the entire Excel file and allows you to access and manipulate
    its sheets and data
    :param name: The parameter "name" is the name of the sheet in the workbook that you want to read
    :return: a numpy array containing the values from the specified sheet in the given workbook.
    """
    # Reading water depth sheet
    sheet = workbook.sheet_by_name(name)
    # Number of written Rows in sheet
    r = sheet.nrows
    # Number of written Columns in sheet
    c = sheet.ncols
    answ = np.zeros([r - 1, c])
    # Reading each cell in excel sheet 'BC'
    for i in range(1, r):
        for j in range(c):
            answ[i - 1, j] = float(sheet.cell_value(i, j))
    return answ
    