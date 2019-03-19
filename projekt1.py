import numpy as np
import matplotlib.pyplot as plt
import sympy
import math
from scipy.constants import g

def ParseLine(line):
    """ Returns cannon position, target position, velocity vector, wind vector,
    timestamps and string containing function that is defining the terrain.
    """
    line = line[:-1]
    line = line.split("; ")
    
    pos_0_str = line[0][1:-1]
    pos_0 = pos_0_str.split(", ")
    for i in range(len(pos_0)):
        pos_0[i] = float(pos_0[i])
    line.pop(0)

    target_str = line[0][1:-1]
    target = target_str.split(", ")
    for i in range(len(target)):
        target[i] = float(target[i])
    line.pop(0)
    
    v_str = line[0][1:-1]
    v = v_str.split(", ")
    for i in range(len(v)):
        v[i] = float(v[i])
    line.pop(0)

    w_str = line[0][1:-1]
    w = w_str.split(", ")
    for i in range(len(w)):
        w[i] = float(w[i])
    line.pop(0)

    t = []
    while len(line) != 1:
        t.append(float(line[0]))
        line.pop(0)
    
    fun_str = line.pop(0)

    return pos_0, target, v, w, t, fun_str

def StrToCoefficients(fun_str):
    """ Returns the array of coefficients of a function saved in a str.
        Coefficients are in ascending order (from lowest to highest power).
    """
    arr = fun_str.split(' ')

    max_power = int(arr[0][-1])
    power = max_power
    coefficients = [1 for i in range(max_power + 1)]

    for i in range(len(arr)):
        if 'x' in arr[i] and arr[i][-1] == 'x':
            arr[i] += '^1'

    i = 0
    while i < len(arr):
        if arr[i] != '+' and arr[i] != '-':
            if power == int(arr[i][-1]) and 'x' in arr[i]:
                coefficients[max_power - power] *= float(arr[i][:-3])
                power -= 1
                i += 1
            elif power == 0:
                coefficients[max_power] = float(arr[i])
                i += 1
            else:
                coefficients[max_power - power] = 0
                power -= 1
        elif arr[i] == '-':
            coefficients[max_power - power] *= -1
            i += 1
        else:
            i += 1

    coefficients.reverse()
    return(coefficients)

def PolyCoefficients(x, coeffs):
    """ Returns a polynomial for ``x`` values for the ``coeffs`` provided and
    a equation that sympy is capable of solving

    The coefficients must be in ascending order (``x**0`` to ``x**o``).
    """
    y = 0
    z = sympy.symbols('z', real=True)
    equation = 0
    for i in range(len(coeffs)):
        y += coeffs[i]*x**i
        equation += coeffs[i]*z**i

    return y, equation

def GetEquationOfMotion(v, w, pos0):
    """ Returns a equation of motion for the bullet that sympy
    is capable of solving.
    v - velocity vector
    w - wind vector
    pos0 - starting positiion
    """
    vel_x = v[0] + w[0]
    vel_y = v[1] + w[1]

    z = sympy.symbols('z', real=True)
    equation = (vel_y/vel_x) * (z - pos0[0]) - ( g/ (2*(vel_x**2)) ) *\
    (z - pos0[0])**2 + pos0[1]

    return equation

def GetYOfMotion(v, w, x, pos0):
    """ Returns y values for the bullet motion
    v - velocity vector
    w - wind vector
    pos0 - starting positiion
    """
    vel_x = v[0] + w[0]
    vel_y = v[1] + w[1]
    position = [x - pos0[0], pos0[1]]

    y = (vel_y/vel_x) * position[0] - ( g/ (2*(vel_x**2)) ) *\
    position[0]**2 + position[1]

    return y

def GetImpactPoint(terrain_eq, motion_eq, cannon_position):
    """ Calculates and returns the x coordinate of impact point
    (calculates points of intersections and ignores the cannon location.
    Returns the first point of impact with terrain)
    """
    z = sympy.symbols("z", real=True)
    equation = terrain_eq - motion_eq

    coordinates = sympy.solve(equation, z)

    for i in range(len(coordinates) - 1):
        if coordinates[i] == cannon_position:
            coordinates.pop(i)

    coordinates.sort()

    return float(coordinates[0])

def RunCalculations(pos_0, target, v, w, t, fun_str, out_name):
    """ Calculates everything that's necessary to plot all the plots,
    returns coordinates of the impact point
    """
    # create x values, calculate y values and equations for further calculations
    x_terrain = np.linspace(pos_0[0]-0.03, target[0]+0.03, 10000)
    y_terrain, terrain_equation = PolyCoefficients(x_terrain, 
                                                   StrToCoefficients(fun_str))

    motion_equation = GetEquationOfMotion(v, w, pos_0)
    impact_point_x = GetImpactPoint(terrain_equation, 
                                    motion_equation, pos_0[0])
    x_shot = np.linspace(pos_0[0], impact_point_x, 10000)
    y_shot = GetYOfMotion(v, w, x_shot, pos_0)

    # mark bullet positions after time
    z = sympy.symbols("z", real=True)
    mark = []
    for i in t:
        mark.append(i)
        value = motion_equation.subs(z, i)
        if i > impact_point_x:
            mark.append(motion_equation.subs(z, impact_point_x))
            mark[-2] = impact_point_x
        else:
            mark.append(value)

    y_times = [mark[i] for i in range(len(mark)) if i%2]
    x_times = [mark[i] for i in range(len(mark)) if not i%2]
    plt.plot(x_times, y_times, "ro", label='Timestamps', zorder=10)
    
    # mark positions of cannon, target and hit point
    plt.plot(pos_0[0], pos_0[1], "m*", label="Cannon", zorder=10)
    plt.plot(target[0], target[1], "k*", label="Target", zorder=10)
    plt.plot(impact_point_x, motion_equation.subs(z, impact_point_x), 'ko', 
             label="Impact point", markersize=12, zorder=5)

    # plot trajectory of the bullet
    plt.plot(x_shot, y_shot, 'r', label='Trajectory')

    # plot the terrain
    plt.plot(x_terrain, y_terrain, 'g')
    plt.fill_between(x_terrain, y_terrain, color='green', label='terrain')

    # restrict the area included
    plt.ylim(bottom=0)
    plt.xlim(pos_0[0]-0.03, target[0]+0.03)

    # show the legend
    ax = plt.gca()
    ax.legend()

    # save the masterpiece
    plt.savefig(out_name)

    # show the masterpiece
    plt.show()
    
    return [impact_point_x, motion_equation.subs(z, impact_point_x)]

def CalculateMaxHeight(v_0, w):
    """ Calculate and return maximum reached height
    """
    max_height = (v_0[1] + w[1])**2 / (2*g)
    return max_height

def CalculateYVelocityInTime(v_y0, timestamp):
    """ Calculate and return Y velocity at a given moment in time
    """
    v = v_y0 - timestamp*g
    return v

def GetVelocityInTime(v_0, w, timestamp):
    """ Calculate and return velocity vector at a given moment in time
    """
    v = [v_0[0] + w[0], CalculateYVelocityInTime(v_0[1], timestamp)]
    return v;

def DetermineHit(impact_point, target):
    """ Determine wheter the target has been hit or not.
    Returns 1 if yes, 0 if not
    """
    distance = math.hypot(target[0] - impact_point[0], 
                          target[1] - impact_point[1])
    
    if distance <= 0.05:
        return 1
    else:
        return 0

def ParseOutput(v_0, w, t, impact_point, target):
    """ Prepare and return the perfect line for output file
    """
    output_str = "(" + str(impact_point[0]) + ", " +\
    str(impact_point[1]) + "); "
    
    max_height = CalculateMaxHeight(v_0, w)
    
    output_str += str(max_height) + "; "
    
    for i in t:
        v_timestamp = GetVelocityInTime(v_0, w, i)
        output_str += "[" + str(v_timestamp[0]) + ", " +\
        str(v_timestamp[1]) + "]; "
    
    hit = DetermineHit(impact_point, target)
    
    output_str += str(hit) + "\n"
    
    return output_str

if __name__ == "__main__":
    with open("input.txt", "r") as file_input:
        with open("output.txt", "a+") as file_output:
            i = 1
            for line in file_input:
                plot_output_file_name = str(i) + ".png"
    
                pos_0, target, v_0, w, t, fun_str = ParseLine(line)
    
                impact_point = RunCalculations(pos_0, target, v_0, w, t, fun_str, 
                                               plot_output_file_name)
                i += 1
            
                output_line = ParseOutput(v_0, w, t, impact_point, target)
                file_output.write(output_line)
    