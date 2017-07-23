"""
A jupyter notebook to play with quintic polynomial trajectory generation
"""
#%%
import numpy as np
from numpy.linalg import inv

#%%


def jmt(start, end, time):
    """
    The function receives a start and an end configuration
    and the time required to achieve the end configuration,
    and return a 5th order trajectory polynomial that satisfies
    start and end conditions
    start = (s[0] , u[0],  acc[0])
    end =  (s[t=time] , u[t=time],  acc[t=time])
    """
    # the first 3 coefficients are from s(t=0) derivatives
    c_1, c_2, c_3 = start[0], start[1], start[2] / 2

    matrix_a = np.array([
        [time**3,           time**4,       time**5],
        [3 * time**2,   4 * time**3,   5 * time**4],
        [6 * time,     12 * time**2,  20 * time**3]
    ])

    vector_b = np.array(
        [end[0] - (start[0] + start[1] * time + 0.5 * start[2] * time**2),
         end[1] - (start[1] + start[2] * time),
         end[2] - start[2]])

    x_r = inv(matrix_a).dot(vector_b)  # solve AX = B system

    c_4, c_5, c_6 = x_r[0], x_r[1], x_r[2]

    return [c_1, c_2, c_3, c_4, c_5, c_6]


#%%
def close_enough(poly, target_poly, eps=0.01):
    """
    Returns True if the coefficients of the two polynomials have the following
    a. they are the same number
    b. they are close enough pointwise.
    """
    if type(target_poly) != type(poly):
        target_poly = list(target_poly)
    if len(poly) != len(target_poly):
        raise Exception("Incorrect number of terms")
    for term, target_term in zip(poly, target_poly):
        if abs(term - target_term) > eps:
            print("At least one of the terms is not close enough")
            return False
    return True

#%%


def passed(test: bool) -> str:
    """Convert True or False , to passed or failed"""
    return 'passed' if test else 'failed'


#%%
#  test case : 1
START = [0, 10, 0]
END = [10, 10, 0]
TIME = 1
TRAJECTORY = [0.0, 10.0, 0.0, 0.0, -0.0, 0.0]
SOLUTION = jmt(START, END, TIME)
print('Test 1 :', passed(close_enough(SOLUTION, TRAJECTORY)))


#%%
#  test case : 2
START = [0, 10, 0]
END = [20, 15, 20]
TIME = 2
TRAJECTORY = [0.0, 10.0, 0.0, 0.0, -0.625, 0.3125]
SOLUTION = jmt(START, END, TIME)
print('Test 3 :', passed(close_enough(SOLUTION, TRAJECTORY)))


#%%
#  test case : 3
START = [5, 10, 2]
END = [-30, -20, -4]
TIME = 5
TRAJECTORY = [5.0, 10.0, 1.0, -3.0, 0.64, -0.0432]
SOLUTION = jmt(START, END, TIME)
print('Test 3 :', passed(close_enough(SOLUTION, TRAJECTORY)))
