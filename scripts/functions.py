import numpy as np
from copy import copy

pi = np.pi


def jacobian(model, q, nlink, rbdl, delta=0.0001):

    ndof = model.dof_count
    # Alocacion de memoria para el Jacobiano
    J = np.zeros((3, ndof))
    # Posicion inicial (con q)
    x = rbdl.CalcBodyToBaseCoordinates(model, q, nlink, np.zeros(3))
    # Iteracion para la derivada de cada columna
    for i in xrange(ndof):
        dq = copy(q);
        # Incremento en el angulo i-esimo
        dq[i] = dq[i]+delta
        # Transformacion homogenea luego del incremento (q+dq)
        dx = rbdl.CalcBodyToBaseCoordinates(model, dq, nlink, np.zeros(3))
        # Aproximacion al Jacobiano de posicion usando diferencias finitas
        J[0:3,i] = (dx-x)/delta
    return J
