"""
Generates a 3D Seashell figure.

This script generates a 3D plot of a seashell, using the mathematical
model found in Jorge Picado's 2009 paper: 
SEASHELLS: THE PLAINNESS AND BEAUTY OF THEIR MATHEMATICAL DESCRIPTION
https://www.mat.uc.pt/~picado/conchas/eng/article.pdf#toolbar=0

Author: Petra Fisher
Date: August 2024
"""

import numpy as np
import matplotlib.pyplot as plt

def main(D, A, alpha, beta, mu, omega, phi, a, b, L, P, W1, W2, N):
    """
    Generates and plots the 3D sea shell, using the model's
    equations and matplotlib. Please see the parameters section
    below for descriptions of the 14 parameters.
    
    Parameters:
    -----------
    D : int
        Direction of coiling (1 for dextral, -1 for sinistral).
    A : float
        Size of the spiral aperture.
    alpha : float
        Equiangular angle of the spiral.
    beta : float
        Enlarging angle of the spiral.
    mu : float
        Rotation around the horizontal axis.
    omega : float
        Rotation around the OZ axis and horizontal axis.
    phi : float
        Rotation around a vector orthogonal to the ellipse plane.
    a : float
        Semimajor axis of the generating curve (ellipse).
    b : float
        Semiminor axis of the generating curve (ellipse).
    L : float
        Height of each nodule.
    P : float
        Position angle of the nodule on the generating curve.
    W1 : float
        Length of each nodule along the generating curve.
    W2 : float
        Length of each nodule along the helico-spiral.
    N : int
        Number of nodules per complete revolution.
    """
    theta_array = np.linspace(0, 2*np.pi, 100)
    s_array = np.linspace(0, 2*np.pi, 100)

    nodule_array = np.zeros((N, len(theta_array)))
    if L != 0 and N != 0 and W1 != 0 and W2 != 0:
        for n in range(N):
            for t_index, t in enumerate(theta_array):
                nodule_array[n][t_index] = (2*np.pi/(n+1))*((n+1)*t/(2*np.pi))\
                                            - int((n+1)*t/(2* np.pi))

    x_array = np.zeros((len(s_array), len(theta_array)))
    y_array = np.zeros((len(s_array), len(theta_array)))
    z_array = np.zeros((len(s_array), len(theta_array)))

    for s_index, s in enumerate(s_array):
        for t_index, t in enumerate(theta_array):
            if L != 0 and N != 0 and W1 != 0 and W2 != 0:
                k = L*np.exp(-((2*(s-P)/W1)**2) - ((2*nodule_array[0][t_index]/W2)**2))
            else:
                k = 0
            h = np.sqrt(1 / ((np.cos(s) / a) ** 2 + (np.sin(s) / b) ** 2)) + k
            x_array[s_index][t_index] = D*(A*np.sin(beta)*np.cos(t)
                                            + h*(np.cos(s+phi)*np.sin(t+omega)
                                                + np.sin(mu)*np.sin(s+phi)*np.cos(t+omega)
                                            ))*np.exp(t/np.tan(alpha))
            y_array[s_index][t_index] = (A*np.sin(beta)*np.sin(t)
                                            + h*(np.cos(s+phi)*np.sin(t+omega)
                                                 + np.sin(mu)*np.sin(s+phi)*np.cos(t+omega)
                                            ))* np.exp(t/np.tan(alpha))
            z_array[s_index][t_index] = (-A*np.cos(beta) + h*np.cos(mu)*np.sin(s+phi)) \
                                        * np.exp(t/np.tan(alpha))

    plt.figure()
    ax = plt.axes(projection="3d")
    ax.plot_surface(x_array, y_array, z_array, cmap='viridis', edgecolor='none')
    ax.set_title("Seashell")
    plt.show()


if __name__ == "__main__":
    main(D=1, A=100, alpha=95*np.pi/180, beta=25*np.pi/180,
         mu=0*np.pi/180, omega=0*np.pi/180, phi=360*np.pi/180,
          a=10, b=20, L=0, P=0, W1=0, W2=0, N=0)
