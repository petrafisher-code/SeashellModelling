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

def main(D, A, alpha, beta, mu, omega, phi, a, b, L, P, W1, W2, N, image):
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
    image: string
        Name of the image to use for the surface colour (JPG).
    """
    # Define the ranges of s and theta, as seen in (Picado, 2009)
    theta_array = np.linspace(-4 * np.pi, 4 * np.pi, 200)
    s_array = np.linspace(np.radians(-270), np.radians(90), 50)

    # Calculate nodule array if appropriate
    nodule_array = np.zeros(len(theta_array))
    if L != 0 and N != 0 and W1 != 0 and W2 != 0:
        for t_index, t in enumerate(theta_array):
            nodule_array[t_index] = (2*np.pi/N)*((N*t)/(2*np.pi) - np.floor((N*t)/(2*np.pi)))

    # Define x, y, and z arrays
    x_array = np.zeros((len(s_array), len(theta_array)))
    y_array = np.zeros((len(s_array), len(theta_array)))
    z_array = np.zeros((len(s_array), len(theta_array)))

    # Iterate through s and theta array values, calculating x, y, and z.
    for s_index, s in enumerate(s_array):
        for t_index, t in enumerate(theta_array):
            if L != 0 and N != 0 and W1 != 0 and W2 != 0:
                k = L*np.exp(-(2*(s-P)/W1)**2 - (2*nodule_array[t_index]/W2)**2)
            else:
                k = 0
            h = 1/(np.sqrt((np.cos(s) / a) ** 2 + (np.sin(s) / b) ** 2)) + k
            e = np.exp(t / np.tan(alpha))
            x_array[s_index][t_index] = D*(A*np.sin(beta)*np.cos(t+omega)
                                            + h*(np.cos(s+phi)*np.cos(t+omega)
                                                - np.sin(mu)*np.sin(s+phi)*np.sin(t+omega)))*e
            y_array[s_index][t_index] = (A*np.sin(beta)*np.sin(t+omega)
                                            + h*(np.cos(s+phi)*np.sin(t+omega)
                                                +np.sin(mu)*np.sin(s+phi)*np.cos(t+omega)))*e
            z_array[s_index][t_index] = (-A*np.cos(beta) + h*np.cos(mu)*np.sin(s+phi))*e

    # Calculate the surface overlays from the specified JPG image
    img = plt.imread(f'surface_images/{image}.jpg')
    img = img / 255.0
    s_normalized = (s_array - s_array.min()) / (s_array.max() - s_array.min())
    t_normalized = (theta_array - theta_array.min()) / (theta_array.max() - theta_array.min())
    s_idx = (s_normalized * (img.shape[0] - 1)).astype(int)
    t_idx = (t_normalized * (img.shape[1] - 1)).astype(int)
    facecolors = img[s_idx][:, t_idx]

    # Create figure (plotting two subplots)
    fig = plt.figure(figsize=(12, 6))

    # Black wireframe image
    ax1 = fig.add_subplot(122, projection='3d')
    ax1.plot_wireframe(x_array, y_array, z_array, color='black')
    ax1.view_init(elev=30, azim=-60)
    ax1.axis('off')

    # Colored surface plot with image overlay
    ax2 = fig.add_subplot(121, projection='3d')
    ax2.plot_surface(x_array, y_array, z_array, rstride=1, cstride=1, facecolors=facecolors,
                     edgecolor='none')
    ax2.view_init(elev=30, azim=-60)

    plt.show()

if __name__ == "__main__":
    # Hamish's example:
    main(D=1, A=100, alpha=np.radians(95), beta=np.radians(25),
         mu=np.radians(0), omega=np.radians(0), phi=np.radians(180),
          a=10, b=20, L=0, P=0, W1=0, W2=0, N=0, image="pink_and_purple")

    # PAUA
    main(D=1, A=9, alpha=np.radians(45), beta=np.radians(50),
         mu=np.radians(1), omega=np.radians(-30), phi=np.radians(-10),
         a=8, b=10, L=0, P=0, W1=0, W2=0, N=0, image="paua")

    # SHELL
    main(D=1, A=80, alpha=np.radians(86), beta=np.radians(10),
         mu=np.radians(5), omega=np.radians(0), phi=np.radians(-45),
          a=20, b=20, L=10, P=np.radians(40), W1=180, W2=0.4, N=10, image="shell")
