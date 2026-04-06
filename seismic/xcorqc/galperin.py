import numpy as np

# Galperin angle
alpha_deg = 35.2643897
alpha_rad = np.radians(alpha_deg)

# forward: UVW -> ZNE and the inverse: ZNE->UVW
UVW2ZNE = np.array([
    [np.cos(alpha_rad), np.cos(alpha_rad), np.cos(alpha_rad)],
    [np.sin(alpha_rad), -0.5 * np.sin(alpha_rad), -0.5 * np.sin(alpha_rad)],
    [0, np.sqrt(3) / 2 * np.sin(alpha_rad), -np.sqrt(3) / 2 * np.sin(alpha_rad)]
])
ZNE2UVW = np.linalg.inv(UVW2ZNE)

def galperin_transform(A, B, C, direction=1):
    """
    Transform between Galperin (U,V,W) and standard (Z,N,E).
    Parameters
    ----------
    A, B, C : array-like
        Input components:
        - if direction=1  → (U, V, W)
        - if direction=-1 → (Z, N, E)

    direction : int
        +1 : Galperin → ZNE (forward)
        -1 : ZNE → Galperin (inverse)
    Returns
    -------
    tuple
        Transformed components:
        - (Z, N, E) if direction=1
        - (U, V, W) if direction=-1
    """

    if direction not in (1, -1):
        raise ValueError("direction must be +1 (forward) or -1 (inverse)")
    # end if

    # Choose matrix based on direction
    if direction == 1:
        M = UVW2ZNE
    else:
        M = ZNE2UVW
    # end if

    # Stack inputs
    X = np.array([A, B, C])

    # Transform
    Y = M @ X

    return tuple(Y)
# end func

import numpy as np

def galperin_transform(E, N, Z, direction=1):
    """
    Transform between Galperin (U,V,W) and standard (Z,N,E).

    Convention:
        U ↔ E
        V ↔ N
        W ↔ Z

    Parameters
    ----------
    E, N, Z : array-like
        Input components:
        - if direction=1  → (U, V, W) mapped as (E, N, Z)
        - if direction=-1 → (Z, N, E)

    direction : int
        +1 : Galperin → ZNE
        -1 : ZNE → Galperin

    Returns
    -------
    tuple
        - (Z, N, E) if direction=1
        - (E, N, Z) if direction=-1  (i.e., U, V, W)
    """

    if direction not in (1, -1):
        raise ValueError("direction must be +1 or -1")
    # end if

    if direction == 1:
        # Input is (U,V,W) but passed as (E,N,Z)
        ts = np.array([E, N, Z])
        Z_out, N_out, E_out = UVW2ZNE @ ts
        return Z_out, N_out, E_out
    else:
        # Input is (Z,N,E)
        ts = np.array([E, N, Z])  # careful: reordering
        U, V, W = ZNE2UVW @ ts

        # Map back: U→E, V→N, W→Z
        return U, V, W
    # end if
# end func