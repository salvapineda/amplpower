import numpy as np


def find_min_max(a, b, c, xmin, xmax, ymin, ymax, zmin, zmax):
    """
    Find min and max of f(x,y,z) = a x² + b x y cos(z) + c x y sin(z)
    within given bounds.

    Args:
        a, b, c: Coefficients (non-zero)
        xmin, xmax: Bounds for x (xmin > 0)
        ymin, ymax: Bounds for y (ymin > 0)
        zmin, zmax: Bounds for z (within -π/2 to π/2)

    Returns:
        (min_val, max_val)
    """

    # Define the function
    def f(x, y, z):
        return a * x**2 + b * x * y * np.cos(z) + c * x * y * np.sin(z)

    # Generate all boundary combinations
    boundaries = [(x, y, z) for x in (xmin, xmax) for y in (ymin, ymax) for z in (zmin, zmax)]

    # Add special z candidates (from KKT conditions)
    special_z = []
    if c != 0:
        z1 = np.arctan2(-b, c)  # From ∂f/∂y = 0
        if zmin <= z1 <= zmax:
            special_z.append(z1)
    if b != 0:
        z2 = np.arctan2(c, b)  # From ∂f/∂z = 0
        if zmin <= z2 <= zmax:
            special_z.append(z2)

    # Add points with special z values
    for z in special_z:
        for x in (xmin, xmax):
            for y in (ymin, ymax):
                boundaries.append((x, y, z))
        # Check interior x if a < 0
        if a < 0:
            for y in (ymin, ymax):
                x = -(b * y * np.cos(z) + c * y * np.sin(z)) / (2 * a)
                if xmin < x < xmax:
                    boundaries.append((x, y, z))

    # Evaluate all candidates
    min_val, max_val = float("inf"), float("-inf")

    for point in boundaries:
        val = f(*point)
        if val < min_val:
            min_val = val
        if val > max_val:
            max_val = val

    return min_val, max_val
