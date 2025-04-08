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

    # Make a list of candidates
    candidates = []
    y_critical = [ymin, ymax]
    z_critical = [zmin, zmax]
    for zc in [-np.pi, -np.pi / 2, 0, np.pi / 2, np.pi, np.arctan2(c, b), np.arctan2(-c, -b), np.arctan2(-b, c), np.arctan2(b, -c)]:
        if zmin <= zc <= zmax:
            z_critical += [zc]
    for y in y_critical:
        for z in z_critical:
            x_critical = [xmin, xmax]
            if a != 0:
                for xc in [-(b * y * np.cos(z) + c * y * np.sin(z)) / (2 * a), (b * y * np.cos(z) + c * y * np.sin(z)) / (2 * a)]:
                    if xmin < xc < xmax:
                        x_critical += [xc]
            for x in x_critical:
                candidates.append((x, y, z))

    # Evaluate all candidates
    min_val, max_val = float("inf"), float("-inf")
    for point in candidates:
        val = f(*point)
        if val < min_val:
            min_val = val
        if val > max_val:
            max_val = val
    return min_val, max_val
