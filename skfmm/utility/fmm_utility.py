import numpy as np

def contour_lengths_2d(X, Y, phi, contour_value=0.0):
    try:
        import pylab as pl
    except ImportError:
        raise Exception("This function requires matplotlib")

    contour_set = pl.contour(X, Y, phi, [contour_value])
    contour_lengths = []
    for path in contour_set.collections[0].get_paths():
        local_length = 0.0
        path_list = list(path.iter_segments())
        for i in range(len(path_list)-1):
            p0, p1 = path_list[i][0], path_list[i+1][0]
            local_length += (np.linalg.norm(p0-p1)).sum()
        contour_lengths.append(local_length)

    return np.array(contour_lengths)

def total_contour_length_2d(X, Y, phi, contour_value=0):
    return contour_lengths_2d(X, Y, phi, contour_value).sum()

def contour_area_2d(phi, dx, contour_value=0.0):
    phi = np.asfarray(phi)
    dx = float(dx)
    dxn = dx ** len(phi.shape)
    return (phi < contour_value).sum() * dxn

def normal_direction(phi):
    vectors = np.array(np.gradient(phi))
    norms = np.apply_along_axis(np.linalg.norm, 0, vectors)
    nv = vectors / norms

def curvature():
    raise NotImplemented("this function is not implimented yet!")
