import numpy as np

def twoPiBound(heading : float) -> float:
    if heading < 0:
        heading = heading + 2 * np.pi

    return heading