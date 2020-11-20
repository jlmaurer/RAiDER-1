import numpy as np

from geopy.distance import geodesic

def geodesicDistance(X):
    ''' 
    Calculate the pair-wise differences between lat/lon points on an ellipsoid 
    
    Parameters
    ----------
    X - a N x (2,3) matrix that represents either (lat, lon) coordinates or 
         (lat, lon, heights) coordinates on the WGS84 ellipsoid

    Returns
    -------
    Y - a condensed distance matrix Y. For each i and j (where i < j < m), where
         m is the number of original obesrvations. the metric dist(X[i], X[j] is 
         computed and stored in entry ij
    '''
    [N, D] = X.shape
    Y = np.zeros((N, N))

    for i1 in range(N):
        for i2 in range(N):
            if (i1 < i2) and (i2 < N):
                Y[i1,i2] = geodesic(X[i1,:2], X[i2,:2]).km

    return Y[~(Y==0)]
    

