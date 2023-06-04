import random
import numpy as np
import numba as nb

from time import time

from numpy.polynomial import Polynomial

from . import svd

def ransac_sequential(
    a,b,
    niter=100, 
    accuracy=0.1, 
    trial_frac=20, 
    inlier_frac=0.8,
):
    """ RANSAC https://en.wikipedia.org/wiki/Random_sample_consensus
    niter – maximum number of iterations allowed in the algorithm
    accuracy – threshold value to determine when a data point fits a model
    trial_frac – fraction of trial data points for fitting the model
    inlier_frac – fraction of close data points required
    """
    npnts = a.shape[0]
    besterr = np.inf
    bestfit = None
    bestinliers = None
    allids = list(range(npnts))
    maxinliers = 0
    thresh = int(np.floor(npnts*inlier_frac))
    nbatch = int(np.floor(npnts*trial_frac))
    curtime = time()
    for kk in range(niter):
        print('iter',kk,'maxinliers',maxinliers,'thresh',thresh,'besterr',besterr)
        maybeinliers = np.random.choice(allids,nbatch)
        a_,b_ = create_fit_matrices(a,b,maybeinliers)
        maybemodel,scales = svd._fit_x(a_,b_)
        b_calc = calculate_function(a,maybemodel)
        alsoinliers = ( np.abs(b_calc-b) < accuracy ).squeeze()
        if np.sum(alsoinliers)>maxinliers:
            maxinliers = np.sum(alsoinliers)
        if np.sum(alsoinliers) >= thresh:
            a_,b_ = create_fit_matrices(a,b,alsoinliers)
            b_calc_ = calculate_function(a_,maybemodel)
            bettermodel,scales = svd._fit_x(a_,b_)
            b_calc = calculate_function(a_,bettermodel*scales)
            thiserr = np.max(np.abs(b_calc-b_)) # option 1
            #thiserr = np.sum((b_calc-b_)**2) # option 2
            if thiserr < besterr:
                bestfit = bettermodel
                besterr = thiserr
                bestinliers = alsoinliers
            print(
                'iter>>>',kk,
                'minerr>>',np.min(np.abs(b_calc-b_)),
                'maxerr>>',np.max(np.abs(b_calc-b_)),
                'npoints>>>',len(b_calc),
                'len(maybeinliers)>>>',len(maybeinliers),
                'np.sum(alsoinliers)>>>',np.sum(alsoinliers),
                'thresh>>>',thresh,
                'besterr>>>',besterr,
                'thiserr>>>',thiserr,
            )
    bestinliers = np.where(bestinliers)[0]
    outliers = sorted(set(allids)-set(bestinliers))
    outliers = np.array(outliers)
    print('fit completed in %d iterations'%(kk+1))
    print('npnts>>>',npnts)
    print('maxinliers,thresh>>>',maxinliers,thresh)
    print('besterr>>>',besterr)
    print('%f sec. elapsed'%(time()-curtime))
    return bestfit,outliers,maxinliers,thresh

#@nb.njit(cache=True)
def create_fit_matrices(a,b,index):
    """ Create matrices for the linear SVD fit using 
        the active row and column indexes """
    #print('create_fit_matrices',index,index.shape)
    a_ = a[index,:]
    b_ = b[index,:]
    return a_,b_

#@nb.njit(cache=True)
def calculate_function(a,p):
    """ calculate linear function with matrix 'a' on params """
    return a@p

#####################
##### OLD CODE

def ransac_parallel_DEL(x, y, 
    order=3, 
    n=20, 
    k=100, 
    t=0.1, 
    #d=100, 
    f=0.8
):
    """ Thanks https://en.wikipedia.org/wiki/Random_sample_consensus
  
    n – minimum number of data points required to fit the model
    k – maximum number of iterations allowed in the algorithm
    t – threshold value to determine when a data point fits a model
    d – number of close data points required to assert that a model fits well to data
    f – fraction of close data points required
    """
    x = np.array(x)
    y = np.array(y)
    besterr = np.inf
    bestfit = None
    bestinliers = None
    indexes = list(range(len(x)))
    maxinliers = 0
    thresh = int(np.floor(len(x)*f))
    curtime = time()
    for kk in range(k):
        maybeinliers = random.sample(indexes,n)
        maybemodel = polyfit(x[maybeinliers], y[maybeinliers], order)
        alsoinliers = np.abs(np.polyval(maybemodel, x)-y) < t
        if sum(alsoinliers)>maxinliers:
            maxinliers = sum(alsoinliers)
        if sum(alsoinliers) >= thresh:
            bettermodel = polyfit(x[alsoinliers], y[alsoinliers], order)
            #thiserr = np.sum(np.abs(np.polyval(bettermodel, x[alsoinliers])-y[alsoinliers]))
            thiserr = np.max(np.abs(np.polyval(bettermodel, x[alsoinliers])-y[alsoinliers]))
            if thiserr < besterr:
                bestfit = bettermodel
                besterr = thiserr
                bestinliers = alsoinliers
    print('fit completed in %d iterations'%(kk+1))
    print('maxinliers,thresh>>>',maxinliers,thresh)
    print('len(x)>>>',len(x))
    print('%f sec. elapsed'%(time()-curtime))
    bestinliers = np.where(bestinliers)[0]
    outliers = sorted(set(indexes)-set(bestinliers))
    outliers = np.array(outliers)
    return bestfit,outliers,maxinliers,thresh