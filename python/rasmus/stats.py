# python libs
from math import *
import cmath
import random
import os

# rasmus libs
from rasmus import util
from rasmus import algorithms







def mean(vals):
    return sum(vals) / float(len(vals))

def median(vals):
    lenvals = len(vals)
    sortvals = util.sort(vals)
    
    if lenvals % 2 == 0:
        return (sortvals[lenvals / 2] + sortvals[lenvals / 2 - 1]) / 2.0
    else:
        return sortvals[lenvals / 2]

def mode(vals):
    top = 0
    topkey = None
    for key, val in util.histDict(vals).iteritems():
        if val > top:
            top = val
            topkey = key
    return topkey


def msqerr(vals1, vals2):
    """Mean squared error"""
    
    assert len(vals1) == len(vals2), "lists are not the same length"
    
    
    return mean([(vals1[i] - vals2[i]) ** 2 
                 for i in xrange(len(vals1))])
    
    

def variance(vals):
    """Variance"""
    u = mean(vals)
    return sum(map(lambda x: (x - u)**2, vals)) / float(len(vals)-1)

def sdev(vals):
    """Standard deviation"""
    return sqrt(variance(vals))

def serror(vals):
    """Stanadrd error"""
    return sdev(vals) / sqrt(len(vals))

def covariance(lst1, lst2):
    """Covariance"""
    m1 = mean(lst1)
    m2 = mean(lst2)
    tot = 0.0
    for i in xrange(len(lst1)):
        tot += (lst1[i] - m1) * (lst2[i] - m2)    
    return tot / (len(lst1)-1)


def covmatrix(mat):
    """Covariance Matrix"""
    size = len(mat)
    
    return util.list2matrix(map(lambda (i,j): cov(mat[i], mat[j]), 
                            util.range2(size, size)),
                            size, size)

def corrmatrix(mat):
    """Correlation Matrix"""
    size = len(mat)
    
    return util.list2matrix(map(lambda (i,j): corr(mat[i], mat[j]), 
                            util.range2(size, size)),
                            size, size)


def corr(lst1, lst2):
    """Pearson's Correlation"""
    num = covariance(lst1, lst2)
    denom = float(sdev(lst1) * sdev(lst2))
    if denom != 0:
        return num / denom
    else:
        return 1e1000


def qqnorm(data, plot=None):
    """Quantile-quantile plot"""
    
    data2 = util.sort(data)
    norm = [random.normalvariate(0, 1) for x in range(len(data2))]
    norm.sort()
    
    if plot == None:
        return util.plot(data2, norm)
    else:
        plot.plot(data2, norm)
        return plot



def fitLine(xlist, ylist):
    """2D regression"""
    
    xysum = 0
    xxsum = 0
    n = len(xlist)        
    for i in range(n):
        xysum += xlist[i] * ylist[i]
        xxsum += xlist[i] * xlist[i]
    avgx = mean(xlist)
    avgy = mean(ylist)

    if (xxsum - n*avgx*avgx) == 0:
        slope = 1e10
    else:
        slope = (xysum - n*avgx*avgy) / float(xxsum - n*avgx*avgx)

    inter = (avgy*xxsum - avgx*xysum) / float(xxsum - n*avgx*avgx)

    return (slope, inter)


def fitLineError(xlist, ylist, slope, inter):
    error = 0
    n = len(xlist)
    
    for i in range(n):
        error += ((xlist[i]*slope + inter) - ylist[i]) ** 2
    return error / n


def percentile(vals, perc, rounding=-1, sort=True):
    """Give the value at a percentile
       
       rounding -- round down if -1 or round up for 1
    """
    
    if sort:
        vals2 = util.sort(vals)
    else:
        vals2 = vals
    n = len(vals2)
    if rounding == -1:
        return vals2[util.clamp(int(perc * n), 0, n-1)]
    elif rounding == 1:
        return vals2[util.clamp(int(ceil(perc * n)), 0, n-1)]
    else:
        raise Exception("rounding must be 1 or -1")


def smooth(vals, radius):
    """
    return an averaging of vals using a radius
    
    Note: not implemented as fast as possible
    runtime: O(len(vals) * radius)
    """
    
    vals2 = []
    vlen = len(vals)
    
    for i in xrange(vlen):
        radius2 = min(i, vlen - i - 1, radius)
        vals2.append(mean(vals[i-radius2:i+radius2+1]))
    
    return vals2


def smooth2(x, y, xradius):
    """
    return an averaging of x and y using xradius
    
    x must be sorted least to greatest
    """

    vlen = len(x)
    assert vlen == len(y)
    
    # simple case
    if vlen == 0:
        return [], []
    
    x2 = []
    y2 = []
    
    start = min(x)
    end = max(x)
    xtot = x[0]
    ytot = y[0]
    
    low = 0
    high = 0
    
    for i in xrange(vlen):
        xi = x[i]
    
        xradius2 = min(xi - start, end - xi, xradius)
    
        # move window
        while x[low] < xi - xradius2:
            xtot -= x[low]
            ytot -= y[low]
            low += 1
        while x[high] < xi + xradius2:
            high += 1
            xtot += x[high]
            ytot += y[high]
        
        denom = float(high - low + 1)
        x2.append(xtot / denom)
        y2.append(ytot / denom)
    
    return x2, y2



def smooth_old(x, radius):
    """
    return an averaging of vals using a radius
    
    Note: not implemented as fast as possible
    runtime: O(len(vals) * radius)
    """
    
    vlen = len(x)
    
    # simple case
    if vlen == 0:
        return []
    
    x2 = []
    
    tot = x[0]
    
    low = 0
    high = 0
    
    for i in range(vlen):
        xi = x[i]
    
        xradius2 = min(i, vlen - i - 1, xradius)
    
        # move window
        while x[low] < xi - xradius2:
            xtot -= x[low]
            ytot -= y[low]
            low += 1
        while x[high] < xi + xradius2:
            high += 1
            xtot += x[high]
            ytot += y[high]
        
        denom = float(high - low + 1)
        x2.append(xtot / denom)
        y2.append(ytot / denom)
    
    return x2, y2


def factorial(x, k=1):
    """Simple implementation of factorial"""
    
    n = 1
    for i in xrange(int(k)+1, int(x)+1):
        n *= i
    return n


def choose(n, k):
    return factorial(n, n - k) / factorial(k)



def sample(weights):
    """
    Randomly choose an int between 0 and len(probs)-1 using
    the weights stored in list probs.
    
    item i will be chosen with probability weights[i]/sum(weights)
    """
    
    probs = util.oneNorm(weights)
    
    cdf = [0]
    for i in range(1, len(probs)):
        cdf.append(cdf[-1] + probs[i-1])
    
    pick = random.random()
    
    low,top = algorithms.binsearch(cdf, pick)
    
    assert low != None
    
    return low
    

    
def chyper(m, n, M, N, report=0):
    '''
    calculates cumulative probability based on
    hypergeometric distribution
    over/under/both (report = 0/1/2)
    (uses /seq/compbio02/software-Linux/misc/chyper)
    '''

    assert( (type(m) == type(n) == type(M) == type(N) == int)
            and m <= n and m <= M and n <= N)

    command = "chyper %d %d %d %d 2>/dev/null" % (m, n, M, N)
    stream = os.popen(command)
    val = stream.read()
    if val == '':
        raise Exception("error in chyper")
    else:
        val = val.strip()
        vals = map(float, val.split(' ')[4:6])
        
    if report == 0:
        #p-val for over-repr.
        return vals[0]
    elif report == 1:
        #p-val for under-repr.
        return vals[1]
    elif report == 2:
        #tuple (over, under)
        return vals
    else:
        raise "unknown option"


def rhyper(m, n, M, N, report=0):
    from rpy import r

    '''
    calculates cumulative probability based on
    hypergeometric distribution
    over/under/both (report = 0/1/2)
    (uses R through RPy)
    '''
    
    assert( (type(m) == type(n) == type(M) == type(N) == int)
            and m <= n and m <= M and n <= N)
    
    
    if report == 0:
        #p-val for over-repr.
        return r.phyper(m-1, M, N-M, n, lower_tail=False)
    elif report == 1:
        #p-val for under-repr.
        return r.phyper(m, M, N-M, n)
    elif report == 2:
        #tuple (over, under)
        return r.phyper(m-1, M, N-M, n, lower_tail=False), r.phyper(m, M, N-M, n)
    else:
        raise "unknown option"

    
#=============================================================================
# Distributions
#

def bionomialPdf(n, params):
    p, k = params
    return choose(n, k) * (p ** k) * ((1-p) ** (n - k))

def gaussianPdf(x, params):
    return 1/sqrt(2*pi) * exp(- x**2 / 2.0)

def normalPdf(x, params):
    mu, sigma = params
    return 1.0/(sigma * sqrt(2.0*pi)) * exp(- (x - mu)**2 / (2.0 * sigma**2))

def normalCdf(x, params):
    mu, sigma = params
    return (1 + erf((x - mu)/(sigma * sqrt(2)))) / 2.0

def logNormalPdf(x, params):
    """mu and sigma are the mean and standard deviation of the 
       variable's logarithm"""
    
    mu, sigma = params
    return 1/(x * sigma * sqrt(2*pi)) * \
           exp(- (log(x) - mu)**2 / (2.0 * sigma**2))

def logNormalCdf(x, params):
    """mu and sigma are the mean and standard deviation of the 
       variable's logarithm"""
    
    mu, sigma = params
    return (1 + erf((log(x) - mu)/(sigma * sqrt(2)))) / 2.0


def poissonPdf(x, params):
    p = params[0]
    if x > 0 and p > 0:
        return exp(p) * p **x / x
    else:
        return 0.0


def gammaPdf(x, params):
    alpha, beta = params
    if x <= 0 or alpha <= 0 or beta <= 0:
        return 0.0
    else:
        return (exp(-x * beta) * (x ** (alpha - 1)) * (beta ** alpha)) / \
           gamma(alpha)

def gammaPdf2(x, params):
    alpha, beta = params
    if x <= 0 or alpha <= 0 or beta <= 0:
        return 0.0
    else:
        return exp(- x * beta + (alpha - 1)*log(x) + alpha * log(beta) -
                   gammaln(alpha))


def gammaCdf(x, params):
    alpha, beta = params
    if x <= 0:
        return 0
    else:
        return gammainc(alpha, x * beta) / gamma(alpha)


def betaPdf2(x, params):
    """A simpler implementation of beta distribution but will overflow
       for values of alpha and beta near 100
    """
    
    alpha, beta = params
    if 0 < x < 1 and alpha > 0 and beta > 0:
        return gamma(alpha + beta) / (gamma(alpha)*gamma(beta)) * \
               x ** (alpha-1) * (1-x)**(beta-1)
    else:
        return 0.0

def betaPdf(x, params):
    alpha, beta = params
    
    if 0 < x < 1 and alpha > 0 and beta > 0:
        return e**(gammaln(alpha + beta) - (gammaln(alpha) + gammaln(beta)) + \
                   (alpha-1) * log(x) +  (beta-1) * log(1-x))
    else:
        return 0.0
    


def betaPdf3(x, params):
    alpha, beta = map(int, params)
    if 0 < x < 1 and alpha > 0 and beta > 0:
        n = min(alpha-1, beta-1)
        m = max(alpha-1, beta-1)
        
        prod1 = 1
        for i in range(1,n+1):
            prod1 *= ((n+i)*x*(1-x))/i
        
        prod2 = 1
        if alpha > beta:
            for i in range(n+1, m+1):
                prod2 *= ((n+i)*x)/i
        else:
            for i in range(n+1, m+1):
                prod2 *= ((n+i)*(1-x))/i
        
        return prod1 * prod2 * (alpha + beta - 1)
    else:
        return 0.0


def gamma(x):
    """
    Lanczos approximation to the gamma function. 
    
    found on http://www.rskey.org/gamma.htm   
    """
    
    ret = 1.000000000190015 + \
          76.18009172947146 / (x + 1) + \
          -86.50532032941677 / (x + 2) + \
          24.01409824083091 / (x + 3) + \
          -1.231739572450155 / (x + 4) + \
          1.208650973866179e-3 / (x + 5) + \
          -5.395239384953e-6 / (x + 6)
    
    return ret * sqrt(2*pi)/x * (x + 5.5)**(x+.5) * exp(-x-5.5)



def gammaln(xx):
    """
    From numerical alogrithms in C

    float gammln(float xx)
    Returns the value ln[(xx)] for xx > 0.
    {
        Internal arithmetic will be done in double precision, a nicety that you can omit if five-figure
        accuracy is good enough.
        double x,y,tmp,ser;
        static double cof[6]={76.18009172947146,-86.50532032941677,
             24.01409824083091,-1.231739572450155,
             0.1208650973866179e-2,-0.5395239384953e-5};
        int j;
        y=x=xx;
        tmp=x+5.5;
        tmp -= (x+0.5)*log(tmp);
        ser=1.000000000190015;
        for (j=0;j<=5;j++) ser += cof[j]/++y;
        return -tmp+log(2.5066282746310005*ser/x);
    }
    """

    cof = [76.18009172947146,-86.50532032941677,
         24.01409824083091,-1.231739572450155,
         0.1208650973866179e-2,-0.5395239384953e-5]
    
    y = x = xx
    tmp = x + 5.5
    tmp -= (x + 0.5) * log(tmp)
    ser = 1.000000000190015
    
    for j in range(6):
        y += 1
        ser += cof[j] / y
    
    return - tmp + log(2.5066282746310005 * ser / x)
    



GAMMA_INCOMP_ACCURACY = 1000
def gammainc(a, x):
    # found on http://www.rskey.org/gamma.htm
    
    ret = 0
    term = 1.0/x
    for n in xrange(GAMMA_INCOMP_ACCURACY):
        term *= x/(a+n)
        ret += term
        if term < .0001:
            break
    return x**a * exp(-x) * ret


def erf(x):
    # http://www.theorie.physik.uni-muenchen.de/~serge/erf-approx.pdf
    
    a = 8/(3*pi) * (pi - 3)/(4 - pi)
    axx = a * x * x
    
    if x >= 0:
        return sqrt(1 - exp(-x*x * (4.0/pi + axx)/(1 + axx)))
    else:
        return - sqrt(1 - exp(-x*x * (4.0/pi + axx)/(1 + axx)))
    


def chiSquare(rows, expected=None, nparams=0):
    # ex: rows = [[1,2,3],[1,4,5]]
    assert(util.equal(map(len,rows)))

    if 0 in map(sum,rows): return 0,1.0
    cols = zip(* rows)
    if 0 in map(sum,cols): return 0,1.0

    if not expected:
        expected = make_expected(rows)

    chisq = 0
    for obss,exps in zip(rows,expected):
        for obs, exp in zip(obss, exps):
            chisq += ((obs-exp)**2)/exp

    df = max(len(rows)-1, 1)*max(len(rows[0])-1, 1) - nparams

    p = chi_square_lookup(chisq,df)

    return chisq,p


def make_expected(rows):
    rowtotals = map(sum, rows)
    coltotals = map(sum, zip(* rows))
    grandtotal = float(sum(rowtotals))

    expected = []
    for row,rowtotal in zip(rows,rowtotals):
        expected_row = []
        for obs, coltotal in zip(row, coltotals):
            exp = rowtotal * coltotal / grandtotal
            expected_row.append(exp)
        expected.append(expected_row)
    return expected


def chiSquareFit(xbins, ybins, func, nsamples, nparams, minsamples=5):
    sizes = [xbins[i+1] - xbins[i] for i in xrange(len(xbins)-1)]
    sizes.append(sizes[-1])
    
    # only focus on bins that are large enough
    counts = [ybins[i] * sizes[i] * nsamples for i in xrange(len(xbins)-1)]
    
    expected = []
    for i in xrange(len(xbins)-1):
        expected.append((func(xbins[i]) + func(xbins[i+1]))/2.0 * 
                         sizes[i] * nsamples)
        
    # ensure we have enough expected samples in each bin
    ind = util.find(util.gefunc(minsamples), expected)
    counts = util.mget(counts, ind)
    expected = util.mget(expected, ind)
    
    if len(counts) == 0:
        return [0, 1], counts, expected
    else:
        return chiSquare([counts], [expected], nparams), counts, expected


chi_square_table = {
    1: [1.64, 2.71, 3.84, 5.02, 6.64, 10.83],
    2: [3.22, 4.61, 5.99, 7.38, 9.21, 13.82],
    3: [4.64, 6.25, 7.82, 9.35, 11.34, 16.27],
    4: [5.99, 7.78, 9.49, 11.14, 13.28, 18.47],
    5: [7.29, 9.24, 11.07, 12.83, 15.09, 20.52],
    6: [8.56, 10.64, 12.59, 14.45, 16.81, 22.46],
    7: [9.80, 12.02, 14.07, 16.01, 18.48, 24.32],
    8: [11.03, 13.36, 15.51, 17.53, 20.09, 26.12],
    9: [12.24, 14.68, 16.92, 19.02, 21.67, 27.88],
    10: [13.44, 15.99, 18.31, 20.48, 23.21, 29.59],
    11: [14.63, 17.28, 19.68, 21.92, 24.72, 31.26],
    12: [15.81, 18.55, 21.03, 23.34, 26.22, 32.91],
    13: [16.98, 19.81, 22.36, 24.74, 27.69, 34.53],
    14: [18.15, 21.06, 23.68, 26.12, 29.14, 36.12],
    15: [19.31, 22.31, 25.00, 27.49, 30.58, 37.70],
    16: [20.47, 23.54, 26.30, 28.85, 32.00, 39.25],
    17: [21.61, 24.77, 27.59, 30.19, 33.41, 40.79],
    18: [22.76, 25.99, 28.87, 31.53, 34.81, 42.31],
    19: [23.90, 27.20, 30.14, 32.85, 36.19, 43.82],
    20: [25.04, 28.41, 31.41, 34.17, 37.57, 45.31],
    21: [26.17, 29.62, 32.67, 35.48, 38.93, 46.80],
    22: [27.30, 30.81, 33.92, 36.78, 40.29, 48.27],
    23: [28.43, 32.01, 35.17, 38.08, 41.64, 49.73],
    24: [29.55, 33.20, 36.42, 39.36, 42.98, 51.18],
    25: [30.68, 34.38, 37.65, 40.65, 44.31, 52.62],
    26: [31.79, 35.56, 38.89, 41.92, 45.64, 54.05],
    27: [32.91, 36.74, 40.11, 43.19, 46.96, 55.48],
    28: [34.03, 37.92, 41.34, 44.46, 48.28, 56.89],
    29: [35.14, 39.09, 42.56, 45.72, 49.59, 58.30],
    30: [36.25, 40.26, 43.77, 46.98, 50.89, 59.70]
}


def chi_square_lookup(value, df):
    
    ps = [0.20, 0.10, 0.05, 0.025, 0.01, 0.001]
    
    if df <= 0:
        return 1.0    
    
    row = chi_square_table[min(df, 30)]

    for i in range(0,len(row)):
        if row[i] >= value:
            i = i-1
            break
    
    if i == -1: return 1
    else: return ps[i]


def ttest(lst1, lst2):
    sdevdist = sqrt(var(lst1)/len(lst1) + var(lst2)/len(lst2))
    t = abs(mean(lst1) - mean(lst2)) / sdevdist
    df = len(lst2) + len(lst2) - 2
    
"""
t-table

 	0.1  	0.05  	0.01  	0.001
1 	6.31 	12.71 	63.66 	636.62
2 	2.92 	4.30 	9.93 	31.60
3 	2.35 	3.18 	5.84 	12.92
4 	2.13 	2.78 	4.60 	8.61
5 	2.02 	2.57 	4.03 	6.87
6 	1.94 	2.45 	3.71 	5.96
7 	1.89 	2.37 	3.50 	5.41
8 	1.86 	2.31 	3.36 	5.04
9 	1.83 	2.26 	3.25 	4.78
10 	1.81 	2.23 	3.17 	4.59
11 	1.80 	2.20 	3.11 	4.44
12 	1.78 	2.18 	3.06 	4.32
13 	1.77 	2.16 	3.01 	4.22
14 	1.76 	2.14 	2.98 	4.14
15 	1.75 	2.13 	2.95 	4.07
16 	1.75 	2.12 	2.92 	4.02
17 	1.74 	2.11 	2.90 	3.97
18 	1.73 	2.10 	2.88 	3.92
19 	1.73 	2.09 	2.86 	3.88
20 	1.72 	2.09 	2.85 	3.85
21 	1.72 	2.08 	2.83 	3.82
22 	1.72 	2.07 	2.82 	3.79
23 	1.71 	2.07 	2.82 	3.77
24 	1.71 	2.06 	2.80 	3.75
25 	1.71 	2.06 	2.79 	3.73
26 	1.71 	2.06 	2.78 	3.71
27 	1.70 	2.05 	2.77 	3.69
28 	1.70 	2.05 	2.76 	3.67
29 	1.70 	2.05 	2.76 	3.66
30 	1.70 	2.04 	2.75 	3.65
40 	1.68 	2.02 	2.70 	3.55
60 	1.67 	2.00 	2.66 	3.46
120 1.66 	1.98 	2.62 	3.37  
"""    

"""
r	90%	95%	97.5%	99.5%
1	3.07766	6.31371	12.7062	63.656
2	1.88562	2.91999	4.30265	9.92482
3	1.63774	2.35336	3.18243	5.84089
4	1.53321	2.13185	2.77644	4.60393
5	1.47588	2.01505	2.57058	4.03212
10	1.37218	1.81246	2.22814	3.16922
30	1.31042	1.69726	2.04227	2.74999
100	1.29007	1.66023	1.98397	2.62589
infty	1.28156	1.64487	1.95999	2.57584
"""


def spearman(vec1, vec2):
    """Spearman's rank test"""
    
    assert len(vec1) == len(vec2), "vec1 and vec2 are not the same length"
    
    n = len(vec1)
    rank1 = util.sortInd(vec1)
    rank2 = util.sortInd(vec2)
    
    R = sum((vec1[i] - vec2[i])**2 for i in xrange(n))
    
    Z = (6*R - n*(n*n - 1)) / (n*(n + 1) * sqrt(n - 1))
    
    return Z
    


# input:
#   xdata, ydata  - data to fit
#   func          - a function of the form f(x, params)
#
def fitCurve(xdata, ydata, func, paramsInit):   
    import scipy

    y = scipy.array(ydata)
    p0 = scipy.array(paramsInit)
    
    def error(params):
        y2 = scipy.array(map(lambda x: func(x, params), xdata))
        return y - y2

    params, msg = scipy.optimize.leastsq(error, p0)
    
    resid = error(params)
    
    return list(params), sum(resid*resid)

    
def fitDistrib(func, paramsInit, data, start, end, step, perc=1.0):
    xdata, ydata = util.distrib(data, low=start, width=step)
    ydata = [i / perc for i in ydata]
    xdata = util.histbins(xdata)
    params, resid = fitCurve(xdata, ydata, func, paramsInit)
    return params, resid
    

def plotfuncFit(func, paramsInit, xdata, ydata, start, end, step, plot = None,
                **options):
    if not plot:
        plot = util.Gnuplot()
    
    options.setdefault('style', 'boxes')
    
    params, resid = fitCurve(xdata, ydata, func, paramsInit)
    plot.plot(util.histbins(xdata), ydata, **options)
    plot.plotfunc(lambda x: func(x, params), start, end, step)
    
    return plot, params, resid
    

def plotdistribFit(func, paramsInit, data, start, end, step, plot = None,
                   **options):
    xdata, ydata = util.distrib(data, low=start, width=step)
    return plotfuncFit(func, paramsInit, xdata, ydata, start, end, step/10, plot,
                       **options)


    


def solveCubic(a, b, c, real=True):
    """solves x^3 + ax^2 + bx + x = 0 for x"""
    
    p = b - a*a / 3.0
    q = c + (2*a*a*a - 9*a*b) / 27.0
    
    # special case: avoids division by zero later on
    if p == q == 0:
        return [- a / 3.0]
    
    # 
    # u = (q/2 +- sqrt(q^2/4 + p^3/27))^(1/3)
    #
    
    # complex math is used to find complex roots
    sqrteqn = cmath.sqrt(q*q/4.0 + p*p*p/27.0)
    
    # find fist cube root
    u1 = (q/2.0 + sqrteqn)**(1/3.0)
    
    # special case: avoids division by zero later on
    if u1 == 0:
        u1 = (q/2.0 - sqrteqn)**(1/3.0)
    
    # find other two cube roots
    u2 = u1 * complex(-.5, -sqrt(3)/2)
    u3 = u1 * complex(-.5, sqrt(3)/2)
    
    # finds roots of cubic polynomial
    root1 = p / (3*u1) - u1 - a / 3.0
    root2 = p / (3*u2) - u2 - a / 3.0
    root3 = p / (3*u3) - u3 - a / 3.0
    
    if real:
        return [root1.real, root2.real, root3.real]
    else:
        return [root1, root2, root3]
    
    
    
    
    
    
    
    
    
    

