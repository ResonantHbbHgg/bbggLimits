import os
import math

# This script is provided by Olivier
# It reads the grid points of Non-Res parameter scan and provides various useful tools

def loadMapping_():
    """
    Load the mapping file
    """
    mapping_file = os.path.join(os.environ['CMSSW_BASE'], 'src', 'HiggsAnalysis/bbggLimits/data', 'list_all_translation_1507.txt')

    # 324 is dummy but it's the SM => we replace with results from BM we have
    dummy_points = [910, 985, 990]

    data = []
    with open(mapping_file) as f:
        index = 0
        for line in f.readlines():
            if index in dummy_points:
                index += 1
                continue

            tokens = line.split()
            row = {
                    'point': index,
                    'lambda': float(tokens[0]),
                    'yt': float(tokens[1]),
                    'c2': float(tokens[2]),
                    'cg': float(tokens[3]),
                    'c2g': float(tokens[4])
                    }

            data.append(row)
            index += 1

    # Remove duplicate
    seen = set()
    new_data = []
    for d in data:
        p = d['point']
        del d['point']
        t = tuple(d.items())
        d['point'] = p

        new_data.append(d)
        if t not in seen:
            seen.add(t)
            # new_data.append(d)

    return sorted(new_data, key=lambda t: t['point'])


def getPointFromParameters(l, yt, c2, cg, c2g, mapping = None):
    """
    Returns the point matching the 5 parameters, or None if there's none
    """
    if mapping==None:
        mapping = loadMapping_()

    point = [x['point'] for x in mapping if x['lambda'] == l and x['yt'] == yt and x['c2'] == c2 and x['cg'] == cg and x['c2g'] == c2g]

    if len(point) == 0:
        return None
    elif len(point) == 1:
        return point[0]
    else:
        print  'WARNING: More than 1 point mapping to a given set of parameters: %s' % (', '.join([str(p) for p in point]))
        print '\t We will take the first one'
        # raise Exception('More than 1 point mapping to a given set of parameters: %s' % (', '.join([str(p) for p in point])))
        return point[0]

def getParametersFromPoint(point, mapping = None, as_dict = False):
    """
    Returns a tuple of lambda, yt, c2, cg and c2g for a given point
    """
    if mapping==None:
        mapping = loadMapping_()

    if point < 0:
        return None

    point = [x for x in mapping if x['point'] == point]

    if len(point) == 0:
        return None
    elif len(point) != 1:
        raise Exception('More than 1 point mapping to a given set of parameters: %s' % (', '.join([str(p['point']) for p in point])))
    else:
        p = point[0]
        if as_dict:
            return p
        else:
            return (p['lambda'], p['yt'], p['c2'], p['cg'], p['c2g'])

def getLambdaScanPoints(yt=1, c2=0, cg=0, c2g=0):
    """
    Return a list of point suitable for the lambda scan
    """

    points = [x for x in mapping if x['yt'] == yt and x['c2'] == c2 and x['cg'] == cg and x['c2g'] == c2g]

    return [x['point'] for x in sorted(points, key=lambda t: t['lambda'])]

def getYtScanPoints(klambda=1, c2=0, cg=0, c2g=0, mapping=None):
    """
    Return a list of point suitable for the yt scan
    """
    if mapping==None:
        mapping = loadMapping_()


    points = [x for x in mapping if x['lambda'] == klambda and x['c2'] == c2 and x['cg'] == cg and x['c2g'] == c2g]

    return [x['point'] for x in sorted(points, key=lambda t: t['yt'])]

def getC2ScanPoints(klambda=1, yt=1, cg=0, c2g=0, mapping = None):
    """
    Return a list of point suitable for the c2 scan
    """
    if mapping==None:
        mapping = loadMapping_()

    points = [x for x in mapping if x['lambda'] == klambda and x['yt'] == yt and x['cg'] == cg and x['c2g'] == c2g]

    return [x['point'] for x in sorted(points, key=lambda t: t['c2'])]

def getCgScanPoints(klambda=1, yt=1, c2=0, c2g=0, mapping=None):
    """
    Return a list of point suitable for the cg scan
    """
    if mapping==None:
        mapping = loadMapping_()

    points = [x for x in mapping if x['lambda'] == klambda and x['yt'] == yt and x['c2'] == c2 and x['c2g'] == c2g]

    return [x['point'] for x in sorted(points, key=lambda t: t['cg'])]

def getC2gScanPoints(klambda=1, yt=1, c2=0, cg=0, mapping=None):
    """
    Return a list of point suitable for the c2g scan
    """
    if mapping==None:
        mapping = loadMapping_()

    points = [x for x in mapping if x['lambda'] == klambda and x['yt'] == yt and x['c2'] == c2 and x['cg'] == cg]

    return [x['point'] for x in sorted(points, key=lambda t: t['c2g'])]

def getPoints(f, mapping=None, *f_args):
    """
    Returns a list of points passing the filter f
    """
    try:
        filt = f_args[0]
    except:
        filt = None
    print 'Filter for points is:', filt
    if mapping==None:
        mapping = loadMapping_()

    points = [x['point'] for x in mapping if f(x, filt)]

    return points

def getCrossSectionForParameters(l, yt, c2, cg, c2g):
    """
    Get the theoretical cross-section for a given set of parameters in fb
    """

    params = (l, yt, c2, cg, c2g)

    A = [2.09078, 10.1517, 0.282307, 0.101205, 1.33191, -8.51168, -1.37309, 2.82636,
         1.45767, -4.91761, -0.675197, 1.86189, 0.321422, -0.836276, -0.568156]

    # From https://github.com/cms-hh/Plotting/blob/nonResonant/nonResonant/5Dfunction.py#L38
    def f(kl, kt, c2, cg, c2g):
      return A[0]*kt**4 + A[1]*c2**2 + (A[2]*kt**2 + A[3]*cg**2)*kl**2 + A[4]*c2g**2 + \
          (A[5]*c2 + A[6]*kt*kl)*kt**2 + (A[7]*kt*kl + A[8]*cg*kl)*c2 + A[9]*c2*c2g + \
          (A[10]*cg*kl + A[11]*c2g)*kt**2+ (A[12]*kl*cg + A[13]*c2g)*kt*kl + A[14]*cg*c2g*kl

    # From https://twiki.cern.ch/twiki/bin/view/LHCPhysics/LHCHXSWGHH#Current_recommendations_for_di_H
    hh_sm_xs = 33.45

    # Uncertainty
    hh_sm_xs_up = hh_sm_xs * (1 + math.sqrt(4.3**2 + 5**2 + 2.3**2 + 2.1**2) / 100)
    hh_sm_xs_down = hh_sm_xs * (1 - math.sqrt(6**2 + 5**2 + 2.3**2 + 2.1**2) / 100)

    res = f(*params)
    return res * hh_sm_xs, res * hh_sm_xs_down, res * hh_sm_xs_up

def getCrossSectionForPoint(point, mapping=None):
    """
    Get the theoretical cross-section for a given point in fb
    """
    return getCrossSectionForParameters(*getParametersFromPoint(point, mapping))

if __name__ == "__main__":

  # Load the mapping file from the `data` folder
  mapping = loadMapping_()
