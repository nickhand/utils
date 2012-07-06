#
# A module for handling catalogs of sources
# 
#

import numpy
import pickle, copy
from flipper import *

defaultCols = { \
        'ra'   : { 'desc':'Right Ascension of the object'   , 'order': 0, 'type': float, 'fmt' : '%10.4f', 'default': -999.}, \
        'dec'  : { 'desc':'Declination of the object'       , 'order': 1, 'type': float, 'fmt' : '%10.4f', 'default': -999.}, \
        's/n'  : { 'desc':'signal to noise of the detection', 'order': 2, 'type': float, 'fmt' : '%10.4f', 'default': -999.}
        }

def compare_by (fieldname): 
    def compare_two_dicts (a, b): 
        return cmp(a[fieldname], b[fieldname]) 
    return compare_two_dicts 

def convertRADecDegreesToSexagesimal(ra, dec):

    if ra < 0:
        ra += 360.
    ra_hr = int(ra/15.)
    minsec = numpy.mod(ra,15.)/15.*60
    ra_min = int(minsec)
    ra_sec = numpy.mod(minsec,1.)*60

    dec_deg = abs(int(dec))
    minsec = numpy.mod(abs(dec),1.)*60
    dec_min = int(minsec)
    dec_sec = numpy.mod(minsec,1.)*60

    return (ra_hr, ra_min, ra_sec), (dec_deg, dec_min, dec_sec)
    
def convertSexagesimalToRaDecDegrees(ra, dec):
    ra_hr, ra_min, ra_sec= ra.split(':')

    sign = ra_hr[0]
    if sign == '-':
        sign_ra = -1.
        ra_hr = ra_hr[1:]
    else:
        sign_ra = 1.

    ra = sign_ra * 15. * ( float(ra_hr) + float(ra_min)/60. + float(ra_sec)/3600. )

    dec_deg, dec_min, dec_sec= dec.split(':')

    sign = dec_deg[0]
    if sign == '-':
        sign_dec = -1.
        dec_deg  = dec_deg[1:]
    else:
        sign_dec = 1.

    dec =  sign_dec * (float(dec_deg) + float(dec_min)/60. + float(dec_sec)/3600.)

    return ra, dec

class catalog( list ):

    def __init__( self, cols = None):
        if cols == None:
            self.cols = defaultCols
        else:
            self.cols = cols
        self.ncol = len(self.cols.keys())
        self.sep = ' '

    def selectByID(self, id):
        k = 0
        for row in self:
            if row['id'] == id:
                return k
            k += 1
        
        return None

    def nRows(self):
        i = 0
        for row in self:
            i += 1
        return i

    def partition(self, colName, N, reverse=True):
        self.sortBy(colName)
        if reverse:
            self.reverse()
        out = []
        size = self.nRows()/N
        for i in range(0, self.nRows(), size):
            if len(out) == N:
                for row in self[i:i+size]:
                    out[-1].append(copy.copy(row))
            else:
                sub = catalog(cols = self.cols)
                for row in self[i:i+size]:
                    sub.append(copy.copy(row))
            
                out.append(sub)

        return out

    def select( self, colName, func ):
        sub = catalog(cols = self.cols)
        for row in self:
            if colName != None:
                if func(row[colName]):
                    sub.append(copy.copy(row))
            else:
                if func(row):
                    sub.append(copy.copy(row))
        return sub

    def shave(self, minIndex, maxIndex, key, reverse=True):

        colNames = self.cols
        catNew = catalog(cols = colNames)
        cnt = 0

        self.sortBy(key)
        if (reverse):
            self.reverse()
        for row in self:
            if cnt >= minIndex and cnt < maxIndex:
                catNew.addRow(row)
            if cnt >= maxIndex:
                break

            cnt += 1
            
        return catNew
        
    def appendCatalog( self, cat ):

        colNames = self.cols
        catNew = catalog(cols = colNames)
        for row in self:
            newrow = row
            catNew.addRow(newrow)
        for row in cat:
            newrow = row
            catNew.addRow(newrow)

        return catNew 

    def addRow( self, row ):
        newRow = {}
        inputKeys = row.keys()
        for k in self.cols.keys():
            if k in inputKeys:
                newRow[k]=row[k]
            else:
                newRow[k]=self.cols[k]['default']
        self.append(newRow)

    def addCol( self, name, desc, order, dtype, fmt, default ):

        cnames = self.cols.keys()
        if name in cnames:
            raise ValueError("Name %s already in column list" % name)
        for n in cnames:
            if self.cols[n]['order'] >=  order:
                self.cols[n]['order'] += 1
        
        self.cols[name] = {'desc': desc   , 'order': order, 'type': dtype, 'fmt' : fmt, 'default': default}
        self.ncol = len(self.cols.keys())
       


    def sortBy( self, col ):
        self.sort(compare_by(col))
        

    def __str__( self ):
        string = ""
        colNames = self.cols.keys()
        for colName in self.cols.keys():
            colNames[self.cols[colName]['order']] = colName
        for row in self:
            for k in colNames:
                if self.cols[k]['default'] == row[k] and 'defaultFmt' in self.cols[k].keys():
                    string += self.cols[k]['defaultFmt'] % row[k]
                    string += self.sep
                    continue
                try:
                    string += self.cols[k]['fmt'] % row[k]
                except:
                    print "Failed to format col %s with this value %s." % (k, str(row[k])), self.cols[k]['default'] == row[k], self.cols[k]['default']
                    raise
                string += self.sep
            string += '\n'
        return string

    def colNamesSorted( self ):
        colNames = self.cols.keys()
        for colName in self.cols.keys():
            colNames[self.cols[colName]['order']] = colName
        return colNames

    def colsString( self, comment = False ):
        colNames = self.colNamesSorted()
        cs = ""
        for i in range(self.ncol):
            if comment:
                cs += "# %d : %s - %s\n" % (i, colNames[i], self.cols[colNames[i]]['desc'])
            else:
                cs += "%d : %s - %s\n" % (i, colNames[i], self.cols[colNames[i]]['desc'])
        return cs

    def writeRegionFile( self, filename, radius = .1, raCol = 'ra', decCol = 'dec' ):
        string = ""
        f = file( filename, 'w' )
        header = "global color=black font='helvetica 10 normal' select=1 highlite=1 edit=1 move=1 delete=1 include=1 fixed=0 source\n" 
        f.write(header)
        f.write("wcs\n")
        for row in self:
            if ('s/n' in row.keys() and row['s/n'] >= 5) or ('snr' in row.keys() and row['snr'] >= 5):
                f.write("circle( %f, %f, .1 )#color=red\n" % (row[raCol], row[decCol]))
            else:
                f.write("circle( %f, %f, .1 )\n" % (row[raCol], row[decCol]))
        f.close()

    def read( self, filename, binary = False ):
        if binary:
            pass

        f = file(filename)
        lines = f.readlines()
        for line in lines:
            if line[0] == '#' or len(line) <= 1:
                continue
            fields = line.split()
            if fields[0][0] == '#' or len(fields) == 0:
                continue
            self.append({})
            for colName in self.cols.keys():
                col = self.cols[colName]
                if col['order'] < 0: 
                    continue
                if fields[col['order']] == col['default']:
                    self[-1][colName] = col['default']
                else:
                    try:
                        self[-1][colName] = col['type'](fields[col['order']])
                    except:
                        self[-1][colName] = eval(col['type'])(fields[col['order']])
        f.close()

    def write( self, filename ):
        f = file(filename, 'w')
        pickle.dump( self, f )
        f.close()

    def writeASCII( self, filename, header=True, mode = 'w' ):
        f = file( filename, mode ) 
        if header:
            f.write(self.colsString(comment=True))
        f.write(self.__str__())
        fd = flipperDict.flipperDict()
        fd['cols'] = self.cols
        fd.writeToFile("%s.dict" % filename)

    def arrayFromCol(self, colName):
        c = []
        for row in self:
            c.append(row[colName])
        return numpy.array(c)

    #XXX under construction
    def plotCol1VsCol2(self, col1, col2, err1 = None, err2 = None ):
        cat = catalog.readFromPickleFile("merged.cat")
        rf = cat.arrayFromCol("Recovered_val")
        cf = cat.arrayFromCol("InputRadio_S_148")
        cf /=1000.
        re = cat.arrayFromCol("Recovered_err")
        # pylab.errorbar(cf[numpy.where(cf!=-.999)], rf[numpy.where(cf!=-.999)], yerr=re[numpy.where(cf!=-.999)], fmt='.')
        pylab.errorbar(cf[numpy.where(cf!=-.999)], rf[numpy.where(cf!=-.999)], yerr=None, fmt='.')
        pylab.gca().set_yscale('log')
        pylab.gca().set_xscale('log')
        pylab.plot([.0001,10],[.0001,10])
        utils.saveAndShow()

    def findNearestMatch( self, ra, dec, rad = None ):
        """
        Given a row/source, find the closest row/source in this catalog
        @return distance (deg), closest row
        """
        closestRow = self[0]
        distToClosestRow = distanceBetweenSources(closestRow, ra, dec )
        rowsWithinRadius = []
        if rad != None and distToClosestRow < rad:
            rowsWithinRadius.append([distToClosestRow, closestRow])
        for myrow in self[1:]:
            d = distanceBetweenSources(myrow, ra, dec)
            if d < distToClosestRow:
                distToClosestRow = d
                closestRow = myrow
            if rad != None and d < rad:
                rowsWithinRadius.append([d,myrow])
        return distToClosestRow, closestRow, rowsWithinRadius

def readFromASCIIFile( filename, colFile= None, binary = False ):
    if binary:
        pass

    if colFile == None:
        colFile = "%s.dict" % filename
    colDict = flipperDict.flipperDict()
    colDict.readFromFile(colFile)
    cat = catalog(cols = colDict['cols'])
    cat.read(filename)
    return cat

def readFromPickleFile( filename ):
    f = file(filename)
    cat = pickle.load(f)
    f.close()
    return cat

fromClusterCols = { \
        'ra'   : { 'desc':'Right Ascension of the object'   , 'order': 0, 'type': float, 'fmt' : '%10.3f', 'default': -999.}, \
        'dec'  : { 'desc':'Declination of the object'       , 'order': 1, 'type': float, 'fmt' : '%10.3f', 'default': -999.}, \
        'val'  : { 'desc':'val of the detection',           'order': 2, 'type': float, 'fmt' : '%10.3e', 'default': -999.}, \
        'err'  : { 'desc':'1-sigma error of the detection', 'order': 3, 'type': float, 'fmt' : '%10.3e', 'default': -999.}, \
        's/n'  : { 'desc':'signal to noise of the detection', 'order': 4, 'type': float, 'fmt' : '%10.3f', 'default': -999.}, \
        'npix' : { 'desc':'Number of pixels in cluster'     , 'order': 5, 'type': int,   'fmt' : '%10d  ', 'default': -999 }  \
        }
        
        

def catalogFromGroupList( clusterList ):
    cat = catalog(cols = fromClusterCols)
    for cl in clusterList:
        err = cl.valMax*(1./cl.s_nMax)
        cat.addRow( {'ra':cl.ra, 'dec':cl.dec, 's/n':cl.s_nMax, 'err':err, 'val':cl.valMax, 'npix':len(cl.pixels)} )
    return cat

def mergeCatalogs( cat1, name1, cat2, name2, dist=1.5, cat1Only = False ):
    """
    dist = distance between objects at which point you declare a match in arcmin
    """
    newCols = { }
    colNames1 = cat1.cols.keys()
    colNames2 = cat2.cols.keys()
    for colName in colNames1:
        if name1 != None:
            newCols["%s_%s" % (name1, colName)] = copy.deepcopy(cat1.cols[colName] )
        else:
            newCols[colName] = cat1.cols[colName] 
    for colName in colNames2:
        newCols["%s_%s" % (name2, colName)] = copy.deepcopy(cat2.cols[colName])
        newCols["%s_%s" % (name2, colName)]['order'] += len(colNames1)
    newCols['matched'] = {'desc':'match between catalogs','order': len(colNames1)+len(colNames2), 'type': float, 'fmt' : ' %s ', 'default': "N"}
    cat = catalog(cols = newCols)
    double2 = []
    for row1 in cat1:
        row = {}
        for colName in colNames1:
            if name1 != None:
                row["%s_%s" % (name1, colName)] = row1[colName]
            else:
                row[colName] = row1[colName]
        matches = False
        for row2 in cat2:
            cosdec = numpy.cos((row1['dec'] + row2['dec'])/2 * numpy.pi/180)
            xdiff = (row1['ra']-row2['ra'])*cosdec
            ydiff = row1['dec'] - row2['dec']
            d = numpy.sqrt(xdiff**2 + ydiff**2)
            if d < dist/60.:
                for colName in colNames2:
                    row["%s_%s" % (name2, colName)] = row2[colName]
                double2.append(row2)
                matches = True
                row['matched'] = 'Y'
        if not matches:
            for colName in colNames2:
                row["%s_%s" % (name2, colName)] = cat2.cols[colName]['default']
            row['matched'] = 'N'
        cat.addRow(row)
    if not cat1Only:
        for row2 in cat2:
            if row2 in double2:
                continue
            row = {}
            for colName in colNames1:
                row["%s_%s" % (name1, colName)] = cat1.cols[colName]['default']
            for colName in colNames2:
                row["%s_%s" % (name2, colName)] = row2[colName]
            row['matched'] = 'N'
            cat.addRow(row)
    return cat

def getDecimalRADecFromRow( row ):
    if 'ra' in row.keys():
        return row['ra'], row['dec']
    elif 'ra_s' in row.keys():
        return convertSexagesimalToRaDecDegrees( row['ra_s'], row['dec_s'] )
    else:
        raise ValueError('No key ra or ra_s in row: %s' % str(row))

def distanceBetweenSources(row1, ra2, dec2):
    """
    Compute distance in degrees between to catalog entries
    """
    ra1, dec1 = getDecimalRADecFromRow( row1 )
    cosdec = numpy.cos(numpy.pi/180. * (dec1+dec2)/2)
    dra = (ra1-ra2)*cosdec
    ddec = dec1-dec2
    return (dra**2+ddec**2)**0.5


def read(catname):
    try:
        cat = readFromPickleFile( catname )
    except:
        cat = readFromASCIIFile( catname )
    return cat
