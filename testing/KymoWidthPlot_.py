from ij import IJ;
from ij.gui import ProfilePlot
from ij.gui import Plot
from ij.measure import CurveFitter
from java.util import Arrays
from math import *

def main():
    imp = IJ.getImage();
    dim=imp.getDimensions();
    width = dim[0];
    height=dim[1];
    FWHM= list();
    intensity = list();
    for ii in range(0,height):
    	f,iVal,_ = widthresult=fitProfile(imp,width,ii);
    	FWHM.append(f);
        intensity.append(iVal);
    	
    print(FWHM);
    print(intensity)

def fitProfile(imp,width,plotpos_y):
    IJ.makeLine(0, plotpos_y+0.5, width, plotpos_y+0.5)
    pp = ProfilePlot(imp);
    y =pp.getProfile();
    #pp.createWindow();
    x = range(0,len(y));
    xavg = sum(x)/len(x);

    initguess = [min(y), max(y), xavg,initwidth,1]
    # The formula for a supergaussian
    superGaussian = "y = a + b * exp(-1*(pow ((x-c)*(x-c), e ) / pow( 2*d*d , e ) ) )"
    cf = CurveFitter(x,y);
    cf.doCustomFit(superGaussian,initguess,False)
    p=cf.getParams();
    xC = p[2];
    par_d = p[3];
    par_e = p[4];
    FWHM = 2*sqrt(2)*par_d*pow(log(2),1/(2*par_e));
    x0 = xC-FWHM/2;
    x1 = xC + FWHM/2;
    widthResult=[xC,x0,x1];
    intensity = max(y);
    return FWHM, intensity, widthResult;

main()
