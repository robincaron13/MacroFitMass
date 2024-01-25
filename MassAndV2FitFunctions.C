
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TMath.h"
#include "TString.h"

#include "TCanvas.h"
#include "TPaveText.h"
#include "TH1D.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TArrow.h"

#include <Riostream.h>

//_______________________________________________________________________________
// Inv Mass signal functions
//
//_______________________________________________________________________________
Double_t CrystalBall(Double_t *x,Double_t *par) {
    
    Double_t t = (x[0]-par[1])/par[2];
    if (par[3] < 0) t = -t;
    Double_t absAlpha = fabs((Double_t)par[3]);
    if (t >= -absAlpha) {
        return par[0]*(exp(-0.5*t*t));
    }
    else {
        Double_t a =  TMath::Power(par[4]/absAlpha,par[4])*exp(-0.5*absAlpha*absAlpha);
        Double_t b = par[4]/absAlpha - absAlpha;
        return par[0]*(a/TMath::Power(b - t, par[4]));
    }
}



Bool_t reject;
//_______________________________________________________________________________
Double_t CrystalBallExtended(Double_t *x,Double_t *par){
    
    Double_t t = (x[0]-par[1])/par[2];
    if (par[3] < 0) t = -t;
    
    Double_t absAlpha = fabs((Double_t)par[3]);
    Double_t absAlpha2 = fabs((Double_t)par[5]);
    
    if (t >= -absAlpha && t < absAlpha2){ // gaussian core
        return par[0]*(exp(-0.5*t*t));
    }
    
    if (t < -absAlpha){ //left tail
        Double_t a =  TMath::Power(par[4]/absAlpha,par[4])*exp(-0.5*absAlpha*absAlpha);
        Double_t b = par[4]/absAlpha - absAlpha;
        return par[0]*(a/TMath::Power(b - t, par[4]));
    }
    if (t >= absAlpha2){ //right tail
        Double_t c =  TMath::Power(par[6]/absAlpha2,par[6])*exp(-0.5*absAlpha2*absAlpha2);
        Double_t d = par[6]/absAlpha2 - absAlpha2;
        return  par[0]*(c/TMath::Power(d + t, par[6]));
    }
    return 0. ;
}

//_______________________________________________________________________________
Double_t FuncJpsi_NA60(Double_t *x, Double_t *par){
    
    Double_t s1=par[1];
//  Double_t s2=par[2];
    Double_t t1=par[3];
    Double_t t2=par[7];
    Double_t s00=par[2];
    Double_t t0=0.;
    Double_t t=(x[0]-s1)/s00;
    
    if (t >= t1 && t < t2) {
        t0 = 1;
    } else if (t < t1) {
        t0 = 1+TMath::Power(par[4]*(t1-t),par[5]-(par[6]*TMath::Sqrt(t1-t)));
    } else if (t >= t2) {
        t0 = 1+TMath::Power(par[8]*(t-t2),par[9]-(par[10]*TMath::Sqrt(t-t2)));
    }
    Double_t FitJPsi = par[0]*TMath::Exp(-0.5*((t/t0)*(t/t0)));
    
    return FitJPsi;
}

//_______________________________________________________________________________
Double_t MyGauss(Double_t *x,Double_t *par){
    Double_t val = TMath::Power( ((x[0]-par[1])/par[2]),2) ;
    return par[0]*( exp(-0.5 * val));
}


// New fit functions added
//
//
//
//
//
Double_t TwoCrystalBallExtended(Double_t *x,Double_t *par)
{
    Double_t mUpsi1 =  9.46030;  //
    Double_t mUpsi2 =  10.02326; //
    Double_t sig_fac16r2_1 = 1.0019;
    Double_t sig_fac16s2_1 = 1.03;
    Double_t sig_fac2_1 = sig_fac16s2_1;

    Double_t sum = 0;            //
    Double_t t = (x[0]-par[1])/par[2];
    if (par[3] < 0) t = -t;
    Double_t absAlpha = fabs((Double_t)par[3]);
    Double_t absAlpha2 = fabs((Double_t)par[5]);
    if (t < -absAlpha){ //left tail
        Double_t a =  TMath::Power(par[4]/absAlpha,par[4])*exp(-0.5*absAlpha*absAlpha);
        Double_t b = par[4]/absAlpha - absAlpha;
        sum += par[0]*(a/TMath::Power(b - t, par[4]));
    }
    else if (t >= -absAlpha && t < absAlpha2){ // gaussian core
        sum += par[0]*(exp(-0.5*t*t));
    }
    else if (t >= absAlpha2){ //right tail
        Double_t c =  TMath::Power(par[6]/absAlpha2,par[6])*exp(-0.5*absAlpha2*absAlpha2);
        Double_t d = par[6]/absAlpha2 - absAlpha2;
        sum += par[0]*(c/TMath::Power(d + t, par[6]));
    }
    else
        sum += 0;
    
    t = (x[0]-(par[1]*mUpsi2/mUpsi1))/(par[2]*sig_fac2_1);
    if (par[3] < 0) t = -t;
    
    absAlpha = fabs((Double_t)par[3]);
    absAlpha2 = fabs((Double_t)par[5]);
    
    if (t < -absAlpha) //left tail
    {
        Double_t a =  TMath::Power(par[4]/absAlpha,par[4])*exp(-0.5*absAlpha*absAlpha);
        Double_t b = par[4]/absAlpha - absAlpha;
        
        sum += par[7]*(a/TMath::Power(b - t, par[4]));
    }
    else if (t >= -absAlpha && t < absAlpha2) // gaussian core
    {
        sum += par[7]*(exp(-0.5*t*t));
    }
    else if (t >= absAlpha2) //right tail
    {
        Double_t c =  TMath::Power(par[6]/absAlpha2,par[6])*exp(-0.5*absAlpha2*absAlpha2);
        Double_t d = par[6]/absAlpha2 - absAlpha2;
        sum += par[7]*(c/TMath::Power(d + t, par[6]));
    } else
        sum += 0;
    
    return sum ;
}


Double_t ThreeIndCrystalBallExtended(Double_t *x,Double_t *par)
{
    return CrystalBallExtended(x,par) + CrystalBallExtended(x,&par[7]) + CrystalBallExtended(x,&par[14]);
}

Double_t TwoIndCrystalBallExtended(Double_t *x,Double_t *par)
{
    return CrystalBallExtended(x,par) + CrystalBallExtended(x,&par[7]);
}

Double_t ThreeCrystalBallExtended(Double_t *x,Double_t *par)
{
    
    Double_t mUpsi1 =  9.46030;
    Double_t mUpsi2 = 10.02326;
    Double_t mUpsi3 = 10.3552;
    
    Double_t sig_facPbPb2over1 = 1.03;
    Double_t sig_facPbPb3over1 = 1.06;
    
    Double_t sig_fac16r2_1 = 1.0019;
    Double_t sig_fac16r3_1 = 1.0304;
    
    Double_t sig_fac16s2_1 = 1.03501;
    Double_t sig_fac16s3_1 = 1.06615;
    Double_t sig_fac2_1 = sig_fac16s2_1;
    Double_t sig_fac3_1 = sig_fac16s3_1;

    Double_t sum = 0;
    Double_t t = (x[0]-par[1])/par[2];
    if (par[3] < 0) t = -t;
    
    Double_t absAlpha = fabs((Double_t)par[3]);
    Double_t absAlpha2 = fabs((Double_t)par[5]);
    
    if (t < -absAlpha) //left tail
    {
        Double_t a =  TMath::Power(par[4]/absAlpha,par[4])*exp(-0.5*absAlpha*absAlpha);
        Double_t b = par[4]/absAlpha - absAlpha;
        sum += par[0]*(a/TMath::Power(b - t, par[4]));
    } else if (t >= -absAlpha && t < absAlpha2) // gaussian core
    {
        sum += par[0]*(exp(-0.5*t*t));
    } else if (t >= absAlpha2) //right
    {
        Double_t c =  TMath::Power(par[6]/absAlpha2,par[6])*exp(-0.5*absAlpha2*absAlpha2);
        Double_t d = par[6]/absAlpha2 - absAlpha2;
        sum += par[0]*(c/TMath::Power(d + t, par[6]));
    } else
        sum += 0;
    
    // STATE 2S
    t = (x[0]-(par[1]*mUpsi2/mUpsi1))/(par[2]*sig_facPbPb2over1);
    if (par[3] < 0) t = -t;
    
    absAlpha = fabs((Double_t)par[3]);
    absAlpha2 = fabs((Double_t)par[5]);
    
    if (t < -absAlpha) //left tail
    {
        Double_t a =  TMath::Power(par[4]/absAlpha,par[4])*exp(-0.5*absAlpha*absAlpha);
        Double_t b = par[4]/absAlpha - absAlpha;
        sum += par[7]*(a/TMath::Power(b - t, par[4]));
    }
    else if (t >= -absAlpha && t < absAlpha2) // gaussian core
    {
        sum += par[7]*(exp(-0.5*t*t));
    }
    else if (t >= absAlpha2) //right tail
    {
        Double_t c =  TMath::Power(par[6]/absAlpha2,par[6])*exp(-0.5*absAlpha2*absAlpha2);
        Double_t d = par[6]/absAlpha2 - absAlpha2;
        sum += par[7]*(c/TMath::Power(d + t, par[6]));
    } else
        sum += 0;
    
    // STATE 3S
    t = (x[0]-(par[1]*mUpsi3/mUpsi1))/(par[2]*sig_facPbPb3over1);
    if (par[3] < 0) t = -t;
    
    absAlpha = fabs((Double_t)par[3]);
    absAlpha2 = fabs((Double_t)par[5]);
    
    if (t < -absAlpha) //left
    {
        Double_t a =  TMath::Power(par[4]/absAlpha,par[4])*exp(-0.5*absAlpha*absAlpha);
        Double_t b = par[4]/absAlpha - absAlpha;
        
        sum += par[8]*(a/TMath::Power(b - t, par[4]));
    }
    else if (t >= -absAlpha && t < absAlpha2) // gaussian core
    {
        sum += par[8]*(exp(-0.5*t*t));
    }
    else if (t >= absAlpha2) //right tail
    {
        Double_t c =  TMath::Power(par[6]/absAlpha2,par[6])*exp(-0.5*absAlpha2*absAlpha2);
        Double_t d = par[6]/absAlpha2 - absAlpha2;
        
        sum += par[8]*(c/TMath::Power(d + t, par[6]));
    } else
        sum += 0;
    
    return sum ;
}


Double_t RandomPhi(Double_t*x, Double_t *par){
    return (1 + 2*par[0]*cos(2*x[0]) );
}

Double_t TwoNA60(Double_t *x,Double_t *par)
{
    
    Double_t mUpsi1 =  9.46030;
    Double_t mUpsi2 = 10.02326;
    Double_t mUpsi3 = 10.3552;
    Double_t sig_fac16r2_1 = 1.04759;
    Double_t sig_fac16r3_1 = 1.07586;
    Double_t sig_fac16s2_1 = 1.008576;
    Double_t sig_fac16s3_1 = 1.07999;
    Double_t sig_fac2_1 = sig_fac16r2_1;
    Double_t sig_fac3_1 = sig_fac16r3_1;
    
    Double_t s1=par[1];
    //  Double_t s2=par[2];
    Double_t t1=par[3];
    Double_t t2=par[7];
    Double_t s00=par[2];
    Double_t t0=0.;
    Double_t t=(x[0]-s1)/s00;
    
    Double_t sum = 0;
    if (t >= t1 && t < t2) {
        t0 = 1;
    } else if (t < t1) {
        t0 = 1+TMath::Power(par[4]*(t1-t),par[5]-(par[6]*TMath::Sqrt(t1-t)));
    } else if (t >= t2) {
        t0 = 1+TMath::Power(par[8]*(t-t2),par[9]-(par[10]*TMath::Sqrt(t-t2)));
    }
    sum = par[0]*TMath::Exp(-0.5*((t/t0)*(t/t0)));
    
    
    t = (x[0]-(s1*mUpsi2/mUpsi1))/(s00*sig_fac2_1);
    if ( t >= t1 && t <= t2) {
        t0 = 1;
    } else if (t < t1) {
        t0 = 1+TMath::Power(par[4]*(t1-t),par[5]-(par[6]*TMath::Sqrt(t1-t)));
    }  else if (t >= t2) {
        t0 = 1+TMath::Power(par[8]*(t-t2),par[9]-(par[10]*TMath::Sqrt(t-t2)));
    }
    sum += par[11]*TMath::Exp(-0.5*((t/t0)*(t/t0)));
    
    return sum ;
}

Double_t ThreeNA60(Double_t *x,Double_t *par)
{
    
    Double_t mUpsi1 =  9.46030;
    Double_t mUpsi2 = 10.02326;
    Double_t mUpsi3 = 10.3552;
    
    Double_t sig_facPbPb2over1 = 1.03;
    Double_t sig_facPbPb3over1 = 1.06;
    
    Double_t sig_fac16r2_1 = 1.04759;
    Double_t sig_fac16r3_1 = 1.07586;
    Double_t sig_fac16s2_1 = 1.008576;
    Double_t sig_fac16s3_1 = 1.07999;
    Double_t sig_fac2_1 = sig_fac16s2_1;
    Double_t sig_fac3_1 = sig_fac16s3_1;
    
    Double_t s1=par[1];
    //  Double_t s2=par[2];
    Double_t t1=par[3];
    Double_t t2=par[7];
    Double_t s00=par[2];
    Double_t t0=0.;
    Double_t t=(x[0]-s1)/s00;
    
    Double_t sum = 0;
    if (t >= t1 && t < t2) {
        t0 = 1;
    } else if (t < t1) {
        t0 = 1+TMath::Power(par[4]*(t1-t),par[5]-(par[6]*TMath::Sqrt(t1-t)));
    } else if (t >= t2) {
        t0 = 1+TMath::Power(par[8]*(t-t2),par[9]-(par[10]*TMath::Sqrt(t-t2)));
    }
    sum = par[0]*TMath::Exp(-0.5*((t/t0)*(t/t0)));
    
    
    
    t = (x[0]-(s1*mUpsi2/mUpsi1))/(s00*sig_facPbPb2over1);
    if ( t >= t1 && t <= t2) {
        t0 = 1;
    } else if (t < t1) {
        t0 = 1+TMath::Power(par[4]*(t1-t),par[5]-(par[6]*TMath::Sqrt(t1-t)));
    } else if (t >= t2) {
        t0 = 1+TMath::Power(par[8]*(t-t2),par[9]-(par[10]*TMath::Sqrt(t-t2)));
    }
    sum += par[11]*TMath::Exp(-0.5*((t/t0)*(t/t0)));
    
    
    
    t = (x[0]-(s1*mUpsi3/mUpsi1))/(s00*sig_facPbPb3over1);
    if ( t>= t1 && t <= t2) {
        t0 = 1;
    } else if (t < t1) {
        t0 = 1+TMath::Power(par[4]*(t1-t),par[5]-(par[6]*TMath::Sqrt(t1-t)));
    } else if (t >= t2) {
        t0 = 1+TMath::Power(par[8]*(t-t2),par[9]-(par[10]*TMath::Sqrt(t-t2)));
    }
    sum += par[12]*TMath::Exp(-0.5*((t/t0)*(t/t0)));
    
    return sum ;
}


//_______________________________________________________________________________
//end of the new part of the fit code
//
//
//_______________________________________________________________________________
// Inv Mass background functions
//
//_______________________________________________________________________________
Double_t ExpBeforePeak(Double_t *x,Double_t *par){
    return par[0]*(exp(x[0]*par[1]));
}
//_______________________________________________________________________________
Double_t ExpAfterPeak(Double_t *x,Double_t *par){
    return par[0]*(exp(x[0]*par[1]));
}
//_______________________________________________________________________________
Double_t backFcn(Double_t *x,Double_t *par){
    if ( (reject && x[0] > 9.2 && x[0] < 10.4 ) | (reject && x[0] > 2.8 && x[0] < 3.4 ) ) {
        TF1::RejectPoint();
        return 0;
    }
    return ExpBeforePeak(x,par) + ExpAfterPeak(x,&par[2]);}
//_______________________________________________________________________________
Double_t backFcnPol(Double_t *x,Double_t *par){
    return par[0] + par[1]*x[0] + par[2]*x[0]*x[0] + par[3]*x[0]*x[0]*x[0];
}



//_______________________________________________________________________________
Double_t backFcnVWGaus(Double_t *x,Double_t *par){
    
    Double_t sigma = par[2] + par[3]*((x[0]-par[1])/par[1]);
    
    Double_t FitBck = par[0]*TMath::Exp(-(x[0]-par[1])*(x[0]-par[1])/(2.*sigma*sigma));
    
    if ( (reject && x[0] > 9.2 && x[0] < 10.4 ) | (reject && x[0] > 2.8 && x[0] < 3.4 ) ) {
        TF1::RejectPoint();
        return 0;
    }
    return FitBck;
}



TF1 * fitFunc;  // fit function pointer

const int NPAR = 2; // number of function parameters;

//____________________________________________________________________
double f(double * x, double * p) {
    // function used to fit the data
    return p[1]*TMath::Sin( p[0] * x[0] );
}

//____________________________________________________________________
double df_dPar(double * x, double * p) {
    // derivative of the function w.r..t parameters
    // use calculated derivatives from TF1::GradientPar
    
    double grad[NPAR];
    // p is used to specify for which parameter the derivative is computed
    int ipar = int(p[0] );
    assert (ipar >=0 && ipar < NPAR );
    
    assert(fitFunc);
    fitFunc->GradientPar(x, grad);
    
    return grad[ipar];
}

//____________________________________________________________________
double IntegralError(int npar, double * c, double * errPar,
                     double * covMatrix = 0) {
    // calculate the error on the integral given the parameter error and
    // the integrals of the gradient functions c[]
    
    double err2 = 0;
    for (int i = 0; i < npar; ++i) {
        if (covMatrix == 0) { // assume error are uncorrelated
            err2 += c[i] * c[i] * errPar[i] * errPar[i];
        } else {
            double s = 0;
            for (int j = 0; j < npar; ++j) {
                s += covMatrix[i*npar + j] * c[j];
            }
            err2 += c[i] * s;
        }
    }
    
    return TMath::Sqrt(err2);
}








//_______________________________________________________________________________
Double_t backFcnQVWGaus(Double_t *x,Double_t *par){
    Double_t sigma = par[2]+par[3]*((x[0]-par[1])/par[1])+par[4]*((x[0]-par[1])/par[1])*((x[0]-par[1])/par[1]);
    Double_t FitBck = par[0]*TMath::Exp(-(x[0]-par[1])*(x[0]-par[1])/(2.*sigma*sigma));
    return FitBck > 0 ? FitBck : 0;
}
//_______________________________________________________________________________
//Double_t backFcnPolR(Double_t *x,Double_t *par){
//	Double_t xm = x[0]-2.2;
//	Double_t np = 1e-2;
//	Double_t dx = 1e-1;
//	Double_t pol2 = 1 + par[1]*xm + par[2]*xm*xm;
//	pol2 = TMath::Max(pol2,0.);
//	Double_t pol3 = 1 + par[3]*xm + par[4]*xm*xm + par[5]*xm*xm*xm;
//	Double_t r = 0;
//	if (TMath::Abs(pol3)>np) {
//		r = par[0]*pol2/pol3;
//	} else {
//		Int_t ht = 3;
//		Int_t nit = 0;
//		r = 0;
//		xm -= (ht+1)*dx;
//		for (Int_t it=ht; it>=-ht; it--) {
//			xm += dx;
//			if (it==0) continue;
//			pol2 = 1 + par[1]*xm + par[2]*xm*xm;
//			pol3 = 1 + par[3]*xm + par[4]*xm*xm + par[5]*xm*xm*xm;
//			if (TMath::Abs(pol3)>np) {
//				nit++;
//				r += par[0]*pol2/pol3;
//			}
//		}
//		r /= nit;
//	}
//
//	return TMath::Max(r,1e-4);
//	//	return r;
//}
//_______________________________________________________________________________
Double_t backFcnPolR(Double_t *x,Double_t *par){
    Double_t xm = x[0]-2.2;
    Double_t np = 1e-1;
    Double_t nn = 1e-1;
    Double_t pol2 = 10 + par[1]*xm + par[2]*xm*xm;
    pol2 = TMath::Max(pol2,0.);
    //	Double_t pol3 = TMath::Abs(1 + par[3]*xm + par[4]*xm*xm + par[5]*xm*xm*xm);
    Double_t pol3 = TMath::Max(10 + par[3]*xm + par[4]*xm*xm + par[5]*xm*xm*xm,nn);
    
    //	Double_t r = 0.;
    //	if (TMath::Abs(pol3)<np) {
    //		pol3 = TMath::Power(TMath::Abs(pol3),1-(np-TMath::Abs(pol3))/np);
    ////		pol3 = TMath::Power(TMath::Abs(pol3),TMath::Abs(pol3));
    //	}
    //	pol3 = TMath::Sqrt(TMath::Power(pol3,2));
    //	pol3 = TMath::Sqrt(TMath::Abs(pol3));
    //	pol3 = TMath::Max(pol3,np);
    //	pol3 = TMath::Power(pol3,2);
    //	pol3 = TMath::Max(pol3,np);
    //	if (pol3<np) {
    //		return par[0];
    //	}
    //	pol3 += (pol3<np) ? (np-pol3) : 0.;
    Double_t r = par[0]*pol2/pol3;
    //	Double_t r = par[0]*TMath::Min(pol2/pol3,2.);
    //	cout << xm << " " << par[0] << " " << par[1] << " " << par[2] << " " << par[3] << " " << par[4] << " " << par[5] << " " << pol2 << " " << pol3 << " " << r << endl;
    //	return TMath::Max(TMath::Min(r,1.5*par[0]),0.);
    //	Double_t rr = r;
    //	if (r>1.5*par[0]) {
    //		rr = 2.0*par[0]-((r-2.0*par[0])*(r-2.0*par[0])/(0.5*par[0]));
    //	}
    //	rr = TMath::Min(rr,r);
    //	cout << xm << " " << pol2 << " " << pol3 << " " << r << " " << rr << " " << par[0] << endl;
    return TMath::Max(r,1e-4);
    //	return r;
}
//_______________________________________________________________________________
//Double_t backFcnPolR(Double_t *x,Double_t *par){
//	Double_t xm = x[0]-2.2;
//	Double_t np = 1e-1;
//	Double_t nn = 1e-3;
//	Double_t pol2 = 1 + par[1]*xm + par[2]*xm*xm;
//	pol2 = TMath::Max(pol2,0.);
//	//	Double_t pol3 = TMath::Abs(1 + par[3]*xm + par[4]*xm*xm + par[5]*xm*xm*xm);
//	Double_t pol3 = TMath::Max(1 + par[3]*xm + par[4]*xm*xm + par[5]*xm*xm*xm,nn);
//
//	//	Double_t r = 0.;
//	if (TMath::Abs(pol3)<np) {
//		pol3 = TMath::Power(TMath::Abs(pol3),1-(np-TMath::Abs(pol3))/np);
//		//		pol3 = TMath::Power(TMath::Abs(pol3),TMath::Abs(pol3));
//	}
//	//	pol3 = TMath::Sqrt(TMath::Power(pol3,2));
//	//	pol3 = TMath::Sqrt(TMath::Abs(pol3));
//	//	pol3 = TMath::Max(pol3,np);
//	//	pol3 = TMath::Power(pol3,2);
//	//	pol3 = TMath::Max(pol3,np);
//	//	if (pol3<np) {
//	//		return par[0];
//	//	}
//	//	pol3 += (pol3<np) ? (np-pol3) : 0.;
//	Double_t r = par[0]*pol2/pol3;
//	//	Double_t r = par[0]*TMath::Min(pol2/pol3,2.);
//	//	cout << xm << " " << par[0] << " " << par[1] << " " << par[2] << " " << par[3] << " " << par[4] << " " << par[5] << " " << pol2 << " " << pol3 << " " << r << endl;
//	//	return TMath::Max(TMath::Min(r,1.5*par[0]),0.);
//	//	Double_t rr = r;
//	//	if (r>1.5*par[0]) {
//	//		rr = 2.0*par[0]-((r-2.0*par[0])*(r-2.0*par[0])/(0.5*par[0]));
//	//	}
//	//	rr = TMath::Min(rr,r);
//	//	cout << xm << " " << pol2 << " " << pol3 << " " << r << " " << rr << " " << par[0] << endl;
//	return TMath::Max(r,1e-4);
//	//	return r;
//}
//_______________________________________________________________________________
Double_t backFcnCheb1(Double_t *x,Double_t *par){
    Double_t xmin = 2;
    Double_t xmax = 5;
    Double_t xx = (2.0 * x[0] - xmin - xmax)/(xmax-xmin);
    const int order = 5;
    Double_t fT[order+1] = {0};
    if (order == 1) return par[0];
    if (order == 2) return par[0] + xx*par[1];
    //build the polynomials
    fT[0] = 1;
    fT[1] = xx;
    for (int i=1; i<order; ++i)
        fT[i+1] = 2*xx*fT[i] - fT[i-1];
    
    Double_t sum;
    sum = par[0]*fT[0];
    //  double sum = fT[0];
    for (int i=1; i<= order; ++i)
        sum += par[i]*fT[i];
    
    //  if (sum < 1 && sum > 0) ? TMath::Sqrt(sum):sum;
    //  return sum;
    return TMath::Max(sum,0.0001);
    //  return par[0]*sum;
}

//_______________________________________________________________________________
Double_t backFcnChebJpsi(Double_t *x,Double_t *par){
    Double_t xmin = 2;
    Double_t xmax = 5;
    Double_t xx = (2.0 * x[0] - xmin - xmax)/(xmax-xmin);
    const int order = 3;
    Double_t fT[order+1] = {0};
    if (order == 1) return par[0];
    if (order == 2) return par[0] + xx*par[1];
    //build the polynomials
    fT[0] = 1;
    fT[1] = xx;
    for (int i=1; i<order; ++i)
        fT[i+1] = 2*xx*fT[i] - fT[i-1];
    
    Double_t sum;
    sum = par[0]*fT[0];
    //  double sum = fT[0];
    for (int i=1; i<= order; ++i)
        sum += par[i]*fT[i];
    
    //  if (sum < 1 && sum > 0) ? TMath::Sqrt(sum):sum;
    //  return sum;
    return TMath::Max(sum,0.0001);
    //  return par[0]*sum;
}

Double_t backFcnCheb(Double_t *x,Double_t *par){
    Double_t xmin = 7;
    Double_t xmax = 13;
    Double_t xx = (2.0 * x[0] - xmin - xmax)/(xmax-xmin);
    const int order = 5;
    Double_t fT[order+1] = {0};
    if (order == 1) return par[0];
    if (order == 2) return par[0] + xx*par[1];
    //build the polynomials
    fT[0] = 1;
    fT[1] = xx;
    for (int i=1; i<order; ++i)
        fT[i+1] = 2*xx*fT[i] - fT[i-1];
    
    Double_t sum;
    sum = par[0]*fT[0];
    //  double sum = fT[0];
    for (int i=1; i<= order; ++i)
        sum += par[i]*fT[i];
    
    if (reject && x[0] > 9.2 && x[0] < 10.4) {
        TF1::RejectPoint();
        return 0;
    }
    //  if (sum < 1 && sum > 0) ? TMath::Sqrt(sum):sum;
    //  return sum;
    return TMath::Max(sum,0.0001);
    //  return par[0]*sum;
}
//_______________________________________________________________________________
// Inv Mass signal + background functions
//
//_______________________________________________________________________________
Double_t TotalG(Double_t *x,Double_t *par){
    return ExpBeforePeak(x,par) + ExpAfterPeak(x,&par[2]) + MyGauss(x,&par[4]);
}

//_______________________________________________________________________________
Double_t TotalCB(Double_t *x,Double_t *par){
    return ExpBeforePeak(x,par) + ExpAfterPeak(x,&par[2]) + CrystalBall(x,&par[4]);
}
//_______________________________________________________________________________
Double_t TotalCBExt(Double_t *x,Double_t *par){
    return ExpBeforePeak(x,par) + ExpAfterPeak(x,&par[2]) + CrystalBallExtended(x,&par[4]);
}
//_______________________________________________________________________________
Double_t TotalTwoCBExt(Double_t *x,Double_t *par){
    return backFcn(x,par) + TwoCrystalBallExtended(x,&par[4]);
}
//_______________________________________________________________________________
Double_t TotalTwoIndCBExtFcn(Double_t *x,Double_t *par){
    return backFcn(x,par) + TwoIndCrystalBallExtended(x,&par[4]);
}
//_______________________________________________________________________________
Double_t TotalTwoNA60Exp(Double_t *x,Double_t *par){
    return ExpBeforePeak(x,par) + ExpAfterPeak(x,&par[2]) + TwoNA60(x,&par[4]);
}
//_______________________________________________________________________________
Double_t TotalTwoIndCBExt(Double_t *x,Double_t *par){
    return ExpBeforePeak(x,par) + ExpAfterPeak(x,&par[2]) + TwoIndCrystalBallExtended(x,&par[4]);
}
//_______________________________________________________________________________
Double_t TotalThreeCBExt(Double_t *x,Double_t *par){
    return backFcn(x,par) + ThreeCrystalBallExtended(x,&par[4]);
}
//_______________________________________________________________________________
Double_t TotalThreeNA60Exp(Double_t *x,Double_t *par){
    return ExpBeforePeak(x,par) + ExpAfterPeak(x,&par[2]) + ThreeNA60(x,&par[4]);
}
//_______________________________________________________________________________
Double_t TotalThreeIndCBExt(Double_t *x,Double_t *par){
    return ExpBeforePeak(x,par) + ExpAfterPeak(x,&par[2]) + ThreeIndCrystalBallExtended(x,&par[4]);
}
//_______________________________________________________________________________
Double_t TotalCBPol(Double_t *x,Double_t *par){
    return backFcnPol(x,par) + CrystalBall(x,&par[4]);
}
//_______________________________________________________________________________
Double_t TotalCBVWGaus(Double_t *x,Double_t *par){
    return backFcnVWGaus(x,par) + CrystalBall(x,&par[4]);
}
//_______________________________________________________________________________
Double_t TotalCBExtPol(Double_t *x,Double_t *par){
    return backFcnPol(x,par) + CrystalBallExtended(x,&par[4]);
}
//_______________________________________________________________________________
Double_t TotalCBExtVWGaus(Double_t *x,Double_t *par){
    return backFcnVWGaus(x,par) + CrystalBallExtended(x,&par[4]);
}
//_______________________________________________________________________________
Double_t TotalTwoIndCBExtVWGaus(Double_t *x,Double_t *par){
    return backFcnVWGaus(x,par) + TwoIndCrystalBallExtended(x,&par[4]);
}
//_______________________________________________________________________________
Double_t TotalTwoIndCBExtCheb(Double_t *x,Double_t *par){
    return backFcnChebJpsi(x,par) + TwoIndCrystalBallExtended(x,&par[4]);
}

//_______________________________________________________________________________
Double_t TotalTwoCBExtVWGaus(Double_t *x,Double_t *par){
    return backFcnVWGaus(x,par) + TwoCrystalBallExtended(x,&par[4]);
}
//_______________________________________________________________________________
Double_t TotalTwoNA60VWGaus(Double_t *x,Double_t *par){
    return backFcnVWGaus(x,par) + TwoNA60(x,&par[4]);
}
//_______________________________________________________________________________
Double_t TotalThreeCBExtVWGaus(Double_t *x,Double_t *par){
    return backFcnVWGaus(x,par) + ThreeCrystalBallExtended(x,&par[4]);
}
Double_t TotalThreeIndCBExtVWGaus(Double_t *x,Double_t *par){
    return backFcnVWGaus(x,par) + ThreeIndCrystalBallExtended(x,&par[4]);
}
Double_t TotalThreeCBExtVWGaus2(Double_t *x,Double_t *par){
    return backFcnVWGaus(x,par) + ThreeCrystalBallExtended(x,&par[7]);
}
Double_t TotalThreeIndCBExtVWGaus2(Double_t *x,Double_t *par){
    return backFcnVWGaus(x,par) + ThreeIndCrystalBallExtended(x,&par[7]);
}

//_______________________________________________________________________________
Double_t TotalThreeNA60VWGaus(Double_t *x,Double_t *par){
    return backFcnVWGaus(x,par) + ThreeNA60(x,&par[4]);
}
Double_t TotalThreeNA60VWGaus2(Double_t *x,Double_t *par){
    return backFcnVWGaus(x,par) + ThreeNA60(x,&par[7]);
}






//_______________________________________________________________________________
Double_t TotalCBQVWGaus(Double_t *x,Double_t *par){
    return backFcnQVWGaus(x,par) + CrystalBall(x,&par[5]);
}
//_______________________________________________________________________________
Double_t TotalCBExtQVWGaus(Double_t *x,Double_t *par){
    return backFcnQVWGaus(x,par) + CrystalBallExtended(x,&par[5]);
}
//_______________________________________________________________________________
Double_t TotalNA60QVWGaus(Double_t *x,Double_t *par){
    return backFcnQVWGaus(x,par) + FuncJpsi_NA60(x,&par[5]);
}
//_______________________________________________________________________________
Double_t TotalCBExtPolR(Double_t *x,Double_t *par){
    return backFcnPolR(x,par) + CrystalBallExtended(x,&par[6]);
}
//_______________________________________________________________________________
Double_t TotalNA60PolR(Double_t *x,Double_t *par){
    return backFcnPolR(x,par) + FuncJpsi_NA60(x,&par[6]);
}
//_______________________________________________________________________________
Double_t TotalCBExtCheb(Double_t *x,Double_t *par){
    return backFcnCheb(x,par) + CrystalBallExtended(x,&par[7]);
}
//_______________________________________________________________________________
Double_t TotalTwoCBExtCheb(Double_t *x,Double_t *par){
    return backFcnCheb(x,par) + TwoCrystalBallExtended(x,&par[6]);
}


//_______________________________________________________________________________
Double_t TotalTwoNA60Cheb(Double_t *x,Double_t *par){
    return backFcnCheb(x,par) + TwoNA60(x,&par[6]);
}

//_______________________________________________________________________________
Double_t TotalThreeCBExtCheb(Double_t *x, Double_t *par){
    return backFcnCheb(x,par) + ThreeCrystalBallExtended(x,&par[6]);
}

Double_t TotalThreeCBExtCheb2(Double_t *x, Double_t *par){
    return backFcnCheb(x,par) + ThreeIndCrystalBallExtended(x,&par[7]);
}

Double_t TotalThreeIndCB2Exp(Double_t *x, Double_t *par){
    return backFcn(x,par) + ThreeIndCrystalBallExtended(x,&par[4]);
}


Double_t TotalThreeCBExtChebJpsi(Double_t *x, Double_t *par){
    return backFcnCheb1(x,par) + ThreeIndCrystalBallExtended(x,&par[6]);
}
//_______________________________________________________________________________
Double_t TotalThreeNA60Cheb(Double_t *x, Double_t *par){
    return backFcnCheb(x,par) + ThreeNA60(x,&par[6]);
}

//_______________________________________________________________________________
Double_t TotalNA60Cheb(Double_t *x,Double_t *par){
    return backFcnCheb(x,par) + FuncJpsi_NA60(x,&par[7]);
}



//_______________________________________________________________________________
Double_t MyDoubleGauss(Double_t *x,Double_t *par){
    return MyGauss(x,par) + MyGauss(x,&par[3]);
}




// v2 vs. Mass background functions
//_______________________________________________________________________________
Double_t v2BackgroundExp(Double_t *x, Double_t*par){
    return exp(par[0] + par[1]*x[0]); //par[0]*x[0] + par[1];
}
//_______________________________________________________________________________
//Double_t v2Background(Double_t *x, Double_t*par){
//  return par[0]*x[0]*x[0] + par[1]*x[0] + par[2];
//}
// v2 background functions pol2 historical
//_______________________________________________________________________________
Double_t v2Background(Double_t *x, Double_t*par){
    return par[0]*((x[0]-3.1)*(x[0]-3.1)) + par[1]*(x[0]-3.1) + par[2];
}
// v2 background functions pol2
//_______________________________________________________________________________
Double_t v2BackgroundPol2(Double_t *x, Double_t*par){
    return par[0] + par[1]*(x[0]-3.1) + par[2]*((x[0]-3.1)*(x[0]-3.1));
}
// v2 background functions pol3
//_______________________________________________________________________________
Double_t v2BackgroundPol3(Double_t *x, Double_t*par){
    return par[0] + par[1]*(x[0]-3.1) + par[2]*((x[0]-3.1)*(x[0]-3.1)) + par[3]*((x[0]-3.1)*(x[0]-3.1)*(x[0]-3.1));
}
Double_t v2BackgroundPol3Y(Double_t *x, Double_t*par){
    return par[0] + par[1]*(x[0]-9.46) + par[2]*((x[0]-9.46)*(x[0]-9.46)) + par[3]*((x[0]-9.46)*(x[0]-9.46)*(x[0]-9.46))+ par[4]*((x[0]-9.46)*(x[0]-9.46)*(x[0]-9.46)*(x[0]-9.46));
}
// v2 background functions pol4
//_______________________________________________________________________________
Double_t v2BackgroundPol4(Double_t *x, Double_t*par){
    return par[0] + par[1]*(x[0]-3.1) + par[2]*((x[0]-3.1)*(x[0]-3.1)) + par[3]*((x[0]-3.1)*(x[0]-3.1)*(x[0]-3.1)) + par[4]*((x[0]-3.1)*(x[0]-3.1)*(x[0]-3.1)*(x[0]-3.1));
}
Double_t v2BackgroundPol4Y(Double_t *x, Double_t*par){
    /*
    if (reject && x[0] > 9.4 ) {
        TF1::RejectPoint();
        return 0;
    }
     */
    //return par[0] + par[1]*(x[0]-9.46) + par[2]*((x[0]-9.46)*(x[0]-9.46)) + par[3]*((x[0]-9.46)*(x[0]-9.46)*(x[0]-9.46))+ par[4]*((x[0]-9.46)*(x[0]-9.46)*(x[0]-9.46)*(x[0]-9.46));
    return par[0] + par[1]*x[0] + par[2]*x[0]*x[0] + par[3]*x[0]*x[0]*x[0] + par[4]*x[0]*x[0]*x[0]*x[0];
}
//_______________________________________________________________________________
Double_t v2BackgroundCheb(Double_t *x,Double_t *par){
    //Double_t xmin = 2.2;
    //Double_t xmax = 4.7;
    Double_t xmin = 8;
    Double_t xmax = 11;
    double xx = (2.0 * x[0]-xmin-xmax)/(xmax-xmin);
    const int order = 4;
    Double_t fT[order+1] = {0};
    if (order == 1) return par[0];
    if (order == 2) return par[0] + xx*par[1];
    // build the polynomials
    fT[0] = 1;
    fT[1] = xx;
    for (int i = 1; i< order; ++i) {
        fT[i+1] =  2 *xx * fT[i] - fT[i-1];
    }
    double sum = par[0]*fT[0];
    for (int i = 1; i<= order; ++i) {
        sum += par[i] * fT[i];
    }
    if (reject && x[0] > 9.1 && x[0] < 9.8) {
        TF1::RejectPoint();
        return 0;
    }
    return sum;
}

//_______________________________________________________________________________
Double_t v2BackgroundChebJpsi(Double_t *x,Double_t *par){
    Double_t xmin = 2.2;
    Double_t xmax = 4.7;
    //Double_t xmin = 8;
    //Double_t xmax = 11;
    double xx = (2.0 * x[0]-xmin-xmax)/(xmax-xmin);
    const int order = 4;
    Double_t fT[order+1] = {0};
    if (order == 1) return par[0];
    if (order == 2) return par[0] + xx*par[1];
    // build the polynomials
    fT[0] = 1;
    fT[1] = xx;
    for (int i = 1; i< order; ++i) {
        fT[i+1] =  2 *xx * fT[i] - fT[i-1];
    }
    double sum = par[0]*fT[0];
    for (int i = 1; i<= order; ++i) {
        sum += par[i] * fT[i];
    }
    if (reject && x[0] > 2.9 && x[0] < 3.3) {
        TF1::RejectPoint();
        return 0;
    }
    return sum;
}

//_______________________________________________________________________________
Double_t ResolutionPol6(Double_t *x, Double_t*par){
    return par[4]*(x[0]-par[5])*(x[0]-par[5])*(x[0]-par[5])*(x[0]-par[5])*(x[0]-par[5])*(x[0]-par[5]) + par[3]*(x[0]-par[5])*(x[0]-par[5])*(x[0]-par[5])*(x[0]-par[5]) + par[2]*(x[0]-par[5])*(x[0]-par[5]) + par[1]*(x[0]-par[5]) + par[0];
}

//_______________________________________________________________________________
Double_t ResolutionPol2(Double_t *x, Double_t*par){
    return  par[2]*(x[0]-par[3])*(x[0]-par[3]) + par[1]*(x[0]-par[3]) + par[0];
}
//_______________________________________________________________________________
Double_t ResolutionPol1(Double_t *x, Double_t*par){
    return  par[1]*x[0] + par[0];
}
//_______________________________________________________________________________
Double_t ResolutionPol0(Double_t *x, Double_t*par){
    return   par[0];
}


Double_t v2BackgroundChebY(Double_t *x,Double_t *par){
    Double_t xmin = 8;
    Double_t xmax = 11;
    double xx = (2.0 * x[0]-xmin-xmax)/(xmax-xmin);
    const int order = 4;
    Double_t fT[order+1] = {0};
    if (order == 1) return par[0];
    if (order == 2) return par[0] + xx*par[1];
    // build the polynomials
    fT[0] = 1;
    fT[1] = xx;
    for (int i = 1; i< order; ++i) {
        fT[i+1] =  2 *xx * fT[i] - fT[i-1];
    }
    double sum = par[0]*fT[0];
    for (int i = 1; i<= order; ++i) {
        sum += par[i] * fT[i];
    }
    if (reject && x[0] > 9.2 && x[0] < 10.4) {
        TF1::RejectPoint();
        return 0;
    }
    return sum;
}
// v2 background pol*exp
//_______________________________________________________________________________
Double_t v2BackgroundPolExp(Double_t *x, Double_t *par){
    //	return (par[2] + par[3]*(x[0]-3.1) + par[4]*(x[0]-3.1)*(x[0]-3.1))*exp(par[0]+par[1]*(x[0]-3.1)); //par[0]*x[0] + par[1];
    return par[0]*(1./exp(par[1]) + par[3]*(x[0]-3.1) + par[4]*(x[0]-3.1)*(x[0]-3.1))*exp(par[1]+par[2]*(x[0]-3.1)); //par[0]*x[0] + par[1];
}
Double_t v2BackgroundPolExpY(Double_t *x, Double_t *par){
    return par[0]*(1./exp(par[1]) + par[3]*(x[0]-9.5) + par[4]*(x[0]-9.5)*(x[0]-9.5) )*exp(par[1]+par[2]*(x[0]-9.5));
}


Double_t v2BackgroundPol3Exp(Double_t *x, Double_t *par){
    if (reject && x[0] > 9.2 && x[0] < 10.4) {
        TF1::RejectPoint();
        return 0;
    }
    return exp(par[0] + par[1]*x[0])*(par[2] + par[3]*x[0] + par[4]*x[0]*x[0] );
}

////_______________________________________________________________________________
//Double_t v2BackgroundPolExp(Double_t *x, Double_t *par){
//	return (par[2] + par[3]*(x[0]-3.1) + par[4]*(x[0]-3.1)*(x[0]-3.1) + par[5]*(x[0]-3.1)*(x[0]-3.1)*(x[0]-3.1))*exp(par[0]+par[1]*(x[0]-3.1)); //par[0]*x[0] + par[1];
//}

//_______________________________________________________________________________
//Double_t v2BackgroundPolExp(Double_t *x, Double_t*par){
//	return (par[2] + par[3]*(x[0]) + par[4]*(x[0])*(x[0]) + par[5]*(x[0])*(x[0])*(x[0]))*exp(par[0]+par[1]*(x[0])); //par[0]*x[0] + par[1];
//}
//_______________________________________________________________________________
//Double_t v2BackgroundPolExp(Double_t *x, Double_t*par){
//	return (par[3] + par[4]*(x[0]) + par[5]*(x[0])*(x[0]))*exp(par[0]+par[1]*(x[0])+par[2]*(x[0])*(x[0])); //par[0]*x[0] + par[1];
//}
//_______________________________________________________________________________
// v2 vs. Mass signal + background functions
//
//_______________________________________________________________________________
Double_t v2vsMassCBVWG(Double_t*x, Double_t *par){ // bkg 4 + sig 5 + flow bkg 3 + flow sig 1= 13
    return (backFcnVWGaus(x, par)*v2Background(x, &par[9]) + par[12]*CrystalBall(x, &par[4]))/(backFcnVWGaus(x, par) + CrystalBall(x, &par[4]));
}
//_______________________________________________________________________________
Double_t v2vsMassCBPol(Double_t*x, Double_t *par){ // bkg 4 + sig 5 + flow bkg 3 + flow sig 1= 13
    return (backFcnPol(x, par)*v2Background(x, &par[9]) + par[12]*CrystalBall(x, &par[4])) /(backFcnPol(x, par) + CrystalBall(x, &par[4]));
}
//_______________________________________________________________________________
Double_t v2vsMassCB2VWG(Double_t*x, Double_t *par){ // bkg 4 + sig 7 + flow bkg 3 + flow sig 1= 15
    return (backFcnVWGaus(x, par)*v2Background(x, &par[11]) + par[14]*CrystalBallExtended(x, &par[4]))/(backFcnVWGaus(x, par) + CrystalBallExtended(x, &par[4]));
}
//_______________________________________________________________________________
Double_t v2vsMassCB2Pol(Double_t*x, Double_t *par){ // bkg 4 + sig 7 + flow bkg 3 + flow sig 1= 15
    return (backFcnPol(x, par)*v2Background(x, &par[11]) + par[14]*CrystalBallExtended(x, &par[4])) /(backFcnPol(x, par) + CrystalBallExtended(x, &par[4]));
}
//_______________________________________________________________________________
Double_t v2vsMassCB2QVWG(Double_t*x, Double_t *par){ // bkg 5 + sig 7 + flow bkg 3 + flow sig 1 = 16
    return (backFcnQVWGaus(x, par)*v2Background(x, &par[12]) + par[15]*CrystalBallExtended(x, &par[5]))/(backFcnQVWGaus(x, par) + CrystalBallExtended(x, &par[5]));
}
//_______________________________________________________________________________
Double_t v2vsMassCB2PolR(Double_t*x, Double_t *par){ // bkg 6 + sig 7 + flow bkg 3 + flow sig 1 = 17
    return (backFcnPolR(x, par)*v2Background(x, &par[13]) + par[16]*CrystalBallExtended(x, &par[6])) /(backFcnPolR(x, par) + CrystalBallExtended(x, &par[6]));
}
//_______________________________________________________________________________
Double_t v2vsMassNA60QVWG(Double_t*x, Double_t *par){ // bkg 5 + sig 11 + flow bkg 3 + flow sig 1 = 20
    return (backFcnQVWGaus(x, par)*v2Background(x, &par[16]) + par[19]*FuncJpsi_NA60(x, &par[5]))/(backFcnQVWGaus(x, par) + FuncJpsi_NA60(x, &par[5]));
}
//_______________________________________________________________________________
Double_t v2vsMassNA60PolR(Double_t*x, Double_t *par){ // bkg 6 + sig 7 + flow bkg 3 + flow sig 1 = 17
    return (backFcnPolR(x, par)*v2Background(x, &par[17]) + par[20]*FuncJpsi_NA60(x, &par[6])) /(backFcnPolR(x, par) + FuncJpsi_NA60(x, &par[6]));
}

//_______________________________________________________________________________
Double_t v2vsMassCB2QVWGExp(Double_t*x, Double_t *par){ // bkg 5 + sig 7 + flow bkg 3 + flow sig 1 = 16
    return (backFcnQVWGaus(x, par)*v2BackgroundExp(x, &par[12]) + par[15]*CrystalBallExtended(x, &par[5]))/(backFcnQVWGaus(x, par) + CrystalBallExtended(x, &par[5]));
}
//_______________________________________________________________________________
Double_t v2vsMassCB2PolRExp(Double_t*x, Double_t *par){ // bkg 6 + sig 7 + flow bkg 3 + flow sig 1 = 17
    return (backFcnPolR(x, par)*v2BackgroundExp(x, &par[13]) + par[16]*CrystalBallExtended(x, &par[6])) /(backFcnPolR(x, par) + CrystalBallExtended(x, &par[6]));
}
//_______________________________________________________________________________
Double_t v2vsMassCB2ChebExp(Double_t*x, Double_t *par){ // bkg 7 + sig 7 + flow bkg 3 + flow sig 1 = 18
    return (backFcnCheb(x, par)*v2BackgroundExp(x, &par[14]) + par[17]*CrystalBallExtended(x, &par[7]))/(backFcnCheb(x, par) + CrystalBallExtended(x, &par[7]));
}
//_______________________________________________________________________________
Double_t v2vsMassNA60QVWGExp(Double_t*x, Double_t *par){ // bkg 5 + sig 11 + flow bkg 3 + flow sig 1 = 20
    return (backFcnQVWGaus(x, par)*v2BackgroundExp(x, &par[16]) + par[19]*FuncJpsi_NA60(x, &par[5]))/(backFcnQVWGaus(x, par) + FuncJpsi_NA60(x, &par[5]));
}
//_______________________________________________________________________________
Double_t v2vsMassNA60PolRExp(Double_t*x, Double_t *par){ // bkg 6 + sig 7 + flow bkg 3 + flow sig 1 = 17
    return (backFcnPolR(x, par)*v2BackgroundExp(x, &par[17]) + par[20]*FuncJpsi_NA60(x, &par[6])) /(backFcnPolR(x, par) + FuncJpsi_NA60(x, &par[6]));
}
//_______________________________________________________________________________
Double_t v2vsMassNA60ChebExp(Double_t*x, Double_t *par){ // bkg 7 + sig 11 + flow bkg 3 + flow sig 1 = 22
    return (backFcnCheb(x, par)*v2BackgroundExp(x, &par[18]) + par[21]*FuncJpsi_NA60(x, &par[7]))/(backFcnCheb(x, par) + FuncJpsi_NA60(x, &par[7]));
}

//_______________________________________________________________________________
Double_t v2vsMassCB2QVWGPol2(Double_t*x, Double_t *par){ // bkg 5 + sig 7 + flow bkg 3 + flow sig 1 = 16
    return (backFcnQVWGaus(x, par)*v2BackgroundPol2(x, &par[12]) + par[15]*CrystalBallExtended(x, &par[5]))/(backFcnQVWGaus(x, par) + CrystalBallExtended(x, &par[5]));
}
//_______________________________________________________________________________
Double_t v2vsMassCB2PolRPol2(Double_t*x, Double_t *par){ // bkg 6 + sig 7 + flow bkg 3 + flow sig 1 = 17
    return (backFcnPolR(x, par)*v2BackgroundPol2(x, &par[13]) + par[16]*CrystalBallExtended(x, &par[6])) /(backFcnPolR(x, par) + CrystalBallExtended(x, &par[6]));
}
//_______________________________________________________________________________
Double_t v2vsMassCB2ChebPol2(Double_t*x, Double_t *par){ // bkg 7 + sig 7 + flow bkg 3 + flow sig 1 = 18
    return (backFcnCheb(x, par)*v2BackgroundPol2(x, &par[14]) + par[17]*CrystalBallExtended(x, &par[7]))/(backFcnCheb(x, par) + CrystalBallExtended(x, &par[7]));
}
//_______________________________________________________________________________
Double_t v2vsMassNA60QVWGPol2(Double_t*x, Double_t *par){ // bkg 5 + sig 11 + flow bkg 3 + flow sig 1 = 20
    return (backFcnQVWGaus(x, par)*v2BackgroundPol2(x, &par[16]) + par[19]*FuncJpsi_NA60(x, &par[5]))/(backFcnQVWGaus(x, par) + FuncJpsi_NA60(x, &par[5]));
}
//_______________________________________________________________________________
Double_t v2vsMassNA60PolRPol2(Double_t*x, Double_t *par){ // bkg 6 + sig 7 + flow bkg 3 + flow sig 1 = 17
    return (backFcnPolR(x, par)*v2BackgroundPol2(x, &par[17]) + par[20]*FuncJpsi_NA60(x, &par[6])) /(backFcnPolR(x, par) + FuncJpsi_NA60(x, &par[6]));
}
//_______________________________________________________________________________
Double_t v2vsMassNA60ChebPol2(Double_t*x, Double_t *par){ // bkg 7 + sig 11 + flow bkg 3 + flow sig 1 = 22
    return (backFcnCheb(x, par)*v2BackgroundPol2(x, &par[18]) + par[21]*FuncJpsi_NA60(x, &par[7]))/(backFcnCheb(x, par) + FuncJpsi_NA60(x, &par[7]));
}
//_______________________________________________________________________________
Double_t v2vsMassCB2QVWGPol3(Double_t*x, Double_t *par){ // bkg 5 + sig 7 + flow bkg 4 + flow sig 1 = 17
    return (backFcnQVWGaus(x, par)*v2BackgroundPol3(x, &par[12]) + par[16]*CrystalBallExtended(x, &par[5]))/(backFcnQVWGaus(x, par) + CrystalBallExtended(x, &par[5]));
}

//_______________________________________________________________________________
Double_t v2vsMassCB2PolRPol3(Double_t*x, Double_t *par){ // bkg 6 + sig 7 + flow bkg 4 + flow sig 1 = 18
    return (backFcnPolR(x, par)*v2BackgroundPol3(x, &par[13]) + par[17]*CrystalBallExtended(x, &par[6])) /(backFcnPolR(x, par) + CrystalBallExtended(x, &par[6]));
}




Double_t HarmonicFlowComponents(Double_t*x, Double_t *par){
    return par[0]*( 1 + 2*( par[1]*cos(2*(x[0])) ))  ;
}




//_______________________________________________________________________________
Double_t v2vsMassCB2ChebPol3(Double_t*x, Double_t *par){ // bkg 6 + sig 7 + flow bkg 4 + flow sig 1 = 18
    return (backFcnCheb1(x, par)*v2BackgroundPol3(x, &par[13]) + par[17]*CrystalBallExtended(x, &par[6]))/(backFcnCheb1(x, par) + CrystalBallExtended(x, &par[6]));
}
//_______________________________________________________________________________
Double_t v2vsMassCB2ChebPol3Y(Double_t*x, Double_t *par){ // bkg 6 + sig 7 + flow bkg 4 + flow sig 1 = 18
    return (backFcnCheb(x, par)*v2BackgroundPol3Y(x, &par[13]) + par[17]*CrystalBallExtended(x, &par[6]))/(backFcnCheb(x, par) + CrystalBallExtended(x, &par[6]));
}



//_______________________________________________________________________________
Double_t v2vsMassCB2Pol3Y(Double_t*x, Double_t *par){ // bkg 4 + sig 7 + flow bkg 4 + flow sig 1 = 16
    return (backFcn(x, par)*v2BackgroundPol3Y(x, &par[11]) + par[15]*CrystalBallExtended(x, &par[4]))/(backFcn(x, par) + CrystalBallExtended(x, &par[4]));
}

//_______________________________________________________________________________
Double_t v2vsMassCB2Pol4(Double_t*x, Double_t *par){ // bkg 4 + sig 7 + flow bkg 4 + flow sig 1 = 16
    return (backFcn(x, par)*v2BackgroundPol4(x, &par[11]) + par[16]*CrystalBallExtended(x, &par[4]))/(backFcn(x, par) + CrystalBallExtended(x, &par[4]));
}

//_______________________________________________________________________________
Double_t v2vsMassCB2Pol4_PeakPol(Double_t*x, Double_t *par){ // bkg 4 + sig 7 + flow bkg 4 + flow sig 1 = 16
    return (backFcn(x, par)*v2BackgroundPol4(x, &par[11]) + ( par[16] + par[17]*(x[0]-par[5]) + par[18]*(x[0]-par[5])*(x[0]-par[5])  )*CrystalBallExtended(x, &par[4]))/(backFcn(x, par) + CrystalBallExtended(x, &par[4]));
}

//_______________________________________________________________________________
Double_t v2vsMassCB2ChebPol4_PeakPol(Double_t*x, Double_t *par){ // bkg 4 + sig 7 + flow bkg 4 + flow sig 1 = 16
    return (backFcnChebJpsi(x, par)*v2BackgroundPol4(x, &par[11]) + ( par[16]   )*CrystalBallExtended(x, &par[4]))/(backFcnChebJpsi(x, par) + CrystalBallExtended(x, &par[4]));
}


//_______________________________________________________________________________
Double_t v2vsMassCB2Pol4Y(Double_t*x, Double_t *par){ // bkg 4 + sig 7 + flow bkg 4 + flow sig 1 = 16
    return (backFcn(x, par)*v2BackgroundPol4Y(x, &par[11]) + par[16]*CrystalBallExtended(x, &par[4]))/(backFcn(x, par) + CrystalBallExtended(x, &par[4]));
}

Double_t v2vsMassTwoCB2Pol3Y(Double_t*x, Double_t *par){ // bkg 4 + sig 8 + flow bkg 4 + flow sig 1 = 17
    return (backFcn(x, par)*v2BackgroundPol3Y(x, &par[12]) + par[16]*TwoCrystalBallExtended(x, &par[4]))/(backFcn(x, par) + TwoCrystalBallExtended(x, &par[4]));
}


//_______________________________________________________________________________
Double_t v2vsMassCB2GausPol3Y(Double_t*x, Double_t *par){ // bkg 4 + sig 7 + flow bkg 4 + flow sig 1 = 16
    return (backFcnVWGaus(x, par)*v2BackgroundPol3Y(x, &par[11]) + par[15]*CrystalBallExtended(x, &par[4]))/(backFcnVWGaus(x, par) + CrystalBallExtended(x, &par[4]));
}

//_______________________________________________________________________________
Double_t v2vsMassCB2GausPol4(Double_t*x, Double_t *par){ // bkg 4 + sig 7 + flow bkg 4 + flow sig 1 = 16
    return (backFcnVWGaus(x, par)*v2BackgroundPol4(x, &par[11]) + par[16]*CrystalBallExtended(x, &par[4]))/(backFcnVWGaus(x, par) + CrystalBallExtended(x, &par[4]));
}

//_______________________________________________________________________________
Double_t v2vsMassCB2GausPol4_PeakPol(Double_t*x, Double_t *par){ // bkg 4 + sig 7 + flow bkg 4 + flow sig 1 = 16
    return (backFcnVWGaus(x, par)*v2BackgroundPol4(x, &par[11]) + ( par[16]  )*CrystalBallExtended(x, &par[4]))/(backFcnVWGaus(x, par) + CrystalBallExtended(x, &par[4]));
}

//_______________________________________________________________________________
Double_t v2vsMassCB2GausPol4Y(Double_t*x, Double_t *par){ // bkg 4 + sig 7 + flow bkg 4 + flow sig 1 = 16
    
    return (backFcnVWGaus(x, par)*v2BackgroundPol4Y(x, &par[11]) + par[16]*CrystalBallExtended(x, &par[4]))/(backFcnVWGaus(x, par) + CrystalBallExtended(x, &par[4]));
}

//_______________________________________________________________________________
Double_t v2vsMassCB2GausPol4Y_2S(Double_t*x, Double_t *par){ // bkg 4 + sig 7 + flow bkg 4 + flow sig 1 = 16
    
    return (backFcnVWGaus(x, par)*v2BackgroundPol4Y(x, &par[18]) + (par[23]*CrystalBallExtended(x, &par[4]) + par[24]*CrystalBallExtended(x, &par[11]) ) )/(backFcnVWGaus(x, par) + CrystalBallExtended(x, &par[4]) + CrystalBallExtended(x, &par[11]) );
}
//_______________________________________________________________________________
Double_t v2vsMassCB2Pol4Y_2S(Double_t*x, Double_t *par){ // bkg 4 + sig 7 + flow bkg 4 + flow sig 1 = 16
    
    return (backFcn(x, par)*v2BackgroundPol4Y(x, &par[18]) + (par[23]*CrystalBallExtended(x, &par[4]) + par[24]*CrystalBallExtended(x, &par[11]) ) )/(backFcn(x, par) + CrystalBallExtended(x, &par[4]) + CrystalBallExtended(x, &par[11]) );
}



Double_t v2vsMassTwoCB2GausPol3Y(Double_t*x, Double_t *par){ // bkg 4 + sig 8 + flow bkg 4 + flow sig 1 = 17
    return (backFcnVWGaus(x, par)*v2BackgroundPol3Y(x, &par[12]) + par[16]*TwoCrystalBallExtended(x, &par[4]))/(backFcnVWGaus(x, par) + TwoCrystalBallExtended(x, &par[4]));
}


//_______________________________________________________________________________
Double_t v2vsMassAlpha(Double_t*x, Double_t *par){ // bkg 6 + sig 7  = 13
    return (CrystalBallExtended(x, &par[6]))/(backFcnCheb1(x, par) + CrystalBallExtended(x, &par[6]));
}
Double_t v2vsMassAlphaCheb(Double_t*x, Double_t *par){ // bkg 6 + sig 7  = 13
    return (CrystalBallExtended(x, &par[6]))/(backFcnCheb(x, par) + CrystalBallExtended(x, &par[6]));
}
Double_t v2vsMassAlphaGaus(Double_t*x, Double_t *par){ // bkg 6 + sig 7  = 13
    return (CrystalBallExtended(x, &par[4]))/(backFcnVWGaus(x, par) + CrystalBallExtended(x, &par[4]));
}
Double_t v2vsMassAlphaExpo(Double_t*x, Double_t *par){ // bkg 6 + sig 7  = 13
    return (CrystalBallExtended(x, &par[4]))/(backFcn(x, par) + CrystalBallExtended(x, &par[4]));
}
Double_t v2vsMassThreeAlphaVWG(Double_t*x, Double_t *par){ // bkg 4 + sig 9  = 13
    return (ThreeCrystalBallExtended(x, &par[4]))/(backFcnVWGaus(x, par) + ThreeCrystalBallExtended(x, &par[4]));
}

Double_t v2vsMassAlphaVWG(Double_t*x, Double_t *par){ // bkg 7 + sig 12  = 18
    return (CrystalBallExtended(x, &par[4]))/(backFcnVWGaus(x, par) + CrystalBallExtended(x, &par[4]));
}

Double_t v2vsMassThreeIndAlphaVWG(Double_t*x, Double_t *par){ // bkg 4 + sig 21  = 25
    return (ThreeIndCrystalBallExtended(x, &par[4]))/(backFcnVWGaus(x, par) + ThreeIndCrystalBallExtended(x, &par[4]));
}


Double_t v2vsMassThreeNA60AlphaVWG(Double_t*x, Double_t *par){ // bkg 4 + sig 13  = 17
    return (ThreeNA60(x, &par[4]))/(backFcnVWGaus(x, par) + ThreeNA60(x, &par[4]));
}
//_______________________________________________________________________________
Double_t v2vsMassSignificanceCheb(Double_t*x, Double_t *par){ // bkg 7 + sig 7  = 14
    return (CrystalBallExtended(x, &par[6]))/ TMath::Sqrt( backFcnCheb(x, par) + CrystalBallExtended(x, &par[6]) );
}
Double_t SignalBackgroundCheb(Double_t*x, Double_t *par){ // bkg 7 + sig 7  = 14
    return (CrystalBallExtended(x, &par[6]))/ backFcnCheb(x, par) ;
}

//_______________________________________________________________________________
Double_t v2vsMassSignificanceGaus(Double_t*x, Double_t *par){ // bkg 4 + sig 7  = 11
    return (CrystalBallExtended(x, &par[4]))/ TMath::Sqrt( backFcnVWGaus(x, par) + CrystalBallExtended(x, &par[4]) );
}
Double_t SignalBackgroundGaus(Double_t*x, Double_t *par){ // bkg 4 + sig 7  = 11
    return (CrystalBallExtended(x, &par[4]))/ backFcnVWGaus(x, par) ;
}

//_______________________________________________________________________________
Double_t v2vsMassSignificanceExpo(Double_t*x, Double_t *par){ // bkg 4 + sig 7  = 11
    return (CrystalBallExtended(x, &par[4]))/ TMath::Sqrt( backFcn(x, par) + CrystalBallExtended(x, &par[4]) );
}
Double_t SignalBackgroundExpo(Double_t*x, Double_t *par){ // bkg 4 + sig 7  = 11
    return (CrystalBallExtended(x, &par[4]))/ backFcn(x, par) ;
}




//_______________________________________________________________________________
Double_t v2vsMassNA60QVWGPol3(Double_t*x, Double_t *par){ // bkg 5 + sig 11 + flow bkg 4 + flow sig 1 = 21
    return (backFcnQVWGaus(x, par)*v2BackgroundPol3(x, &par[16]) + par[20]*FuncJpsi_NA60(x, &par[5]))/(backFcnQVWGaus(x, par) + FuncJpsi_NA60(x, &par[5]));
}
//_______________________________________________________________________________
Double_t v2vsMassNA60PolRPol3(Double_t*x, Double_t *par){ // bkg 6 + sig 11 + flow bkg 4 + flow sig 1 = 22
    return (backFcnPolR(x, par)*v2BackgroundPol3(x, &par[17]) + par[21]*FuncJpsi_NA60(x, &par[6])) /(backFcnPolR(x, par) + FuncJpsi_NA60(x, &par[6]));
}
//_______________________________________________________________________________
Double_t v2vsMassNA60ChebPol3(Double_t*x, Double_t *par){ // bkg 7 + sig 11 + flow bkg 4 + flow sig 1 = 23
    return (backFcnCheb(x, par)*v2BackgroundPol3(x, &par[18]) + par[22]*FuncJpsi_NA60(x, &par[7]))/(backFcnCheb(x, par) + FuncJpsi_NA60(x, &par[7]));
}
//_______________________________________________________________________________
Double_t v2vsMassCB2QVWGPol4(Double_t*x, Double_t *par){ // bkg 5 + sig 7 + flow bkg 5 + flow sig 1 = 18
    return (backFcnQVWGaus(x, par)*v2BackgroundPol4(x, &par[12]) + par[17]*CrystalBallExtended(x, &par[5]))/(backFcnQVWGaus(x, par) + CrystalBallExtended(x, &par[5]));
}
//_______________________________________________________________________________
Double_t v2vsMassCB2PolRPol4(Double_t*x, Double_t *par){ // bkg 6 + sig 7 + flow bkg 5 + flow sig 1 = 19
    return (backFcnPolR(x, par)*v2BackgroundPol4(x, &par[13]) + par[18]*CrystalBallExtended(x, &par[6])) /(backFcnPolR(x, par) + CrystalBallExtended(x, &par[6]));
}
//_______________________________________________________________________________
Double_t v2vsMassCB2ChebPol4(Double_t*x, Double_t *par){ // bkg 7 + sig 7 + flow bkg 5 + flow sig 1 = 20
    return (backFcnCheb(x, par)*v2BackgroundPol4(x, &par[14]) + par[19]*CrystalBallExtended(x, &par[7]))/(backFcnCheb(x, par) + CrystalBallExtended(x, &par[7]));
}
//_______________________________________________________________________________
Double_t v2vsMassNA60QVWGPol4(Double_t*x, Double_t *par){ // bkg 5 + sig 11 + flow bkg 5 + flow sig 1 = 22
    return (backFcnQVWGaus(x, par)*v2BackgroundPol4(x, &par[16]) + par[21]*FuncJpsi_NA60(x, &par[5]))/(backFcnQVWGaus(x, par) + FuncJpsi_NA60(x, &par[5]));
}
//_______________________________________________________________________________
Double_t v2vsMassNA60PolRPol4(Double_t*x, Double_t *par){ // bkg 6 + sig 11 + flow bkg 5 + flow sig 1 = 23
    return (backFcnPolR(x, par)*v2BackgroundPol4(x, &par[17]) + par[22]*FuncJpsi_NA60(x, &par[6])) /(backFcnPolR(x, par) + FuncJpsi_NA60(x, &par[6]));
}
//_______________________________________________________________________________
Double_t v2vsMassNA60ChebPol4(Double_t*x, Double_t *par){ // bkg 7 + sig 11 + flow bkg 5 + flow sig 1 = 24
    return (backFcnCheb(x, par)*v2BackgroundPol4(x, &par[18]) + par[23]*FuncJpsi_NA60(x, &par[7]))/(backFcnCheb(x, par) + FuncJpsi_NA60(x, &par[7]));
}
//_______________________________________________________________________________
Double_t v2vsMassCB2QVWGCheb(Double_t*x, Double_t *par){ // bkg 5 + sig 7 + flow bkg 5 + flow sig 1 = 18
    return (backFcnQVWGaus(x, par)*v2BackgroundCheb(x, &par[12]) + par[17]*CrystalBallExtended(x, &par[5]))/(backFcnQVWGaus(x, par) + CrystalBallExtended(x, &par[5]));
}




//_______________________________________________________________________________
Double_t v2vsMassCB2VWGChebJpsi(Double_t*x, Double_t *par){
    return (backFcnVWGaus(x, par)*v2BackgroundChebJpsi(x, &par[11]) + par[16]*CrystalBallExtended(x, &par[4]))/(backFcnVWGaus(x, par) + CrystalBallExtended(x, &par[4]));
}

//_______________________________________________________________________________
Double_t v2vsMassCB2VWGChebJpsi_PeakPol(Double_t*x, Double_t *par){
    return (backFcnVWGaus(x, par)*v2BackgroundChebJpsi(x, &par[11]) + ( par[16]   )*CrystalBallExtended(x, &par[4]))/(backFcnVWGaus(x, par) + CrystalBallExtended(x, &par[4]));
}

//_______________________________________________________________________________
Double_t v2vsMassCB2ChebJpsi(Double_t*x, Double_t *par){
    return (backFcn(x, par)*v2BackgroundChebJpsi(x, &par[11]) + par[16]*CrystalBallExtended(x, &par[4]))/(backFcn(x, par) + CrystalBallExtended(x, &par[4]));
}

//_______________________________________________________________________________
Double_t v2vsMassCB2ChebJpsi_PeakPol(Double_t*x, Double_t *par){
    return (backFcn(x, par)*v2BackgroundChebJpsi(x, &par[11]) + ( par[16] + par[17]*(x[0]-par[5]) + par[18]*(x[0]-par[5])*(x[0]-par[5]) )*CrystalBallExtended(x, &par[4]))/(backFcn(x, par) + CrystalBallExtended(x, &par[4]));
}

//_______________________________________________________________________________
Double_t v2vsMassCB2ChebChebJpsi_PeakPol(Double_t*x, Double_t *par){
    return (backFcnChebJpsi(x, par)*v2BackgroundChebJpsi(x, &par[11]) + ( par[16]  )*CrystalBallExtended(x, &par[4]))/(backFcnChebJpsi(x, par) + CrystalBallExtended(x, &par[4]));
}

//_______________________________________________________________________________
Double_t v2vsMassCB2PolRCheb(Double_t*x, Double_t *par){ // bkg 6 + sig 7 + flow bkg 5 + flow sig 1 = 19
    return (backFcnPolR(x, par)*v2BackgroundCheb(x, &par[13]) + par[18]*CrystalBallExtended(x, &par[6])) /(backFcnPolR(x, par) + CrystalBallExtended(x, &par[6]));
}
//_______________________________________________________________________________
Double_t v2vsMassCB2ChebCheb(Double_t*x, Double_t *par){ // bkg 7 + sig 7 + flow bkg 5 + flow sig 1 = 20
    return (backFcnCheb(x, par)*v2BackgroundCheb(x, &par[14]) + par[19]*CrystalBallExtended(x, &par[7]))/(backFcnCheb(x, par) + CrystalBallExtended(x, &par[7]));
}
//_______________________________________________________________________________
Double_t v2vsMassNA60QVWGCheb(Double_t*x, Double_t *par){ // bkg 5 + sig 11 + flow bkg 5 + flow sig 1 = 22
    return (backFcnQVWGaus(x, par)*v2BackgroundCheb(x, &par[16]) + par[21]*FuncJpsi_NA60(x, &par[5]))/(backFcnQVWGaus(x, par) + FuncJpsi_NA60(x, &par[5]));
}





Double_t v2vsMassThreeNA60ChebPol3(Double_t*x, Double_t *par){ // bkg 6 + sig 13 + flow bkg 4 + flow sig 1 = 24
    return (backFcnCheb(x, par)*v2BackgroundPol3Y(x, &par[19]) + par[23]*ThreeNA60(x, &par[6]))/(backFcnCheb(x, par) + ThreeNA60(x, &par[6]));
}

Double_t v2vsMassThreeNA60Cheb(Double_t*x, Double_t *par){ // bkg 4 + sig 13 + flow bkg 5 + flow sig 1 = 23
    return (backFcn(x, par)*v2BackgroundChebY(x, &par[17]) + par[22]*ThreeNA60(x, &par[4]))/(backFcn(x, par) + ThreeNA60(x, &par[4]));
}
Double_t v2vsMassCB2ExpCheb(Double_t*x, Double_t *par){ // bkg 4 + sig 7 + flow bkg 5 + flow sig 1 = 17
    return (backFcn(x, par)*v2BackgroundChebY(x, &par[11]) + par[16]*CrystalBallExtended(x, &par[4]))/(backFcn(x, par) + CrystalBallExtended(x, &par[4]));
}

Double_t v2vsMassCB2ExpPol3(Double_t*x, Double_t *par){ // bkg 4 + sig 7 + flow bkg 5 + flow sig 1 = 17
    return (backFcn(x, par)*v2BackgroundPol3Exp(x, &par[11]) + par[16]*CrystalBallExtended(x, &par[4]))/(backFcn(x, par) + CrystalBallExtended(x, &par[4]));
}
Double_t v2vsMassCB2VWGPol3(Double_t*x, Double_t *par){ // bkg 4 + sig 7 + flow bkg 5 + flow sig 1 = 17
    return (backFcnVWGaus(x, par)*v2BackgroundPol3Exp(x, &par[11]) + par[16]*CrystalBallExtended(x, &par[4]))/(backFcnVWGaus(x, par) + CrystalBallExtended(x, &par[4]));
}

Double_t v2vsMassThreeCB2VWGChebY(Double_t*x, Double_t *par){ // bkg 4 + sig 9 + flow bkg 5 + flow sig 1 = 19
    return (backFcnVWGaus(x, par)*v2BackgroundChebY(x, &par[13]) + par[18]*ThreeCrystalBallExtended(x, &par[4]))/(backFcnVWGaus(x, par) + ThreeCrystalBallExtended(x, &par[4]));
}

Double_t v2vsMassThreeCB2ExpoCheb(Double_t*x, Double_t *par){ // bkg 4 + sig 7 + flow bkg 5 + flow sig 1 = 17
    return (backFcn(x, par)*v2BackgroundChebY(x, &par[11]) + par[16]*CrystalBallExtended(x, &par[4]))/(backFcn(x, par) + CrystalBallExtended(x, &par[4]));
}

Double_t v2vsMassTwoCB2ExpoCheb(Double_t*x, Double_t *par){ // bkg 4 + sig 8 + flow bkg 5 + flow sig 1 = 18
    return (backFcn(x, par)*v2BackgroundChebY(x, &par[12]) + par[17]*TwoCrystalBallExtended(x, &par[4]))/(backFcn(x, par) + TwoCrystalBallExtended(x, &par[4]));
}

Double_t v2vsMassCB2GausCheb(Double_t*x, Double_t *par){ // bkg 4 + sig 7 + flow bkg 5 + flow sig 1 = 17
    return (backFcnVWGaus(x, par)*v2BackgroundChebY(x, &par[11]) + par[16]*CrystalBallExtended(x, &par[4]))/(backFcnVWGaus(x, par) + CrystalBallExtended(x, &par[4]));
}


Double_t v2vsMassCB2VWGCheb(Double_t*x, Double_t *par){ // bkg 4 + sig 7 + flow bkg 5 + flow sig 1 = 17
    return (backFcnVWGaus(x, par)*v2BackgroundCheb(x, &par[11]) + par[16]*CrystalBallExtended(x, &par[4]))/(backFcnVWGaus(x, par) + CrystalBallExtended(x, &par[4]));
}


Double_t v2vsMassTwoCB2VWGCheb(Double_t*x, Double_t *par){ // bkg 4 + sig 8 + flow bkg 5 + flow sig 1 = 18
    return (backFcnVWGaus(x, par)*v2BackgroundCheb(x, &par[12]) + par[17]*TwoCrystalBallExtended(x, &par[4]))/(backFcnVWGaus(x, par) + TwoCrystalBallExtended(x, &par[4]));
}
Double_t v2vsMassThreeIndCB2VWGCheb(Double_t*x, Double_t *par){ // bkg 4 + sig 21 + flow bkg 5 + flow sig 1 = 31
    return (backFcnVWGaus(x, par)*v2BackgroundCheb(x, &par[25]) + par[30]*ThreeIndCrystalBallExtended(x, &par[4]))/(backFcnVWGaus(x, par) + ThreeIndCrystalBallExtended(x, &par[4]));
}
Double_t v2vsMassThreeCB2VWGCheb(Double_t*x, Double_t *par){ // bkg 4 + sig 9 + flow bkg 5 + flow sig 1 = 19
    return (backFcnVWGaus(x, par)*v2BackgroundCheb(x, &par[13]) + par[18]*ThreeCrystalBallExtended(x, &par[4]))/(backFcnVWGaus(x, par) + ThreeCrystalBallExtended(x, &par[4]));
}
Double_t v2vsMassThreeIndCB2ExpCheb(Double_t*x, Double_t *par){ // bkg 4 + sig 21 + flow bkg 5 + flow sig 1 = 31
    return (backFcn(x, par)*v2BackgroundCheb(x, &par[25]) + par[30]*ThreeIndCrystalBallExtended(x, &par[4]))/(backFcn(x, par) + ThreeIndCrystalBallExtended(x, &par[4]));
}
Double_t v2vsMassThreeCB2ExpCheb(Double_t*x, Double_t *par){ // bkg 4 + sig 21 + flow bkg 5 + flow sig 1 = 31
    return (backFcn(x, par)*v2BackgroundCheb(x, &par[13]) + par[18]*ThreeCrystalBallExtended(x, &par[4]))/(backFcn(x, par) + ThreeCrystalBallExtended(x, &par[4]));
}




Double_t v2vsMassCB2ExpCheb_1S_2S(Double_t*x, Double_t *par){ // 25
    return (backFcn(x, par)*v2BackgroundCheb(x, &par[18]) + (par[23]*CrystalBallExtended(x, &par[4]) + par[24]*CrystalBallExtended(x, &par[11])  ))/(backFcn(x, par) + CrystalBallExtended(x, &par[4]) + CrystalBallExtended(x, &par[11]) );
}

Double_t v2vsMassCB2VWGCheb_1S_2S(Double_t*x, Double_t *par){ // 25
    return (backFcnVWGaus(x, par)*v2BackgroundCheb(x, &par[18]) + (par[23]*CrystalBallExtended(x, &par[4]) + par[24]*CrystalBallExtended(x, &par[11])  ))/(backFcnVWGaus(x, par) + CrystalBallExtended(x, &par[4]) + CrystalBallExtended(x, &par[11]) );
    
}



//_______________________________________________________________________________
Double_t v2vsMassNA60PolRCheb(Double_t*x, Double_t *par){ // bkg 6 + sig 11 + flow bkg 5 + flow sig 1 = 23
    return (backFcnPolR(x, par)*v2BackgroundCheb(x, &par[17]) + par[22]*FuncJpsi_NA60(x, &par[6])) /(backFcnPolR(x, par) + FuncJpsi_NA60(x, &par[6]));
}
//_______________________________________________________________________________
Double_t v2vsMassNA60ChebCheb(Double_t*x, Double_t *par){ // bkg 7 + sig 11 + flow bkg 5 + flow sig 1 = 24
    return (backFcnCheb(x, par)*v2BackgroundCheb(x, &par[18]) + par[23]*FuncJpsi_NA60(x, &par[7]))/(backFcnCheb(x, par) + FuncJpsi_NA60(x, &par[7]));
}

//_______________________________________________________________________________
Double_t v2vsMassCB2QVWGPolExp(Double_t*x, Double_t *par){ // bkg 5 + sig 7 + flow bkg 5 + flow sig 1 = 18
    return (backFcnQVWGaus(x, par)*v2BackgroundPolExp(x, &par[12]) + par[17]*CrystalBallExtended(x, &par[5]))/(backFcnQVWGaus(x, par) + CrystalBallExtended(x, &par[5]));
}
//_______________________________________________________________________________
Double_t v2vsMassThreeCB2QVWGPolExp(Double_t*x, Double_t *par){ // bkg 5 + sig 21 + flow bkg 5 + flow sig 1 = 31
    return (backFcnQVWGaus(x, par)*v2BackgroundPolExp(x, &par[25]) + par[30]*ThreeIndCrystalBallExtended(x, &par[5]))/(backFcnQVWGaus(x, par) + ThreeIndCrystalBallExtended(x, &par[5]));
}


//_______________________________________________________________________________
Double_t v2vsMassCB2PolExp(Double_t*x, Double_t *par){ // bkg 4 + sig 21 + flow bkg 5 + flow sig 1 = 30
    return (backFcn(x, par)*v2BackgroundPolExp(x, &par[24]) + par[29]*ThreeIndCrystalBallExtended(x, &par[4]))/(backFcn(x, par) + ThreeIndCrystalBallExtended(x, &par[4]));
}

//_______________________________________________________________________________
Double_t v2vsMassCB2PolRPolExp(Double_t*x, Double_t *par){ // bkg 6 + sig 7 + flow bkg 5 + flow sig 1 = 19
    return (backFcnPolR(x, par)*v2BackgroundPolExp(x, &par[13]) + par[18]*CrystalBallExtended(x, &par[6])) /(backFcnPolR(x, par) + CrystalBallExtended(x, &par[6]));
}
//_______________________________________________________________________________
Double_t v2vsMassCB2ChebPolExp(Double_t*x, Double_t *par){ // bkg 7 + sig 7 + flow bkg 5 + flow sig 1 = 20
    return (backFcnCheb(x, par)*v2BackgroundPolExp(x, &par[14]) + par[19]*CrystalBallExtended(x, &par[7]))/(backFcnCheb(x, par) + CrystalBallExtended(x, &par[7]));
}
//_______________________________________________________________________________
Double_t v2vsMassNA60QVWGPolExp(Double_t*x, Double_t *par){ // bkg 5 + sig 11 + flow bkg 5 + flow sig 1 = 22
    return (backFcnQVWGaus(x, par)*v2BackgroundPolExp(x, &par[16]) + par[21]*FuncJpsi_NA60(x, &par[5]))/(backFcnQVWGaus(x, par) + FuncJpsi_NA60(x, &par[5]));
}
//_______________________________________________________________________________
Double_t v2vsMassNA60PolRPolExp(Double_t*x, Double_t *par){ // bkg 6 + sig 11 + flow bkg 5 + flow sig 1 = 23
    return (backFcnPolR(x, par)*v2BackgroundPolExp(x, &par[17]) + par[22]*FuncJpsi_NA60(x, &par[6])) /(backFcnPolR(x, par) + FuncJpsi_NA60(x, &par[6]));
}
//_______________________________________________________________________________
Double_t v2vsMassNA60ChebPolExp(Double_t*x, Double_t *par){ // bkg 7 + sig 11 + flow bkg 5 + flow sig 1 = 24
    return (backFcnCheb(x, par)*v2BackgroundPolExp(x, &par[18]) + par[23]*FuncJpsi_NA60(x, &par[7]))/(backFcnCheb(x, par) + FuncJpsi_NA60(x, &par[7]));
}

//_______________________________________________________________________________
Double_t v2vsPtPol2Po1(Double_t *x, Double_t*par){
    
    //par[1]=4;
    
    if(x[0]<par[1]){
        return ( MyGauss(x, &par[0]) )*x[0] + ( par[3]*x[0] ) + par[4]*0;

    }
    else{
        return ( MyGauss(x, &par[0]) )*x[0] + ( par[3]*x[0] ) + par[4]*(x[0]-par[1]);
        
    }
    
}

//_______________________________________________________________________________
Double_t v2vsCentralityPol2(Double_t *x, Double_t*par){
    return ( par[0] + par[1]*(x[0]-par[3]) + par[2]*(x[0]-par[3])*(x[0]-par[3]) ) ;
}
//_______________________________________________________________________________
Double_t v2vsCentralityPol1(Double_t *x, Double_t*par){
    return par[0] + par[1]*x[0] ;
}
//_______________________________________________________________________________
Double_t v2vsCentralityPol0(Double_t *x, Double_t*par){
    return par[0]  ;
}
//_______________________________________________________________________________
Double_t PtvsCentralityInverseP(Double_t *x, Double_t*par){
    Double_t a = TMath::Power(par[3]*x[0],-1);
    Double_t b = TMath::Power(par[2]+a,-1);

    //return (  par[0]/(1 + par[1]*inverse_p )  );
    return par[0]*b+par[1];

}


// pt vs. Mass background functions
//_______________________________________________________________________________
Double_t ptBackground(Double_t *x, Double_t*par){
    return par[0]*x[0]*x[0] + par[1]*x[0] + par[2];
}

// pt vs. Mass background functions
//_______________________________________________________________________________
Double_t ptBackgroundPol3(Double_t *x, Double_t*par){
    return par[0] + par[1]*(x[0]-3.1) + par[2]*((x[0]-3.1)*(x[0]-3.1)) + par[3]*((x[0]-3.1)*(x[0]-3.1)*(x[0]-3.1));
}

//_______________________________________________________________________________
// pT vs. Mass signal + background functions
//
//_______________________________________________________________________________
Double_t ptvsMassCBVWG(Double_t*x, Double_t *par){ // bkg 4 + sig 5 + flow bkg 3 + flow sig 1= 13
    return (backFcnVWGaus(x, par)*ptBackground(x, &par[9])  + par[12]*CrystalBall(x, &par[4]))/(backFcnVWGaus(x, par) + CrystalBall(x, &par[4]));
}
//_______________________________________________________________________________
Double_t ptvsMassCBPol(Double_t*x, Double_t *par){ // bkg 4 + sig 5 + flow bkg 3 + flow sig 1= 13
    return (backFcnPol(x, par)*ptBackground(x, &par[9])  + par[12]*CrystalBall(x, &par[4])) /(backFcnPol(x, par) + CrystalBall(x, &par[4]));
}
//_______________________________________________________________________________
Double_t ptvsMassCB2VWG(Double_t*x, Double_t *par){ // bkg 4 + sig 7 + flow bkg 3 + flow sig 1= 15
    return (backFcnVWGaus(x, par)*ptBackground(x, &par[11])  + par[14]*CrystalBallExtended(x, &par[4]))/(backFcnVWGaus(x, par) + CrystalBallExtended(x, &par[4]));
}
//_______________________________________________________________________________
Double_t ptvsMassCB2Pol(Double_t*x, Double_t *par){ // bkg 4 + sig 7 + flow bkg 3 + flow sig 1= 15
    return (backFcnPol(x, par)*ptBackground(x, &par[11])  + par[14]*CrystalBallExtended(x, &par[4])) /(backFcnPol(x, par) + CrystalBallExtended(x, &par[4]));
}
//_______________________________________________________________________________
Double_t ptvsMassCB2QVWG(Double_t*x, Double_t *par){ // bkg 5 + sig 7 + flow bkg 3 + flow sig 1 = 16
    return (backFcnQVWGaus(x, par)*ptBackground(x, &par[12]) + par[15]*CrystalBallExtended(x, &par[5]))/(backFcnQVWGaus(x, par) + CrystalBallExtended(x, &par[5]));
}
//_______________________________________________________________________________
Double_t ptvsMassCB2PolR(Double_t*x, Double_t *par){ // bkg 6 + sig 7 + flow bkg 3 + flow sig 1 = 17
    return (backFcnPolR(x, par)*ptBackground(x, &par[13]) + par[16]*CrystalBallExtended(x, &par[6])) /(backFcnPolR(x, par) + CrystalBallExtended(x, &par[6]));
}
//_______________________________________________________________________________
Double_t ptvsMassNA60QVWG(Double_t*x, Double_t *par){ // bkg 5 + sig 11 + flow bkg 3 + flow sig 1 = 20
    return (backFcnQVWGaus(x, par)*ptBackground(x, &par[16]) + par[19]*FuncJpsi_NA60(x, &par[5]))/(backFcnQVWGaus(x, par) + FuncJpsi_NA60(x, &par[5]));
}
//_______________________________________________________________________________
Double_t ptvsMassNA60PolR(Double_t*x, Double_t *par){ // bkg 6 + sig 7 + flow bkg 3 + flow sig 1 = 17
    return (backFcnPolR(x, par)*ptBackground(x, &par[17]) + par[20]*FuncJpsi_NA60(x, &par[6])) /(backFcnPolR(x, par) + FuncJpsi_NA60(x, &par[6]));
}
//_______________________________________________________________________________
Double_t ptvsMassCB2QVWGPol3(Double_t*x, Double_t *par){ // bkg 5 + sig 7 + pt bkg 4 + pt sig 1 = 16
    return (backFcnQVWGaus(x, par)*ptBackgroundPol3(x, &par[12]) + par[16]*CrystalBallExtended(x, &par[5]))/(backFcnQVWGaus(x, par) + CrystalBallExtended(x, &par[5]));
}
//_______________________________________________________________________________
Double_t ptvsMassCB2PolRPol3(Double_t*x, Double_t *par){ // bkg 6 + sig 7 + flow bkg 4 + flow sig 1 = 17
    return (backFcnPolR(x, par)*ptBackgroundPol3(x, &par[13]) + par[17]*CrystalBallExtended(x, &par[6])) /(backFcnPolR(x, par) + CrystalBallExtended(x, &par[6]));
}
//_______________________________________________________________________________
Double_t ptvsMassCB2ChebPol3(Double_t*x, Double_t *par){ // bkg 7 + sig 7 + flow bkg 4 + flow sig 1 = 18
    return (backFcnCheb(x, par)*ptBackgroundPol3(x, &par[14]) + par[18]*CrystalBallExtended(x, &par[7])) /(backFcnCheb(x, par) + CrystalBallExtended(x, &par[7]));
}
//_______________________________________________________________________________
Double_t ptvsMassNA60QVWGPol3(Double_t*x, Double_t *par){ // bkg 5 + sig 11 + flow bkg 4 + flow sig 1 = 20
    return (backFcnQVWGaus(x, par)*ptBackgroundPol3(x, &par[16]) + par[20]*FuncJpsi_NA60(x, &par[5]))/(backFcnQVWGaus(x, par) + FuncJpsi_NA60(x, &par[5]));
}
//_______________________________________________________________________________
Double_t ptvsMassNA60PolRPol3(Double_t*x, Double_t *par){ // bkg 6 + sig 11 + flow bkg 4 + flow sig 1 = 21
    return (backFcnPolR(x, par)*ptBackgroundPol3(x, &par[17]) + par[21]*FuncJpsi_NA60(x, &par[6])) /(backFcnPolR(x, par) + FuncJpsi_NA60(x, &par[6]));
}
//_______________________________________________________________________________
Double_t ptvsMassNA60ChebPol3(Double_t*x, Double_t *par){ // bkg 7 + sig 11 + flow bkg 4 + flow sig 1 = 22
    return (backFcnCheb(x, par)*ptBackgroundPol3(x, &par[18]) + par[22]*FuncJpsi_NA60(x, &par[7])) /(backFcnCheb(x, par) + FuncJpsi_NA60(x, &par[7]));
}

// y vs. Mass background functions
//_______________________________________________________________________________
Double_t yBackground(Double_t *x, Double_t*par){
    return par[0]*x[0]*x[0] + par[1]*x[0] + par[2];
}

// y vs. Mass background functions
//_______________________________________________________________________________
Double_t yBackgroundPol3(Double_t *x, Double_t*par){
    return par[0] + par[1]*(x[0]-3.1) + par[2]*((x[0]-3.1)*(x[0]-3.1)) + par[3]*((x[0]-3.1)*(x[0]-3.1)*(x[0]-3.1));
}

//_______________________________________________________________________________
// y vs. Mass signal + background functions
//
//_______________________________________________________________________________
Double_t yvsMassCBVWG(Double_t*x, Double_t *par){ // bkg 4 + sig 5 + flow bkg 3 + flow sig 1= 13
    return (backFcnVWGaus(x, par)*yBackground(x, &par[9])  + par[12]*CrystalBall(x, &par[4]))/(backFcnVWGaus(x, par) + CrystalBall(x, &par[4]));
}
//_______________________________________________________________________________
Double_t yvsMassCBPol(Double_t*x, Double_t *par){ // bkg 4 + sig 5 + flow bkg 3 + flow sig 1= 13
    return (backFcnPol(x, par)*yBackground(x, &par[9])  + par[12]*CrystalBall(x, &par[4])) /(backFcnPol(x, par) + CrystalBall(x, &par[4]));
}
//_______________________________________________________________________________
Double_t yvsMassCB2VWG(Double_t*x, Double_t *par){ // bkg 4 + sig 7 + flow bkg 3 + flow sig 1= 15
    return (backFcnVWGaus(x, par)*yBackground(x, &par[11])  + par[14]*CrystalBallExtended(x, &par[4]))/(backFcnVWGaus(x, par) + CrystalBallExtended(x, &par[4]));
}
//_______________________________________________________________________________
Double_t yvsMassCB2Pol(Double_t*x, Double_t *par){ // bkg 4 + sig 7 + flow bkg 3 + flow sig 1= 15
    return (backFcnPol(x, par)*yBackground(x, &par[11])  + par[14]*CrystalBallExtended(x, &par[4])) /(backFcnPol(x, par) + CrystalBallExtended(x, &par[4]));
}
//_______________________________________________________________________________
Double_t yvsMassCB2QVWG(Double_t*x, Double_t *par){ // bkg 5 + sig 7 + flow bkg 3 + flow sig 1 = 16
    return (backFcnQVWGaus(x, par)*yBackground(x, &par[12]) + par[15]*CrystalBallExtended(x, &par[5]))/(backFcnQVWGaus(x, par) + CrystalBallExtended(x, &par[5]));
}
//_______________________________________________________________________________
Double_t yvsMassCB2PolR(Double_t*x, Double_t *par){ // bkg 6 + sig 7 + flow bkg 3 + flow sig 1 = 17
    return (backFcnPolR(x, par)*yBackground(x, &par[13]) + par[16]*CrystalBallExtended(x, &par[6])) /(backFcnPolR(x, par) + CrystalBallExtended(x, &par[6]));
}
//_______________________________________________________________________________
Double_t yvsMassNA60QVWG(Double_t*x, Double_t *par){ // bkg 5 + sig 11 + flow bkg 3 + flow sig 1 = 20
    return (backFcnQVWGaus(x, par)*yBackground(x, &par[16]) + par[19]*FuncJpsi_NA60(x, &par[5]))/(backFcnQVWGaus(x, par) + FuncJpsi_NA60(x, &par[5]));
}
//_______________________________________________________________________________
Double_t yvsMassNA60PolR(Double_t*x, Double_t *par){ // bkg 6 + sig 7 + flow bkg 3 + flow sig 1 = 17
    return (backFcnPolR(x, par)*yBackground(x, &par[17]) + par[20]*FuncJpsi_NA60(x, &par[6])) /(backFcnPolR(x, par) + FuncJpsi_NA60(x, &par[6]));
}
//_______________________________________________________________________________
Double_t yvsMassCB2QVWGPol3(Double_t*x, Double_t *par){ // bkg 5 + sig 7 + y bkg 4 + y sig 1 = 16
    return (backFcnQVWGaus(x, par)*yBackgroundPol3(x, &par[12]) + par[16]*CrystalBallExtended(x, &par[5]))/(backFcnQVWGaus(x, par) + CrystalBallExtended(x, &par[5]));
}
//_______________________________________________________________________________
Double_t yvsMassCB2PolRPol3(Double_t*x, Double_t *par){ // bkg 6 + sig 7 + flow bkg 4 + flow sig 1 = 17
    return (backFcnPolR(x, par)*yBackgroundPol3(x, &par[13]) + par[17]*CrystalBallExtended(x, &par[6])) /(backFcnPolR(x, par) + CrystalBallExtended(x, &par[6]));
}
//_______________________________________________________________________________
Double_t yvsMassCB2ChebPol3(Double_t*x, Double_t *par){ // bkg 7 + sig 7 + flow bkg 4 + flow sig 1 = 18
    return (backFcnCheb(x, par)*yBackgroundPol3(x, &par[14]) + par[18]*CrystalBallExtended(x, &par[7])) /(backFcnCheb(x, par) + CrystalBallExtended(x, &par[7]));
}
//_______________________________________________________________________________
Double_t yvsMassNA60QVWGPol3(Double_t*x, Double_t *par){ // bkg 5 + sig 11 + flow bkg 4 + flow sig 1 = 20
    return (backFcnQVWGaus(x, par)*yBackgroundPol3(x, &par[16]) + par[20]*FuncJpsi_NA60(x, &par[5]))/(backFcnQVWGaus(x, par) + FuncJpsi_NA60(x, &par[5]));
}
//_______________________________________________________________________________
Double_t yvsMassNA60PolRPol3(Double_t*x, Double_t *par){ // bkg 6 + sig 11 + flow bkg 4 + flow sig 1 = 21
    return (backFcnPolR(x, par)*yBackgroundPol3(x, &par[17]) + par[21]*FuncJpsi_NA60(x, &par[6])) /(backFcnPolR(x, par) + FuncJpsi_NA60(x, &par[6]));
}
//_______________________________________________________________________________
Double_t yvsMassNA60ChebPol3(Double_t*x, Double_t *par){ // bkg 7 + sig 11 + flow bkg 4 + flow sig 1 = 22
    return (backFcnCheb(x, par)*yBackgroundPol3(x, &par[18]) + par[22]*FuncJpsi_NA60(x, &par[7])) /(backFcnCheb(x, par) + FuncJpsi_NA60(x, &par[7]));
}

////_______________________________________________________________________________
//Double_t fline(Double_t *x, Double_t *par)
//{
//	if (reject && x[0] > rejectLow  && x[0] < rejectHigh) {
//		TF1::RejectPoint();
//		return 0;
//	}
//	if (expORpol==0) return (par[0]*(exp(x[0]*par[1])) + par[2]*(exp(x[0]*par[3])));    // two exponential functions
//	if (expORpol==1) return par[0] + par[1]*x[0] + par[2]*x[0]*x[0] + par[3]*x[0]*x[0]*x[0];  // pol
//	if (expORpol==2) return backFcnVWGaus(x, par); // vwg
//	if (expORpol==3) return backFcnQVWGaus(x, par); // vwg
//	if (expORpol==4) return backFcnPolR(x, par); // vwg
//	return 0.;
//}






