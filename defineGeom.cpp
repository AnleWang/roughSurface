#include "defineGeom.h"

int main(){
  initParams();
  addFractal();
  dumpFractal();
  if (dumpSpecSmooth) 
    dumpSpectra();
  if (dumpACF)
    dumpHHACF();
}

void initParams(){
  //set default values
  lengthX = 1.; lengthY = 1.;
  nx = 256; ny = 256;
  randSeed = 4579;

  //set default value of rough surface
  lambdaS = 0.01, lambdaR = 0.5, Hurst = 0.8;
  stepHeight = 0.1;
  rmsValueSet = 1.;

  //output options
  dumpStep = 0; dumpSpecSmooth = 1; dumpSpecStep = 0;
  dumpACF = 0;

  ifstream input("params.in");
  if (!input.is_open()) cerr << "#Using default values..." << endl;
  else{
    double param;
    std::string ROL;
    std::size_t NIS = std::string::npos; //NIS = NOT IN STRING
    while(!input.eof()){
      input >> param; getline(input,ROL);
      if      (ROL.find("# lengthX")    !=NIS) lengthX = param;
      else if (ROL.find("# lengthY")    !=NIS) lengthY = param;
      else if (ROL.find("# nx")         !=NIS) nx = param;
      else if (ROL.find("# ny")         !=NIS) ny = param;
      else if (ROL.find("# randSeed")   !=NIS) randSeed = param;
      else if (ROL.find("# lambdaShort")!=NIS) lambdaS = param;
      else if (ROL.find("# lambdaRoll") !=NIS) lambdaR = param;
      else if (ROL.find("# Hurst")      !=NIS) Hurst = param;
      else if (ROL.find("# rmsValueSet")!=NIS) rmsValueSet = param;
      else if (ROL.find("# stepHeight") !=NIS) stepHeight = param;
      else if (ROL.find("# dumpStep")   !=NIS) dumpStep = param;
      else if (ROL.find("# dumpSpecSmooth")!=NIS) dumpSpecSmooth = param;
      else if (ROL.find("# dumpSpecStep")!=NIS)dumpSpecStep = param;
      else if (ROL.find("# dumpACF")    !=NIS) dumpACF = param;
      else continue;
    }
  }
  input.close();

  //sanity check of parameters
  if (((nx&(nx-1))!=0) || ((ny&(ny-1))!=0))
  {cerr << "### nx and/or ny is not a power of 2. Caution!!!\n\n";}
  if ((lengthX<=0)||(lengthY<=0))
  {cerr << "### lengthX and/or lengthY not poistive... "; exit(0);}

  if (2*lambdaS > lengthX) cerr << "### 2*lambdaS>lengthX. Caution!\n";
  if (2*lambdaS > lengthY) cerr << "### 2*lambdaS>lengthY. Caution!\n";
  nqxMax = lengthX/lambdaS;
  nqyMax = lengthY/lambdaS;
  if (2*nqxMax>nx) cerr << "### 2*qxMax>nx. Caution!\n";
  if (2*nqyMax>ny) cerr << "### 2*qyMax>ny. Caution!\n";
  Lint nqx = nqxMax > nx ? nx : nqxMax;
  Lint nqy = nqyMax > ny ? ny : nqyMax;

  //write the params out
  ofstream output("params.out");
  output << lengthX     << "\t\t# lengthX" << endl;
  output << lengthY     << "\t\t# lengthY" << endl;
  output << nx          << "\t\t# nx"      << endl;
  output << ny          << "\t\t# ny"      << endl;
  output << randSeed    << "\t\t# randSeed"<< endl;
  output << lambdaS     << "\t\t# lambdaShort"<< endl;
  output << lambdaR     << "\t\t# lambdaRoll"<< endl;
  output << Hurst       << "\t\t# Hurst"   << endl;
  output << rmsValueSet << "\t\t# rmsValueSet" << endl;
  output << stepHeight  << "\t\t# stepHeight" << endl;
  output << dumpStep    << "\t\t# dumpStep" << endl;
  output << dumpSpecSmooth << "\t\t# dumpSpecSmooth" << endl;
  output << dumpSpecStep<< "\t\t# dumpSpecStep" << endl;
  output << dumpACF     << "\t\t# dumpACF" << endl;
  output.close();

  //derive the essential params
  dqx = 2*PI/lengthX; dqx2 = dqx*dqx;
  dqy = 2*PI/lengthY; dqy2 = dqy*dqy;
  srand(randSeed);
  nxH = nx/2; nyHP1 = (ny/2)+1;
  sizeCompl = nx*nyHP1;
  sizeReal = 2*sizeCompl; //????
  nxny = nx*ny;
 
  dx = lengthX / nx; dqx2 = TWOPI / lengthX; dqx2 *= dqx2;
  dy = lengthY / ny; dqy2 = TWOPI / lengthY; dqy2 *= dqy2; 

  equilPos  = (double *) fftw_malloc(sizeReal*sizeof(double));
  equilStep = (double *) fftw_malloc(sizeReal*sizeof(double));

  for (Lint k = 0; k < sizeReal; ++k) {
  equilPos[k]  = 0.;
  equilStep[k] = 0.;
  }
}

void addFractal(){
  //allocate memory for height spectrum and define fftw_plan
  Complex* specF = (Complex *) fftw_malloc(sizeCompl*sizeof(Complex));
  double*  specR = (double  *) fftw_malloc(sizeReal*sizeof(double));
  fftw_plan specF2R = 
  fftw_plan_dft_c2r_2d (nx, ny, (fftw_complex*) specF, specR,FFTW_ESTIMATE);
  
  double qRoll = TWOPI/lambdaR;
  double qMax  = TWOPI/lambdaS;
  double msHeight = 0, msGrad = 0, msCurv = 0;
  for (Lint k = 0; k < sizeCompl; ++k){
    double q = getQ(k);
    double qRed = q/qRoll;
    if ( (q>qMax) || (q==0) ) continue;
    double prefac = sqrtSpec(qRed, Hurst);

    double phase =  rrand()*TWOPI;
    specF[k] = prefac * exp( Complex(0,1) * phase);
    get_iqxiqy(k);
    double weight = 2; // weightQ not usable, unless elaDim>0 !!!
    if ( (iqy==0) || (iqy==ny/2) ) {
      weight = 1;
      if ( (iqx==0) || (iqx==nxH) ) {
        // Fourier transform is purely real
        phase = (phase>PI) ? -1 : 1;
        specF[k] = prefac * phase;
      } else if (iqx>nxH) {
        // complex conjugate has already been drawn
        Lint kConj = getLinCIndex(nx-iqx,iqy);
        specF[k]  = conj(specF[kConj]);
    } }
    double q2 = q*q;
    double specLoc = specF[k].real()*specF[k].real()
                   + specF[k].imag()*specF[k].imag();
    specLoc *= weight; msHeight += specLoc;
    specLoc *= q2; msGrad += specLoc;
    specLoc *= q2; msCurv += specLoc;    
  }
  
  int fAlertCont = 0;
  if ( (2*nqxMax>nx) || (2*nqyMax>ny) ){
    fAlertCont = 1;
    double msHeightCont, msGradCont, msCurvCont;
    continuumAnalysis(&msHeightCont, &msGradCont, &msCurvCont);
    msHeight = msHeightCont;
    msGrad   = msGradCont;
    msCurv   = msCurvCont;
  }
 
 /* 
  double rmsGrad = sqrt(msGrad);
  for (Lint k = 0; k < sizeCompl; ++k){ 
    specF[k] /= rmsGrad;
  }
  */
   
  double rmsHeight = sqrt(msHeight);
  rmsHeight *= rmsValueSet*rmsValueSet; // modified here
  for (Lint k = 0; k < sizeCompl; ++k){ 
    specF[k] /= rmsHeight;
  }
  
  // transform to real space
  fftw_execute(specF2R);

  double heightMin = specR[0], heightMax = specR[0], meanHeight = specR[0];
  for (Lint k = 1; k < sizeReal; ++k) {
    heightMin = specR[k] < heightMin ? specR[k] : heightMin;
    heightMax = specR[k] > heightMax ? specR[k] : heightMax;
  }
  double heightDiff = heightMax - heightMin;

  // set lowest value of height to zero, only if ID==0
  //heightMin = heightMax;

  for (Lint k = 0; k < sizeReal; ++k) equilPos[k] += specR[k]-heightMin;
  // dump information on sheet
  ofstream output("params.out",ofstream::app);
  output << rmsValueSet << "\t# rmsHeight Fourier\n";
  output  << sqrt(msGrad/msHeight)   << "\t\t# rmsGrad Fourier\n";
  output << sqrt(msCurv/msHeight)   << "\t\t# rmsCurv Fourier\n";
  output << heightDiff << "\t# absolute height\n";
  output.close();

  // free memory
  fftw_free(specF);
  fftw_free(specR);

}

double getQ(Lint k){
  get_iqxiqy(k);
  int jqx = abs(iqx-nx);
  jqx = (iqx<jqx)? iqx : jqx;
  double q2 = jqx*jqx*dqx2 + iqy*iqy*dqy2;
  return (sqrt(q2));
}

double sqrtSpec(double qRed, double Hurst){
  //double y = pow( 1./sqrt(1.+qRed*qRed) , 1+Hurst);
  double y = pow(1./qRed, 1+Hurst);
  return(y);
}

void continuumAnalysis(double *msHeight, double *msGrad, double *msCurv){
  if (lengthX!=lengthY) cerr << "#Impossible to continuum analysis..." << endl;
  *msHeight = 0; *msGrad = 0; *msCurv = 0;

  double stepsPerTwo = 50;
  double resol = pow(2.,1./stepsPerTwo);
  double q0 = 0.6*TWOPI/lengthX; // q0 should be fully contained in integral
  double qR = TWOPI/lambdaR;

  double qOld = TWOPI/lambdaS;
  while (qOld>q0) {
    double qNew = qOld/resol;
    double qMean = (qOld+qNew)/2;
    double qMean2 = qMean*qMean;
    double dQ = qOld-qNew;
    // add to the mean squares
    double sqrtS = sqrtSpec(qMean/qR, Hurst);
    double localSpec = sqrtS*sqrtS*TWOPI*qMean*dQ;
    *msHeight += localSpec; localSpec *= qMean2;
    *msGrad   += localSpec; localSpec *= qMean2;
    *msCurv   += localSpec;
    // renew q's
    qOld = qNew;
  }
}

void dumpFractal(){
  ofstream dumP("equilPos.dat");
  dumP << "#" << nx << "\t" << ny << "\n\n";
  for (Lint ix=0; ix<=nx; ++ix) {
    for (Lint iy=0; iy<=ny; ++iy) {
      Lint k = getLinRIndex(ix%nx, iy%ny);
      dumP << ix*dx << "\t" << iy*dy << "\t" << equilPos[k] << "\n";
    } dumP << endl;
  } dumP.close();
  
  ofstream dumpExp("exp.dat"); 
  dumpExp << nx << endl;
  dumpExp << ny << endl;
  dumpExp << dx << endl;
  dumpExp << 5.0 << endl;
  
  for (Lint ix = 0; ix < nx; ++ix){
    for (Lint iy = 0; iy < ny; ++iy){
      Lint k = getLinRIndex(ix%nx, iy%ny);
      dumpExp << equilPos[k] << endl;
    }
    dumpExp << endl;
  }
  dumpExp.close();

}

void dumpSpectra(){
  //calculate the height-height-ACF
  vector<vector<double> > equilPos3;
  equilPos3.resize(nx);
  for (ptrdiff_t ix = 0; ix < nx; ix++){equilPos3[ix].resize(ny);}

  for (ptrdiff_t ix = 0; ix < nx; ++ix){
      for (ptrdiff_t iy = 0; iy < ny; ++iy){
        ptrdiff_t k = getLinRIndex(ix,iy);
        equilPos3[ix][iy] = equilPos[k];
      }
  }

    double selfACF[nx];
    for (int ix = 0; ix < nx; ix++){selfACF[ix] = 0.;}

    ofstream outSelfACF("selfACF.dat");
    int nWidth = nx;
    for (ptrdiff_t iWidth = 0; iWidth < nWidth; ++iWidth){
      double deltaR = dx*iWidth;
      int neighN = iWidth;
      for (ptrdiff_t ix = 0; ix < nx; ++ix){
        for (ptrdiff_t iy = 0; iy < ny; ++iy){
          int jxR = (ix+neighN+nx)%nx;
          int jyT = (iy+neighN+ny)%ny;
          double rDummyR = equilPos3[ix][iy]*equilPos3[jxR][iy];
          double rDummyT = equilPos3[ix][iy]*equilPos3[ix][jyT];
          selfACF[iWidth] += rDummyR+rDummyT;
        }
      }
         selfACF[iWidth] /= 2.*nx*ny;
         outSelfACF << deltaR << "\t" << selfACF[iWidth] << endl;
    }
     outSelfACF.close();

     int nR = nWidth;
     //fourier transform
     fftw_complex *in, *out;
     fftw_plan planForward;

     in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nR);
     out= (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nR);

     for (int iR = 0; iR < nR; ++iR){
       in[iR][0] = selfACF[iR];
       in[iR][1] = 0.;
     }

     planForward = fftw_plan_dft_1d(nR,in, out,FFTW_FORWARD,FFTW_ESTIMATE);
     fftw_execute(planForward);

     vector<complex<double> > fourierPos;
     fourierPos.resize(nR);

     for (int qx = 0; qx < nR; qx++){
       fourierPos[qx].real(out[qx][0]);
       fourierPos[qx].imag(out[qx][1]);
     }

     fftw_destroy_plan(planForward) ;
     fftw_free(in);
     fftw_free(out);

     ofstream outSpectrum("spectrum.dat");
     for (int qx = 0; qx < nR; qx++){
       outSpectrum << qx*2.*PI/nx << "\t" << real(fourierPos[qx])*PI/double(qx) << endl;
     }
     outSpectrum.close();
}

void dumpHHACF(){
  vector<vector <double> > equilRealPos;
  equilRealPos.resize(nx);
  for (Lint ix = 0; ix < nx; ++ix){equilRealPos[ix].resize(ny);}
  for (Lint ix = 0; ix < nx; ++ix){
    for (Lint iy = 0; iy < ny; ++iy){
      Lint k = getLinRIndex(ix,iy);
      equilRealPos[ix][iy] = equilPos[k];
    }
  }

  double corrReal[nx/2];

  //calculate the height-difference-autocorrelation function
  int dR = 0;
  int iT = 0;
  double ddR = 1.2;
  int idR = (int) ddR;
  while ((dR+idR) <= nx/2){
    for (int iX = 0; iX < nx; ++iX){
      double rD1 = 0;
      for (int iY = 0; iY < ny; ++iY){
        //draw random angle
	int nSample = 2;
	for (int iSample = 0; iSample < nSample; ++iSample){
	  double randAngle = 2.*PI*rand()*1.0/RAND_MAX;
	  double rDX = dR*cos(randAngle) + rand()*1.0/RAND_MAX;
          double rDY = dR*sin(randAngle) + rand()*1.0/RAND_MAX;
           
          int jX = ((int)(iX+rDX)+2*nx) % nx;
          int jY = ((int)(iY+rDY)+2*ny) % ny;
	  
	  double dHeight = equilRealPos[iX][iY] - equilRealPos[jX][jY];
	  rD1 += dHeight * dHeight / nSample;
	}
      }
      corrReal[iT] += rD1;
    }
  corrReal[iT] /= (nx*2.*ny);

  iT += 1;
  dR += idR;
  ddR *= 1.2;
  idR = ddR;
  }
  
  //dump out the correlation of real surface
  ofstream corr("corr.Real.dat");
  double slope = 1, curvature = 1;
  dR = 0;
  iT = 0;
  ddR = 1.2;
  idR = (int) ddR;
  while ((dR+idR) <= nx/2){
    double dxTrue = idR * dx;
    corr << dx*dR << "\t" << corrReal[iT] << endl;
    iT += 1;
    dR += idR;
    ddR *= 1.2;
    idR = ddR;
  }
  corr.close();

}
