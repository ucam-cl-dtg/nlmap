#include <NLMAP/LinearisedLS.hh>
//#include <NLMAP/gaussianelimination.hh>
#include <NLMAP/GaussJordanEliminator.hh>
#include <iostream>

LinearisedLS::LinearisedLS(
			   const REAL* x,
			   const REAL* y,
			   const REAL* z,
			   const REAL* d,
			   const REAL* sigma,
			   const int n) : mA(0), mB(0) {

  // calculate the A matrix - 3 columns, nC2 rows
  mDim = (int)(n*(n-1)/2.0);
  if (mDim<=0.0) throw NLMAPException();

  mA = new REAL*[mDim];
  for (int i=0; i<mDim; i++) {
    mA[i] = new REAL[3];
  }

  mB = new REAL[mDim];

  // Now go through each pair of two
  int count=0;
  for (int p1=0; p1<n; p1++) {
    for (int p2=0; p2<p1; p2++) {
      mA[count][0]=2.0*(x[p1]-x[p2]);
      mA[count][1]=2.0*(y[p1]-y[p2]);
      mA[count][2]=2.0*(z[p1]-z[p2]);
      mB[count]=d[p2]*d[p2] - x[p2]*x[p2]
	- y[p2]*y[p2] - z[p2]*z[p2] -
	d[p1]*d[p1] + x[p1]*x[p1] +
	y[p1]*y[p1] + z[p1]*z[p1];
      count++;
    }
  }

  if (count!=mDim) throw NLMAPException();
	
}

LinearisedLS::~LinearisedLS() {
  if (mA) {
    for (int i=0; i<mDim; i++) {
      delete[] mA[i];
    }
    delete[] mA;
  }
  
  if (mB) {
    delete[] mB;
  }

}


XYZData LinearisedLS::GetPosition() {
  // Transpose A, AT
  REAL AT[3][mDim];
  REAL **InvATA = new REAL*[3];
  for (int i=0; i<3; i++) {
    InvATA[i]=new REAL[3];
  }
  

  for (int i=0; i<mDim; i++) {
    for (int j=0; j<3; j++) {
      AT[j][i]=mA[i][j];
    }
  }

  
  // Make ATA and invert it
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      InvATA[i][j] = 0.0;
      for (int k=0; k<mDim; k++) {
	InvATA[i][j] += AT[i][k]*AT[j][k];
      }
    }
  }
  REAL **y = new REAL*[3]; // Dummy data
  for (int i=0; i<3; i++) y[i] = new REAL[1];
  y[0][0]=1.0;
  y[1][0]=1.0;
  y[2][0]=1.0;

  GaussJordanEliminator::Eliminate(InvATA,3,y,3);

  for (int i=0; i<3; i++) delete[] y[i];
  delete[] y;
  


  // Multiply out  P=(InvATA)ATB
  REAL ap[3],p[3];

  for (int r=0; r<3; r++) {
    ap[r]=0.0;
    for (int i=0; i<mDim; i++) {
      ap[r]+=mB[i]*AT[r][i];
    }
  }

  for (int i=0; i<3; i++) {
    p[i]=0;
    for (int j=0; j<3; j++) {
      p[i]+=InvATA[i][j]*ap[j];
    }
  }

  for (int i=0; i<3; i++) {
    delete[] InvATA[i];
  }
  delete[] InvATA;

  XYZData pd;
  pd.x = p[0];
  pd.y = p[1];
  pd.z = p[2];
  return pd;
}
