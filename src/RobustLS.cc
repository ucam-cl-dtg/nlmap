#include <NLMAP/RobustLS.hh>
#include <NLMAP/GaussJordanEliminator.hh>

#include <iostream>
#include <cmath>

#define CONST 1.5


RobustLS::RobustLS(
		   REAL *x,
		   REAL *y,
		   REAL *z,
		   REAL *d,
		   REAL *s,
		   int n) : 
  mX(x), mY(y), mZ(z), mN(n) {


// calculate the A matrix - 3 columns, nC2 rows
  mDim = (int)(n*(n-1)/2.0);
  if (mDim<=0.0) throw NLMAPException();

  mA = new REAL*[mDim];
  for (int i=0; i<mDim; i++) {
    mA[i] = new REAL[3];
  }

  mB = new REAL[mDim];
  mS = new REAL[mDim];

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

      mS[count] = 2.0*sqrt(s[p1]*s[p1]*d[p1]*d[p1] +
			   s[p2]*s[p2]*d[p2]*d[p2] );

      count++;
    }
  }

  if (count!=mDim) throw NLMAPException();

}




RobustLS::~RobustLS() {
   if (mA) {
    for (int i=0; i<mDim; i++) {
      delete[] mA[i];
    }
    delete[] mA;
  }
  
  if (mB) {
    delete[] mB;
  }

  if (mH) {
    for (int i=0; i<mDim; i++) delete[] mH[i];
    delete[] mH;
  }
}





void RobustLS::LeastSquares() {
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
  


  // Get (InvATA)AT
  REAL mult[3][mDim];
  for (int i=0; i<mDim; i++) {
    mult[0][i]=0.0;
    mult[1][i]=0.0;
    mult[2][i]=0.0;
    for (int j=0; j<3; j++) {
      mult[0][i] += InvATA[0][j]*AT[j][i];
      mult[1][i] += InvATA[1][j]*AT[j][i];
      mult[2][i] += InvATA[2][j]*AT[j][i];
    }
  }
  

  // Calculate the point
  for (int j=0; j<3; j++) {
    mP[j]=0.0;
    for (int i=0; i<mDim; i++) 
      mP[j]+= mult[j][i] * mB[i];
  }

 //   // Store the hat matrix mH = mA*mult diagonals
//    mH = new REAL[mDim];
//    for (int k=0; k<mDim; k++) {
//      mH[k]=0.0;
//      for (int q=0;q<3; q++) {
//        mH[k] += mult[q][k] * mA[k][q];
//      }
//    }


  mH = new REAL*[mDim];
  for (int i=0; i<mDim; i++) mH[i]=new REAL[mDim];

  for (int r=0; r<mDim; r++) {
    for (int c=0; c<mDim; c++) {
      mH[r][c] = 0.0;
      for (int i=0; i<3; i++) 
	mH[r][c] += mA[r][i] * mult[i][c];
    }
  }

  for (int i=0; i<3; i++) {
    delete[] InvATA[i];
  }
  delete[] InvATA;
}



XYZData RobustLS::GetPosition(const int maxit) {
 
  int iteration=0;
  // while not converged
  REAL lastErr=1.0, errChange=1.0;

  while (iteration<maxit && errChange>0.00001) {
    LeastSquares();

    // Figure out the residuals
    REAL bhat[mDim];
    REAL sumsq=0.0;
    for (int i=0;i<mDim; i++) {
      bhat[i] = 0.0;
      for (int j=0; j<mDim; j++) {
	bhat[i] += mH[i][j] * mB[j];
      }
    }

    for (int i=0; i<mDim; i++) {
      REAL res=mB[i]-bhat[i];
      sumsq+=res*res;

      if (res < -CONST*mS[i]) {
	mB[i]=res-CONST*mS[i];
      }
      if (res > CONST*mS[i]) {
	mB[i] = res + CONST*mS[i];
      }
    }

    // Calculate the new variances
    for (int i=0; i<mDim; i++) {
      mS[i]=sqrt((1.0-mH[i][i])*sumsq/(mDim-3));
    }

    errChange = fabs(sumsq/(mDim-3)-lastErr)/lastErr;
    lastErr = sumsq/(mDim-3);

    std::cout << "Iteration: " << iteration << " " << mP[0] << " " << mP[1] 
	      << " " << mP[2] << " " << lastErr << std::endl;

    iteration++;
  }

  XYZData pd;
  pd.x = mP[0];
  pd.y = mP[1];
  pd.z = mP[2];
  return pd;
}
