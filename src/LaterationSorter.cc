/*
  Copyright (C) 2004 Robert K. Harle

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

  Email: robert.harle@achilles.org
  Email: rkh23@cantab.net
*/

#include "NLMAP/LaterationSorter.hh"
#include "NLMAP/NLMAPExceptions.hh"

#include <cmath>

//----------------------------------------
// Evaluate the studentized residuals for
// the current parameters and return the largest
// POSITIVE residual
//----------------------------------------
int LaterationSorter::GetMaxResidualIndex() {
  LaterationData *mLatData = dynamic_cast<LaterationData *>(mData);

  int n = mLatData->GetInputDataSize();

  // Calculate derivates of residuals
  REAL res_der[n][3];
  REAL res[n];
  for (int i=0;i<n;i++) {
    REAL d = sqrt( (mFunc->GetParams()[0]-mLatData->GetX(i))*(mFunc->GetParams()[0]-mLatData->GetX(i)) +
		   (mFunc->GetParams()[1]-mLatData->GetY(i))*(mFunc->GetParams()[1]-mLatData->GetY(i)) +
		   (mFunc->GetParams()[2]-mLatData->GetZ(i))*(mFunc->GetParams()[2]-mLatData->GetZ(i)) );
      res[i]=mLatData->GetMeasurement(i)-d;
      res_der[i][0] = (mFunc->GetParams()[0]-mLatData->GetX(i))/d;
      res_der[i][1] = (mFunc->GetParams()[1]-mLatData->GetY(i))/d;
      res_der[i][2] = (mFunc->GetParams()[2]-mLatData->GetZ(i))/d;
  }

  // Multiply this by its transpose to make a 3x3 matrix
  REAL mat[3][3];
  for (int j=0; j<3;j++) {
    for (int i=0;i<3; i++) {
      mat[i][j] = 0.0;
      for (int k=0; k<n; k++) {
	mat[i][j]+=res_der[k][i]*res_der[k][j];
      }
    }
  }

  // Invert this matrix (closed form for 3x3)
  REAL det = mat[0][0] * (mat[1][1]*mat[2][2]-mat[1][2]*mat[2][1]) -
    mat[0][1] * (mat[1][0]*mat[2][2]-mat[1][2]*mat[2][0]) -+
    mat[0][2] * (mat[1][0]*mat[2][1]-mat[1][1]*mat[2][0]);

  if (det==0) throw SingularMatrix();

  REAL mat_inv[3][3];
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      int minrow,maxrow,mincol, maxcol;
      if (i==0) {minrow=1;maxrow=2;}
      else if (i==1) {minrow=0;maxrow=2;}
      else {minrow=0;maxrow=1;}
      if (j==0) {mincol=1;maxcol=2;}
      else if (j==1) {mincol=0;maxcol=2;}
      else {mincol==0;maxcol=1;}

      mat_inv[i][j]=(mat[minrow][mincol]*mat[maxrow][maxcol] -
		     mat[minrow][maxcol]*mat[maxrow][mincol])/det;
      if ((i+j)%2==1) mat_inv[i][j]*=-1;
    }
  }


  REAL h1[n][3];
  for (int c=0; c<n; c++) {
    h1[c][0]=res_der[c][0]*mat_inv[0][0] + res_der[c][1]*mat_inv[1][0] + res_der[c][2]*mat_inv[2][0];
    h1[c][1]=res_der[c][0]*mat_inv[0][1] + res_der[c][1]*mat_inv[1][1] + res_der[c][2]*mat_inv[2][1];
    h1[c][2]=res_der[c][0]*mat_inv[0][2] + res_der[c][1]*mat_inv[1][2] + res_der[c][2]*mat_inv[2][2];
  }

  REAL hat_diag[n];
  for (int c=0; c<n; c++) {
    hat_diag[c]= h1[c][0]*res_der[c][0] +
      h1[c][1]*res_der[c][1] +
      h1[c][2]*res_der[c][2];
  }

  REAL maxr=0.0;
  int maxidx=-1;
  // Now studentize the residuals
  for (int c=0; c<n; c++) {
    if (hat_diag[c]==1.0) throw ResidualSorterException();
    REAL sr=0.0;

    // Externally studentize if we can
    if (n>4) {
      REAL sum=0.0;
      for (int i=0; i<n; i++) {
	if (i!=c) sum+=res[i]*res[i];
      }
      sum/=(n-3-1);

      if (sum==0.0) throw ResidualSorterException();
      sr=res[c]/(sqrt(sum)*sqrt(1.0-hat_diag[c]));
    }

    // Internally studentize
    else {
      if (mFunc->GetError()==-1.0) throw ResidualSorterException();
      sr = res[c]/(mFunc->GetError()*sqrt(1.0-hat_diag[c]));
    }

    if (maxidx==-1 || sr>maxr) {
      maxidx=c;
      maxr=sr;
    }
  }
  return maxidx;

};
