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
#ifndef NLM_EXCEPTIONS_HH
#define NLM_EXCEPTIONS_HH

#include <exception>
#include <cstring>

//--------------------------------
// Base exception class for NLMaP
//--------------------------------
class NLMAPException : public std::exception {
public:
  virtual char * what() { return "Exception thrown in NLMAP"; }
};




//--------------------------------
// Exceptions thrown by NonLinearModel
// and accompanying classes
//--------------------------------

// General failure to converge
class FailedToConverge : public NLMAPException {
public:
  char * what() { return "Nonlinear model failed to converge"; }
};

// Model failed to converge fast enough
class MaxIterations : public NLMAPException {
public:
  char * what() { return "Nonlinear model took too many iterations"; }
};
 

// Encountered a singular matrix
class SingularMatrix : public NLMAPException {
public:
  char * what() { return "Singular matrix encountered"; }
};

// General error in the FitFunction implementation
class FitFunctionException : public NLMAPException {
public:
  char * what() { return "Exception in FitFunction implementation"; }
};

class IndexOutOfBounds : public NLMAPException {
public:
  char * what() { return "Attempt to access non-existent data"; }
};


//--------------------------------
// Exceptions thrown by IterativeModel
// and accompanying classes
//--------------------------------

// General error in the ResidualSorter implementation
class ResidualSorterException : public NLMAPException {
public:
  char * what() { return "Exception when trying to evaluate residuals"; }
};



class ModelingError : public NLMAPException {
public:
  ModelingError(char *f){strcpy(msg,f);}
  ModelingError() {strcpy(msg,"Error in iterative modeling algorithm");}
  char * what() { return msg;}
private: 
  char  msg[512];
};


class ModelingFailure : public NLMAPException {
public:
  char * what() { return "Unable to model data to sufficient accuracy"; }
};

class InvalidData : public NLMAPException {
public:
  char * what() { return "Attempt to add invalid data"; }
};


#endif
