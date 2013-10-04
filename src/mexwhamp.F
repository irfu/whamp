#include "fintrf.h"
C======================================================================
#if 0
C     
C     matsq.F
C     .F file needs to be preprocessed to generate .for equivalent
C     
#endif
C     
C     matsq.f
C
C     squares the input matrix
      
C     This is a MEX-file for MATLAB.
C     Copyright 1984-2011 The MathWorks, Inc.
C     $Revision: 1.13.2.5 $
C======================================================================
C     Gateway routine
      subroutine mexFunction(nlhs, plhs, nrhs, prhs)
      use comin
      use comoutput
 
C     Declarations
      implicit none

C     mexFunction arguments:
      mwPointer plhs(*), prhs(*)
      integer nlhs, nrhs
 
C     Function declarations:
      mwPointer mxCreateDoubleMatrix, mxGetPr
!      mwPointer, external ::  mxCreateDoubleMatrix, mxGetPr
      mwPointer, external :: mxIsComplex, mxIsNumeric
      mwPointer, external ::  mxGetM, mxGetN
 
C     Pointers to input/output mxArrays:
      mwPointer x_pr, y_pr
 
C     Array information:
      mwPointer m, n
      mwSize size
 
C     Arguments for computational routine:
C     Maximum size = numel
      integer*4 numel
      parameter(numel = 1000)
      real*8  x(numel), y(numel)
 
C-----------------------------------------------------------------------
C     Check for proper number of arguments. 
      if(nrhs .ne. 1) then
         call mexErrMsgIdAndTxt ('MATLAB:matsq:nInput',
     +                           'One input required.')
      elseif(nlhs .ne. 1) then
         call mexErrMsgIdAndTxt ('MATLAB:matsq:nOutput',
     +                           'One output required.')
      endif

C      x_pr = mxGetPr(prhs(1))
      call input_size(1)
      call mxCopyPtrToReal8(mxGetPr(prhs(1)),DN,size)



C     Create matrix for the return argument.
C      plhs(1) = mxCreateDoubleMatrix(m,n,0)
C      plhs(1) = mxCreateDoubleMatrix(1,10,0)
      plhs(1) = mxCreateDoubleMatrix(1,10,1)
C      x_pr = mxGetPr(prhs(1))
      y_pr = mxGetPr(plhs(1))
C      call mxCopyPtrToReal8(x_pr,x,size)

C     Set input using comin

C     Call the computational subroutine.
      !call matsq(y, x, m, n)

C     Set output using comoutput


C     Load the data into y_pr, which is the output to MATLAB.
      call mxCopyReal8ToPtr(DN,y_pr,10)     

      return

      contains
      subroutine input_size(ind)
      integer ind
C     Get the size of the input array.
      m = mxGetM(prhs(ind))
      n = mxGetN(prhs(ind))
      size = m*n

C     Column * row must be smaller than numel
      if(size .gt. numel) then
         call mexErrMsgIdAndTxt ('MATLAB:matsq:mSize',
     +                           'Row * column must be <= 1000.')
      endif
      
C     Check that the array is numeric (not strings).
      if(mxIsNumeric(prhs(ind)) .eq. 0) then
         call mexErrMsgIdAndTxt ('MATLAB:matsq:NonNumeric',
     +                           'Input must be a numeric array.')
      endif

      end !input_size
      
      end

C-----------------------------------------------------------------------
C     Set whamp input


