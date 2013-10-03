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
 
C     Declarations
      implicit none

C     mexFunction arguments:
      mwPointer plhs(*), prhs(*)
      integer nlhs, nrhs
 
C     Function declarations:
      mwPointer mxCreateDoubleMatrix, mxGetPr
      integer mxIsNumeric
      mwPointer mxGetM, mxGetN
 
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

C     Get the size of the input array.
      m = mxGetM(prhs(1))
      n = mxGetN(prhs(1))
      size = m*n

C     Column * row must be smaller than numel
      if(size .gt. numel) then
         call mexErrMsgIdAndTxt ('MATLAB:matsq:mSize',
     +                           'Row * column must be <= 1000.')
      endif
      
C     Check that the array is numeric (not strings).
      if(mxIsNumeric(prhs(1)) .eq. 0) then
         call mexErrMsgIdAndTxt ('MATLAB:matsq:NonNumeric',
     +                           'Input must be a numeric array.')
      endif

C     Create matrix for the return argument.
      plhs(1) = mxCreateDoubleMatrix(m,n,0)
      x_pr = mxGetPr(prhs(1))
      y_pr = mxGetPr(plhs(1))
      call mxCopyPtrToReal8(x_pr,x,size)

C     Call the computational subroutine.
      call matsq(y, x, m, n)

C     Load the data into y_pr, which is the output to MATLAB.
      call mxCopyReal8ToPtr(y,y_pr,size)     

      return
      end

C-----------------------------------------------------------------------
C     Computational routine

      subroutine matsq(y, x, m, n)
      mwSize m, n, i, j
      real*8 x(m,n), y(m,n)
C     
      do 20 i=1,m
         do 10 j=1,n
            y(i,j)= x(i,j)**2
 10      continue
 20   continue
      return
      end
