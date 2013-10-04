#include "fintrf.h"
C======================================================================
#if 0
C     
#endif
C     
C======================================================================
C     Gateway routine
      subroutine mexFunction(nlhs, plhs, nrhs, prhs)
      use comin
      use comoutput
 
C     Declarations
      implicit none
      integer ind ! index for input and output
      integer ii  ! index in loops

C     mexFunction arguments:
      mwPointer plhs(*), prhs(*)
      integer nlhs, nrhs
 
C     Function declarations:
      mwPointer mxCreateDoubleMatrix, mxGetPr
      mwPointer, external :: mxIsComplex, mxIsNumeric
      mwPointer, external :: mxGetM, mxGetN
 
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
      if(nrhs .ne. 11) then
         call mexErrMsgIdAndTxt ('MATLAB:matsq:nInput',
     +                           'One input required.')
      elseif(nlhs .ne. 1) then
         call mexErrMsgIdAndTxt ('MATLAB:matsq:nOutput',
     +                           'One output required.')
      endif

C      x_pr = mxGetPr(prhs(1))
      ind = 1
      call input_size(ind)
      call mxCopyPtrToReal8(mxGetPr(prhs(ind)), XC ,size)
      ind = 2
      call input_size(ind)
      call mxCopyPtrToReal8(mxGetPr(prhs(ind)), PZL ,size)
      ind = 3
      call input_size(ind)
      call mxCopyPtrToReal8(mxGetPr(prhs(ind)), zfirst ,size)
      ind = 4
      call input_size(ind)
      call mxCopyPtrToReal8(mxGetPr(prhs(ind)), DN ,size)
      ind = 5
      call input_size(ind)
      call mxCopyPtrToReal8(mxGetPr(prhs(ind)), TA ,size)
      ind = 6
      call input_size(ind)
      call mxCopyPtrToReal8(mxGetPr(prhs(ind)), DD ,size)
      ind = 7
      call input_size(ind)
      call mxCopyPtrToReal8(mxGetPr(prhs(ind)), AA ,size)
      ind = 8
      call input_size(ind)
      call mxCopyPtrToReal8(mxGetPr(prhs(ind)), ASS ,size)
      ind = 9
      call input_size(ind)
      call mxCopyPtrToReal8(mxGetPr(prhs(ind)), VD ,size)
      ind = 10
      call input_size(ind)
      call mxCopyPtrToReal8(mxGetPr(prhs(ind)), PM ,size)
      ind = 11
      call input_size(ind)
      call mxCopyPtrToReal8(mxGetPr(prhs(ind)), ZM ,size)


C     Call the computational subroutine.
      call whamp_engine

C     Set output using comoutput


C     Create matrix for the return argument.
C      plhs(1) = mxCreateDoubleMatrix(m,n,0)
      ind = 1
      plhs(1) = mxCreateDoubleMatrix(1,kperpSize,0)
      y_pr = mxGetPr(plhs(ind))
      call mxCopyReal8ToPtr(kperpOUT,y_pr,kperpSize)     
      plhs(1) = mxCreateDoubleMatrix(1,10,1)


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
C      if(size .gt. numel) then
C         call mexErrMsgIdAndTxt ('MATLAB:matsq:mSize',
C     +                           'Row * column must be <= 1000.')
C      endif
      
C     Check that the array is numeric (not strings).
C      if(mxIsNumeric(prhs(ind)) .eq. 0) then
C         call mexErrMsgIdAndTxt ('MATLAB:matsq:NonNumeric',
C     +                           'Input must be a numeric array.')
C      endif
      end !input_size
      
      end
