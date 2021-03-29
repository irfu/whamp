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
      integer,external ::  mxClassIDFromClassName
 
C     Function declarations:
      mwPointer mxCreateNumericMatrix
      mwPointer mxCreateDoubleMatrix, mxGetPr, mxGetPi
      mwPointer, external :: mxIsComplex, mxIsNumeric
      mwPointer, external :: mxGetM, mxGetN
 
C     Pointers to input/output mxArrays:
      mwPointer x_pr       ! input
      mwPointer y_pr, y_pi ! output
 
C     Array information:
      mwSize m, n
      mwSize size
 
C-----------------------------------------------------------------------
C     Check for proper number of arguments. 
      if(nrhs .ne. 13) then
         call mexErrMsgIdAndTxt ('MATLAB:matsq:nInput',
     +                           'One input required.')
      elseif(nlhs .ne. 21) then
         call mexErrMsgIdAndTxt ('MATLAB:matsq:nOutput',
     +                           'One output required.')
      endif

C      x_pr = mxGetPr(prhs(1))
      ind = 1
      call input_size(ind)
      x_pr = mxGetPr(prhs(ind))
      call mxCopyPtrToReal8(x_pr, XC, size)
      ind = 2
      call input_size(ind)
      x_pr = mxGetPr(prhs(ind))
      call mxCopyPtrToReal8(x_pr, PZL, size)
      ind = 3
      call input_size(ind)
      x_pr = mxGetPr(prhs(ind))
      call mxCopyPtrToReal8(x_pr, cycleZFirst, size)
      ind = 4
      call input_size(ind)
      x_pr = mxGetPr(prhs(ind))
      call mxCopyPtrToReal8(x_pr, DN, size)
      ind = 5
      call input_size(ind)
      x_pr = mxGetPr(prhs(ind))
      call mxCopyPtrToReal8(x_pr, TA, size)
      ind = 6
      call input_size(ind)
      x_pr = mxGetPr(prhs(ind))
      call mxCopyPtrToReal8(x_pr, DD, size)
      ind = 7
      call input_size(ind)
      x_pr = mxGetPr(prhs(ind))
      call mxCopyPtrToReal8(x_pr, AA, size)
      ind = 8
      call input_size(ind)
      x_pr = mxGetPr(prhs(ind))
      call mxCopyPtrToReal8(x_pr, ASS, size)
      ind = 9
      call input_size(ind)
      x_pr = mxGetPr(prhs(ind))
      call mxCopyPtrToReal8(x_pr, VD, size)
      ind = 10
      call input_size(ind)
      x_pr = mxGetPr(prhs(ind))
      call mxCopyPtrToReal8(x_pr, PM, size)
      ind = 11
      call input_size(ind)
      x_pr = mxGetPr(prhs(ind))
      call mxCopyPtrToReal8(x_pr, ZM, size)
      ind = 12
      call input_size(ind)
      x_pr = mxGetPr(prhs(ind))
      call mxCopyPtrToReal8(x_pr, XOI, size)
      ind = 13
      call input_size(ind)
      x_pr = mxGetPr(prhs(ind))
      call mxCopyPtrToInteger4(x_pr,maxIterations,size)

C     Call the computational subroutine.
      call whamp_engine

C     Set output using comoutput

C     Create matrix for the return argument.
C     plhs(1) = mxCreateDoubleMatrix(m,n,0)
      size = kperpSize * kparSize
      ind = 1
      plhs(ind) = mxCreateDoubleMatrix(1,kperpSize,0)
      y_pr = mxGetPr(plhs(ind))
      call mxCopyReal8ToPtr(kperpOut,y_pr,kperpSize)
      ind = 2
      plhs(ind) = mxCreateDoubleMatrix(1,kparSize,0)
      y_pr = mxGetPr(plhs(ind))
      call mxCopyReal8ToPtr(kparOut,y_pr,kparSize)
      ind = 3
      plhs(ind) = mxCreateDoubleMatrix(kperpSize,kparSize,1)
      y_pr = mxGetPr(plhs(ind))
      y_pi = mxGetPi(plhs(ind))
      call mxCopyComplex16ToPtr(fOUT,y_pr,y_pi,size)
      ind = 4
      plhs(ind) = mxCreateDoubleMatrix(kperpSize,kparSize,1)
      y_pr = mxGetPr(plhs(ind))
      y_pi = mxGetPi(plhs(ind))
      call mxCopyComplex16ToPtr(ExOUT,y_pr,y_pi,size)
      ind = 5
      plhs(ind) = mxCreateDoubleMatrix(kperpSize,kparSize,1)
      y_pr = mxGetPr(plhs(ind))
      y_pi = mxGetPi(plhs(ind))
      call mxCopyComplex16ToPtr(EyOUT,y_pr,y_pi,size)
      ind = 6
      plhs(ind) = mxCreateDoubleMatrix(kperpSize,kparSize,1)
      y_pr = mxGetPr(plhs(ind))
      y_pi = mxGetPi(plhs(ind))
      call mxCopyComplex16ToPtr(EzOUT,y_pr,y_pi,size)
      ind = 7
      plhs(ind) = mxCreateDoubleMatrix(kperpSize,kparSize,1)
      y_pr = mxGetPr(plhs(ind))
      y_pi = mxGetPi(plhs(ind))
      call mxCopyComplex16ToPtr(BxOUT,y_pr,y_pi,size)
      ind = 8
      plhs(ind) = mxCreateDoubleMatrix(kperpSize,kparSize,1)
      y_pr = mxGetPr(plhs(ind))
      y_pi = mxGetPi(plhs(ind))
      call mxCopyComplex16ToPtr(ByOUT,y_pr,y_pi,size)
      ind = 9
      plhs(ind) = mxCreateDoubleMatrix(kperpSize,kparSize,1)
      y_pr = mxGetPr(plhs(ind))
      y_pi = mxGetPi(plhs(ind))
      call mxCopyComplex16ToPtr(BzOUT,y_pr,y_pi,size)
      ind = 10
      plhs(ind) = mxCreateDoubleMatrix(kperpSize,kparSize,1)
      y_pr = mxGetPr(plhs(ind))
      y_pi = mxGetPi(plhs(ind))
      call mxCopyComplex16ToPtr(SxOUT,y_pr,y_pi,size)
      ind = 11
      plhs(ind) = mxCreateDoubleMatrix(kperpSize,kparSize,1)
      y_pr = mxGetPr(plhs(ind))
      y_pi = mxGetPi(plhs(ind))
      call mxCopyComplex16ToPtr(SyOUT,y_pr,y_pi,size)
      ind = 12
      plhs(ind) = mxCreateDoubleMatrix(kperpSize,kparSize,1)
      y_pr = mxGetPr(plhs(ind))
      y_pi = mxGetPi(plhs(ind))
      call mxCopyComplex16ToPtr(SzOUT,y_pr,y_pi,size)
      ind = 13
      plhs(ind) = mxCreateDoubleMatrix(kperpSize,kparSize,0)
      y_pr = mxGetPr(plhs(ind))
      call mxCopyReal8ToPtr(EBOUT,y_pr,size)
      ind = 14
      plhs(ind) = mxCreateDoubleMatrix(kperpSize,kparSize,0)
      y_pr = mxGetPr(plhs(ind))
      call mxCopyReal8ToPtr(VGPOUT,y_pr,size)
      ind = 15
      plhs(ind) = mxCreateDoubleMatrix(kperpSize,kparSize,0)
      y_pr = mxGetPr(plhs(ind))
      call mxCopyReal8ToPtr(VGZOUT,y_pr,size)
      ind = 16
      plhs(ind) = mxCreateDoubleMatrix(kperpSize,kparSize,0)
      y_pr = mxGetPr(plhs(ind))
      call mxCopyReal8ToPtr(SGPOUT,y_pr,size)
      ind = 17
      plhs(ind) = mxCreateDoubleMatrix(kperpSize,kparSize,0)
      y_pr = mxGetPr(plhs(ind))
      call mxCopyReal8ToPtr(SGZOUT,y_pr,size)
      ind = 18
      plhs(ind) = mxCreateDoubleMatrix(kperpSize,kparSize,0)
      y_pr = mxGetPr(plhs(ind))
      call mxCopyReal8ToPtr(uOUT,y_pr,size)

      ind = 19
      plhs(ind) = mxCreateNumericMatrix(kperpSize,kparSize,
     +      mxClassIDFromClassName('int32'),0)
      y_pr = mxGetPr(plhs(ind))
      call mxCopyInteger4ToPtr(flagSolutionFoundOUT,y_pr,size)
      ind = 20
      plhs(ind) = mxCreateNumericMatrix(kperpSize,kparSize,
     +      mxClassIDFromClassName('int32'),0)
      y_pr = mxGetPr(plhs(ind))
      call mxCopyInteger4ToPtr(flagTooHeavilyDampedOUT,y_pr,size)
      ind = 21
      plhs(ind) = mxCreateNumericMatrix(kperpSize,kparSize,
     +      mxClassIDFromClassName('int32'),0)
      y_pr = mxGetPr(plhs(ind))
      call mxCopyInteger4ToPtr(flagNoConvergenceOUT,y_pr,size)

      return

      contains
      subroutine input_size(ind)
      integer ind
C     Get the size of the input array.
      m = mxGetM(prhs(ind))
      n = mxGetN(prhs(ind))
      size = m*n

      end !input_size
      
      end
