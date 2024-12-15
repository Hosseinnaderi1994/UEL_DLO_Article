************************************************************************
! This User Element subroutine has been written by Hossein Naderi (naderise@msu.edu)
! one can contact me via hosseinnaderi6270949@gmail.com.
! This UEL has been written for the article"Constitutive modeling of... 
! diffusion-limited oxidation coupled with a large deformation theory for polymer degradation".
      
! User element for DLO coupled with large elastic deformation theory.  
!This is for plane strain, axisymetric element .
!
! nodal variables are the displacements and the oxygen concentration.
! 
! This subroutine is for the following element types: UPE4 & UAX4
!  
! MASTER ELEMENT :
!      
!              A eta (=xi_2)
!  4-node      |
!   quad       |Face 3
!        4-----------3
!        |     |     |
!        |     |     |
!  Face 4|     ------|---> xi (=xi_1)
!        |           | Face2
!        |           |
!        1-----------2
!          Face 1
!
***********************************************************************
! User element statement in the input file (set U? values as needed FOR plain strain or axisymmetric element):
!
!  2D elements
!  *User Element,Nodes=4,Type=U?,Iproperties=2,Properties=8,Coordinates=2,Variables=?,Unsymm
!  1,2,11
!
!********************************************************************************************
!     Global SDV's (used for visualization): extent of oxidations, Von Mises stress, Shear modulus
!
!     Local SDV's (used for the solution procedure)
!       j = 0
!       do k = 1,nIntPt
!          svars(1+j) = P_sci ---- P_sci at integ pt k
!          svars(2+j) = P_ref ---- P_ref at integ pt k     
!          j = j + nlSdv
!       end loop over k
!
!     In the input file, set 'User output variables'= number of global SDV's
!
!     In the input file, set 'ngSdv'= number of global SDV's
!
!     In the input file, set 'nlSdv'= number of local SDV's
!
!     In the input file, set 'varibles'=(nlSdv*nIntPt)
!*********************************************************************************************
!
!     Material Properties Vector
** G     Kbulk   gama  chi   Vmol      Rgas 
!0.5e9,  1.0e6, -0.01, 0.1, 0.14675, 8.3145
! theta, nlSdv ngSdv
!  373,    2,   4
!     --------------------------------------------------------------
!     G    = props(1)     ! Shear modulus
!     Kbulk= props(2)  ! Bulk modulus
!     gama = props(3)  ! Shrinkage parameter
!     chi  = props(4)  ! Chi parameter
!     Vmol = props(5)  ! Volume of a mole of fluid particles
!     Rgas = props(6)  ! Universal gas constant
!     nlSdv= jprops(1) ! Number of local sdv's per integ pt
!     ngSdv= jprops(2) ! Number of global sdv's per integ pt
!
!***********************************************************************

      module global

      ! This module is used to transfer SDV's from the UEL to the UVARM so that SDV's can be visualized on a dummy mesh
      !
      !  globalSdv(X,Y,Z):
      !   X - element pointer
      !   Y - integration point pointer
      !   Z - SDV pointer
      !
      !  numElem:
      !   Total number of elements in the real mesh, the dummy mesh needs to have the same number of elements, and 
      !   the dummy mesh needs to have the same number of integ points.  You must set that parameter value here.
      !
      !  ElemOffset:
      !   Offset between element numbers on the real mesh and dummy mesh.  That is set in the input file, and 
      !    that value must be set here the same.

      integer numElem,ElemOffset,err

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Set the number of UEL elements used here
      parameter(numElem=110)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Set the offset here for UVARM plotting, must match input file!
      parameter(ElemOffset=1000)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real*8, allocatable :: globalSdv(:,:,:)
   
      end module global

***********************************************************************

      SUBROUTINE UVARM(UVAR,DIRECT,T,TIME,DTIME,CMNAME,ORNAME,
     1 NUVARM,NOEL,NPT,LAYER,KSPT,KSTEP,KINC,NDI,NSHR,COORD,
     2 JMAC,JMATYP,MATLAYO,LACCFLA)

      ! This subroutine is used to transfer SDV's from the UEL onto the dummy mesh for viewing.  Note that an offset of
      !  ElemOffset is used between the real mesh and the dummy mesh. 
     
      use global

      include 'ABA_PARAM.INC'

      CHARACTER*80 CMNAME,ORNAME
      CHARACTER*3 FLGRAY(15)
      DIMENSION UVAR(NUVARM),DIRECT(3,3),T(3,3),TIME(2)
      DIMENSION ARRAY(15),JARRAY(15),JMAC(*),JMATYP(*),COORD(*)
C      uvars for visualization(!P_sci,!P_ref)
        uvar(1) = globalSdv(noel-ElemOffset,npt,1) 
        uvar(2) = globalSdv(noel-ElemOffset,npt,2)
        uvar(3) = globalSdv(noel-ElemOffset,npt,3) 
        uvar(4) = globalSdv(noel-ElemOffset,npt,4) 
      return
      end subroutine uvarm

****************************************************************************

      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     +     PROPS,NPROPS,coords,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,
     +     TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,
     +     PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,
     +     NJPROP,PERIOD)

      use global
*
      IMPLICIT NONE
     
*     VARIABLES DEFINED IN UEL, PASSED BACK TO ABAQUS
*
      REAL(8) :: RHS,AMATRX,SVARS,ENERGY
*
*     VARIABLES PASSED INTO UEL 
*
      REAL(8) :: PROPS,coords,Uall,DUall,Vel,Accn,TIME,
     1  DTIME,PARAMS,ADLMAG,PREDEF,DDLMAG,PNEWDT,PERIOD
      INTEGER :: NDOFEL,NRHS,NSVARS,NPROPS,MCRD,NNODE,JTYPE,KSTEP,KINC,
     1  JELEM,NDLOAD,JDLTYP,NPREDF,LFLAGS,MLVARX,MDLOAD,JPROPS,NJPROP
*
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     1 SVARS(*),ENERGY(8),coords(MCRD,NNODE),Uall(NDOFEL),
     2 DUall(MLVARX,*),Vel(NDOFEL),Accn(NDOFEL),TIME(2),PARAMS(*),
     3 JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     4 PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)

      integer lenJobName,lenOutDir,nDim,nInt,nIntS
      character*256 jobName,outDir,fileName
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      parameter(nInt=1)  ! number of volume integration pionts
      parameter(nIntS=1) ! number of surface integration points
      ! ********************************************************************************
      !*********************************************************************************
      ! Perform initial checks
      ! Open the debug/error message file
      !
      call getJobName(jobName,lenJobName)
      call getOutDir(outDir,lenOutDir)
      fileName = outDir(1:lenOutDir)//'\aaMSGS_'//
     +     jobName(1:lenJobName)//'.dat'
      open(unit=80,file=fileName,status='unknown')
      ! Check the procedure type, this should be a coupled
      !  temperature displacement or pore pressure displacement
      !  which are any of the following (64,65,72,73)
      !
      if((lflags(1).eq.64).or.(lflags(1).eq.65).or.
     +     (lflags(1).eq.72).or.(lflags(1).eq.73)) then
         !
         ! all is good
         !
      else
         write(*,*) 'Abaqus does not have the right procedure'
         write(*,*) 'go back and chekc the procedure type'
         write(*,*) 'lflags(1)=',lflags(1)
         write(80,*) 'Abaqus does not have the right procedure'
         write(80,*) 'go back and chekc the procedure type'
         write(80,*) 'lflags(1)=',lflags(1)
         call xit
      endif
      ! Make sure Abaqus knows you are doing a large
      !  deformation problem, I think this only matters
      !  when it comes to output in viewer
      !
      if(lflags(2).eq.0) then
         !
         ! lflags(2)=0 -> small disp.
         ! lflags(2)=1 -> large disp.
         !
         write(*,*) 'Abaqus thinks you are doing'
         write(*,*) 'a small displacement analysis'
         write(*,*) 'go in and set nlgeom=yes'
         write(80,*) 'Abaqus thinks you are doing'
         write(80,*) 'a small displacement analysis'
         write(80,*) 'go in and set nlgeom=yes'
         call xit
      endif
      ! Check to see if you are doing a general
      !  step or a linear purturbation step
      !
      if(lflags(4).eq.1) then
         !
         ! lflags(4)=0 -> general step
         ! lflags(4)=1 -> linear purturbation step
         !
         write(*,*) 'Abaqus thinks you are doing'
         write(*,*) 'a linear purturbation step'
         write(80,*) 'Abaqus thinks you are doing'
         write(80,*) 'a linear purturbation step'
         call xit         
      endif
      ! Do nothing if a ``dummy'' step
      !
      if(dtime.eq.0.0) return
      !
      ! Done with initial checks
!************************************************************************************************
      ! Call the paricular element to perform the analysis
      !
      if(jtype.eq.1) then
         !
         ! This is a plane strain analysis
         !
         nDim = 2
         call UPE4(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     +        PROPS,NPROPS,coords,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,
     +        TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,
     +        PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,
     +        NJPROP,PERIOD,
     +        nDim,nInt,nIntS)
         !
         !
      elseif(jtype.eq.2) then
         !
         ! This is an axisymmetric analysis
         !
         nDim = 2
         call UAX4(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     +        PROPS,NPROPS,coords,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,
     +        TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,
     +        PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,
     +        NJPROP,PERIOD,
     +        nDim,nInt,nIntS)
         !
      else
         !
         ! We have a problem...
         !
         write(*,*) 'Element type not supported, jtype=',jtype
         write(80,*) 'Element type not supported, jtype=',jtype
         call xit
         !
      endif
      !
      ! Done with this element, RHS and AMATRX already returned
      !  as output from the specific element routine called
      return
      end subroutine uel

************************************************************************
!----------------------plain strain element-----------------------------
************************************************************************         

      subroutine UPE4(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     +     PROPS,NPROPS,coords,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,
     +     TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,
     +     PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,
     +     NJPROP,PERIOD,
     +     nDim,nInt,nIntS)

      use global
*
      IMPLICIT NONE
*
*     VARIABLES DEFINED IN UEL, PASSED BACK TO ABAQUS
*
      REAL(8) :: RHS,AMATRX,SVARS,ENERGY
*
*     VARIABLES PASSED INTO UEL 
*
      REAL(8) :: PROPS,coords,Uall,DUall,Vel,Accn,TIME,
     1  DTIME,PARAMS,ADLMAG,PREDEF,DDLMAG,PNEWDT,PERIOD
      INTEGER :: NDOFEL,NRHS,NSVARS,NPROPS,MCRD,NNODE,JTYPE,KSTEP,KINC,
     1  JELEM,NDLOAD,JDLTYP,NPREDF,LFLAGS,MLVARX,MDLOAD,JPROPS,NJPROP
*
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     1 SVARS(*),ENERGY(8),coords(MCRD,NNODE),Uall(NDOFEL),
     2 DUall(MLVARX,*),Vel(NDOFEL),Accn(NDOFEL),TIME(2),PARAMS(*),
     3 JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     4 PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)

      real*8 u(nNode,2),du(nNode,ndofel),thetaNew(nNode)
      real*8 thetaOld(nNode),dtheta(nNode),CR(nNode)
      real*8 CROld(nNode),dCR(nNode),uNew(nNode,ndofel)
      real*8 uOld(nNode,ndofel),u_t(nNode,ndofel),v(nNode,2)
      real*8 coordsC(mcrd,nNode),Le,Gshear

      integer i,j,k,l,m,n,nIntPt,nDim,intpt,pOrder,a1,b1,a11,b11,face
      integer nInt,ii,jj,a12,b12,pe,nSdv,stat,q,nIntV,nIntPtV,p,ngSdv
      integer nlSdv,kk,lenJobName,lenOutDir,faceFlag,nIntS

      real*8 Iden(3,3),Ru(2*nNode,1),Rc(nNode,1),body(3)
      real*8 Kuu(2*nNode,2*nNode),Kcc(nNode,nNode),sh0(nNode),detMapJ0
      real*8 dshxi(nNode,2),dsh0(nNode,2),dshC0(nNode,2),detMapJ0C
      real*8 Fc_tau(3,3),Fc_t(3,3),detFc_tau,detFc_t,w(nInt)
      real*8 Psci_t,Pref_t,CR_tau,CR_t,CRn,CRno,dCRdt,dCRdX(2,1)
      real*8 sh(nNode),detMapJ,dsh(nNode,2),detMapJC
      real*8 dshC(nNode,2),mu_tau,F_tau(3,3)
      real*8 F_t(3,3),detF_tau,xi(nInt,2),detF,TR_tau(3,3),T_tau(3,3)
      real*8 SVM,SpTanMod(3,3,3,3),Rsci,Rref,dCdotdC,Dox
      real*8 dDoxdC,dRscidC,dRrefdC,Psci_tau,Pref_tau
      real*8 Smat(3,1),Bmat(3,2*nNode),BodyForceRes(2*nNode,1),Qmat(4,4)
      real*8 Gmat(4,2*nNode),G0mat(4,2*nNode),Amat(4,4),wS(nIntS)
      real*8 xLocal(nIntS),yLocal(nIntS),Kcu(nNode,2*nNode),detF_t
      real*8 Kuc(2*nNode,nNode),Nvec(1,nNode),AmatUC(3,1),TanFac
      real*8 SpUCMod(3,3),SpCUMod(3,3,3),SpCUModFac(3,3),AmatCU(2,4)

      real*8 zero,one,two,half,Pi,three,third
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,Pi=3.141592653d0,
     +     three=3.d0,third=1.d0/3.d0)


      ! Get element parameters
      nlSdv = jprops(1) ! number of local sdv's per integ point
      ngSdv = jprops(2) ! number of global sdv's per integ point
      ! Allocate memory for the globalSdv's used for viewing  results on the dummy mesh
      !
      if(.not.allocated(globalSdv)) then
         !
         ! allocate memory for the globalSdv's
         ! numElem needs to be set in the MODULE nInt needs to be set in the UEL
         !
         stat=0
         allocate(globalSdv(numElem,nInt,ngSdv),stat=err)
         if(stat.ne.0) then
            write(*,*) '//////////////////////////////////////////////'
            write(*,*) 'error when allocating globalSdv'
            write(*,*) '//////////////////////////////////////////////'
            write(*,*) '   stat=',stat
            write(*,*) '  ngSdv=',ngSdv
            write(*,*) '   nInt=',nInt
            write(*,*) 'numElem=',numElem
            write(*,*) '  nNode=',nNode
            write(*,*) 'lbound(globalSdv)',lbound(globalSdv)
            write(*,*) 'ubound(globalSdv)',ubound(globalSdv)
            write(*,*) '//////////////////////////////////////////////'
            call xit
         endif
         write(*,*) '-------------------------------------------------'
         write(*,*) '----------- globalSDV ALLOCATED -----------------'
         write(*,*) '-------------------------------------------------'
         write(*,*) '---------- YOU PUT NUMBER OF ELEMENTS -----------'
         write(*,*) '---------- numElem=',numElem
         write(*,*) '---------- UPE4 ELEMENTS ------------------------'
         write(*,*) '-------------------------------------------------'
         write(*,*) '---------- YOU PUT NUMBER OF POINTS -------------'
         write(*,*) '---------- nInt =',nInt
         write(*,*) '---------- nIntS=',nIntS
         write(*,*) '-------------------------------------------------'
         write(*,*) '---------- YOU PUT NUMBER OF SDVs ---------------'
         write(*,*) '---------- ngSdv=',ngSdv,nlSdv
         write(*,*) '-------------------------------------------------'
      endif
      ! Identity tensor
      call onem(Iden)
      ! Initialize the residual and tangent matrices to zero.
      Ru  = zero
      Rc = zero
      Kuu = zero
      Kcc = zero
      Kuc = zero
      Kcu = zero
      Energy = zero
      ! Body forces
      body(1:3) = zero
      ! Obtain nodal displacements and oxygen concentration
      !
      k = 0
      do i=1,nNode
         do j=1,nDim
            k = k + 1
            u(i,j) = Uall(k)
            du(i,j) = DUall(k,1)
            uOld(i,j) = u(i,j) - du(i,j)
         enddo
         k = k + 1
         CR(i) = Uall(k)
         dCR(i) = DUall(k,1)
         CROld(i) = CR(i) - dCR(i)
      enddo
      ! Obtain current nodal coordinates
      !
      do i=1,nNode
         do j=1,nDim
            coordsC(j,i) = coords(j,i) + u(i,j)
         enddo
      enddo 
      ! Impose any time-stepping changes on the increments of oxygen concentration or displacement if you wish
      !
      ! oxygen concentration increment
      !
      do i=1,nNode
         if(dabs(dCR(i)).gt.1.0E-1) then
            pnewdt = 0.5
            return
         endif
      enddo
      !
      ! displacement increment, based on element diagonal
      !
      Le = dsqrt(((coordsC(1,1)-coordsC(1,3))**two) + 
     +     ((coordsC(2,1)-coordsC(2,3))**two))
      !
      do i=1,nNode
         do j=1,nDim
            if(dabs(du(i,j)).gt.10.0*Le) then
               pnewdt = 0.5
               return
            endif
         enddo
      enddo
!print results for the first increment in log file
      if(JELEM.eq.1) then
      write(*,*) '#####################################################'
      write(*,*) '*********************start a new iteration'
      write(*,*) '#####################################################'
      write(*,*) '//////////////////////////////////////////////'
      write(*,*) '*********** KINC,dtime',KINC,dtime
       write(*,*) '*********** jtype,Time',jtype,time(1)
      endif
      if((kinc.le.4).and.(kinc.ge.0)) then
      write(*,*) '*********************start a new element'
      write(*,*) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
      write(*,*) '---------- KINC,JELEM,time(1)',KINC,JELEM,time(1)
      write(*,*) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
      write(*,*) '---------- ui',u(1,1),u(2,1),u(3,1),u(4,1)
      write(*,*) '---------- uj',u(1,2),u(2,2),u(3,2),u(4,2)      
      write(*,*) '---------- CR',CR
      write(*,*) '---------- coordsC',coordsC
      endif

      ! Get the deformation gradient for use in the 'F-bar' method.
      ! Obtain shape functions and their local gradients at the element centriod, that means xi=eta=zeta=0.0, and nIntPt=1
      !
      if(nNode.eq.4) then
         call calcShape2DLinear(1,zero,1,sh0,dshxi)
      else
         write(*,*) 'Incorrect number of nodes: nNode.ne.4'
         write(80,*) 'Incorrect number of nodes: nNode.ne.4'
         call xit
      endif
      ! Map shape functions from local to global reference coordinate system
      !
      if(mcrd.eq.2) then
         call mapShape2Da(nNode,dshxi,coords,dsh0,detMapJ0,stat)
         if(stat.eq.0) then
            pnewdt = 0.5
            return
         endif
      elseif(mcrd.eq.3) then
         call mapShape2D(nNode,dshxi,coords,dsh0,detMapJ0,stat)
         if(stat.eq.0) then
            pnewdt = 0.5
            return
         endif
      else
         ! big problem
         write(*,*) 'Unexpected error, mcrd=',mcrd
         write(80,*) 'Unexpected error, mcrd=',mcrd
         call xit
      endif
      ! Map shape functions from local to global current coordinate system
      !
      if(mcrd.eq.2) then
         call mapShape2Da(nNode,dshxi,coordsC,dshC0,detMapJ0C,stat)
         if(stat.eq.0) then
            pnewdt = 0.5
            return
         endif
      elseif(mcrd.eq.3) then
         call mapShape2D(nNode,dshxi,coordsC,dshC0,detMapJ0C,stat)
         if(stat.eq.0) then
            pnewdt = 0.5
            return
         endif
      else
         ! big problem
         write(*,*) 'Unexpected error, mcrd=',mcrd
         write(80,*) 'Unexpected error, mcrd=',mcrd
         call xit
      endif
      ! Calculate the deformation gradient at the element centriod at the the begining and end of the increment for use in 
      !  the `F-bar' method. `Tau' represents the end of the increment and `t' the previous increment.
      !
      Fc_tau = Iden
      Fc_t = Iden
      do i=1,nDim
         do j=1,nDim
            do k=1,nNode
               Fc_tau(i,j) = Fc_tau(i,j) + dsh0(k,j)*u(k,i)
               Fc_t(i,j) = Fc_t(i,j) + dsh0(k,j)*uOld(k,i)
            enddo
         enddo
      enddo
      !
      ! modify for plane-strain
      !
      Fc_tau(3,3) = one
      Fc_t(3,3) = one
      !
      ! 2D plane-strain implementation detF
      !
      detFc_t = Fc_t(1,1)*Fc_t(2,2) - Fc_t(1,2)*Fc_t(2,1)
      detFc_tau = Fc_tau(1,1)*Fc_tau(2,2) - Fc_tau(1,2)*Fc_tau(2,1)
      ! With the deformation gradient known at the element centriod we are now able to implement the `F-bar' method later
      !
      !***********************************************************************************************************
      ! Begin the loop over body integration points
      !***********************************************************************************************************
      ! Obtain integration point local coordinates and weights
      !
      if(nNode.eq.4) then
         ! gauss integration for a rectangular element
         if(nInt.eq.4) then
            call xint2D4pt(xi,w,nIntPt) ! 4-pt integration, nInt=4 above
         elseif(nInt.eq.1) then
            call xint2D1pt(xi,w,nIntPt) ! 1-pt integration, nInt=1 above
         else
            write(*,*) 'Invalid number of int points, nInt=',nInt
            write(80,*) 'Invalid number of int points, nInt=',nInt
            call xit
         endif
      else
         write(*,*) 'Incorrect number of nodes: nNode.ne.4'
         write(80,*) 'Incorrect number of nodes: nNode.ne.4'
         call xit
      endif
      
      !******************************************************************************************************
      ! Loop over integration points(calculating the rhs and AMATRX numerically)
      jj = 0 ! jj is used for tracking the state variables
      do intpt=1,nIntPt
         ! Obtain state variables from previous increment
         !
         if((kinc.le.1).and.(kstep.eq.1)) then
            ! this is the first increment, of the first step give initial conditions
             Psci_t=zero
             Pref_t=zero
            else
            ! this is not the first increment, read old values
             Psci_t=svars(1+jj)
             Pref_t  =svars(2+jj)
         endif
 
         ! Obtain shape functions and their local gradients
         !
         if(nNode.eq.4) then
            call calcShape2DLinear(nIntPt,xi,intpt,sh,dshxi)
         else
            write(*,*) 'Incorrect number of nodes: nNode.ne.4'
            write(80,*) 'Incorrect number of nodes: nNode.ne.4'
            call xit
         endif
         

         ! Map shape functions from local to global reference coordinate system
         ! dN/dxi=(jacobian)*dN/etai...
         if(mcrd.eq.2) then
            call mapShape2Da(nNode,dshxi,coords,dsh,detMapJ,stat)
            if(stat.eq.0) then
               pnewdt = 0.5
               return
            endif
         elseif(mcrd.eq.3) then
            call mapShape2D(nNode,dshxi,coords,dsh,detMapJ,stat)
            if(stat.eq.0) then
               pnewdt = 0.5
               return
            endif
         else
            ! big problem
            write(*,*) 'Unexpected error, mcrd=',mcrd
            write(80,*) 'Unexpected error, mcrd=',mcrd
            call xit
         endif
         ! Map shape functions from local to global current coordinate system
         !
         if(mcrd.eq.2) then
            call mapShape2Da(nNode,dshxi,coordsC,dshC,detMapJC,stat)
            if(stat.eq.0) then
               pnewdt = 0.5
               return
            endif
         elseif(mcrd.eq.3) then
            call mapShape2D(nNode,dshxi,coordsC,dshC,detMapJC,stat)
            if(stat.eq.0) then
               pnewdt = 0.5
               return
            endif
         else
            ! big problem
            write(*,*) 'Unexpected error, mcrd=',mcrd
            write(80,*) 'Unexpected error, mcrd=',mcrd
            call xit
      endif
      !!write the result into .log file
      ! if((kinc.lt.1).and.(JELEM.eq.1)) then
      !  write(*,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
      !  write(*,*) 'integration point'
      !  write(*,*) '---------- Iden',Iden
      !  write(*,*) '---------- sh0,dshxi',sh0,dshxi
      !  write(*,*) '---------- xi,w,nIntPt',xi,w,nIntPt
      !  write(*,*) '---------- sh,dshC',sh,dshC
      !endif
         ! Obtain the oxygen concentration and its derivative's at this intPt at the begining and end of the incrment
         CR_tau = zero
         CR_t = zero
         dCRdt = zero
         dCRdX = zero
      !
      !imposing constraint on oxygen concentration(the o2 concentration should be less than environment o2 concentration)
      ! ambient o2 concentration is considered 2.8 mol/m^3
         do k=1,nNode
             CRn=CR(k)
             CRno=CROld(k)
            CR_tau = CR_tau + CRn*sh(k)
            CR_t   = CR_t + CRno*sh(k)
            do i=1,nDim
               dCRdX(i,1) = dCRdX(i,1) + CRn*dshC(k,i)
            enddo
        enddo
     
          
         dCRdt = (CR_tau - CR_t)/dtime
 
         ! Obtain, and modify the deformation gradient at this integration point.  Modify the deformation gradienet for use in the `F-bar'
         !  method.  Also, take care of plane-strain or axisymetric
         !
         F_tau = Iden
         F_t = Iden
         do i=1,nDim
            do j=1,nDim
               do k=1,nNode
                  F_tau(i,j) = F_tau(i,j) + dsh(k,j)*u(k,i)
                  F_t(i,j) = F_t(i,j) + dsh(k,j)*uOld(k,i)
               enddo
            enddo
         enddo
         !
         ! modify F(3,3) for plane-strain 
         !
         F_tau(3,3) = one
         F_t(3,3) = one
         !
         !
         if((nNode.eq.4).and.(nInt.eq.4)) then
            !
            !  2D plane-strain implementation
            !
            detF_t = F_t(1,1)*F_t(2,2) - F_t(1,2)*F_t(2,1)
            detF_tau = F_tau(1,1)*F_tau(2,2) - F_tau(1,2)*F_tau(2,1)
            do i=1,nDim
               do j=1,nDim
                  F_tau(i,j) =((detFc_tau/detF_tau)**half)*F_tau(i,j)
                  F_t(i,j) = ((detFc_t/detF_t)**half)*F_t(i,j)
               enddo
            enddo
         endif
         call mdet(F_tau,detF)
 !print results for the first increment in log file
      if((kinc.le.4).and.(kinc.ge.0)) then
      write(*,*) '*******************************************'
      write(*,*) '****************start loop for each intgpt'
      write(*,*) '---------- intpt',intpt
      write(*,*) '---------- dCRdt,dCRdX,',dCRdt,dCRdX
      endif
         !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
         !
         ! Perform the constitutive time integration at this integ. point 
         !
         call integ(props,nprops,dtime,kstep,kinc,JELEM,
     +     F_tau,F_t,CR_tau,CR_t,Psci_t,Pref_t,  
     +     T_tau,SVM,SpTanMod,Rsci,Rref,dCdotdC,Dox,
     +     dDoxdC,dRscidC,dRrefdC,Psci_tau,Pref_tau,Gshear,
     +     SpUCMod,SpCUModfac)
         !
         !****************************************************************
      !    if((kinc.le.4).and.(kstep.eq.1)) then
      !write(*,*) '*******************************************'
      !write(*,*) '****************After integ subroutine'
      !write(*,*) '----------Rsci,Rref,dCdotdC,Dox',Rsci,Rref,dCdotdC,Dox
      !write(*,*) '--------dDoxdC,dRscidC,dRrefdC',dDoxdC,dRscidC,dRrefdC
      !write(*,*) '---------- Psci_tau,Pref_tau,',Psci_tau,Pref_tau
      !endif
         !
         !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
         ! Save the state variables at this integ point at the end of the increment
         !
             svars(1+jj)= Psci_tau
             svars(2+jj)=Pref_tau  
             jj = jj + nlSdv ! setup for the next intPt
         ! Save the state variables at this integ point in the global array used for plotting field output
         !
         globalSdv(jelem,intPt,1) = Psci_tau   ! sci extent of oxidation
         globalSdv(jelem,intPt,2) = Pref_tau   ! ref extent of oxidation
         globalSdv(jelem,intPt,3) = SVM        ! Von Mises stress
         globalSdv(jelem,intPt,4) = Gshear     ! Shear modulus in intpt
C ***********************************************************************************************************
         ! Compute/update the displacement residual vector
         !
         Smat(1,1) = T_tau(1,1)
         Smat(2,1) = T_tau(2,2)
         Smat(3,1) = T_tau(1,2)
         !
         Bmat = zero
         do kk=1,nNode
            Bmat(1,1+nDim*(kk-1)) = dshC(kk,1)
            Bmat(2,2+nDim*(kk-1)) = dshC(kk,2)
            Bmat(3,1+nDim*(kk-1)) = dshC(kk,2)
            Bmat(3,2+nDim*(kk-1)) = dshC(kk,1)
         enddo
         !
         BodyForceRes = zero
         do kk=1,nNode
            BodyForceRes(1+nDim*(kk-1),1) = sh(kk)*body(1)
            BodyForceRes(2+nDim*(kk-1),1) = sh(kk)*body(2)
         enddo
         !*******************************************************************************************************
         !RHS 
         Ru = Ru + detmapJC*w(intpt)*
     +        (
     +        -matmul(transpose(Bmat),Smat)
     +        + BodyForceRes
     +        )
         ! Compute/update the oxygen concentration residual vector
         !
         do kk=1,nNode
            Nvec(1,kk) = sh(kk)
         enddo
         !
         Rc = Rc + detmapJC*w(intpt)*
     +        (
     +        -transpose(Nvec)*(dCRdt/detF) - Dox*matmul(dshC,dCRdX)+
     +        -transpose(Nvec)*(Rsci)+ transpose(Nvec)*(Rref) 
     +        )
         !************************************************************************************************************
         !AMATRX
         ! Compute/update the displacement tangent matrix
         !
         Gmat = zero
         do kk=1,nNode
            Gmat(1,1+nDim*(kk-1)) = dshC(kk,1)
            Gmat(2,2+nDim*(kk-1)) = dshC(kk,1)
            Gmat(3,1+nDim*(kk-1)) = dshC(kk,2)
            Gmat(4,2+nDim*(kk-1)) = dshC(kk,2)
         enddo

         G0mat = zero
         do kk=1,nNode
            G0mat(1,1+nDim*(kk-1)) = dshC0(kk,1)
            G0mat(2,2+nDim*(kk-1)) = dshC0(kk,1)
            G0mat(3,1+nDim*(kk-1)) = dshC0(kk,2)
            G0mat(4,2+nDim*(kk-1)) = dshC0(kk,2)
         enddo

         Amat = zero
         Amat(1,1) = SpTanMod(1,1,1,1)
         Amat(1,2) = SpTanMod(1,1,2,1)
         Amat(1,3) = SpTanMod(1,1,1,2)
         Amat(1,4) = SpTanMod(1,1,2,2)
         Amat(2,1) = SpTanMod(2,1,1,1)
         Amat(2,2) = SpTanMod(2,1,2,1)
         Amat(2,3) = SpTanMod(2,1,1,2)
         Amat(2,4) = SpTanMod(2,1,2,2)
         Amat(3,1) = SpTanMod(1,2,1,1)
         Amat(3,2) = SpTanMod(1,2,2,1)
         Amat(3,3) = SpTanMod(1,2,1,2)
         Amat(3,4) = SpTanMod(1,2,2,2)
         Amat(4,1) = SpTanMod(2,2,1,1)
         Amat(4,2) = SpTanMod(2,2,2,1)
         Amat(4,3) = SpTanMod(2,2,1,2)
         Amat(4,4) = SpTanMod(2,2,2,2)


         Qmat = zero
         Qmat(1,1) = half*(Amat(1,1)+Amat(1,4)) - half*T_tau(1,1)
         Qmat(2,1) = half*(Amat(2,1)+Amat(2,4)) - half*T_tau(1,2)
         Qmat(3,1) = half*(Amat(3,1)+Amat(3,4)) - half*T_tau(1,2)
         Qmat(4,1) = half*(Amat(4,1)+Amat(4,4)) - half*T_tau(2,2)
         Qmat(1,4) = half*(Amat(1,1)+Amat(1,4)) - half*T_tau(1,1)
         Qmat(2,4) = half*(Amat(2,1)+Amat(2,4)) - half*T_tau(1,2)
         Qmat(3,4) = half*(Amat(3,1)+Amat(3,4)) - half*T_tau(1,2)
         Qmat(4,4) = half*(Amat(4,1)+Amat(4,4)) - half*T_tau(2,2)
            

         if((nNode.eq.4).and.(nInt.eq.4)) then
            ! This is the tangent using the F-bar method with the  4 node fully integrated linear element
            Kuu = Kuu + detMapJC*w(intpt)*
     +           (
     +            matmul(matmul(transpose(Gmat),Amat),Gmat)
     +           + matmul(transpose(Gmat),matmul(Qmat,(G0mat-Gmat)))
     +           )
         else
            !
            ! This is the tangent not using the F-bar method with all other elements
            Kuu = Kuu + detMapJC*w(intpt)*
     +           (
     +           matmul(matmul(transpose(Gmat),Amat),Gmat)
     +           )
         endif
         !
         ! Compute/update the oxygen concentration tangent matrix
         Kcc = Kcc + detmapJC*w(intPt)*
     +        (
     +        +(dCdotdC/detF)*matmul(transpose(Nvec),Nvec)
     +        + Dox*matmul(dshC,transpose(dshC))
     +        + dDoxdC*matmul(matmul(dshC,dCRdX),Nvec)
     +        + (dRscidC)*matmul(transpose(Nvec),Nvec)
     +        - (dRrefdC)*matmul(transpose(Nvec),Nvec)
     +        )
         
         ! Compute/update the oxygen concentration - displacement tangent matrix
         !
         SpCUMod = zero
         do i=1,nDim
            do k=1,nDim
               do l=1,nDim
                  SpCUMod(i,k,l) = SpCUMod(i,k,l)
     +                 + dCRdX(k,1)*SpCUModFac(i,l)
               enddo
            enddo
         enddo
         !
         AmatCU = zero
         AmatCU(1,1) = SpCUMod(1,1,1)
         AmatCU(1,2) = SpCUMod(1,2,1)
         AmatCU(1,3) = SpCUMod(1,1,2)
         AmatCU(1,4) = SpCUMod(1,2,2)
         AmatCU(2,1) = SpCUMod(2,1,1)
         AmatCU(2,2) = SpCUMod(2,2,1)
         AmatCU(2,3) = SpCUMod(2,1,2)
         AmatCU(2,4) = SpCUMod(2,2,2)
         !
         Kcu = Kcu - detMapJC*w(intpt)*
     +        (
     +        matmul(matmul(dshC,AmatCU),Gmat)
     +        )


         ! Compute/update the displacement - oxygen concentration tangent matrix
         AmatUC = zero
         AmatUC(1,1) = SpUCMod(1,1)
         AmatUC(2,1) = SpUCMod(2,2)
         AmatUC(3,1) = SpUCMod(1,2)
         !
         Kuc = Kuc + detMapJC*w(intpt)*
     +        (
     +        matmul(matmul(transpose(Bmat),AmatUC),Nvec)
     +        )

      enddo
      
      !
      ! End the loop over body integration points
      !----------------------------------------------------------------
      ! Return Abaqus the RHS vector and the Stiffness matrix.
      !
      call AssembleElement(nDim,nNode,nDofEl,
     +     Ru,Rc,Kuu,Kuc,Kcu,Kcc,
     +     rhs,amatrx)
      !
      ! End return of RHS and AMATRX
      return
      end subroutine UPE4

************************************************************************
************************************************************************
!----------------------Axisymmetric element-----------------------------
************************************************************************ 
      subroutine UAX4(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     +     PROPS,NPROPS,coords,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,
     +     TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,
     +     PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,
     +     NJPROP,PERIOD,
     +     nDim,nInt,nIntS)

      use global

      IMPLICIT NONE
*
*     VARIABLES DEFINED IN UEL, PASSED BACK TO ABAQUS
*
      REAL(8) :: RHS,AMATRX,SVARS,ENERGY
*
*     VARIABLES PASSED INTO UEL 
*
      REAL(8) :: PROPS,coords,Uall,DUall,Vel,Accn,TIME,
     1  DTIME,PARAMS,ADLMAG,PREDEF,DDLMAG,PNEWDT,PERIOD
      INTEGER :: NDOFEL,NRHS,NSVARS,NPROPS,MCRD,NNODE,JTYPE,KSTEP,KINC,
     1  JELEM,NDLOAD,JDLTYP,NPREDF,LFLAGS,MLVARX,MDLOAD,JPROPS,NJPROP
*
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     1 SVARS(*),ENERGY(8),coords(MCRD,NNODE),Uall(NDOFEL),
     2 DUall(MLVARX,*),Vel(NDOFEL),Accn(NDOFEL),TIME(2),PARAMS(*),
     3 JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     4 PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)

      real*8 u(nNode,2),du(nNode,ndofel),thetaNew(nNode)
      real*8 thetaOld(nNode),dtheta(nNode),CR(nNode)
      real*8 CROld(nNode),dCR(nNode),uNew(nNode,ndofel)
      real*8 uOld(nNode,ndofel),u_t(nNode,ndofel),v(nNode,2)
      real*8 coordsC(mcrd,nNode),Gshear

      integer i,j,k,l,m,n,nIntPt,nDim,intpt,pOrder,a1,b1,a11,b11,face
      integer nInt,ii,jj,a12,b12,pe,nSdv,stat,q,nIntV,nIntPtV,p,ngSdv
      integer nlSdv,kk,faceFlag,nIntS

       real*8 Iden(3,3),CR0,Ru(2*nNode,1),Rc(nNode,1),body(3)
      real*8 Kuu(2*nNode,2*nNode),Kcc(nNode,nNode),sh0(nNode),detMapJ0
      real*8 dshxi(nNode,2),dsh0(nNode,2),dshC0(nNode,2),detMapJ0C,Vmol
      real*8 Fc_tau(3,3),Fc_t(3,3),detFc_tau,detFc_t,w(nInt),Psci0,Pref0
      real*8 Psci_t,Pref_t,F_p1(3,3),CR_tau,CR_t,dCRdt,dCRdX(2,1)
      real*8 xi1_tau,xi3_tau,Co2_tau,CPOOH_tau,CPH_tau,Fn2,CPOOH,CPH
      real*8 sh(nNode),detMapJ,phi_t,dsh(nNode,2),detMapJC,phiLmt,umeror
      real*8 dshC(nNode,2),mu_tau,mu_t,dMUdX(2,1),dMUdt,F_tau(3,3)
      real*8 F_t(3,3),detF_tau,xi(nInt,2),detF,TR_tau(3,3),T_tau(3,3)
      real*8 SVM,SpTanMod(3,3,3,3),Rsci,Rref,dCdotdC,Dox
      real*8 dDoxdC,dRscidC,dRrefdC,Psci_tau,Pref_tau,AR0,ARc,AR_t,AR
      real*8 SmatAX(3,1),BmatAX(3,2*nNode),BodyForceResAX(2*nNode,1)
      real*8 GmatAX(4,2*nNode),G0matAX(4,2*nNode),AmatAX(4,4)
      real*8 QmatAX(4,4),theta0,xi1,xi3
      real*8 xLocal(nIntS),yLocal(nIntS),Kcu(nNode,2*nNode),detF_t
      real*8 Kuc(2*nNode,nNode),Nvec(1,nNode),ResFac,AmatUC(3,1),TanFac
      real*8 SpUCMod(3,3),SpCUMod(3,3,3),SpCUModFac(3,3),AmatCU(2,4)
      

      real*8 zero,one,two,half,Pi,three,third
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,Pi=3.141592653d0,
     +     three=3.d0,third=1.d0/3.d0)


      ! Get element parameters
      nlSdv = jprops(1) ! number of local sdv's per integ point
      ngSdv = jprops(2) ! number of global sdv's per integ point
      ! Allocate memory for the globalSdv's used for viewing  results on the dummy mesh
      !
      if(.not.allocated(globalSdv)) then
         !
         ! allocate memory for the globalSdv's
         ! numElem needs to be set in the MODULE nInt needs to be set in the UEL
         !
         stat=0
         allocate(globalSdv(numElem,nInt,ngSdv),stat=err)
         if(stat.ne.0) then
            write(*,*) '//////////////////////////////////////////////'
            write(*,*) 'error when allocating globalSdv'
            write(*,*) '//////////////////////////////////////////////'
            write(*,*) '   stat=',stat
            write(*,*) '  ngSdv=',ngSdv
            write(*,*) '   nInt=',nInt
            write(*,*) 'numElem=',numElem
            write(*,*) '  nNode=',nNode
            write(*,*) 'lbound(globalSdv)',lbound(globalSdv)
            write(*,*) 'ubound(globalSdv)',ubound(globalSdv)
            write(*,*) '//////////////////////////////////////////////'
            call xit
         endif
         write(*,*) '-------------------------------------------------'
         write(*,*) '----------- globalSDV ALLOCATED -----------------'
         write(*,*) '-------------------------------------------------'
         write(*,*) '----------YOU PUT NUMBER OF ELEMENTS -----------'
         write(*,*) '---------- numElem=',numElem
         write(*,*) '---------- UAX4 ELEMENTS ------------------------'
         write(*,*) '-------------------------------------------------'
         write(*,*) '---------- YOU PUT NUMBER OF POINTS -------------'
         write(*,*) '---------- nInt= ',nInt
         write(*,*) '---------- nIntS=',nIntS
         write(*,*) '-------------------------------------------------'
         write(*,*) '---------- YOU PUT NUMBER OF SDVs ---------------'
         write(*,*) '---------- ngSdv=',ngSdv
         write(*,*) '-------------------------------------------------'
      endif
      ! Identity tensor
      call onem(Iden)
      ! Obtain initial conditions
      theta0 = props(9)
      CR0   = props(10)
      Psci0=zero
      Pref0=zero
      ! Initialize the residual and tangent matrices to zero.
      !
      Ru  = zero
      Rc = zero
      Kuu = zero
      Kcc = zero
      Kuc = zero
      Kcu = zero
      Energy = zero
      ! Body forces
      !
      body(1:3) = zero
      ! Obtain nodal displacements and oxygen concentration
      k = 0
      do i=1,nNode
         do j=1,nDim
            k = k + 1
            u(i,j) = Uall(k)
            du(i,j) = DUall(k,1)
            uOld(i,j) = u(i,j) - du(i,j)
         enddo
         k = k + 1
         CR(i) = Uall(k)
         dCR(i) = DUall(k,1)
         CROld(i) = CR(i) - dCR(i)
       enddo
      ! Obtain current nodal coordinates
      !
      do i=1,nNode
         do j=1,nDim
            coordsC(j,i) = coords(j,i) + u(i,j)
         enddo
      enddo

      !  Get the deformation gradient for use in the  `F-bar' method.
      ! Obtain shape functions and their local gradients at the element centriod, that means xi=eta=zeta=0.0, and nIntPt=1
      if(nNode.eq.4) then
         call calcShape2DLinear(1,zero,1,sh0,dshxi)
      else
         write(*,*) 'Incorrect number of nodes: nNode.ne.4'
         write(80,*) 'Incorrect number of nodes: nNode.ne.4'
         call xit
      endif
      ! Map shape functions from local to global reference coordinate system
      !
      if(mcrd.eq.2) then
         call mapShape2Da(nNode,dshxi,coords,dsh0,detMapJ0,stat)
         if(stat.eq.0) then
            pnewdt = 0.5
            return
         endif
      elseif(mcrd.eq.3) then
         call mapShape2D(nNode,dshxi,coords,dsh0,detMapJ0,stat)
         if(stat.eq.0) then
            pnewdt = 0.5
            return
         endif
      else
         ! big problem
         write(*,*) 'Unexpected error, mcrd=',mcrd
         write(80,*) 'Unexpected error, mcrd=',mcrd
         call xit
      endif
      ! Map shape functions from local to global current coordinate system
      if(mcrd.eq.2) then
         call mapShape2Da(nNode,dshxi,coordsC,dshC0,detMapJ0C,stat)
         if(stat.eq.0) then
            pnewdt = 0.5
            return
         endif
      elseif(mcrd.eq.3) then
         call mapShape2D(nNode,dshxi,coordsC,dshC0,detMapJ0C,stat)
         if(stat.eq.0) then
            pnewdt = 0.5
            return
         endif
      else
         ! big problem
         write(*,*) 'Unexpected error, mcrd=',mcrd
         write(80,*) 'Unexpected error, mcrd=',mcrd
         call xit
      endif
      ! For an axisymmetric problem, find the ``r'' that  shows up in the integrals for axisymmetric
      !  and in the F(3,3), the factors of 2Pi are for integrals i.e., dV = 2 pi r dr dz
      AR0  = zero
      ARc  = zero
      AR_t = zero
      do i=1,nNode
         ! radial coord in ref config at centroid
         AR0  = AR0 + sh0(i)*coords(1,i)
         ! radial coord in current config at centroid
         ARc  = ARc + sh0(i)*(coords(1,i) + u(i,1))
         ! radial coord in current config at centroid in previous step
         AR_t = AR_t + sh0(i)*(coords(1,i) + uOld(i,1))
      enddo
      ! Calculate the deformation gradient at the element centriod at the the begining and end of the increment for use in 
      !  the `F-bar' method. `Tau' represents the end of the increment and `t' the previous increment.
      Fc_tau = Iden
      Fc_t = Iden
      do i=1,nDim
         do j=1,nDim
            do k=1,nNode
               Fc_tau(i,j) = Fc_tau(i,j) + dsh0(k,j)*u(k,i)
               Fc_t(i,j) = Fc_t(i,j) + dsh0(k,j)*uOld(k,i)
            enddo
         enddo
      enddo
      !
      ! modify for axisymmetric
      !
      Fc_tau(3,3) = ARc/AR0
      Fc_t(3,3) = AR_t/AR0
      ! axisymmetric implementation detF
      call mdet(Fc_tau,detFc_tau)
      call mdet(Fc_t,detFc_t)
      !
      ! With the deformation gradient known at the element centriod we are now able to implement the `F-bar' method later
      
      !***********************************************************************************************************
      ! Begin the loop over body integration points
      !***********************************************************************************************************
      ! Obtain integration point local coordinates and weights
      !
      if(nNode.eq.4) then
         !
         ! gauss integration for a rectangular element
         if(nInt.eq.4) then
            call xint2D4pt(xi,w,nIntPt) ! 4-pt integration, nInt=4 above
         elseif(nInt.eq.1) then
            call xint2D1pt(xi,w,nIntPt) ! 1-pt integration, nInt=1 above
         else
            write(*,*) 'Invalid number of int points, nInt=',nInt
            write(80,*) 'Invalid number of int points, nInt=',nInt
            call xit
         endif
      else
         write(*,*) 'Incorrect number of nodes: nNode.ne.4'
         write(80,*) 'Incorrect number of nodes: nNode.ne.4'
         call xit
      endif
      ! Loop over integration points
      !
      jj = 0 ! jj is used for tracking the state variables
      do intpt=1,nIntPt
         ! Obtain state variables from previous increment
         !
         if((kinc.le.1).and.(kstep.eq.1)) then
            ! this is the first increment, of the first step give initial conditions
            Psci_t  = Psci0
            Pref_t  = Pref0
         else
            ! this is not the first increment, read old values
            Psci_t = svars(1+jj)
            Pref_t =svars(2+jj)
         endif
         ! Obtain shape functions and their local gradients
         !
         if(nNode.eq.4) then
            call calcShape2DLinear(nIntPt,xi,intpt,sh,dshxi)
         else
            write(*,*) 'Incorrect number of nodes: nNode.ne.4'
            write(80,*) 'Incorrect number of nodes: nNode.ne.4'
            call xit
         endif
         ! Map shape functions from local to global reference coordinate system
         !
         if(mcrd.eq.2) then
            call mapShape2Da(nNode,dshxi,coords,dsh,detMapJ,stat)
            if(stat.eq.0) then
               pnewdt = 0.5
               return
            endif
         elseif(mcrd.eq.3) then
            call mapShape2D(nNode,dshxi,coords,dsh,detMapJ,stat)
            if(stat.eq.0) then
               pnewdt = 0.5
               return
            endif
         else
            ! big problem
            write(*,*) 'Unexpected error, mcrd=',mcrd
            write(80,*) 'Unexpected error, mcrd=',mcrd
            call xit
         endif


         ! Map shape functions from local to global current coordinate system
         !
         if(mcrd.eq.2) then
            call mapShape2Da(nNode,dshxi,coordsC,dshC,detMapJC,stat)
            if(stat.eq.0) then
               pnewdt = 0.5
               return
            endif
         elseif(mcrd.eq.3) then
            call mapShape2D(nNode,dshxi,coordsC,dshC,detMapJC,stat)
            if(stat.eq.0) then
               pnewdt = 0.5
               return
            endif
         else
            ! big problem
            write(*,*) 'Unexpected error, mcrd=',mcrd
            write(80,*) 'Unexpected error, mcrd=',mcrd
            call xit
         endif
         ! For an axisymmetric problem, find the ``r'' that  shows up in the integrals for axisymmetric
         !  and in the F(3,3), the factors of 2Pi are for integrals  i.e., dV = 2 pi r dr dz
         AR0  = zero
         AR   = zero
         AR_t = zero
         do i=1,nNode
            AR0 = AR0 + sh(i)*coords(1,i)
            AR  = AR  + sh(i)*(coords(1,i) + u(i,1))
            AR_t = AR_t + sh(i)*(coords(1,i) + uOld(i,1))
         enddo
         AR0  = two*Pi*AR0
         AR   = two*Pi*AR
         AR_t = two*Pi*AR_t
        
         ! Obtain the oxygen concentration and its derivative's at this intPt at the begining and end of the incrment
         CR_tau = zero
         CR_t = zero
         dCRdt = zero
         dCRdX = zero
         do k=1,nNode
            CR_tau = CR_tau + CR(k)*sh(k)
            CR_t   = CR_t + CROld(k)*sh(k)
            do i=1,nDim
               dCRdX(i,1) = dCRdX(i,1) + CR(k)*dshC(k,i)
            enddo
         enddo
         dCRdt = (CR_tau - CR_t)/dtime

         ! Obtain, and modify the deformation gradient at this integration point.  Modify the deformation gradienet for use in the `F-bar'
         !  method.  Also, take care of plane-strain or axisymetric
         !
         F_tau = Iden
         F_t = Iden
         do i=1,nDim
            do j=1,nDim
               do k=1,nNode
                  F_tau(i,j) = F_tau(i,j) + dsh(k,j)*u(k,i)
                  F_t(i,j) = F_t(i,j) + dsh(k,j)*uOld(k,i)
               enddo
            enddo
         enddo
         !
         ! modify F(3,3) for axisymetric, give R/R0
         !
         F_tau(3,3) = AR/AR0
         F_t(3,3) = AR_t/AR0
         !
         ! Modify the deformation gradient for the `F-bar' method only when using the 4 node fully integrated linear
         !  element, do not use the `F-bar' method for any other element
         !
         if((nNode.eq.4).and.(nInt.eq.4)) then
            call mdet(F_tau,detF_tau)
            call mdet(F_t,detF_t)
            F_tau = ((detFc_tau/detF_tau)**third)*F_tau
            F_t = ((detFc_tau/detF_tau)**third)*F_t
         endif
         call mdet(F_tau,detF)


         !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
         !
         ! Perform the constitutive time integration at this integ. point 
         !
        call integ(props,nprops,dtime,kstep,kinc,JELEM,
     +     F_tau,F_t,CR_tau,CR_t,Psci_t,Pref_t,  
     +     T_tau,SVM,SpTanMod,Rsci,Rref,dCdotdC,Dox,
     +     dDoxdC,dRscidC,dRrefdC,Psci_tau,Pref_tau,Gshear,
     +     SpUCMod,SpCUModfac)
         !
         !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

         ! Save the state variables at this integ point at the end of the increment
        ! Save the state variables at this integ point at the end of the increment
         !
         svars(1+jj) = Psci_tau
         svars(2+jj) = Pref_tau
         jj = jj + nlSdv ! setup for the next intPt
         ! Save the state variables at this integ point in the global array used for plotting field output
         !
         globalSdv(jelem,intPt,1) = Psci_tau   ! sci extent of oxidation
         globalSdv(jelem,intPt,2) = Pref_tau   ! ref extent of oxidation
         !
         ! Compute/update the displacement residual vector
         !
         SmatAx(1,1) = T_tau(1,1)
         SmatAx(2,1) = T_tau(2,2)
         SmatAx(3,1) = T_tau(1,2)
         SmatAx(4,1) = T_tau(3,3)
         !
         BmatAx = zero
         do kk=1,nNode
            BmatAx(1,1+nDim*(kk-1)) = dshC(kk,1)
            BmatAx(2,2+nDim*(kk-1)) = dshC(kk,2)
            BmatAx(3,1+nDim*(kk-1)) = dshC(kk,2)
            BmatAx(3,2+nDim*(kk-1)) = dshC(kk,1)
            BmatAx(4,1+nDim*(kk-1)) = sh(kk)/(AR/(two*Pi))
         enddo
         !
         BodyForceResAX = zero
         do kk=1,nNode
            BodyForceResAx(1+nDim*(kk-1),1) = sh(kk)*body(1)
            BodyForceResAx(2+nDim*(kk-1),1) = sh(kk)*body(2)
      enddo
         !
      !*******************************************************************************************************
         !RHS 
         Ru = Ru + detmapJC*w(intpt)*AR*
     +        (
     +        -matmul(transpose(BmatAx),SmatAx)
     +        + BodyForceResAx
     +        )
         ! Compute/update the chemical potential residual vector
         !
         do kk=1,nNode
            Nvec(1,kk) = sh(kk)
         enddo
         Rc = Rc + detmapJC*w(intpt)*AR*
     +        (
     +        transpose(Nvec)*dCRdt + Dox*matmul(dshC,dCRdX)+
     +        transpose(Nvec)*Rsci- transpose(Nvec)*Rref 
     +        )
         !************************************************************************************************************
         !AMATRX
         ! Compute/update the displacement tangent matrix
         !
         GmatAx = zero
         do kk=1,nNode
            GmatAx(1,1+nDim*(kk-1)) = dshC(kk,1)
            GmatAx(2,2+nDim*(kk-1)) = dshC(kk,1)
            GmatAx(3,1+nDim*(kk-1)) = dshC(kk,2)
            GmatAx(4,2+nDim*(kk-1)) = dshC(kk,2)
            GmatAx(5,1+nDim*(kk-1)) = sh(kk)/(AR/(two*Pi))
         enddo

         G0matAx = zero
         do kk=1,nNode
            G0matAx(1,1+nDim*(kk-1)) = dshC0(kk,1)
            G0matAx(2,2+nDim*(kk-1)) = dshC0(kk,1)
            G0matAx(3,1+nDim*(kk-1)) = dshC0(kk,2)
            G0matAx(4,2+nDim*(kk-1)) = dshC0(kk,2)
            G0matAX(5,1+nDim*(kk-1)) = sh0(kk)/ARc
         enddo

         AmatAx = zero
         AmatAx(1,1) = SpTanMod(1,1,1,1)
         AmatAx(1,2) = SpTanMod(1,1,2,1)
         AmatAx(1,3) = SpTanMod(1,1,1,2)
         AmatAx(1,4) = SpTanMod(1,1,2,2)
         AmatAx(1,5) = SpTanMod(1,1,3,3)
         AmatAx(2,1) = SpTanMod(2,1,1,1)
         AmatAx(2,2) = SpTanMod(2,1,2,1)
         AmatAx(2,3) = SpTanMod(2,1,1,2)
         AmatAx(2,4) = SpTanMod(2,1,2,2)
         AmatAx(2,5) = SpTanMod(2,1,3,3)
         AmatAx(3,1) = SpTanMod(1,2,1,1)
         AmatAx(3,2) = SpTanMod(1,2,2,1)
         AmatAx(3,3) = SpTanMod(1,2,1,2)
         AmatAx(3,4) = SpTanMod(1,2,2,2)
         AmatAx(3,5) = SpTanMod(1,2,3,3)
         AmatAx(4,1) = SpTanMod(2,2,1,1)
         AmatAx(4,2) = SpTanMod(2,2,2,1)
         AmatAx(4,3) = SpTanMod(2,2,1,2)
         AmatAx(4,4) = SpTanMod(2,2,2,2)
         AmatAx(4,5) = SpTanMod(2,2,3,3)
         AmatAx(5,1) = SpTanMod(3,3,1,1)
         AmatAx(5,2) = SpTanMod(3,3,2,1)
         AmatAx(5,3) = SpTanMod(3,3,1,2)
         AmatAx(5,4) = SpTanMod(3,3,2,2)
         AmatAx(5,5) = SpTanMod(3,3,3,3)

         QmatAx = zero
         QmatAx(1,1) = third*(AmatAx(1,1)+AmatAx(1,4)+AmatAx(1,5)) 
     +        - (two/three)*T_tau(1,1)
         QmatAx(2,1) = third*(AmatAx(2,1)+AmatAx(2,4)+AmatAx(2,5))
     +        - (two/three)*T_tau(1,2)
         QmatAx(3,1) = third*(AmatAx(3,1)+AmatAx(3,4)+AmatAx(3,5))
     +        - (two/three)*T_tau(1,2)
         QmatAx(4,1) = third*(AmatAx(4,1)+AmatAx(4,4)+AmatAx(4,5))
     +        - (two/three)*T_tau(2,2)
         QmatAx(5,1) = third*(AmatAx(5,1)+AmatAx(5,4)+AmatAx(5,5))
     +        - (two/three)*T_tau(3,3)
         QmatAx(1,4) = QmatAx(1,1)
         QmatAx(2,4) = QmatAx(2,1)
         QmatAx(3,4) = QmatAx(3,1)
         QmatAx(4,4) = QmatAx(4,1)
         QmatAx(5,4) = QmatAx(5,1)
         QmatAx(1,5) = QmatAx(1,1)
         QmatAx(2,5) = QmatAx(2,1)
         QmatAx(3,5) = QmatAx(3,1)
         QmatAx(4,5) = QmatAx(4,1)
         QmatAx(5,5) = QmatAx(5,1)
            

         if((nNode.eq.4).and.(nInt.eq.4)) then
            !
            ! This is the tangent using the F-bar method with the 4 node fully integrated linear element
            !
            Kuu = Kuu + detMapJC*w(intpt)*AR*
     +           (
     +           matmul(matmul(transpose(GmatAx),AmatAx),GmatAx)
     +           + matmul(transpose(GmatAx),matmul(QmatAx,
     +           (G0matAx-GmatAx)))
     +           )
         else
            !
            ! This is the tangent NOT using the F-bar method with all other elements
            !
            Kuu = Kuu + detMapJC*w(intpt)*AR*
     +           (
     +           matmul(matmul(transpose(GmatAx),AmatAx),GmatAx)
     +           )
         endif
         ! Compute/update the oxygen cocentration tangent matrix
         !
         
         Kcc = Kcc + detmapJC*AR*w(intPt)*
     +       (-dCdotdC*matmul(transpose(Nvec),Nvec)
     +        - Dox*matmul(dshC,transpose(dshC))
     +        - dDoxdC*matmul(matmul(dshC,dCRdX),Nvec)
     +        - dRscidC*matmul(transpose(Nvec),Nvec)
     +        + dRrefdC*matmul(transpose(Nvec),Nvec)
     +        )
         ! Compute/update the oxygen concentration - displacement tangent matrix.
         !
         SpCUMod = zero
         do i=1,nDim
            do k=1,nDim
               do l=1,nDim
                  SpCUMod(i,k,l) = SpCUMod(i,k,l)
     +                 +  dCRdX(k,1)*SpCUModFac(i,l)
               enddo
            enddo
         enddo
         !
         AmatCU = zero
         AmatCU(1,1) = SpCUMod(1,1,1)
         AmatCU(1,2) = SpCUMod(1,2,1)
         AmatCU(1,3) = SpCUMod(1,1,2)
         AmatCU(1,4) = SpCUMod(1,2,2)
         AmatCU(1,5) = SpCUMod(1,3,3)
         AmatCU(2,1) = SpCUMod(2,1,1)
         AmatCU(2,2) = SpCUMod(2,2,1)
         AmatCU(2,3) = SpCUMod(2,1,2)
         AmatCU(2,4) = SpCUMod(2,2,2)
         AmatCU(2,5) = SpCUMod(2,3,3)
         !
         Kcu = Kcu - detMapJC*w(intpt)*AR*
     +        (
     +        matmul(matmul(dshC,AmatCU),GmatAX)
     +        )


         ! Compute/update the displacement - chemical potential tangent matrix
         !
         AmatUC = zero
         AmatUC(1,1) = SpUCMod(1,1)
         AmatUC(2,1) = SpUCMod(2,2)
         AmatUC(3,1) = SpUCMod(1,2)
         AmatUC(4,1) = SpUCMod(3,3)
         !
         Kuc = Kuc + detMapJC*w(intpt)*AR*
     +        (
     +        matmul(matmul(transpose(BmatAX),AmatUC),Nvec)
     +        )


      enddo
      !
      ! End the loop over body integration points
      !----------------------------------------------------------------
      ! Return Abaqus the RHS vector and the Stiffness matrix.
      !
      call AssembleElement(nDim,nNode,nDofEl,
     +     Ru,Rc,Kuu,Kuc,Kcu,Kcc,
     +     rhs,amatrx)
      !
      ! End return of RHS and AMATRX
      !----------------------------------------------------------------
      return
      end subroutine UAX4

************************************************************************
************************************************************************

      subroutine integ(props,nprops,dtime,kstep,kinc,JELEM,
     +     F_tau,F_t,CR_tau,CR_t,Psci_t,Pref_t,  
     +     T_tau,SVM,SpTanMod,Rsci_tau,Rref_tau,dCdotdC,Dox,
     +     dDoxdC,dRscidC,dRrefdC,Psci_tau,Pref_tau,Gshear,
     +     SpUCMod,SpCUModfac)
      ! This subroutine computes everything required for the time integration of the AMATRX AND RHS.
      !
      ! Inputs:
      !  1) material parameters, props(nprops)
      !  2) time increment, dtime
      !  3) deformation gradient, F_tau(3,3)
      !  4) deformation gradient at previous time step, F_t(3,3)
      !  5) OXYGEN CONCENTRATION, CR_TAU
      !  6) OXYGEN CONCENTRATION at previous time step, CR_t
      !  8) old sci extent of oxidation, Psci_t
      !  10) old ref extent of oxidation, Pref_t
      !
      ! Outputs:
      !  1) Cauchy stress, T_tau(3,3)
      !  2) spatial tangent modulus, SpTanMod(3,3,3,3)
      !  3) chain scision reaction rate, Rsci_tau
      !  4) network reformation reaction rate, Rref_tau
      !  5) derivative of the Cdot with C, dCdotdC
      !  6) diffusivity of the integration point, Dox
      !  7) derivative of the Dox with C, dDoxdC
      !  8) derivative of the Rsci with C, dRscidC
      !  9) derivative of the Rref with C, dRrefdC
      !  10)new sci extent of oxidation, Psci_tau
      !  11) new ref extent of oxidation, Pref_tau
      !  16) displacement - chemical potential modulus terms, SpUCMod
      !  17) chemical potential - displacement modulus terms, SpCUModfac

      implicit none

      integer i,j,k,l,m,n,I1,J1,nprops,stat,KSTEP,KINC,JELEM,nlSdv,ngSdv

      real*8 alpha,beta,Ksci,Kref,Esci,Eref,D0,r0sci,r0ref,lan
      real*8 Ed,Hsci,Href,Co1,co2,vsci,vref,mu0,theta,R
      real*8 Iden(3,3),props(nprops),Gshear,Gushear,Kbulk,gama,chi,Vmol
      real*8 detF,trF,F_t(3,3),dtime,JP,F_tau(3,3),Finv(3,3),FinvT(3,3)
      real*8 B_tau(3,3),trB_tau,C_tau(3,3),trC_tau,Fpinv(3,3),Fe(3,3)
      real*8 Feinv(3,3),detFe,FeinvT(3,3),Ce_tau(3,3),Be_tau(3,3)
      real*8 Ce_iso(3,3),Be_iso(3,3),trCe_iso,Lambda,T_tau(3,3),SVM
      real*8 TR_tau(3,3),dTRdF(3,3,3,3),spTanMod(3,3,3,3),Me(3,3),trMe 
      real*8 CR_tau,CR_t,phi_tau,dCdotdC,Psci_t,Pref_t,mu_tau
      real*8 Fnsci,Fnref,Rsci_tau,Rref_tau,Psci_tau,Pref_tau,dtdC
      real*8 dRscidC,dRrefdC,Dox,dDoxdC,dDoxdPsci,dDoxdPref,dPscidt,dCR
      real*8 dPrefdt,SpUCMod(3,3),SpCUModfac(3,3),dTdJP(3,3),dJPdPsci
       
      real*8 zero,one,two,three,third,half
      parameter(zero=0.d0,one=1.d0,two=2.d0,three=3.d0,third=1.d0/3.d0,
     +     half=1.d0/2.d0)
      !
      !
      !reaction rate functions parameter
      !
      parameter(alpha=4.45,beta=1.1,Ksci=2.84E-9,Kref=2.76E-9,
     + Esci=2.3E4,Eref=2.4E4,D0=1.0E-6,r0sci=20.0,r0ref=20.0,
     + lan=0.7,Ed=3.5E4,Hsci=339.0E4,Href=708.39E4,Co1=0.2,co2=0.4,
     + vsci=0.5,vref=0.5,theta=353.0,mu0=3.88E5)
      !
      ! Identity tensor
      !
      call onem(Iden)
      ! Obtain material properties
      !
       Gushear = props(1)  ! Shear modulus
       Kbulk  = props(2)  ! Bulk modulus
       gama   = props(3)  ! Shrinkage parameter
       chi    = props(4)  ! Chi parameter
       Vmol   = props(5)  ! Volume of a mole of fluid particles
       R      = props(6)  ! Universal gas constant
       !theta= props(8)  ! Absolute temperature
       !
      if((kinc.lt.1).and.(JELEM.eq.1)) then
      write(*,*) '*********************start in integration point'
      write(*,*) 'checking the value of material props and reaction c'
      write(*,*) '---------- G,Kbulk,gama,',Gushear,Kbulk,gama
      write(*,*) '---------- chi,Vmol,R,theta,mu0',chi,Vmol,R,theta,mu0
      write(*,*) '---------- alpha,beta,Ksci,Kref,',alpha,beta,Ksci,Kref
      write(*,*) '---------- Esci,Eref,D0,r0sci',Esci,Eref,D0,r0sci   
      write(*,*) '---------- r0ref,lan,Ed,Hsci,',r0ref,lan,Ed,Hsci
      write(*,*) '---------- Href,Co1,co2,vsr',Href,Co1,co2,vsci,vref
      endif
C **********************************************************************************************
      ! Compute the F/F-1/Fe/C/B
C **********************************************************************************************
      !
      ! Compute the inverse of F_tau, Fe, its determinant, and its transpose
      call matInv3D(F_tau,Finv,detF,stat)
      if(stat.eq.0) then
         write(*,*) 'Problem: detF.lt.zero'
         call xit
      endif
      FinvT=transpose(Finv)
      trF=F_tau(1,1) + F_tau(2,2) + F_tau(3,3)
      !
      !plastic deformation gradient determinant
      !
      JP=exp(three*gama*psci_t)
      !
      ! Compute the left Cauchy-Green tensor and its trace
      !
      B_tau = matmul(F_tau,transpose(F_tau))
      trB_tau = B_tau(1,1) + B_tau(2,2) + B_tau(3,3)
      ! Compute the right Cauchy-Green tensor and its trace
      !
      C_tau = matmul(transpose(F_tau),F_tau)
      trC_tau = C_tau(1,1) + C_tau(2,2) + C_tau(3,3)
      !
      !compute Fe/Ce
        Fpinv= JP**(-one/THREE)*Iden
        Fe=matmul(F_tau,Fpinv)
         call matInv3D(Fe,Feinv,detFe,stat)
          if(stat.eq.0) then
            write(*,*) 'Problem: detF.lt.zero'
          call xit
          endif
      ! !compute the Ce_dis for calculation of ractive forces
        Ce_tau=matmul(transpose(Fe),Fe)
        Be_tau=matmul(Fe,transpose(Fe))
        FeinvT=transpose(Feinv)
        Lambda=sqrt(trC_tau/three)
C ***********************************************************************************************
      ! Compute the Cauchy stress and spatial tangent stiffness
C ***********************************************************************************************
        !
        !The influence of aging on shear modulus
        Gshear=Gushear*(1-co1*Psci_t+co2*Pref_t)
        !Gshear=Gushear
       T_tau = (Gshear*(B_tau-Iden) + Kbulk*dlog(detFe)*Iden)/detF
       !
       ! compue the Von Mises stress from Cauchy stress second order tensor
       !
       SVM=SQRT((ONE/TWO)*(((T_tau(1,1)-T_tau(2,2))**two)+
     +  ((T_tau(2,2)-T_tau(3,3))**two)+((T_tau(3,3)-T_tau(1,1))**two)+
     +  two*three*(T_tau(1,1)**two+T_tau(2,2)**two+T_tau(3,3)**two)))
       ! Compute the 1st Piola stress
       !
      TR_tau = Gshear*(F_tau - FinvT) + Kbulk*dlog(detFe)*FinvT
      ! Compute dTRdF, the so-called material tangent modulus
      !
      dTRdF = zero
      do i=1,3
         do j = 1,3
            do k = 1,3
               do l = 1,3          
                 dTRdF(i,j,k,l) = dTRdF(i,j,k,l)
     +                 + Gshear*Iden(i,k)*Iden(j,l)
     +                 + Gshear*Finv(l,i)*Finv(j,k)
     +                 + Kbulk*Finv(j,i)*Finv(l,k)
     +                 - Kbulk*dlog(detFe)*Finv(l,i)*Finv(j,k)
               enddo
            enddo
         enddo
      enddo
      !
      ! Calculate the so-called spatial tangent modulus, based on the push forward of the material tangent modulus
      !
      SpTanMod = zero
      do i=1,3
         do j=1,3
            do k=1,3
               do l=1,3
                  do m=1,3
                     do n=1,3
                        SpTanMod(i,j,k,l) = SpTanMod(i,j,k,l) + 
     +                       (dTRdF(i,m,k,n)*F_tau(j,m)*F_tau(l,n))/detF
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      
	 !
       !compute the mandel stress
       
        Me=detF*matmul(matmul(transpose(Fe),T_tau),FeinvT)
        trMe=Me(1,1)+Me(2,2)+Me(3,3)
!print results for the first increment in log file
      if((kinc.le.4).and.(kinc.ge.0)) then
      write(*,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
      write(*,*) '******************start in integ subroutine' 
      write(*,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'      
      write(*,*) '---------- detF,detFe,trMe',detF,detFe,trMe
      write(*,*) '---------- F_tau,Jp',F_tau,Jp
      write(*,*) '---------- Fpinv',Fpinv
      write(*,*) '---------- T_tau',T_tau
      endif
**************************************************************************************************************
      ! compute the derivative of oxygen concentration rate with  respect to the oxygen concentration(dCdot/dC)
**************************************************************************************************************     
      ! Compute the perturbation on the oxygen concentration
      !
      dCdotdC = (one)/(dtime)
C ****************************************************************************************************
      ! Compute the reaction rate functions
C *****************************************************************************************************
      !compute chemical potential at the intpt
      if (CR_tau.lt.0) then
          CR_tau=1.0e-15
      endif
      !
      phi_tau =(one/(one+Vmol*CR_tau))
      !
      !We have used phi=0.999 rather than 1.0 to eliminate numerical difficulties
      !with the ln(1-phi)
      if(phi_tau.ge.1) then
          phi_tau=0.9999
      endif
      mu_tau=mu0+R*theta*(dlog(1-phi_tau )+phi_tau +chi*(phi_tau)**two)
      !
      !compute the reactive forces at current time step
     
      Fnsci=(co1*Gushear/two)*
     +        (three*(Lambda**two-one)-two*dlog(detF))+
     +           Hsci*(1-Psci_t)+gama*trMe+mu_tau
      Fnref=(-co2*Gushear/two)*
     +        (three*(Lambda**two-one)-two*dlog(detF))
     +            +Href*(1-Pref_t)+mu_tau
      !
      !Compute the reaction rate functions
      !
      Rsci_tau=(alpha*CR_tau/(1+Beta*CR_tau))*Ksci*exp(-Esci/(R*theta))*
     +           *(1-Psci_t)*Fnsci
      Rref_tau=(alpha*CR_tau/(1+Beta*CR_tau))*Kref*exp(-Eref/(R*theta))*
     +           *(1-Pref_t)*Fnref
     
!print results for the first increment in log file
      if((kinc.le.4).and.(kinc.ge.0)) then
      write(*,*) '**************************************** '   
      write(*,*) '----------reaction rate functions in intgpt '
      write(*,*) '---------- CR_tau ',CR_tau
      write(*,*) '---------- dCdotdC,',dCdotdC
      write(*,*) '---------- phi_tau,',phi_tau   
      write(*,*) '---------- mu_tau,',mu_tau
      write(*,*) '---------- Fnsci,Fnref',Fnsci,Fnref
      write(*,*) '---------- Rsci_tau,Rref_tau,',Rsci_tau,Rref_tau
      endif
 
C ***************************************************************************************************
      ! update the extent of oxidation
C ***************************************************************************************************
      !  new sci extent of oxidation, Psci_tau
      !
      !reaction 2
       Psci_tau = Psci_t+dtime*Vsci*Rsci_tau/R0sci
       if (Psci_tau.gt.one) then
           Psci_tau=one
      endif
      !
      !reaction6
       Pref_tau = Pref_t+dtime*Vref*Rref_tau/R0ref
        if (Pref_tau.gt.one) then 
            Pref_tau=one
      endif

**************************************************************************************************************
      ! compute the derivative of reaction rate functions with  respect to the oxygen concentration(dRscidC, dRrefdC)
**************************************************************************************************************         
      ! derivative of the Rsci with C, dRscidC
      !
      dRscidC=(alpha/(1+beta*CR_tau)**two)*Ksci*exp(-Esci/(R*theta))*
     +           *(1-Psci_t)*Fnsci
      ! derivative of the Rref with C, dRrefdC
      !
      dRrefdC=(alpha/(1+beta*CR_tau)**two)*Kref*exp(-Eref/(R*theta))*
     +           *(1-Pref_t)*Fnref
!print results for the first increment in log file
      if((kinc.le.4).and.(kinc.ge.0)) then
      write(*,*) '**************************************** '      
      write(*,*) '---------- deravatives'
      write(*,*) '---------- Psci_tau, Pref_tau',Psci_tau, Pref_tau
      write(*,*) '---------- dRscidC,dRrefdC,',dRscidC,dRrefdC
      endif
C ************************************************************************************************
      ! compute Dox / dDoxdC
C ************************************************************************************************
      ! diffusivity of the integration point, Dox
      !
       Dox=D0*exp(-Ed/(R*theta))*(1-lan*(Psci_tau+Pref_tau)/two)
      !
      ! derivative of the Dox with respect to C, dDoxdC=dDoxdPsci*dPscidt*dtdC+dDoxdPref*dPrefdt*dtdC
      !
       dDoxdPsci=D0*exp(-Ed/(R*theta))*(-lan/two)
       !
       dDoxdPref=D0*exp(-Ed/(R*theta))*(-lan/two)
       !
       dPscidt=Vsci*Rsci_tau/R0sci
       !
       dPrefdt=Vref*Rref_tau/R0ref
       !
       dCR=CR_tau-CR_t
       if (dCR.le.0) then
          dCR=1.0e-10
      endif
       dtdC=dtime/(dCR)
       dDoxdC=dDoxdPsci*dPscidt*dtdC+dDoxdPref*dPrefdt*dtdC
!print results for the first increment in log file
      if((kinc.le.4).and.(kinc.ge.0)) then
      write(*,*) '**************************************** ' 
      write(*,*) '----------Dox,dPscidt,dPrefdt',Dox,dPscidt,dPrefdt
      write(*,*) '----------dDoxdPref,dDoxdPsci',dDoxdPref,dDoxdPsci
      write(*,*) '----------dDoxdC,dCR',dDoxdC, dCR 
      endif
C ***************************************************************************************************
      ! Compute the CU/UC modulus
C ***************************************************************************************************
      ! Compute the displacement - oxygen concentration modulus (Kuc(dTdC=dTdJP*dJPdPsci*dPscidt*dtdC))
      ! 
      dTdJP= (Kbulk/(detF*JP))*Iden
      dJPdPsci=three*gama*JP
      SpUCMod = dTdJP*dJPdPsci*dPscidt*dtdC
      ! Compute the oxygen concentration - displacement modulus (Kcu)
      !
      SpCUModfac = Dox*Iden
!print results for the first increment in log file
      if((kinc.le.4).and.(kinc.ge.0)) then
      write(*,*) '**************************************** '  
      write(*,*) '---------- Compute the CU/UC modulus'
      write(*,*) '----------JP,dTdJP,dJPdPsci',Jp,dTdJP,dJPdPsci
      write(*,*) '----------SpUCMod',SpUCMod 
      endif
      return
      end subroutine integ

************************************************************************
************************************************************************
************************************************************************
      subroutine AssembleElement(nDim,nNode,ndofel,
     +     Ru,Rc,Kuu,Kuc,Kcu,Kcc,
     +     rhs,amatrx)
      
      !
      ! Subroutine to assemble the local elements residual and tangent
      implicit none

      integer i,j,k,l,m,n,A11,A12,B11,B12,nDim,nNode,nDofEl,nDofN

      real*8 Ru(nDim*nNode,1),Rc(nNode,1),Kuu(nDim*nNode,nDim*nNode)
      real*8 Kcc(nNode,nNode),Kuc(nDim*nNode,nNode),rhs(ndofel,1)
      real*8 Kcu(nNode,nDim*nNode),amatrx(ndofel,ndofel)


      ! Total number of degrees of freedom per node
      !
      nDofN = nDofEl/nNode
      ! init
      !
      rhs(:,1) = 0.d0
      amatrx = 0.d0

      if(nDim.eq.2) then
         !
         ! Assemble the element level residual
         !
         do i=1,nNode
            A11 = nDofN*(i-1)+1
            A12 = nDim*(i-1)+1
            !
            ! displacement
            !
            rhs(A11,1) = Ru(A12,1)
            rhs(A11+1,1) = Ru(A12+1,1)
            !
            ! oxygen concentration
            !
            rhs(A11+2,1) = Rc(i,1)
         enddo
         !
         ! Assemble the element level tangent matrix
         !
         do i=1,nNode
            do j=1,nNode
               A11 = nDofN*(i-1)+1
               A12 = nDim*(i-1)+1
               B11 = nDofN*(j-1)+1
               B12 = nDim*(j-1)+1
               !
               ! displacement
               !
               amatrx(A11,B11) = Kuu(A12,B12)
               amatrx(A11,B11+1) = Kuu(A12,B12+1)
               amatrx(A11+1,B11) = Kuu(A12+1,B12)
               amatrx(A11+1,B11+1) = Kuu(A12+1,B12+1)
               !
               ! oxygen concentration
               !
               amatrx(A11+2,B11+2) = Kcc(i,j)
               !
               ! displacement - oxygen concentration
               !
               amatrx(A11,B11+2) = Kuc(A12,j)
               amatrx(A11+1,B11+2) = Kuc(A12+1,j)
               !
               !oxygen concentration - displacement
               !
               amatrx(A11+2,B11) = Kcu(i,B12)
               amatrx(A11+2,B11+1) = Kcu(i,B12+1)
               !
            enddo
         enddo
      else
         write(*,*) 'How did you get nDim=',nDim
         call xit
      endif

      return
      end subroutine AssembleElement

!****************************************************************************
!     Element subroutines(guass points, shape functions and its deravatives, matrix algebra)
!****************************************************************************

      subroutine xint2D1pt(xi,w,nIntPt)
      !
      ! This subroutine will get the integration point locations and corresponding gauss quadrature weights for 2D elements
      !  using 1 gauss point for integration
      !
      !  xi(nIntPt,2): xi,eta coordinates for the integration pts w(nIntPt):    corresponding integration weights
      !
      implicit none
      !
      integer nIntPt,nDim
      !
      real*8 xi(1,2), w(1)
      ! Initialize
      !
      w = 0.d0
      xi = 0.d0
      ! Number of Gauss points
      !
      nIntPt = 1
      ! Gauss weights
      !
      w = 4.d0
      ! Gauss pt location in master element
      !
      xi(1,1) = 0.d0
      xi(1,2) = 0.d0
      return
      end subroutine xint2D1pt
      
!************************************************************************

      subroutine xint2D4pt(xi,w,nIntPt)
      !
      ! This subroutine will get the integration point locations and corresponding gauss quadrature weights for 2D elements
      !  using 4 gauss points for integration
      !
      !  xi(nIntPt,2): xi,eta coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights
      !
      implicit none
      !
      integer nIntPt,nDim
      !
      real*8 xi(4,2), w(4)
      ! Initialize
      !
      w = 0.d0
      xi = 0.d0
      ! Number of Gauss points
      !
      nIntPt = 4
      ! Gauss weights
      !
      w(1) = 1.d0
      w(2) = 1.d0
      w(3) = 1.d0
      w(4) = 1.d0
      ! Gauss pt locations in master element
      !
      xi(1,1) = -dsqrt(1.d0/3.d0)
      xi(1,2) = -dsqrt(1.d0/3.d0)
      xi(2,1) = dsqrt(1.d0/3.d0)
      xi(2,2) = -dsqrt(1.d0/3.d0)
      xi(3,1) = -dsqrt(1.d0/3.d0)
      xi(3,2) = dsqrt(1.d0/3.d0)
      xi(4,1) = dsqrt(1.d0/3.d0)
      xi(4,2) = dsqrt(1.d0/3.d0)
      return
      end subroutine xint2D4pt

!************************************************************************

      subroutine calcShape2DLinear(nIntPt,xi_int,intpt,sh,dshxi)
      !
      ! Calculate the shape functions and their derivatives at the given integration point in the master element
      !
      !                          eta
      !   4-----------3          |
      !   |           |          |
      !   |           |          |
      !   |           |          |
      !   |           |          |
      !   |           |          O--------- xi
      !   1-----------2        origin at center
      !
      !
      ! sh(i) = shape function of node i at the intpt.
      ! dshxi(i,j) = derivative wrt j direction of shape fn of node i
      !
      implicit none
      !
      integer intpt,nDim,nIntPt
      !
      real*8 xi_int(nIntPt,2),sh(4),dshxi(4,2),xi,eta
      !
      real*8 zero,one,fourth
      parameter(zero=0.d0,one=1.d0,fourth=1.d0/4.d0)
      ! Location in the master element
      !
      xi = xi_int(intpt,1)
      eta = xi_int(intpt,2)
      ! The shape functions
      !
      sh(1) = fourth*(one - xi)*(one - eta)
      sh(2) = fourth*(one + xi)*(one - eta)
      sh(3) = fourth*(one + xi)*(one + eta)
      sh(4) = fourth*(one - xi)*(one + eta) 
      ! The first derivatives
      !
      dshxi(1,1) = -fourth*(one - eta)
      dshxi(1,2) = -fourth*(one - xi)
      dshxi(2,1) = fourth*(one - eta)
      dshxi(2,2) = -fourth*(one + xi)
      dshxi(3,1) = fourth*(one + eta)
      dshxi(3,2) = fourth*(one + xi)
      dshxi(4,1) = -fourth*(one + eta)
      dshxi(4,2) = fourth*(one - xi)
      return
      end subroutine calcShape2DLinear
************************************************************************

      subroutine mapShape2D(nNode,dshxi,coords,dsh,detMapJ,stat)
      !
      !this subroutine Maps derivatives of shape fns from xi-eta-zeta domain to x-y-z domain (in each integration point).
      !
      implicit none
      !
      integer i,j,k,nNode,ieror,stat
      !
      real*8 dshxi(nNode,2),dsh(nNode,2),coords(3,nNode),mapJ(2,2),
     +  mapJ_inv(2,2),detmapJ
      !
      real*8 zero,one,two,half,fourth,eighth
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,fourth=0.25d0,
     +     eighth=1.d0/8.d0)
      
      ! Calculate the mapping Jacobian matrix:
      !
      mapJ = zero
      do i=1,2
        do j=1,2
          do k=1,nNode
              mapJ(i,j) = mapJ(i,j) + dshxi(k,i)*coords(j,k)
          end do
        end do
      end do
      ! Calculate the inverse and the determinant of Jacobian
      !
      call matInv2D(mapJ,mapJ_inv,detMapJ,stat)
      if(stat.eq.0) then
         write(*,*) 'Problem: detF.lt.zero in mapShape2D'
         call xit
      endif
      ! Calculate first derivatives wrt x, y, z
      !
      dsh = transpose(matmul(mapJ_inv,transpose(dshxi)))

      return
      end subroutine mapShape2D

!*************************************************************************

      subroutine mapShape2Da(nNode,dshxi,coords,dsh,detMapJ,stat)
      !
      !this subroutine Maps derivatives of shape fns from xi-eta-zeta domain to x-y-z domain.
      !
      ! This subroutine is exactly the same as the regular mapShape2D with the exception that coords(2,nNode) 
      ! here and coords(3,nNode) in the regular.  I have noticed that a "heat transfer" and  "static" step uses MCRD=2,
      ! but for "coupled-temperature-displacement" you will get MCRD=3, even for a plane analysis.
      !
      implicit none
      !
      integer i,j,k,nNode,ieror,stat
      !
      real*8 dshxi(nNode,2),dsh(nNode,2),coords(2,nNode),mapJ(2,2),
     +  mapJ_inv(2,2),detmapJ
      !
      real*8 zero,one,two,half,fourth,eighth
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,fourth=0.25d0,
     +     eighth=1.d0/8.d0)

      ! Calculate the mapping Jacobian matrix:
      !
      mapJ = zero
      do i=1,2
        do j=1,2
          do k=1,nNode
              mapJ(i,j) = mapJ(i,j) + dshxi(k,i)*coords(j,k)
          end do
        end do
      end do


      ! Calculate the inverse and the determinant of Jacobian
      !
      call matInv2D(mapJ,mapJ_inv,detMapJ,stat)
      if(stat.eq.0) then
         write(*,*) 'Problem: detF.lt.zero in mapShape2Da'
         call xit
      endif


      ! Calculate first derivatives wrt x, y, z
      !
      dsh = transpose(matmul(mapJ_inv,transpose(dshxi)))
      

      return
      end subroutine mapShape2Da
!****************************************************************************
!     matrix algebra subroutines
!****************************************************************************

      subroutine matInv3D(A,A_inv,det_A,istat)
      !
      ! Returns A_inv, the inverse and det_A, the determinant
      ! Note that the det is of the original matrix, not the
      ! inverse
      !
      implicit none
      !
      integer istat
      !
      real*8 A(3,3),A_inv(3,3),det_A,det_A_inv
      istat = 1 
      det_A = A(1,1)*(A(2,2)*A(3,3) - A(3,2)*A(2,3)) -
     +        A(2,1)*(A(1,2)*A(3,3) - A(3,2)*A(1,3)) +
     +        A(3,1)*(A(1,2)*A(2,3) - A(2,2)*A(1,3))
      
      if (det_A .le. 0.d0) then
        write(*,*) 'WARNING: subroutine matInv3D:'
        write(*,*) 'WARNING: det of mat=',det_A
        istat = 0
        return
      end if
          
      det_A_inv = 1.d0/det_A
        
      A_inv(1,1) = det_A_inv*(A(2,2)*A(3,3)-A(3,2)*A(2,3))
      A_inv(1,2) = det_A_inv*(A(3,2)*A(1,3)-A(1,2)*A(3,3))
      A_inv(1,3) = det_A_inv*(A(1,2)*A(2,3)-A(2,2)*A(1,3))
      A_inv(2,1) = det_A_inv*(A(3,1)*A(2,3)-A(2,1)*A(3,3))
      A_inv(2,2) = det_A_inv*(A(1,1)*A(3,3)-A(3,1)*A(1,3))
      A_inv(2,3) = det_A_inv*(A(2,1)*A(1,3)-A(1,1)*A(2,3))
      A_inv(3,1) = det_A_inv*(A(2,1)*A(3,2)-A(3,1)*A(2,2))
      A_inv(3,2) = det_A_inv*(A(3,1)*A(1,2)-A(1,1)*A(3,2))
      A_inv(3,3) = det_A_inv*(A(1,1)*A(2,2)-A(2,1)*A(1,2))
      return
      end subroutine matInv3D
!****************************************************************************

      subroutine matInv2D(A,A_inv,det_A,istat)
      !
      ! Returns A_inv, the inverse, and det_A, the determinant
      ! Note that the det is of the original matrix, not the
      ! inverse
      !
      implicit none
      !
      integer istat
      !
      real*8 A(2,2),A_inv(2,2),det_A,det_A_inv 
      istat = 1
      det_A = A(1,1)*A(2,2) - A(1,2)*A(2,1)
        
      if (det_A .le. 0.d0) then
        write(*,*) 'WARNING: subroutine matInv2D:'
        write(*,*) 'WARNING: det of mat=',det_A
        istat = 0
        return
      end if
            
      det_A_inv = 1.d0/det_A
          
      A_inv(1,1) =  det_A_inv*A(2,2)
      A_inv(1,2) = -det_A_inv*A(1,2)
      A_inv(2,1) = -det_A_inv*A(2,1)
      A_inv(2,2) =  det_A_inv*A(1,1)
      return
      end subroutine matInv2D
!****************************************************************************

      subroutine mdet(A,det)
      !
      ! This subroutine calculates the determinant
      ! of a 3 by 3 matrix [A]
      !
      implicit none
      !
      real*8  A(3,3),det

      det = A(1,1)*A(2,2)*A(3,3) 
     +	  + A(1,2)*A(2,3)*A(3,1)
     +	  + A(1,3)*A(2,1)*A(3,2)
     +	  - A(3,1)*A(2,2)*A(1,3)
     +	  - A(3,2)*A(2,3)*A(1,1)
     +	  - A(3,3)*A(2,1)*A(1,2)
      return
      end subroutine mdet	
!****************************************************************************

      subroutine onem(A)
      !
      ! This subroutine stores the identity matrix in the
      ! 3 by 3 matrix [A]
      !
      implicit none
      !
      integer i,j
      real*8 A(3,3)

      do i=1,3
         do J=1,3
	    if (i .eq. j) then
              A(i,j) = 1.0
            else
              A(i,j) = 0.0
            end if
         end do
      end do
      return
      end subroutine onem
****************************************************************************
