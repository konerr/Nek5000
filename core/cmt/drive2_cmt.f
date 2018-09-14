C> @file drive2_cmt.f mid-level initialization drivers. Not long for this world.
c-----------------------------------------------------------------------
      subroutine nek_cmt_init
      include 'SIZE'
      include 'TOTAL'
      include 'DG'
      include 'CMTDATA'
      if (nio.eq.0) write(6,*)'Set up CMT-Nek'    
      if (toteq.ne.5) then
         if (nio.eq.0) write(6,*)'toteq is low ! toteq = ',toteq
         if (nio.eq.0) write(6,*) 'Reset toteq in SIZE to 5'
         call exitt
      endif
      if (ifrestart) then
         ifheat = .true. ! almost certainly incorrect
      endif
      call setup_cmt_commo
      
c     call setup_cmt_param
      return
      end

!-----------------------------------------------------------------------

      subroutine izero8(a,n)
      integer*8 a(1)
      do i=1,n
         a(i)=0
      enddo
      return
      end
!-----------------------------------------------------------------------
      subroutine setup_cmt_param
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'CMTDATA'
      INCLUDE 'CMTBCDATA'

      real  MixtPerf_R_CpG, MixtPerf_T_DPR, MixtPerf_C_GRT
     >                 ,MixtPerf_Ho_CpTUVW,MixtPerf_Cp_CvR,MixtPerf_R_M
     >                 ,MixtPerf_G_CpR      
      external MixtPerf_R_CpG, MixtPerf_T_DPR, MixtPerf_C_GRT
     >                 ,MixtPerf_Ho_CpTUVW,MixtPerf_Cp_CvR,MixtPerf_R_M
     >                 ,MixtPerf_G_CpR      

      cip_adhoc=10.0
      cvgref     = param(104)
c     gmaref     = param(105)
      molmass    = param(106)
      muref      = param(107)
      coeflambda = param(108)
      suthcoef   = param(109)
      reftemp    = param(110)
      prlam      = param(111)
      pinfty     = param(112)
      rgasref    = MixtPerf_R_M(molmass,dum)
      cpgref     = MixtPerf_Cp_CvR(cvgref,rgasref)
      gmaref     = MixtPerf_G_CpR(cpgref,rgasref) 
! put these in rea file someday
      return
      end
c------------------------------------------------------------------------

      subroutine limiter
! EBDG Stuff. WHERE'S PHI????
      include 'SIZE'
      include 'TOTAL'
      include 'CMTDATA'
      include 'NEKUSE'
      parameter (lxyz=lx1*ly1*lz1)
      common /scrns/ energy(lxyz),scr(lxyz),avstate(toteq),
     > otvar(lx1,ly1,lz1,lelt,1)
      real energy,scr,avstate,otvar
      integer e,eg

      nxyz=lx1*ly1*lz1
      ntot=nxyz*nelt

      epslon=1.0e-13


      rgam=rgasref/(gmaref-1.0)
!      do i=1,ntot
!         rho=max(vtrans(i,1,1,1,irho),1.0e-10)
!!        scr(i,1)=rgam*log(pr(i,1,1,1)/(rho**gmaref))
!         scr(i,1)=log(pr(i,1,1,1)/(rho**gmaref))
!      enddo
!!     call dsop(scr,'MIN',lx1,ly1,lz1)
!      call copy(t(1,1,1,1,4),scr,ntot)
!! elemental entropy minimum
!!     do e=1,nelt
!!        se0(e)=vlmin(scr(1,e),nxyz)
!!     enddo

!      do e=1,nelt
!         call copy(otvar(1,1,1,e,1),u(1,1,1,1,e),nxyz)
!         call copy(otvar(1,1,1,e,5),u(1,1,1,5,e),nxyz)
!      end do
!      rhomin = glmin(otvar(1,1,1,1,1),ntot)
!      emin   = glmin(otvar(1,1,1,1,5),ntot)

      do e=1,nelt

!        rhomin=vlmin(vtrans(1,1,1,e,irho),nxyz)
         rhomin=vlmin(u(1,1,1,1,e),nxyz)

!! positivity-preserving limiter of Zhang and Shu: density
!         rho=vlsc2(bm1(1,1,1,e),u(1,1,1,1,e),nxyz)/volel(e)
!         if (abs(rho-rhomin) .gt. epslon) then
!            theta=min((rho-epslon)/(rho-rhomin+epslon),1.0)
!            do i=1,nxyz
!               uold=u(i,1,1,1,e)
!               u(i,1,1,1,e)=rho+theta*(uold-rho)
!            enddo
!         else
!            theta=1.0
!         endif
!! positivity-preserving limiter of Cheng and Shu: density
!!         rho=vlsc2(bm1(1,1,1,e),u(1,1,1,1,e),nxyz)/volel(e)
!!!         if (abs(rho-rhomin) .gt. epslon) then
!!            do i=1,nxyz
!!               otvar(i,1,1,1,e)=(rho-epslon)/(rho-u(i,1,1,1,e)+epslon)
!!            enddo
!!            theta1=vlmin(otvar(i,1,1,1,e),nxyz)
!!            theta1=min(1.0,theta1)
!!            do i=1,nxyz
!!               uold=u(i,1,1,1,e)
!!               u(i,1,1,1,e)=rho+theta*(uold-rho)
!!            enddo
!!!         else
!!!            theta=1.0
!!!         endif
!! s3 in visit
!         call cfill(t(1,1,1,e,4),theta,nxyz)
!
!! positivity-preserving limiter of Cheng and Shu: internal energy
!         avstate(1)=rho
!         do m=2,toteq
!            avstate(m)=vlsc2(bm1(1,1,1,e),u(1,1,1,m,e),nxyz)/volel(e)
!         enddo
!         eavg=avstate(5)-
!     >              0.5*(avstate(2)**2+avstate(3)**2+avstate(4)**2)/
!     >                   rho
!         eavg=eavg/rho
!
!         call invcol3(vx(1,1,1,e),u(1,1,1,irpu,e),u(1,1,1,irg,e),nxyz)
!         call invcol3(vy(1,1,1,e),u(1,1,1,irpv,e),u(1,1,1,irg,e),nxyz)
!!        if (if3d)
!         call invcol3(vz(1,1,1,e),u(1,1,1,irpw,e),u(1,1,1,irg,e),nxyz)
!! first kinetic energy
!         if (if3d) then
!            call vdot3(scr,
!     >             u(1,1,1,irpu,e),u(1,1,1,irpv,e),u(1,1,1,irpw,e),
!     >             u(1,1,1,irpu,e),u(1,1,1,irpv,e),u(1,1,1,irpw,e),nxyz)
!         else
!            call vdot2(scr,u(1,1,1,irpu,e),u(1,1,1,irpv,e),
!     >                     u(1,1,1,irpu,e),u(1,1,1,irpv,e),nxyz)
!         endif
!         call invcol2(scr,u(1,1,1,irg,e),nxyz)
!         call cmult(scr,0.5,nxyz)
!! then to internal energy
!         call sub3(energy,u(1,1,1,iret,e),scr,nxyz)
!! now mass-specific
!         call invcol2(energy,u(1,1,1,irg,e),nxyz)
!! Compute \theta_x
!         do i=1,nxyz
!            if(energy(i).lt.0.0) then
!               otvar(i,1,1,1,e)=eavg/(eavg-energy(i))
!            else
!               otvar(i,1,1,1,e)=1.0
!            end if
!         enddo
!! Get \theta_i^2
!         theta2=abs(vlmin(otvar(1,1,1,1,e),nxyz))
!         do m=1,toteq
!            do i=1,nxyz
!               uold=u(i,1,1,m,e)
!               u(i,1,1,m,e)=avstate(m)+theta2*(uold-avstate(m))
!            enddo
!         enddo
!
!         call cfill(t(1,1,1,e,5),theta2,nxyz)

! positivity-preserving limiter of Zhang and Shu: internal energy???
!         emin=vlmin(u(1,1,1,5,e),nxyz)
! first kinetic energy
!         if (if3d) then
!            call vdot3(scr,
!     >             u(1,1,1,irpu,e),u(1,1,1,irpv,e),u(1,1,1,irpw,e),
!     >             u(1,1,1,irpu,e),u(1,1,1,irpv,e),u(1,1,1,irpw,e),nxyz)
!         else
!            call vdot2(scr,u(1,1,1,irpu,e),u(1,1,1,irpv,e),
!     >                     u(1,1,1,irpu,e),u(1,1,1,irpv,e),nxyz)
!         endif
!         call invcol2(scr,u(1,1,1,irg,e),nxyz)
!         call cmult(scr,0.5,nxyz)
!         energy_eps=vlmax(scr,nxyz)
!      
!         rhoe=vlsc2(bm1(1,1,1,e),u(1,1,1,iret,e),nxyz)/volel(e)
!         if (abs(rhoe-emin) .gt. energy_eps) then
!            theta=min((rhoe-energy_eps)/(rhoe-emin+epslon),1.0)
!            do i=1,nxyz
!               uold=u(i,1,1,iret,e)
!               u(i,1,1,iret,e)=rhoe+theta*(uold-rhoe)
!            enddo
!         else
!           theta=1.0
!         endif
!         call cfill(t(1,1,1,e,4),theta,nxyz)

         do m=1,toteq
            avstate(m)=vlsc2(bm1(1,1,1,e),u(1,1,1,m,e),nxyz)/volel(e)
         enddo
         rho=avstate(1)
! Entropy-bounded limiter of Lv and Ihme
         e_internal=avstate(5)-
     >              0.5*(avstate(2)**2+avstate(3)**2+avstate(4)**2)/
     >                   rho
         e_internal=e_internal/rho
         eg=gllel(e)
         call cmt_userEOS(1,1,1,eg) ! assigns elm avg to  pres and temp
         do i=1,nxyz
            scr(i)=pr(i,1,1,e)-exp(se0const)*
     >                         (u(i,1,1,1,e)**gmaref)
         enddo
         tau=vlmin(scr,nxyz)
         tau=min(tau,0.0)
         epsebdg(e)=tau/
     >          (tau-(pres-exp(se0const)*rho**gmaref))
         epsebdg(e)=min(epsebdg(e),1.0)
         epsebdg(e)=max(epsebdg(e),0.0)
! diagnostic
         call cfill(t(1,1,1,e,6),epsebdg(e),nxyz)

         do m=1,toteq
            do i=1,nxyz
               uold=u(i,1,1,m,e)
               u(i,1,1,m,e)=uold+epsebdg(e)*(avstate(m)-uold)
            enddo
         enddo

      enddo
      
      return
      end

!-----------------------------------------------------------------------
! JH082118 shock detectors. Test shock detectors in usr file extensively.
!          port them here after they have earned our trust.
!-----------------------------------------------------------------------

      subroutine AVeverywhere(shkdet)
      include 'SIZE'
      include 'TOTAL'
      real shkdet(nelt)
      call rone(shkdet,nelt)
      return
      end

      subroutine limiter_only(shkdet)
! worst ad hoc shock detector ever. read a limiter function from CMTDATA
! epsebdg is the only one we have so far. If it's bigger than some threshold,
! apply AV there too.
      include 'SIZE'
      include 'CMTDATA'
      real shkdet(nelt)
      integer e

      call rzero(shkdet,nelt)
      tol=1.0e-5

      if (.not.time4av) then
         do e=1,nelt
            if (abs(epsebdg(e)) .gt. tol) shkdet(e)=1.0
         enddo
      endif

      emax=glamax(shkdet,nelt)
      if (nio.eq.0) write(6,*) 'max shock detector =',emax

      return
      end
