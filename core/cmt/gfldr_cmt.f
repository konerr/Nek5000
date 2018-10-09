#ifndef NOMPIIO

      subroutine gfldr_cmt(sourcefld)
c
c     generic field file reader
c     reads sourcefld and interpolates all avaiable fields
c     onto current mesh
c
c     memory requirement: 
c     nelgs*nxs**ldim < np*(4*lelt*lx1**ldim)
c
      include 'SIZE'
      include 'TOTAL'
      include 'RESTART'
      include 'GFLDR'
c      include 'SOLN'
      COMMON /solnconsvar/ U(LX1,LY1,LZ1,TOTEQ,lelt) 

      character sourcefld*(*)

      common /scrcg/  pm1(lx1*ly1*lz1,lelv)
      common /nekmpi/ nidd,npp,nekcomm,nekgroup,nekreal

      COMMON /SCRNS/      OTVAR(LX1,LY1,LZ1,lelt,5)
      real                OTVAR

      character*1   hdr(iHeaderSize)

      integer*8 dtmp8

      logical ifbswp, if_byte_swap_test
      real*4 bytetest
      
      integer e

      etime_t = dnekclock_sync()
      if(nio.eq.0) write(6,*) 'call gfldr ',trim(sourcefld) 

      ! open source field file
      ierr = 0
      if(nid.eq.0) then
        open (90,file=sourcefld,status='old',err=100)
        close(90)
        goto 101
 100    ierr = 1
 101  endif
      call err_chk(ierr,' Cannot open source fld file!$')
      call byte_open_mpi(sourcefld,fldh_gfldr,.true.,ierr)

      ! read and parse header
      call byte_read_mpi(hdr,iHeaderSize/4,0,fldh_gfldr,ierr)
      call byte_read_mpi(bytetest,1,0,fldh_gfldr,ierr)

      call mfi_parse_hdr(hdr,ierr)
      call err_chk(ierr,' Invalid header!$')
      ifbswp = if_byte_swap_test(bytetest,ierr)
      call err_chk(ierr,' Invalid endian tag!$')

      nelgs   = nelgr
      nxs     = nxr
      nys     = nyr
      nzs     = nzr
      if(nzs.gt.1) then 
        ldims = 3
      else
        ldims = 2
      endif
      if (ifgtim) time = timer

      ! distribute elements across all ranks
      nels = nelgs/np
      do i = 0,mod(nelgs,np)-1
         if(i.eq.nid) nels = nels + 1
      enddo
      nxyzs      = nxs*nys*nzs
      dtmp8      = nels
      ntots_b    = dtmp8*nxyzs*wdsizr
      rankoff_b  = igl_running_sum(nels) - dtmp8
      rankoff_b  = rankoff_b*nxyzs*wdsizr  
      dtmp8      = nelgs
      nSizeFld_b = dtmp8*nxyzs*wdsizr
      noff0_b    = iHeaderSize + iSize + iSize*dtmp8

      ! do some checks
      if(ldims.ne.ldim) 
     $ call exitti('ldim of source does not match target!$',0)
      if(ntots_b/wdsize .gt. ltots) then
        dtmp8 = nelgs
        lelt_req = dtmp8*nxs*nys*nzs / (np*ltots/lelt)
        lelt_req = lelt_req + 1
        if(nio.eq.0) write(6,*)
     $   'ABORT: buffer too small, increase lelt > ', lelt_req
        call exitt
      endif

      ifldpos = 0
      if(ifgetxr) then
        ! read source mesh coordinates
        call gfldr_getxyz_cmt(xm1s,ym1s,zm1s,ifbswp)
        ifldpos = ldim
      else
        call exitti('source does not contain a mesh!$',0)
      endif

      if(if_full_pres) then
        call exitti('no support for if_full_pres!$',0)
      endif

      if(nelt.ne.nelv) then
        call exitti('no support for conj/HT!$',0)
      endif

      ! initialize interpolation tool using source mesh
      nxf   = 2*nxs
      nyf   = 2*nys
      nzf   = 2*nzs
      nhash = nxs*nys*nzs 
      nmax  = 256

      call fgslib_findpts_setup(inth_gfldr,nekcomm,np,ldim,
     &                          xm1s,ym1s,zm1s,nxs,nys,nzs,
     &                          nels,nxf,nyf,nzf,bb_t,
     &                          nhash,nhash,nmax,tol)

      n = lx1*ly1*lz1
      do e=1,nelt
         call copy(otvar(1,1,1,e,1),u(1,1,1,1,e),n)
         call copy(otvar(1,1,1,e,2),u(1,1,1,2,e),n)
         call copy(otvar(1,1,1,e,3),u(1,1,1,3,e),n)
         call copy(otvar(1,1,1,e,4),u(1,1,1,4,e),n)
         call copy(otvar(1,1,1,e,5),u(1,1,1,5,e),n)
      enddo

      ! read source fields and interpolate
      if(ifgetur) then
c        if(nid.eq.0 .and. loglevel.gt.2) write(6,*) 'reading vel'
        if(nid.eq.0) write(6,*) 'reading vel'
c        call gfldr_getfld(vx,vy,vz,ldim,ifldpos+1,ifbswp)
        call gfldr_getfld_cmt(otvar(1,1,1,1,2),
     >                        otvar(1,1,1,1,3),
     >                        otvar(1,1,1,1,4),ldim,ifldpos+1,ifbswp)
c        call gfldr_getfld_cmt(u(1,1,1,2,1),
c     >                        u(1,1,1,3,1),
c     >                        u(1,1,1,4,1),ldim,ifldpos+1,ifbswp)
        ifldpos = ifldpos + ldim
      endif
      if(ifgetpr) then
c        if(nid.eq.0 .and. loglevel.gt.2) write(6,*) 'reading pr'
        if(nid.eq.0) write(6,*) 'reading pr'
c        call gfldr_getfld(pm1,dum,dum,1,ifldpos+1,ifbswp)
        call gfldr_getfld_cmt(otvar(1,1,1,1,1),dum,dum,
     >                    1,ifldpos+1,ifbswp)
c        call gfldr_getfld_cmt(u(1,1,1,1,1),dum,dum,
c     >                    1,ifldpos+1,ifbswp)
        ifldpos = ifldpos + 1
        if (ifaxis) call axis_interp_ic(pm1)
        call map_pm1_to_pr(pm1,1)
      endif
      if(ifgettr .and. ifheat) then
c        if(nid.eq.0 .and. loglevel.gt.2) write(6,*) 'reading temp'
        if(nid.eq.0) write(6,*) 'reading temp'
c        call gfldr_getfld(t(1,1,1,1,1),dum,dum,1,ifldpos+1,ifbswp)
        call gfldr_getfld_cmt(otvar(1,1,1,1,5),dum,dum,
     >                    1,ifldpos+1,ifbswp)
c        call gfldr_getfld_cmt(u(1,1,1,5,1),dum,dum,
c     >                    1,ifldpos+1,ifbswp)
        ifldpos = ifldpos + 1
      endif
      do i = 1,ldimt-1
         if(ifgtpsr(i)) then
           if(nid.eq.0 .and. loglevel.gt.2) 
     $       write(6,*) 'reading scalar',i
           call gfldr_getfld_cmt(t(1,1,1,1,i+1),dum,dum,
     >                           1,ifldpos+1,ifbswp) 
           ifldpos = ifldpos + 1
         endif
      enddo

      call byte_close_mpi(fldh_gfldr,ierr)
      call fgslib_findpts_free(inth_gfldr)

      do e=1,nelt
         call copy(u(1,1,1,1,e),otvar(1,1,1,e,1),n)
         call copy(u(1,1,1,2,e),otvar(1,1,1,e,2),n)
         call copy(u(1,1,1,3,e),otvar(1,1,1,e,3),n)
         call copy(u(1,1,1,4,e),otvar(1,1,1,e,4),n)
         call copy(u(1,1,1,5,e),otvar(1,1,1,e,5),n)
      enddo

      etime_t = dnekclock_sync() - etime_t
      if(nio.eq.0) write(6,'(A,1(1g8.2),A)')
     &                   ' done :: gfldr  ', etime_t, ' sec'

      return
      end
c-----------------------------------------------------------------------
      subroutine gfldr_getxyz_cmt(xout,yout,zout,ifbswp)

      include 'SIZE'
      include 'GFLDR'
      include 'RESTART'

      real xout(*)
      real yout(*)
      real zout(*)
      logical ifbswp

      integer*8 ioff_b

 
      ioff_b = noff0_b + ldim*rankoff_b
      call byte_set_view(ioff_b,fldh_gfldr)

      nread = ldim*ntots_b/4
      call byte_read_mpi(bufr,nread,-1,fldh_gfldr,ierr)
      if(ifbswp) then
        if(wdsizr.eq.4) call byte_reverse (bufr,nread,ierr)
        if(wdsizr.eq.8) call byte_reverse8(bufr,nread,ierr)
      endif

      call gfldr_buf2vi_cmt (xout,1,bufr,ldim,wdsizr,nels,nxyzs)
      call gfldr_buf2vi_cmt (yout,2,bufr,ldim,wdsizr,nels,nxyzs)
      if(ldim.eq.3)
     $ call gfldr_buf2vi_cmt(zout,3,bufr,ldim,wdsizr,nels,nxyzs)

      return
      end
c-----------------------------------------------------------------------
      subroutine gfldr_getfld_cmt(out1,out2,out3,nldim,ifldpos,ifbswp)

      include 'SIZE'
      include 'GEOM'
      include 'GFLDR'
      include 'RESTART'

      real out1(*)
      real out2(*)
      real out3(*)
      logical ifbswp

      integer*8 ioff_b

      logical ifpts

      integer icalld
      save    icalld
      data    icalld /0/


      ifpts = .false.
      if(icalld.eq.0) then
        ifpts = .true. ! find points
        icalld = 1
      endif

      ! read field data from source fld file
      ioff_b = noff0_b + (ifldpos-1)*nSizeFld_b
      ioff_b = ioff_b  + nldim*rankoff_b
      call byte_set_view(ioff_b,fldh_gfldr)
      nread = nldim*ntots_b/4
      call byte_read_mpi(bufr,nread,-1,fldh_gfldr,ierr)
      if(ifbswp) then
        if(wdsizr.eq.4) call byte_reverse (bufr,nread,ierr)
        if(wdsizr.eq.8) call byte_reverse8(bufr,nread,ierr)
      endif

      ! interpolate onto current mesh
      ntot = lx1*ly1*lz1*nelt
      call gfldr_buf2vi_cmt  (buffld,1,bufr,nldim,wdsizr,nels,nxyzs)
      if (nio.eq.0) write(6,*) 'vi out1'
      call gfldr_intp_cmt    (out1,buffld,ifpts)
      if (nio.eq.0) write(6,*) 'intp out1'
      if(nldim.eq.1) return

      call gfldr_buf2vi_cmt  (buffld,2,bufr,nldim,wdsizr,nels,nxyzs)
      if (nio.eq.0) write(6,*) 'vi out2'
      call gfldr_intp_cmt    (out2,buffld,.false.)
      if (nio.eq.0) write(6,*) 'intp out2'
      if(nldim.eq.2) return

      if(nldim.eq.3) then
        call gfldr_buf2vi_cmt(buffld,3,bufr,nldim,wdsizr,nels,nxyzs)
        if (nio.eq.0) write(6,*) 'vi out3'
        call gfldr_intp_cmt  (out3,buffld,.false.)
        if (nio.eq.0) write(6,*) 'intp out3'
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine gfldr_buf2vi_cmt(vi,index,buf,ldim,wds,nel,nxyz)

      real    vi(*)
      real*4  buf(*)
      integer wds


      do iel = 1,nel
         j = (iel-1)*nxyz
         k = (iel-1)*ldim*nxyz

         if(index.eq.2) k = k+nxyz
         if(index.eq.3) k = k+2*nxyz

         if(wds.eq.4) call copy4r(vi(j+1),buf(k+1)  ,nxyz)
         if(wds.eq.8) call copy  (vi(j+1),buf(2*k+1),nxyz)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine gfldr_intp_cmt(fieldout,fieldin,iffpts)

      include 'SIZE'
      include 'RESTART'
      include 'GEOM'
      include 'GFLDR'

      real    fieldout(*),fieldin(*)
      logical iffpts

      integer*8 i8glsum,nfail,nfail_sum


      nfail = 0
      ntot  = lx1*ly1*lz1*nelt

      toldist = 5e-6
      if(wdsizr.eq.8) toldist = 5e-14

      if(iffpts) then ! locate points (iel,iproc,r,s,t)

        call fgslib_findpts(inth_gfldr,
     &                      grcode,1,
     &                      gproc,1,
     &                      gelid,1,
     &                      grst,ldim,
     &                      gdist,1,
     &                      xm1,1,
     &                      ym1,1,
     &                      zm1,1,ntot)

        do i=1,ntot
           if(grcode(i).eq.1 .and. sqrt(gdist(i)).gt.toldist)
     &       nfail = nfail + 1
           if(grcode(i).eq.2) nfail = nfail + 1
        enddo

        nfail_sum = i8glsum(nfail,1)
        if(nfail_sum.gt.0) then
          if(nio.eq.0) write(6,*)
     &      ' WARNING: Unable to find all mesh points in source fld ',
     &      nfail_sum
        endif

      endif

      ! evaluate inut field at given points
      call fgslib_findpts_eval(inth_gfldr,
     &                         fieldout,1,
     &                         grcode,1,
     &                         gproc,1,
     &                         gelid,1,
     &                         grst,ldim,ntot,
     &                         fieldin)

      return
      end

#endif
