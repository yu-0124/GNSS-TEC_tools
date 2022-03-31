c
      program tomo
c
c     3-d tomography of ionospheric anomalies before earthquakes
c     in Tohoku
c
      implicit real*8 (a-h,o-z)
c
      parameter (nblock_e=10,nblock_n=12,nblock_u=13,maxpara=8000    ! fine version
     1          ,maxstn=2000,maxsat=70,maxdat=15000)
c
      dimension cent(3,-nblock_e:nblock_e,-nblock_n:nblock_n,nblock_u)
     1        ,boundary_e(2*nblock_e+2),boundary_n(2*nblock_n+2)
     2        ,boundary_u(nblock_u+1),jsat(maxdat)
     3        ,jstn(maxdat),xdata(maxdat),ydata(maxdat)
     4        ,b(maxpara,maxpara),binv(maxpara,maxpara),aty(maxpara)
     5        ,z(maxpara),x(maxpara),a(maxpara,maxdat)
     6        ,xyzstn(3,maxstn),xyzsat(3,maxsat),num_sat(maxsat)
     7        ,keisu(3,maxpara)
      character lstn(maxstn)*4,lsit*4,lscr*16,litom*100,lisat*100
     1        ,lires*100
      logical ltikhonov
c
      data cent_lon/135.0/ ! in deg
      data cent_lat/37.0/ ! in deg
      data zero_alt/60.0/ ! in km
c
      data blo_lon/1.0/   ! in deg
      data blo_lat/0.9/   ! in deg
      data blo_alt/60.0/ ! in km
c
      call getarg(1,lscr)
      read(lscr,*) jmode
      call getarg(2,lscr)
      if(jmode.eq.2) then
        read(lscr,*) jcomp
        call getarg(3,lscr)
        read(lscr,*) alt_sample
      else if(jmode.eq.3) then
        read(lscr,*) jcomp
      else 
        read(lscr,*) continuity
        call getarg(3,lscr)
        read(lscr,*) yerr
        call getarg(4,lscr)
        read(lscr,*,end=500) vtec !background VTEC in TECU
        ltikhonov=.true.
  500   continue
        if(ltikhonov) then
         call getarg(5,lscr)
         read(lscr,*)sigma_in_percent  !allowance in Tikhonov constraint in percent of original electron density
        endif        
      endif
c
c-----counting the number of parameters
c
      n_ew=2*nblock_e+1
      n_ns=2*nblock_n+1
      n_ud=nblock_u
      npara=n_ew*n_ns*n_ud
      if(npara.gt.maxpara) stop 99
c
c-----defining the center of the blocks
c
      ip=1
      do ie=-nblock_e,nblock_e
       cent_e=cent_lon + ie*blo_lon
       do in=-nblock_n,nblock_n 
         cent_n=cent_lat + in*blo_lat
         do iu=1,nblock_u
           cent_u=zero_alt + iu*blo_alt
           cent(1,ie,in,iu)=cent_e
           cent(2,ie,in,iu)=cent_n
           cent(3,ie,in,iu)=cent_u
           if(jmode.eq.3) call draw(jcomp,cent(1,ie,in,iu)
     1                          ,blo_lon,blo_lat,blo_alt,0.0d0)
           keisu(1,ip)=ie
           keisu(2,ip)=in
           keisu(3,ip)=iu
           ip=ip+1
         enddo
       enddo
      enddo
      if(jmode.eq.3) stop
c
      ialt=nint((alt_sample-50.)/50.)
      alt_max=zero_alt + (real(ialt)+0.5)*blo_alt
      alt_min=zero_alt + (real(ialt)-0.5)*blo_alt
c
c-----a-priori continuity constraints
c
      if(jmode.eq.0.or.jmode.eq.1) then
c
c     (1) reset to zero
c
        do ipara=1,npara
          aty(ipara)=0.0
          do kpara=1,npara
            b(kpara,ipara)=0.0
          enddo
        enddo 
c
c      (1.1) Tikonov type constraint
       if(ltikhonov)then
 
        do ipara=1,npara
          iu=keisu(3,ipara)
          cent_u=zero_alt+iu*blo_alt
          call chapman(cent_u,vtec,density)
          tikhonov=density * sigma_in_percent / 100.0
          b(ipara,ipara)=1.0/tikhonov**2
        enddo
       endif
c
c     (2) continuity constraints
c
        do ipara=2,npara
          ie1=keisu(1,ipara)
          in1=keisu(2,ipara)
          iu1=keisu(3,ipara)
          do kpara=1,ipara-1
            ie2=keisu(1,kpara)
            in2=keisu(2,kpara)
            iu2=keisu(3,kpara)
            id_e=iabs(ie2-ie1)
            id_n=iabs(in2-in1)
            id_u=iabs(iu2-iu1)
            if(id_e+id_n+id_u.eq.1) then
              b(kpara,ipara)=b(kpara,ipara)-1.0/continuity**2
              b(ipara,kpara)=b(ipara,kpara)-1.0/continuity**2
              b(kpara,kpara)=b(kpara,kpara)+1.0/continuity**2
              b(ipara,ipara)=b(ipara,ipara)+1.0/continuity**2
            endif
          enddo
        enddo 
c
      endif
c
c-----reading site coord
c
      open(1,file='tomo/stn_coor.txt')
      istn=1
 1111 read(1,100,end=1122)lstn(istn),(xyzstn(k,istn),k=1,3)
  100 format(a4,3f14.4)
      istn=istn+1
      goto 1111
 1122 nstn=istn-1
      close(1)
c
c-----reading sat coord
c
      open(1,file='input/tomog2.ipt')
      read(1,'(a100)')lisat
      read(1,'(a100)')litom
      read(1,'(a100)')lires
      close(1)

      open(1,file=lisat)
      isat=1
 2222 read(1,110,end=2233)num_sat(isat),(xyzsat(k,isat),k=1,3)
  110 format(3x,i2,3f16.4)
      isat=isat+1
      goto 2222
 2233 nsat=isat-1
      close(1)
c
c-----reading the STEC_anom
c
      idat=1
      rms=0.0
 3333 read(5,120,end=3344) lsit,ksat,ydata(idat)
  120 format(a4,4x,i2,f10.4)
      rms=rms+ydata(idat)**2
      do istn=1,nstn
        if(lsit.eq.lstn(istn)) then
          jstn(idat)=istn
          goto 3334
        endif
      enddo
      stop 1
 3334 continue
      do isat=1,nsat
        if(ksat.eq.num_sat(isat)) then
          jsat(idat)=isat
          goto 3335
        endif
      enddo
      stop 2
 3335 continue
      if(jmode.eq.2) then
        call ipp(xyzstn(1,istn),xyzsat(1,isat),alt_max,plon1,plat1)
        call ipp(xyzstn(1,istn),xyzsat(1,isat),alt_min,plon2,plat2)
        if(jcomp.eq.1) then
          write(6,'(">")')
          write(6,*)alt_max,plat1
          write(6,*)alt_min,plat2
        else if(jcomp.eq.2) then
          write(6,'(">")')
          write(6,*)plon1,alt_max
          write(6,*)plon2,alt_min
        else if(jcomp.eq.3) then
          write(6,'(">")')
          write(6,*)plon1,plat1
          write(6,*)plon2,plat2
        endif
      endif
      idat=idat+1
      goto 3333
 3344 ndat=idat-1
      if(jmode.eq.2) stop
      rms1=dsqrt(rms/real(ndat))
c
c-----checking for the penetration
c
      do idat=1,ndat
       ksat=jsat(idat)
       kstn=jstn(idat)
       ip=1
       do ie=-nblock_e,nblock_e
        do in=-nblock_n,nblock_n 
          do iu=1,nblock_u
            xlon1=cent(1,ie,in,iu)-0.5*blo_lon
            xlon2=cent(1,ie,in,iu)+0.5*blo_lon
            xlat1=cent(2,ie,in,iu)-0.5*blo_lat
            xlat2=cent(2,ie,in,iu)+0.5*blo_lat
            alt1=cent(3,ie,in,iu)-0.5*blo_alt
            alt2=cent(3,ie,in,iu)+0.5*blo_alt
            call length(xyzstn(1,kstn),xyzsat(1,ksat),
     1           xlon1,xlon2,xlat1,xlat2,alt1,alt2,aleng)
c           write(6,*)ie,in,iu,aleng
            if(jmode.eq.4.and.aleng.gt.0.0) 
     1                   call draw(jcomp,cent(1,ie,in,iu)
     2                    ,blo_lon,blo_lat,blo_alt,aleng)
            a(ip,idat)=aleng/1.e5   ! to make parameters in TECU/100 km
            ip=ip+1
          enddo
        enddo
       enddo
      enddo
      if(jmode.eq.4) stop
c
c------construct normal matrix
c
      do idat=1,ndat
       do k1=1,npara
        do k2=1,npara
          b(k1,k2)=b(k1,k2)+a(k1,idat)*a(k2,idat)/yerr**2
        enddo
        aty(k1)=aty(k1)+a(k1,idat)*ydata(idat)/yerr**2
       enddo
      enddo
c
c------least-squares calculation 
c
      call lsq(b,aty,z,maxpara,npara,x,binv)
c
c------calculated values and residuals
c
      rms=0.0
      open(1,file=lires)
      do idat=1,ndat
       ksat=jsat(idat)
       kstn=jstn(idat)
       ycalc=0.0
       ip=1
       do ie=-nblock_e,nblock_e
        do in=-nblock_n,nblock_n 
          do iu=1,nblock_u
            ycalc=ycalc+a(ip,idat)*x(ip)
            ip=ip+1
          enddo
        enddo
       enddo
       res=ydata(idat)-ycalc
       rms=rms+res**2
       write(1,150)lstn(kstn),num_sat(ksat),ydata(idat),ycalc,res
  150  format(a4,'-PRN',i2,3f10.4)
      enddo
      rms2=dsqrt(rms/real(ndat))
      close(1)
c
c------least-squares calculation and output parameter results
c
      open(1,file=litom)
      write(1,130) nblock_e,nblock_n,nblock_u,cent_lon,cent_lat
     1            ,zero_alt,blo_lon,blo_lat,blo_alt,ndat,rms1,rms2
  130 format('#blocks in e:',i3,1x,'in n:',i3,1x,'in u:',i3,1x
     1      ,'Center (',f6.1,',',f5.1,')',/'Bottom height:'
     2      ,f5.1,' km',1x,'Block size:',f4.1,' x',f4.1,' x' 
     3      ,f6.1,/,'# data:',i4,5x,'rms before:',f6.3,5x
     4      ,'rms after:',f6.3,' TECU',/70('='))
      ip=1
      do ie=-nblock_e,nblock_e
        do in=-nblock_n,nblock_n 
          do iu=1,nblock_u
            xerr=rms2*sqrt(binv(ip,ip))
            write(1,140) ip,ie,in,iu,x(ip),xerr
  140       format(4i4,2f10.4) 
           ip=ip+1
          enddo
        enddo
      enddo
      close(1) 
c
c------termination
c
      stop
      end
c
c====================================================================
c
      subroutine length(xyzsit,xyzsat,xlon1,xlon2,xlat1,xlat2
     1                ,alt1,alt2,aleng)
c
      implicit real*8 (a-h,o-z)
c
      dimension xyzsit(3),xyzsat(3),xyzsip(3),iflg(2,3),xyzpen(3,2)
     1         ,xyz(3)
c
      data erad/6378.d3/
c
      pi=4.0*datan(1.0d0)
      dp=pi/180.0
c
      do k1=1,2
        do k2=1,3
          iflg(k1,k2)=0
        enddo
      enddo
c
      ip=0
c
      xyz(1)=xyzsat(1)-xyzsit(1)
      xyz(2)=xyzsat(2)-xyzsit(2)
      xyz(3)=xyzsat(3)-xyzsit(3)
      sit_height=dsqrt(xyzsit(1)**2 + xyzsit(2)**2 + xyzsit(3)**2)
c     sit_height=erad
c
c-----(1) lower plane
c
      hion=alt1*1.e3
      a=xyz(1)**2+xyz(2)**2+xyz(3)**2
      b=2.0*(xyz(1)*xyzsit(1)+xyz(2)*xyzsit(2)+xyz(3)*xyzsit(3))
      c=sit_height**2-(erad+hion)**2
      epsi=(-b+dsqrt(b**2-4.0*a*c))/2.0/a
      xyzsip(1)=xyzsit(1)+xyz(1)*epsi
      xyzsip(2)=xyzsit(2)+xyz(2)*epsi
      xyzsip(3)=xyzsit(3)+xyz(3)*epsi
c
      alon=datan2(xyzsip(2),xyzsip(1))/dp
      if(alon.lt.0.0) alon=alon+360.0
      alat=datan2(xyzsip(3),dsqrt(xyzsip(1)**2+xyzsip(2)**2))/dp
      if(alon.ge.xlon1.and.alon.lt.xlon2) then
        if(alat.ge.xlat1.and.alat.lt.xlat2) then
          iflg(1,1)=1     ! top surface penetrated
          ip=ip+1
          do k=1,3
            xyzpen(k,ip)=xyzsip(k)
          enddo
        endif
      endif
c
c-----(2) upper plane
c
      hion=alt2*1.e3
      a=xyz(1)**2+xyz(2)**2+xyz(3)**2
      b=2.0*(xyz(1)*xyzsit(1)+xyz(2)*xyzsit(2)+xyz(3)*xyzsit(3))
      c=sit_height**2-(erad+hion)**2
      epsi=(-b+dsqrt(b**2-4.0*a*c))/2.0/a
      xyzsip(1)=xyzsit(1)+xyz(1)*epsi
      xyzsip(2)=xyzsit(2)+xyz(2)*epsi
      xyzsip(3)=xyzsit(3)+xyz(3)*epsi
c
      alon=datan2(xyzsip(2),xyzsip(1))/dp
      if(alon.lt.0.0) alon=alon+360.0
      alat=datan2(xyzsip(3),dsqrt(xyzsip(1)**2+xyzsip(2)**2))/dp
      if(alon.ge.xlon1.and.alon.lt.xlon2) then
        if(alat.ge.xlat1.and.alat.lt.xlat2) then
          iflg(2,1)=1     ! bottom surface penetrated
          ip=ip+1
          do k=1,3
            xyzpen(k,ip)=xyzsip(k)
          enddo
        endif
      endif
c
c-----(3) west plane
c
      eq_lon=xlon1
      a=xyz(1)**2+xyz(2)**2+xyz(3)**2
      distance = dsqrt(a)
      tanphi=dtan(eq_lon*dp)
      b=xyz(1)*tanphi/distance - xyz(2)/distance 
      c= - (tanphi*xyzsit(1)-xyzsit(2))
      epsi=c/b
      xyzsip(1)=xyzsit(1)+xyz(1)*epsi/distance
      xyzsip(2)=xyzsit(2)+xyz(2)*epsi/distance
      xyzsip(3)=xyzsit(3)+xyz(3)*epsi/distance
      gc_height=dsqrt(xyzsip(1)**2 + xyzsip(2)**2 + xyzsip(3)**2)
c
      alat=datan2(xyzsip(3),dsqrt(xyzsip(1)**2+xyzsip(2)**2))/dp
      height=1.e-3*(gc_height - erad)    ! convert to km
      if(alat.ge.xlat1.and.alat.lt.xlat2) then
        if(height.ge.alt1.and.height.lt.alt2) then
          iflg(1,2)=1   ! west surface penetrated
          ip=ip+1
          if(ip.gt.2) goto 5555
          do k=1,3
            xyzpen(k,ip)=xyzsip(k)
          enddo
        endif
      endif
c
c-----(4) east plane
c
      eq_lon=xlon2
      a=xyz(1)**2+xyz(2)**2+xyz(3)**2
      distance = dsqrt(a)
      tanphi=dtan(eq_lon*dp)
      b=xyz(1)*tanphi/distance - xyz(2)/distance 
      c= - (tanphi*xyzsit(1)-xyzsit(2))
      epsi=c/b
      xyzsip(1)=xyzsit(1)+xyz(1)*epsi/distance
      xyzsip(2)=xyzsit(2)+xyz(2)*epsi/distance
      xyzsip(3)=xyzsit(3)+xyz(3)*epsi/distance
      gc_height=dsqrt(xyzsip(1)**2 + xyzsip(2)**2 + xyzsip(3)**2)
c
      alat=datan2(xyzsip(3),dsqrt(xyzsip(1)**2+xyzsip(2)**2))/dp
      height=1.e-3*(gc_height - erad)    ! convert to km
c     write(6,*)1.e-3*gc_height,1.e-3*sit_height,height
      if(alat.ge.xlat1.and.alat.lt.xlat2) then
        if(height.ge.alt1.and.height.lt.alt2) then
          iflg(2,2)=1   ! east surface penetrated
          ip=ip+1
          if(ip.gt.2) goto 5555
          do k=1,3
            xyzpen(k,ip)=xyzsip(k)
          enddo
        endif
      endif
c
c-----(5) south plane
c
      eq_lat=xlat1
      tanphi=dtan(eq_lat*dp)
c
      a=tanphi**2*(xyz(1)**2 + xyz(2)**2) - xyz(3)**2
      b=tanphi**2*(xyz(1)*xyzsit(1) + xyz(2)*xyzsit(2)) 
     1             - xyz(3)*xyzsit(3)
      c=tanphi**2*(xyzsit(1)**2 + xyzsit(2)**2) - xyzsit(3)**2
      epsi=(-b+dsqrt(b**2 - a*c))/a
      xyzsip(3)=xyzsit(3)+xyz(3)*epsi
      if(xyzsip(3)*eq_lat.lt.0.0) then    ! alternative solution is the right one
        epsi=(-b-dsqrt(b**2 - a*c))/a
      endif
c
      xyzsip(1)=xyzsit(1)+xyz(1)*epsi
      xyzsip(2)=xyzsit(2)+xyz(2)*epsi
      xyzsip(3)=xyzsit(3)+xyz(3)*epsi
      gc_height=dsqrt(xyzsip(1)**2 + xyzsip(2)**2 + xyzsip(3)**2)
c
      alon=datan2(xyzsip(2),xyzsip(1))/dp
      alat=datan2(xyzsip(3),dsqrt(xyzsip(1)**2+xyzsip(2)**2))/dp
      if(alon.lt.0.0) alon=alon+360.0
      height=1.e-3*(gc_height - erad)    ! convert to km
      if(alon.ge.xlon1.and.alon.lt.xlon2) then
        if(height.ge.alt1.and.height.lt.alt2) then
          iflg(1,3)=1  ! south surface penetrated
          ip=ip+1
          if(ip.gt.2) goto 5555
          do k=1,3
            xyzpen(k,ip)=xyzsip(k)
          enddo
        endif
      endif
c
c-----(6) north plane
c
      eq_lat=xlat2
      tanphi=dtan(eq_lat*dp)
c
      a=tanphi**2*(xyz(1)**2 + xyz(2)**2) - xyz(3)**2
      b=tanphi**2*(xyz(1)*xyzsit(1) + xyz(2)*xyzsit(2)) 
     1             - xyz(3)*xyzsit(3)
      c=tanphi**2*(xyzsit(1)**2 + xyzsit(2)**2) - xyzsit(3)**2
      epsi=(-b+dsqrt(b**2 - a*c))/a
      xyzsip(3)=xyzsit(3)+xyz(3)*epsi
      if(xyzsip(3)*eq_lat.lt.0.0) then    ! alternative solution is the right one
        epsi=(-b-dsqrt(b**2 - a*c))/a
      endif
      xyzsip(1)=xyzsit(1)+xyz(1)*epsi
      xyzsip(2)=xyzsit(2)+xyz(2)*epsi
      xyzsip(3)=xyzsit(3)+xyz(3)*epsi
      gc_height=dsqrt(xyzsip(1)**2 + xyzsip(2)**2 + xyzsip(3)**2)
c
      alon=datan2(xyzsip(2),xyzsip(1))/dp
      if(alon.lt.0.0) alon=alon+360.0
      height=1.e-3*(gc_height - erad)    ! convert to km
      if(alon.ge.xlon1.and.alon.lt.xlon2) then
        if(height.ge.alt1.and.height.lt.alt2) then
          iflg(2,3)=1  ! north surface penetrated
          ip=ip+1
          if(ip.gt.2) goto 5555
          do k=1,3
            xyzpen(k,ip)=xyzsip(k)
          enddo
        endif
      endif
c
      goto 6666
c
c-----three intersections! take the longest pair
c
 5555 continue
c
      aleng1=0.0
      aleng2=0.0
      aleng3=0.0
      do k=1,3
        aleng1=aleng1 + (xyzpen(k,1)-xyzpen(k,2))**2
        aleng2=aleng2 + (xyzpen(k,1)-xyzsip(k)  )**2
        aleng3=aleng3 + (xyzpen(k,2)-xyzsip(k)  )**2
      enddo
      aleng=dsqrt(aleng1)
      if(aleng.lt.dsqrt(aleng2)) aleng=dsqrt(aleng2)
      if(aleng.lt.dsqrt(aleng3)) aleng=dsqrt(aleng3)
      return
c
c-----(7) calculate penetration length
c
 6666 continue
      aleng=0.0
      if(ip.lt.2) return
      do k=1,3
        aleng=aleng + (xyzpen(k,1)-xyzpen(k,2))**2
      enddo
      aleng=dsqrt(aleng)
c
c-----return
c
      return
      end
c
c====================================================================
c
      subroutine ipp(xyzsit,xyzsat,alt,xlon,xlat)
c
      implicit real*8 (a-h,o-z)
c
      dimension xyzsit(3),xyzsat(3),xyzsip(3),xyz(3)
c
      data erad/6378.d3/
c
      pi=4.0*datan(1.0d0)
      dp=pi/180.0
c
      xyz(1)=xyzsat(1)-xyzsit(1)
      xyz(2)=xyzsat(2)-xyzsit(2)
      xyz(3)=xyzsat(3)-xyzsit(3)
c
      hion=alt*1.e3
      a=xyz(1)**2+xyz(2)**2+xyz(3)**2
      b=2.0*(xyz(1)*xyzsit(1)+xyz(2)*xyzsit(2)+xyz(3)*xyzsit(3))
      c=erad**2-(erad+hion)**2
      epsi=(-b+dsqrt(b**2-4.0*a*c))/2.0/a
      xyzsip(1)=xyzsit(1)+xyz(1)*epsi
      xyzsip(2)=xyzsit(2)+xyz(2)*epsi
      xyzsip(3)=xyzsit(3)+xyz(3)*epsi
c
      xlon=datan2(xyzsip(2),xyzsip(1))/dp
      xlat=datan2(xyzsip(3),dsqrt(xyzsip(1)**2+xyzsip(2)**2))/dp
      if(xlon.lt.0.0) xlon=xlon+360.0
c
      return
      end
c
c=====================================================================
c
      subroutine draw(jcomp,cent,blo_lon,blo_lat,blo_alt,aleng)
c
      implicit real*8 (a-h,o-z)
c
      dimension cent(3)
c
      write(6,'("> -Z",f6.1," km")')1.d-3*aleng
c
      if(jcomp.eq.1) then ! londitudinal wall
        write(6,*) cent(3)-0.5*blo_alt,cent(2)-0.5*blo_lat
        write(6,*) cent(3)+0.5*blo_alt,cent(2)-0.5*blo_lat
        write(6,*) cent(3)+0.5*blo_alt,cent(2)+0.5*blo_lat
        write(6,*) cent(3)-0.5*blo_alt,cent(2)+0.5*blo_lat
        write(6,*) cent(3)-0.5*blo_alt,cent(2)-0.5*blo_lat
      else if(jcomp.eq.2) then ! latitudinal wall
        write(6,*) cent(1)-0.5*blo_lon,cent(3)-0.5*blo_alt
        write(6,*) cent(1)-0.5*blo_lon,cent(3)+0.5*blo_alt
        write(6,*) cent(1)+0.5*blo_lon,cent(3)+0.5*blo_alt
        write(6,*) cent(1)+0.5*blo_lon,cent(3)-0.5*blo_alt
        write(6,*) cent(1)-0.5*blo_lon,cent(3)-0.5*blo_alt
      else if(jcomp.eq.3) then ! horizontal plane
        write(6,*) cent(1)-0.5*blo_lon,cent(2)-0.5*blo_lat
        write(6,*) cent(1)-0.5*blo_lon,cent(2)+0.5*blo_lat
        write(6,*) cent(1)+0.5*blo_lon,cent(2)+0.5*blo_lat
        write(6,*) cent(1)+0.5*blo_lon,cent(2)-0.5*blo_lat
        write(6,*) cent(1)-0.5*blo_lon,cent(2)-0.5*blo_lat
      endif
c
      return
      end
c
c=============================================================
c
      subroutine lsq(b,aty,z,mdim,ndim,x,binv)
c*
      implicit real*8 (a-h,o-z)
c
c declarations
c -------------
c
      dimension b(mdim,mdim),binv(mdim,mdim),aty(mdim),x(mdim)
     1      ,z(mdim)
c
c cholesky separation of normal matrix
c ------------------------------------
c
      call rrsep(b,mdim,ndim,binv)
c
c solve the first equation
c ------------------------
c
      call rsolv2(binv,aty,mdim,ndim,z)
c
c solve the second equation
c -------------------------
c
      call rsolv1(binv,z,mdim,ndim,x)
c
c get inverse matrix of r matrix
c ------------------------------
c
      call rinvs(binv,z,mdim,ndim,b)
c
c square the inverse of r matrix
c ------------------------------
c
      call rsqar(b,mdim,ndim,binv)
c
c return to the main
c ------------------
c
      return
      end
c
c=============================================================
c
      subroutine rrsep(cov,mdim,ndim,r)
c
      implicit real*8 (a-h,o-z)
c
      dimension cov(mdim,mdim),r(mdim,mdim)
c
      r(1,1)=dsqrt(cov(1,1))
      do 1000 k=2,ndim
 1000 r(1,k)=cov(1,k)/r(1,1)
c
      do 1040 k2=2,ndim
c
        do 1010 k3=1,k2-1
 1010   r(k2,k3)=0.0
c
        do 1030 k3=k2,ndim
c
          r(k2,k3)=cov(k2,k3)
          do 1020 k1=1,k2-1
 1020     r(k2,k3)=r(k2,k3)-r(k1,k2)*r(k1,k3)
c
          if(k3.eq.k2) then
            r(k2,k3)=sqrt(r(k2,k3))
          else
            r(k2,k3)=r(k2,k3)/r(k2,k2)
          endif
c
 1030   continue
c
 1040 continue
c
      return
      end
c
c=============================================================
c
      subroutine rsolv2(r,scr1,mdim,ndim,scr2)
c
      implicit real*8 (a-h,o-z)
c
      dimension r(mdim,mdim),scr1(mdim),scr2(mdim)
c
      scr2(1)=scr1(1)/r(1,1)
      do 1010 k2=2,ndim
        scr2(k2)=scr1(k2)
        do 1000 k1=1,k2-1
 1000   scr2(k2)=scr2(k2)-r(k1,k2)*scr2(k1)
        scr2(k2)=scr2(k2)/r(k2,k2)
 1010 continue
c
      return
      end
c
c=============================================================
c
      subroutine rsolv1(r,scr1,mdim,ndim,scr2)
c
      implicit real*8 (a-h,o-z)
c
      dimension r(mdim,mdim),scr1(mdim),scr2(mdim)
c
      scr2(ndim)=scr1(ndim)/r(ndim,ndim)
      do 1010 k2=ndim-1,1,-1
        scr2(k2)=scr1(k2)
        do 1000 k1=k2+1,ndim
 1000   scr2(k2)=scr2(k2)-r(k2,k1)*scr2(k1)
        scr2(k2)=scr2(k2)/r(k2,k2)
 1010 continue
c
      return
      end
c
c=============================================================
c
      subroutine rinvs(r,scr,mdim,ndim,rinv)
c
      implicit real*8 (a-h,o-z)
c
      dimension r(mdim,mdim),rinv(mdim,mdim),scr(mdim)
c
      do 1010 i=1,ndim
        do 1000 k=1,ndim
 1000   scr(k)=0.0
        scr(i)=1.0
c
        call rsolv1(r,scr,mdim,ndim,rinv(1,i))
c
 1010 continue
c
      return
      end
c
c=============================================================
c
      subroutine rsqar(r,mdim,ndim,wt)
c
      implicit real*8 (a-h,o-z)
c
      dimension r(mdim,mdim),wt(mdim,mdim)
c
      do 1030 k2=1,ndim
        do 1020 k1=k2,ndim
          wt(k1,k2)=0.0
          do 1010 k=k1,ndim
 1010     wt(k1,k2)=wt(k1,k2)+r(k1,k)*r(k2,k)
          if(k1.ne.k2) wt(k2,k1)=wt(k1,k2)
 1020   continue
 1030 continue
c
      return
      end
c============================================================
c
      subroutine chapman(altitude,vtec,density)
c
      implicit real*8(a-h,o-z)
c
      data h_max/300.0/
      data h_const/80.0/
c
c-----ionospheric profile
c
      tec=0.0
      do ih=0,1000,10
        xi=(real(ih)-h_max)/h_const
        dens_ion=exp((1.0 - xi - exp(-xi))/2.0)
        tec=tec+dens_ion*10.0d3
      enddo
c
      factor=vtec*1.d16/tec
c
      xi=(altitude-h_max)/h_const
      dens_ion=exp((1.0 - xi - exp(-xi))/2.0)
      density = 1.d-11*factor*dens_ion
c
c------end
c
      return
      end

