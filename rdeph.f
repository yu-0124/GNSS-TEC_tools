c
      program rdeph
c
c     a program to read rinex ephsmeris file and satellite position in Earth-fixed
c     coordinates
c
c     1:IODE  2:Crs  3:delta-n  4:m0
c     5:Cuc   6:e    7:Cus      8:root-a
c     9:Toe  10:Cic 11:Omega   12:Cis
c    13:i0   14:Crc 15:omega   16:OmegaDot
c    17:iDot 18-28: not used
c
      parameter (maxsat=40, maxepch=30, maxrec=400)
      implicit real*4 (a-h,o-z)
c
      real*8 dxyz,oel
      dimension oel(28,maxrec),idsat(maxrec),tiempo(maxrec)
     1         ,i1stsat(maxsat),dxyz(3)
      character lscr*1, lbuf*80, lcmt*20,lfile*64
c 
      tstep=0.05        ! output every 3 minutes
c
      tstart=0.0
      tend=24.0
c
c-----reading header
c
 1111 continue
      read(5,100)lcmt
  100 format(60x,a20)
c
      if(lcmt(1:13).ne.'END OF HEADER') goto 1111
c
c-----reading data
c
      krec=1
 2222 read(5,110,end=3333)lbuf
  110 format(a80)
      read(lbuf,130)idsat(krec),iy,imon,iday,ih,imin,sec
  130 format(i2,5i3,f5.1)
      tiempo(krec)=real(ih)+real(imin)/60.0 + sec/3600.0 
      do irec=1,7
        ii=(irec-1)*4
        if(irec.ne.7) then
         read(5,140)(oel(k,krec),k=ii+1,ii+4)
        else
         read(5,140)(oel(k,krec),k=ii+1,ii+1)
         do k=ii+2,ii+4
           oel(k,krec)=0.0
         enddo
        endif
  140 format(3x,4(d19.12))
      enddo
      krec=krec+1
      goto 2222
 3333 close(1)
      mrec=krec-1 
c      
c-----finding 1st appearence of satellites
c
      do isat=1,maxsat
        i1stsat(isat)=0
      enddo
c
      do krec=1,mrec
        if(i1stsat(idsat(krec)).eq.0) i1stsat(idsat(krec))=krec
      enddo
c
c-----calculating satellite positions
c
      time=0.0
c
 4444 continue
c
      if(time.ge.tstart) then
        do isat=1,maxsat
          jrec=i1stsat(isat)
          if(jrec.ne.0) then
           time0=tiempo(jrec)
           call gtxyz(time,time0,oel(1,jrec),dxyz)
           write(6,150)time,isat,(dxyz(k),k=1,3)
  150      format(f8.2,i3,3(1x,d17.11))
          endif
        enddo
      endif
c
      time=time+tstep
      if(time.gt.tend) goto 9999
      goto 4444
c      
 9999 continue
c
      stop
      end
c
c------------------------------------------------------------------------
c
      subroutine gtxyz(time,time0,ele,dxyz)
c
      implicit real*8 (a-h,o-z)
      real*4 time,time0
c
      dimension dxyz(3),ele(28)
      data GM/3986005.d8/
      data omega_dot_e/7292115.d-11/
c 
c     1:IODE  2:Crs  3:delta-n  4:m0
c     5:Cuc   6:e    7:Cus      8:root-a
c     9:Toe  10:Cic 11:Omega   12:Cis
c    13:i0   14:Crc 15:omega   16:OmegaDot
c    17:iDot 18-28: not used
c
c (1) 
      a=ele(8)**2
c
c (2) 
      dnzero=dsqrt(GM/a**3)
c
c (3)
      tk=dble(time-time0)*60.*60.
c
c (4)
      dn=dnzero+ele(3)
c
c (5)
      dmk=ele(4)+dn*tk
c
c (6)
      call kepler(dmk,ele(6),ek)
c
c (7)
      cosvk=(dcos(ek)-ele(6))/(1.0-ele(6)*dcos(ek))
      sinvk=dsqrt(1.0-ele(6)**2)*dsin(ek)/(1.0-ele(6)*dcos(ek))
      vk=datan2(sinvk,cosvk)
c
c (8)
      phik=vk+ele(15)
c
c (9)
      delta_uk=ele(7)*dsin(2.0*phik)+ele(5)*dcos(2.0*phik)
      uk=phik+delta_uk
c
c (10)
      delta_rk=ele(2)*dsin(2.0*phik)+ele(14)*dcos(2.0*phik)
      rk=a*(1.0-ele(6)*dcos(ek))+delta_rk
c
c (11)
      delta_dik=ele(12)*dsin(2.0*phik)+ele(10)*dcos(2.0*phik)
      dik=ele(13)+delta_dik+ele(17)*tk
c
c (12)
      xdashk=rk*dcos(uk)
      ydashk=rk*dsin(uk)
c
c (13)
      omegak=ele(11)+(ele(16)-omega_dot_e)*tk-ele(9)*omega_dot_e
c
c (14)
      dxyz(1)=xdashk*dcos(omegak)-ydashk*dcos(dik)*dsin(omegak)
      dxyz(2)=xdashk*dsin(omegak)+ydashk*dcos(dik)*dcos(omegak)
      dxyz(3)=ydashk*dsin(dik)
c
      return
      end
c
c------------------------------------------------------------------------
c
      subroutine kepler(dmk,e,ek)
c
      implicit real*8 (a-h,o-z)
      logical lfirst
      data thres/1.e-14/
c
      ek=dmk
      lfirst=.true.
      i=1
c
 1000 continue
c
      diff=ek+e*dsin(ek)-dmk
      if(dabs(diff).lt.thres) goto 9999
      partial=1-e*dcos(ek)
      ek=ek-diff/partial
      i=i+1
      goto 1000
c
 9999 niteration=i-1
c
      return
      end

