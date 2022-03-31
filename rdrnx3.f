c
      program rdrnx3     ! for RINEX version 3.02
c
c     a program to read rinex file and output phase difference between
c     L1 and L2 for particular satellite
c
      parameter (maxsat=70, maxepch=17280,maxslip=30,maxtype=27)    ! 5 sec sampling
      parameter (maxglns=30)
      implicit real*4 (a-h,o-z)
c
      real*8 dphase(2,maxsat),dcode(2,maxsat),bias(maxsat)
     1      ,pdif(maxepch,maxsat)
     1      ,vel_light,f_l1,f_l2,array(2,maxepch),pval
     2      ,sum,cdif(maxepch,maxsat),ddata(maxtype)
      dimension jsat(maxsat),time(maxepch),jslip(maxslip)
     1         ,numsat(maxsat),itable(maxsat),jleng(maxslip+1)
     2         ,numepch(maxsat),mfreq(maxglns)
     3         ,fact_meter2tec_g(maxglns),lli(maxsat),lbad(maxsat)
      character lhead*60,lmark*4,lscr*1, lbuf*300, lcmt*20,lfile*64
     1         ,ltype(maxtype)*3,ltypsat(maxsat)*1,ltime1*16,ltime2*16
     2         ,lsys*1
      logical lfirst,lsameday,lonemore
      data vel_light/299792458d0/
      data f_l1_gps/1575420000.d0/
      data f_l2_gps/1227600000.d0/
      data f_l1_galileo/1575420000.d0/
      data f_l2_galileo/1176450000.d0/
      data f_l1_beidou/1561098000.d0/
      data f_l2_beidou/1207140000.d0/
      data f_l1_glonass0/1602000000.d0/     ! central frequency (k=0)
      data f_l2_glonass0/1246000000.d0/     ! central frequency (k=0)
      data f_l1_step/0.5625d6/
      data f_l2_step/0.4375d6/
      data mfreq/1, -4,  5,  6, 1, -4, 5,  0, -2, -7,
     1           0, -1, -2, -7, 0, -1, 4, -3,  3,  2,
     2           4, -3,  3,  2, 0,  0, 0,  0,  0,  0/
      data fac_m3_over_s2/40.308/
      np_arc_thres = 120       ! minimum number for the second longest arc
c 
c  (1) GPS/QZSS
      fact_meter2tec  = (-1.0e-16/fac_m3_over_s2)
     1               *(f_l1_gps**2)*(f_l2_gps**2)
     2               /(f_l1_gps**2 - f_l2_gps**2)
c
c  (2) GLONASS
      do iglns=1,maxglns
        f_l2_glonass = f_l2_glonass0 + mfreq(iglns)*f_l2_step
        f_l1_glonass = f_l1_glonass0 + mfreq(iglns)*f_l1_step
        fact_meter2tec_g(iglns)= (-1.0e-16/fac_m3_over_s2)
     1     *(f_l1_glonass**2)*(f_l2_glonass**2)
     2     /(f_l1_glonass**2 - f_l2_glonass**2)
      enddo
c
c  (3) Galileo
      fact_meter2tec_galileo  = (-1.0e-16/fac_m3_over_s2)
     1               *(f_l1_galileo**2)*(f_l2_galileo**2)
     2               /(f_l1_galileo**2 - f_l2_galileo**2)
c
c  (4) Beidou
      fact_meter2tec_beidou  = (-1.0e-16/fac_m3_over_s2)
     1               *(f_l1_beidou**2)*(f_l2_beidou**2)
     2               /(f_l1_beidou**2 - f_l2_beidou**2)
c
c
c     fact_meter2tec= 6.05/(-0.65)
c     fact_meter2tec= -9.5195857  ! an updated value from Liming
c
      call getarg(1,lfile)
      call getarg(2,lsys)
      call getarg(3,ltime1)
      call getarg(4,ltime2)
c
      tstart=0.0
      tend= 24.0
      read(ltime1,*,end=999) tstart
      read(ltime2,*) tend
  999 continue
c
      lfirst=.true.
      lsameday=.true.
      do isat=1,maxsat
        bias(isat)=0.0d0
      enddo
c
      do isat=1,maxsat
       do iepch=1,maxepch
         pdif(iepch,isat)=0.0d0
       enddo
      enddo
c
c-----reading header
c
      jl1=0
      jp1=0
      jl2=0
      jp2=0
      open(1,file=lfile)
 1111 continue
      read(1,100)lhead,lcmt
  100 format(a60,a20)
c     write(6,100)lhead,lcmt
      if(lcmt.eq.'SYS / # / OBS TYPES '.and.
     1   lhead(1:1).eq.lsys) then
        backspace 1
        read(1,110)ntype,(ltype(itype),itype=1,13)
  110   format(1x,i5,13(1x,a3))
        if(ntype.gt.13) then
          read(1,120)(ltype(itype),itype=14,ntype)
  120     format(6x,13(1x,a3))
        endif
        do itype=1,ntype
          if(lsys.eq.'G') then    ! GPS
            if(ltype(itype).eq.'L1C') jl1=itype
            if(ltype(itype).eq.'L2W') jl2=itype
            if(ltype(itype).eq.'C1C') jp1=itype
            if(ltype(itype).eq.'C2W') jp2=itype
          else if(lsys.eq.'R') then    ! GLONASS
            if(ltype(itype).eq.'L1P') jl1=itype
            if(ltype(itype).eq.'L2P') jl2=itype
            if(ltype(itype).eq.'C1P') jp1=itype
            if(ltype(itype).eq.'C2P') jp2=itype
          else if(lsys.eq.'J') then    ! QZSS
            if(ltype(itype).eq.'L1X') jl1=itype
            if(ltype(itype).eq.'L2X') jl2=itype
            if(ltype(itype).eq.'C1X') jp1=itype
            if(ltype(itype).eq.'C2X') jp2=itype
            if(ltype(itype).eq.'L1C') jl1ca=itype
            if(ltype(itype).eq.'C1C') jp1ca=itype
            if(ltype(itype).eq.'L2L') jl2ca=itype    ! added to read TONGA data
            if(ltype(itype).eq.'C2L') jp2ca=itype    ! added to read TONGA data
          else if(lsys.eq.'E') then    ! Galileo
            if(ltype(itype).eq.'L1X') jl1=itype
            if(ltype(itype).eq.'L5X') jl2=itype
            if(ltype(itype).eq.'C1X') jp1=itype
            if(ltype(itype).eq.'C5X') jp2=itype
            if(ltype(itype).eq.'L1C') jl1alt=itype    ! added to read TONGA data
            if(ltype(itype).eq.'L5Q') jl2alt=itype    ! added to read TONGA data
            if(ltype(itype).eq.'C1C') jp1alt=itype    ! added to read TONGA data
            if(ltype(itype).eq.'C5Q') jp2alt=itype    ! added to read TONGA data
          else if(lsys.eq.'C') then    ! Beidou
            if(ltype(itype).eq.'L2I') jl1=itype
            if(ltype(itype).eq.'L7I') jl2=itype
            if(ltype(itype).eq.'C2I') jp1=itype
            if(ltype(itype).eq.'C7I') jp2=itype
          endif
        enddo
        if(lsys.eq.'J'.and.jl1.eq.0) jl1=jl1ca
        if(lsys.eq.'J'.and.jp1.eq.0) jp1=jp1ca
        if(lsys.eq.'J'.and.jl2.eq.0) jl2=jl2ca    ! added to read TONGA data
        if(lsys.eq.'J'.and.jp2.eq.0) jp2=jp2ca    ! added to read TONGA data
        if(lsys.eq.'E'.and.jl1.eq.0) jl1=jl1alt    ! added to read TONGA data
        if(lsys.eq.'E'.and.jl2.eq.0) jl2=jl2alt    ! added to read TONGA data
        if(lsys.eq.'E'.and.jp1.eq.0) jp1=jp1alt    ! added to read TONGA data
        if(lsys.eq.'E'.and.jp2.eq.0) jp2=jp2alt    ! added to read TONGA data
c       write(6,*)ntype,jl1,jl2,jp1,jp2
      endif
c
      if(lcmt(1:11).eq.'MARKER NAME') lmark=lhead(1:4)
      if(lcmt(1:13).ne.'END OF HEADER') goto 1111
c
c-----reading data
c
      iepch=1
 2222 read(1,130,end=9999)lbuf          ! line starting with '>'
  130 format(a300)
      if(lbuf(3:3).eq.' ') goto 9999
      read(lbuf,140)ih,im,sec,nsat
  140 format(12x,2i3,f11.7,3x,i3)
      time(iepch)=real(ih)+real(im)/60.0 + sec/3600.0 
c     write(6,*) iepch,ih,im,sec,nsat,time(iepch)
c
      do 1000 isat=1,nsat
        read(1,130)lbuf
c       write(6,'(a80)') lbuf(1:80)
c       write(6,'(a300)') lbuf
c       if(isat.eq.1) stop
        read(lbuf,150) ltypsat(isat),numsat(isat)
     1                  ,(ddata(k),lli(k),k=1,ntype)
  150   format(a1,i2,27(f14.3,i1,1x))
        lbad(isat)=0
c       if(lli(jl1).ne.0.or.lli(jl2).ne.0) then   ! flag="4" is all right with GPS, so only flag="1" is flagged
        if(lli(jl1).eq.1.or.lli(jl2).eq.1) then
          lbad(isat)=1
          goto 1000   ! loss-of-lock indicator non zero
        endif
        dphase(1,isat)=ddata(jl1)
        dphase(2,isat)=ddata(jl2)
        dcode(1,isat)=ddata(jp1)
        dcode(2,isat)=ddata(jp2)
c       write(6,*) iepch,isat,ltypsat(isat),numsat(isat)
c    1           ,(dphase(k,isat),k=1,2),(dcode(k,isat),k=1,2)
 1000 continue
c
      if(lfirst) then
        lfirst=.false.
        ksat=0
        do isat=1,nsat
          itable(isat)=0
          if(ltypsat(isat).eq.lsys) then
            ksat=ksat+1
            jsat(ksat)=numsat(isat)
            itable(isat)=ksat
c           write(6,*) isat,ltypsat(isat),numsat(isat),ksat
          endif
        enddo
        msat=ksat
      else
        do 1010 isat=1,nsat
          if(ltypsat(isat).ne.lsys) then     ! pickup GPS
            itable(isat)=0
            goto 1010
          endif
          do k=1,msat
            if(numsat(isat).eq.jsat(k)) then
              itable(isat)=k
              goto 1010
            endif
          enddo
          msat=msat+1
c         write(6,*)msat,' new SV:',numsat(isat)
          jsat(msat)=numsat(isat)
          itable(isat)=msat
 1010   continue
c
        if(lsameday) then
          if(time(iepch).lt.timeb4) then
            lsameday=.false.
            time(iepch)=time(iepch)+24.0
          endif
        else
          time(iepch)=time(iepch)+24.0
        endif
      endif
c
      do 1020 isat=1,nsat
c
        idsat=itable(isat)
        if(idsat.eq.0) goto 1020      ! pickup only specified GNSS
        if(lbad(isat).ne.0) goto 1020 ! loss-of-lock indicator is non-zero
c
        if(dphase(1,isat).eq.0.0.or.dphase(2,isat).eq.0.0) goto 1020
        if(dcode(1,isat).eq.0.0.or.dcode(2,isat).eq.0.0) goto 1020  
c       write(6,*)isat,itable(isat),(dphase(k,isat),k=1,2)
c    1                              ,(dcode(k,isat),k=1,2)
        if(lsys.eq.'G'.or.lsys.eq.'J') then
          f_l2=f_l2_gps
          f_l1=f_l1_gps
        else if(lsys.eq.'R') then
          iglns=jsat(idsat)
          f_l2=f_l2_glonass0 + mfreq(iglns)*f_l2_step
          f_l1=f_l1_glonass0 + mfreq(iglns)*f_l1_step
        else if(lsys.eq.'E') then
          f_l2=f_l2_galileo
          f_l1=f_l1_galileo
        else if(lsys.eq.'C') then
          f_l2=f_l2_beidou
          f_l1=f_l1_beidou
        endif
c
        if(bias(idsat).eq.0.0d0) then
            bias(idsat)=
     1      (vel_light/f_l2)*dphase(2,isat)
     2     -(vel_light/f_l1)*dphase(1,isat)
            pdif(iepch,idsat)=0.00001
        else
            pdif(iepch,idsat)=
     1      (vel_light/f_l2)*dphase(2,isat)
     1     -(vel_light/f_l1)*dphase(1,isat) - bias(idsat)
        endif
c
        cdif(iepch,idsat)=dcode(1,isat)-dcode(2,isat)
c
c       write(6,*)ih,im,jsat(idsat),pdif(iepch,idsat),cdif(iepch,idsat)
 1020 continue
      timeb4=time(iepch)
      iepch=iepch+1
      goto 2222
c
c-----end of file
c
 9999 nepch=iepch-1
      close(1)
c
c     write(6,'(a1,2i5)') lsys,nepch, msat
c
c-----data screening
c
      do 2000 isat=1,msat
c
c     write(6,'(i2," of",i3)') isat,msat
c     
      lonemore=.false.
c
c (1) time window
c
       jepch=1
       do iepch=1,nepch
        if(time(iepch).ge.tstart.and.time(iepch).lt.tend) then
c         write(6,*)time(iepch),tstart,tend
c         write(6,*)isat,iepch,time(iepch),pdif(iepch,isat)
          if(pdif(iepch,isat).ne.0.0) then
            array(1,jepch)=pdif(iepch,isat)
            array(2,jepch)=cdif(iepch,isat)
            jepch=jepch+1
          endif
        endif
       enddo
       mepch=jepch-1
c      write(6,*)isat,nepch,mepch
c
c (2) cycle-slip detection
c
      islip=1
      do iepch=1,mepch
        if(iepch.eq.1) then
          pval=array(1,iepch)
        else
c         write(6,'(i5,2f8.2)') iepch,pval,array(1,iepch)-pval
c         if(dabs(array(1,iepch)-pval).gt.50.0d0) then
          if(dabs(array(1,iepch)-pval).gt.4.0d0) then     ! try this if no data survived
c         if(dabs(array(1,iepch)-pval).gt.5.0d-1) then    ! under EIA, STEC may change 1 TECU in 30 sec, so be careful
            if(islip.le.maxslip) jslip(islip)=iepch
            islip=islip+1
          endif
          pval=array(1,iepch)
        endif       
      enddo
      nslip=islip-1
      if(nslip.gt.maxslip) goto 2000
c
c  (3) longest arc detection
c
      if(nslip.eq.0) then
        istart=1
        iend=mepch
      else if(nslip.eq.1) then
        ileng1=jslip(1)-1
        ileng2=mepch-jslip(1)+1
        if(ileng1.gt.ileng2) then
          istart=1
          iend=jslip(1)-1
          if(ileng2.gt.np_arc_thres) then
            lonemore=.true.
            istart2=jslip(1)
            iend2=mepch
          endif
        else
          istart=jslip(1)
          iend=mepch
          if(ileng1.ge.np_arc_thres) then
            lonemore=.true.
            istart2=1
            iend2=jslip(1)-1
          endif
        endif
      else
        jleng(1)=jslip(1)
        istart=1
        iend=jslip(1)-1
        maxleng=jleng(1)
        do islip=1,nslip-1
          jleng(islip+1)=jslip(islip+1)-jslip(islip)
          if(jleng(islip+1).gt.maxleng) then
            istart=jslip(islip)
            iend=jslip(islip+1)-1
            maxleng=jleng(islip+1)
          endif
        enddo
        jleng(nslip+1)=mepch-jslip(nslip)+1
        if(jleng(nslip+1).gt.maxleng) then
            istart=jslip(nslip)
            iend=mepch
        endif
      endif
c
      numepch(isat)=iend-istart+1
c
      sum=0.0
      do iepch=istart,iend
        sum=sum+(array(2,iepch) - array(1,iepch))
      enddo
      npoint=iend-istart+1
      if(npoint.ne.0) sum=sum/dble(npoint)
c
      if(lonemore) then
        numepch2=iend2-istart2+1
c       write(6,*) 'One more!',numepch2,istart,iend,istart2,iend2
        sum2=0.0
        do iepch=istart2,iend2
          sum2=sum2+(array(2,iepch) - array(1,iepch))
        enddo
        npoint=iend2-istart2+1
        if(npoint.ne.0) sum2=sum2/dble(npoint)
      endif
c     write(6,*) istart,iend,istart2,iend2
c
c-----output data
c
       if(numepch(isat).le.10) goto 2000
c
       write(6,160) jsat(isat),numepch(isat),lmark
  160  format('> sat#',i2, ' #data:',i4, ' site:',a4)
       if(lsys.eq.'G'.or.lsys.eq.'J') then
         fact=fact_meter2tec
       else if(lsys.eq.'R') then
         iglns=jsat(isat)
         fact=fact_meter2tec_g(iglns)
       else if(lsys.eq.'E') then
         fact=fact_meter2tec_galileo
       else if(lsys.eq.'C') then
         fact=fact_meter2tec_beidou
       endif
c
       kepch=0
       do iepch=1,nepch
        if(time(iepch).ge.tstart.and.time(iepch).lt.tend) then
         if(pdif(iepch,isat).ne.0.0) then
          kepch=kepch+1
          if(kepch.ge.istart.and.kepch.le.iend) then
           write(6,170)time(iepch),fact*(pdif(iepch,isat)+sum)
c          write(6,*)time(iepch),fact*(pdif(iepch,isat)+sum)
  170      format(f8.4,f12.4)
          endif
         endif
        endif
       enddo
c
       if(.not.lonemore) goto 2000
c
c-----second longest arc
c
      if(numepch2.le.10) goto 2000
      kepch=0
      lfirst=.true.
      do iepch=1,nepch
        if(time(iepch).ge.tstart.and.time(iepch).lt.tend) then
          if(pdif(iepch,isat).ne.0.0) then
           kepch=kepch+1
           if(kepch.ge.istart2.and.kepch.le.iend2) then
            if(lfirst) then
              write(6,160) jsat(isat),numepch2,lmark
              lfirst=.false.
            endif
            write(6,170)time(iepch),fact*(pdif(iepch,isat)+sum2)
           endif
          endif
        endif
      enddo
c
 2000 continue
c
      stop
      end
c
