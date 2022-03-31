c
      program rdrnx 
c
c     a program to read rinex file and output phase difference between
c     L1 and L2 for different satellite
c
      parameter (maxsat=70, maxepch=17280,maxslip=10,maxtype=9)    ! 5 sec sampling
      implicit real*4 (a-h,o-z)
c
      real*8 dphase(2),dcode(2),pdif(maxepch,maxsat)
     1      ,vel_light,f_l1,f_l2,array(2,maxepch),pval
     2      ,sum,cdif(maxepch,maxsat),ddata(maxtype)
      dimension jsat(maxsat),time(maxepch),jslip(maxslip)
     1         ,numsat(maxsat),itable(maxsat),jleng(maxslip+1)
     2         ,numepch(maxsat)
      character lhead*60,lmark*4,lscr*1, lbuf*80, lcmt*20,lfile*64
     1         ,ltype(maxtype)*2,ltypsat(maxsat)*1
      logical lfirst,lsameday,lonemore
      data vel_light/299792458d0/
      data f_l1/1575420000.d0/
      data f_l2/1227600000.d0/
c 
      fact_meter2tec= 6.05/(-0.65)
c
      lfirst=.true.
      lsameday=.true.
c
      do isat=1,maxsat
       do iepch=1,maxepch
         pdif(iepch,isat)=0.0d0
       enddo
      enddo
c
      tstart=0.0
      tend=24.0
c
c-----reading header
c
      jc1=0
      jp1=0
 1111 continue
      read(5,100)lhead,lcmt
  100 format(a60,a20)
      if(lcmt(1:13).eq.'# / TYPES OF ') then
        backspace 5
        read(5,105)ntype,(ltype(itype),itype=1,ntype)
  105   format(i6,9(4x,a2))
        do itype=1,ntype
          if(ltype(itype).eq.'L1') jl1=itype
          if(ltype(itype).eq.'L2') jl2=itype
          if(ltype(itype).eq.'C1') jc1=itype
          if(ltype(itype).eq.'P1') jp1=itype
          if(ltype(itype).eq.'P2') jp2=itype
        enddo
      endif
c
      if(jp1.eq.0.and.jc1.ne.0) jp1=jc1
c
      if(lcmt(1:11).eq.'MARKER NAME') lmark=lhead(1:4)
      if(lcmt(1:13).ne.'END OF HEADER') goto 1111
c
c-----reading data
c
      iepch=1
 2222 read(5,110,end=9999)lbuf
  110 format(a80)
      if(lbuf(29:29).eq.'4') then       ! comment
        read(lbuf,120)ncomment
  120   format(29x,i3)
        do icomment=1,ncomment
c         read(1,110)lbuf
          read(5,110)lbuf
        enddo
        goto 2222
      else                             ! regular record
        read(lbuf,130)ih,im,sec,nsat,(ltypsat(k),numsat(k),k=1,12)
  130   format(9x,2i3,f11.7,3x,i3,12(a1,i2))
        if(nsat.gt.12) then
          read(5,110)lbuf
          read(lbuf,135)(ltypsat(k),numsat(k),k=13,nsat)
  135     format(32x,12(a1,i2))
        endif
        time(iepch)=real(ih)+real(im)/60.0 + sec/3600.0 
        if(lfirst) then
          lfirst=.false.
          ksat=0
          do isat=1,nsat
            if(ltypsat(isat).eq.'G') then
              ksat=ksat+1
              jsat(ksat)=numsat(isat)
              itable(isat)=isat
            else      ! ignore GLONASS
              itable(isat)=0
            endif
          enddo
          msat=ksat
        else
          do 3333 isat=1,nsat
            if(ltypsat(isat).ne.'G') then     ! ignore GLONASS
              itable(isat)=0
              goto 3333
            endif
            do k=1,msat
              if(numsat(isat).eq.jsat(k)) then
                itable(isat)=k
                goto 3333
              endif
            enddo
            msat=msat+1
            jsat(msat)=numsat(isat)
            itable(isat)=msat
 3333     continue
c
          if(lsameday) then
            if(time(iepch).lt.timeb4) then
              lsameday=.false.
              time(iepch)=time(iepch)+24.0
            endif
          else
            time(iepch)=time(iepch)+24.0
          endif
c
        endif
        do 1000 isat=1,nsat
          read(5,110)lbuf
          kazu=5
          if(ntype.lt.5) kazu=ntype
          read(lbuf,140,err=7777,end=7777)(ddata(k),k=1,kazu)
  140     format(5(f14.3,2x))
          if(ntype.gt.5) then         ! continuation line for the case ntype>5
             read(5,110)lbuf
             read(lbuf,140,err=7777,end=7777)(ddata(k),k=6,ntype)
          endif
          dphase(1)=ddata(jl1)
          dphase(2)=ddata(jl2)
          dcode(1)=ddata(jp1)
          dcode(2)=ddata(jp2)
          if(dcode(1).eq.0.0) dcode(1)=ddata(jc1)
          if(dphase(1).eq.0.0.or.dphase(2).eq.0.0) goto 7777
          if(dcode(1).eq.0.0.or.dcode(2).eq.0.0) goto 7777
          goto 8888
 7777     continue
          goto 1000
 8888     continue 
          idsat=itable(isat)
          if(idsat.eq.0) goto 1000     ! ignore GLONASS
          pdif(iepch,idsat)=
     1      (vel_light/f_l2)*dphase(2)-(vel_light/f_l1)*dphase(1)
          cdif(iepch,idsat)=dcode(1)-dcode(2)
 1000   continue
        timeb4=time(iepch)
        iepch=iepch+1
        goto 2222
      endif
c
c-----end of file
c
 9999 nepch=iepch-1
c
c-----data screening
c
      do 2000 isat=1,msat
c
      lonemore=.false.
c
c (1) time window
c
       jepch=1
       do iepch=1,nepch
        if(time(iepch).gt.tstart.and.time(iepch).lt.tend) then
          if(pdif(iepch,isat).ne.0.0) then
            array(1,jepch)=pdif(iepch,isat)
            array(2,jepch)=cdif(iepch,isat)
            jepch=jepch+1
          endif
        endif
       enddo
       mepch=jepch-1
c
c (2) cycle-slip detection
c
      islip=1
      do iepch=1,mepch
        if(iepch.eq.1) then
          pval=array(1,iepch)
        else
          if(dabs(array(1,iepch)-pval).gt.5.0d-1) then
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
          if(ileng2.ge.120) then
            lonemore=.true.
            istart2=jslip(1)
            iend2=mepch
          endif
        else
          istart=jslip(1)
          iend=mepch
          if(ileng1.ge.120) then
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
      if(lonemore) numepch2=iend2-istart2+1
c
c-----output data
c
       if(numepch(isat).le.10) goto 2000
c
       kepch=0
       lfirst=.true.
       do iepch=1,nepch
        if(time(iepch).gt.tstart.and.time(iepch).lt.tend) then
         if(pdif(iepch,isat).ne.0.0) then
          kepch=kepch+1
          pval=pdif(iepch,isat)
          if(kepch.ge.istart.and.kepch.le.iend) then
           if(lfirst) then
            bias=fact_meter2tec*pdif(iepch,isat)
            write(6,150) jsat(isat),numepch(isat),lmark,bias
  150       format('> sat#',i2,' #data:',i4,' site:',a4,' bias:',f15.4)
            lfirst=.false.
           endif
           write(6,160)time(iepch),fact_meter2tec*pdif(iepch,isat)-bias
  160      format(f7.3,f12.4)
          endif
         endif
        endif
       enddo
c
       if(.not.lonemore) goto 2000
c
c-----output second longest arc
c
       if(numepch2.le.10) goto 2000
c
       kepch=0
       lfirst=.true.
       do iepch=1,nepch
        if(time(iepch).gt.tstart.and.time(iepch).lt.tend) then
         if(pdif(iepch,isat).ne.0.0) then
          kepch=kepch+1
          pval=pdif(iepch,isat)
          if(kepch.ge.istart2.and.kepch.le.iend2) then
           if(lfirst) then
              bias=fact_meter2tec*pdif(iepch,isat)
              write(6,150) jsat(isat),numepch2,lmark,bias
              lfirst=.false.
           endif
           write(6,160)time(iepch),fact_meter2tec*pdif(iepch,isat)-bias
          endif
         endif
        endif
       enddo
c
c
 2000 continue
c
      stop
      end
