cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! BMD PWdW 2015/20016 fABMACS
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine hellof(istep,xx,ff)
      implicit integer*4 (i-n)
      implicit real*8 (a-h,o-z)      
      real*4 bigd(BMAX,BMAX)
      real*8 dpop(2,BMAX,BMAX),pop(BMAX,BMAX)
      real*8 dela(3,BMAX,BMAX),ddela(BMAX,BMAX)
      real*8 decon(BMAX,BMAX),  widths(3)
      real*8 x0(NPARTS*3),ocog(3) !aly
      real*8 xwrp(NPARTS*3)       !for wrapping it
      real*8 alp,bee,cee,width,dih,DelT
      real*8 bolt,dx,radius,scal,omega
      real*8 point(2,3),ormal(3)
      integer look(BMAX,BMAX),ialp
      integer nbin,imabp,np,np1,np2
      save !and save the stuff

      real ff(NPARTS*3),xx(NPARTS*3)
      real*8 widths24(3), widths2(3)
      real*8 jaco(2,NPARTS*3),cog(3) !with ALY      !aly
      real*8 deid(3,4),qa(12)!for dihedral restraint
      real*8 ave1(2),freest
      real*8 orml,sumd,at1,at2,bang,arg


      if(istep.eq.0)then
         open(13,file="sphpoints")!replaced width/2 as sphere center
         do i=1,1
            do j=1,3
               read(13,*) point(i,j)
            enddo
         enddo
         
         bah=1d0
         bah2 = 1d0
         bah3=0d0
         nbin = BMAX !number of bins 
         np=NPARTS
         open(13,file='Reference')
         do i=1,np
            read(13,*) x0(i*3-(3-1)), x0(i*3-(3-2)), x0(i*3-(3-3))
         enddo
         close(13)
c         do i=1,np!116
c            do j=1,3
c               x0(i*3-(3-j))=x0(i*3-(3-j))/10d0!
c            enddo
c         enddo

         np1=NCV1
         np2=NCV2 

         open(11,file='params.in') !user inputs
         read(11,*) temperature
         read(11,*) bee !b parameter in mABP
         read(11,*) cee !c parameter in mABP (is taken as c*dt)
         read(11,*) alp !a/dx in S_a
         read(11,*) (widths(j),j=1,3) ! rectangular
         read(11,*) radius !sphere restraint
         read(11,*) scal !power for shape
         read(11,*) irest !restart or not
         bolt = temperature/11604.467d0 *96.48d0 !hard code TEMPERATURE!!!!!!
         imabp=0
         close(11) !close it
cWTmetaD         omega = bolt*bee*cee !WTmetaD conversion, and muTmetaD
cWTmetaD         DelT = bolt*bee/(1d0-bee) !ditto

         pi=acos(-1d0) !have some pi
         pi2=pi*2d0    !twice
c         alp2 = alp*alp !save the alpha^2
         dx = CVMAXd0/dble(nbin) !grid width
         tol = alp*dx
         tol2=4d0*dx !not too small...
         ialp = int(alp)
         do j=1,nbin !pre-compute delta_alpha "hill functions"
            xatj = (j-1)*dx
            do i=1,nbin
               xati = (i-1)*dx
               arg2=xati-xatj
               ib1 = int((arg2)/dx)+1      !theres an issue here
               if(ib1.gt.nbin)ib1=nbin
               if(ib1.lt.1)ib1=1     
               ib1=nbin*i-(nbin-j)!ib1
               look(i,j)=ib1
!-----MOLLIIIIIII
               ex = arg2/tol
               aex=ex*ex
               aex=sqrt(aex)
               aex2=aex*aex
               func=0d0
               dermoli=0d0
               if(aex.lt.1d0)then
                  func = exp(-scal/(1d0-aex2))/exp(-scal) !so height is 1 at x=0
                  aa=scal*func/((1d0-aex2)*(1d0-aex2))!should be negative?
                  dermoli=aa*2d0*ex/tol
               endif
               dela(1,i,j)=func
                ddela(i,j) =dermoli
!--------------------------------------------------------------------
!              this is used for the "vanishingly small" alpha,
!     do be careful here...
               ex2 = arg2/tol2
               aex2=ex2*ex2
               aex2=sqrt(aex2)
               aex22=aex2*aex2
               if(aex2.lt.1d0)then
                  func2 = exp(-20d0/(1d0-aex22))/exp(-20d0) !so height is 1 at x=0
               else
                  func2=0d0
               endif
               dela(2,i,j)=func2!exp(-arg2*arg2/2d0**2)  !for mabp stuff and muTmetaD
!----------MOLLLIIIII DD
            enddo
         enddo


         do j=1,nbin
            do i=1,nbin                  
               pop(i,j) = 0d0   !initialize arrays
               if(imabp.eq.0)then
               decon(i,j)=0.1d0 !to avoid log(0)
               else
               decon(i,j)=0d0
               endif
               dpop(1,i,j) = 0d0
               dpop(2,i,j) = 0d0
            enddo
         enddo
         if(irest.eq.1)then
            open(99,file='restartABP')
            do j=1,nbin
               do i=1,nbin
                  read(99,111) pop(i,j), dpop(1,i,j), dpop(2,i,j), 
     .                 decon(i,j)
               enddo
            enddo
            close(99)
         endif
         xwrp=xx
      endif                     !end the initial stuff
!do basic PBC here with a cubic box, careful to hardcode the box

      do k=1,3
         widths2(k)=widths(k)*widths(k)
         widths24(k)=widths2(k)/4d0
      enddo


      do l=1,1
      do i=1,np
         do j=i,i
               do k=1,3
                  xxi=xx(i*3-(3-k))
                  xxj=xwrp(j*3-(3-k))
                  dr=xxi-xxj
                  dr2=dr*dr                  
                  if(dr2.gt.widths24(k))then !somebody need unwrappin, just roll j
                     if(dr.gt.0d0)then
                        xx(i*3-(3-k))=xx(i*3-(3-k))-widths(k)                        
                     elseif(dr.lt.0d0)then
                        xx(i*3-(3-k))=xx(i*3-(3-k))+widths(k)                        
                     endif
                  endif
               enddo
         enddo
      enddo
      enddo!l
      xwrp=xx

      jaco = 0d0 !this will hold the RMSD derivative
      !angl1 and angl2 are the CV-in RMSD
      sum=0d0
      do i=1,np1!np1 particles in here
         do j=1,3
            ip=i
            sum=sum+(xx(ip*3-(3-j))-x0(3*ip-(3-j)))**2
            jaco(1,ip*3-(3-j))=(xx(ip*3-(3-j))-x0(3*ip-(3-j)))
c     .           +jaco(1,ip*3-(3-j))
         enddo
      enddo
      angl1=sqrt(0.01d0+sum/dble(np1)/3d0) !root-mean squared
      sum=0d0
      do i=1,np2
         do j=1,3
            ip=i+np1
            sum=sum+(xx(ip*3-(3-j))-x0(3*ip-(3-j)))**2
            jaco(2,ip*3-(3-j))=(xx(ip*3-(3-j))-x0(3*ip-(3-j)))
c     .           +jaco(2,ip*3-(3-j))
         enddo
      enddo
      angl2=sqrt(0.01d0+sum/dble(np2)/3d0) !root-mean squared
      do i=1,np1
         do j=1,3
            ip=i*3-(3-j)
            jaco(1,ip)=jaco(1,ip)/angl1/3d0/dble(np1)
         enddo
      enddo
      do i=1,np2
         do j=1,3
            ip=(i+np1)*3-(3-j)
            jaco(2,ip)=jaco(2,ip)/angl2/3d0/dble(np2)
         enddo
      enddo

      ibin = int(angl1/dx)+1      
      if(ibin.gt.nbin)ibin=nbin
      if(ibin.lt.1)ibin=1     
      ibin1=ibin
      ibin = int(angl2/dx)+1      
      if(ibin.gt.nbin)ibin=nbin
      if(ibin.lt.1)ibin=1     
      ibin2=ibin

!-------------------------------------------------------------------------
!      UPDATE THE BIAS POTENTIAL HERE ---
!-------------------------------------------------------------------------

cWTmetaD      ave1(1) = -dpop(1,ibin1,ibin2)
cWTmetaD      ave1(2) = -dpop(2,ibin1,ibin2)            
cWTmetaD      s=exp(-pop(ibin1,ibin2)/DelT)*omega

cmABP      denom = 1d0+cee*(1d0-bee)*pop(ibin1,ibin2)            
cmABP      ave1(1) = -cee*bee*bolt*dpop(1,ibin1,ibin2)/denom 
cmABP      denom = 1d0+cee*(1d0-bee)*pop(ibin1,ibin2)      
cmABP      ave1(2) = -cee*bee*bolt*dpop(2,ibin1,ibin2)/denom 
      !update grids

      jst=ibin2-(ialp+1) !11 because used 10 in def of tol
      jst=max(jst,1)
      jnd=ibin2+(ialp+1)
      jnd=min(nbin,jnd)
      do j =jst,jnd!1,nbin!jst,jnd!1,nbin!jst,jnd
         boomc= dela(2,ibin2,j) !this is hill with "vanishingly small" width
         booma=dela(1,ibin2,j)  !this is hill for biasing
         ist=ibin1-(ialp+1)           !11 because used 10 in def of tol
         ist=max(ist,1)
         ind=ibin1+(ialp+1)
         ind=min(nbin,ind)
         do i = ist,ind!1,nbin!ist,ind!1,nbin
            boomcc= dela(2,ibin1,i)!second part of "vanishingly small" width
            boom = dela(1,ibin1,i)*booma ! product of two hills for bias
cWTmetaD            pop(i,j) = pop(i,j) + boom*s
cmABP            pop(i,j) = pop(i,j) + boom
!deconvoluted free energy estimate ------------------------------
            decon(i,j)=decon(i,j)+boomcc*boomc
!----------------------------------------------------------------
cWTmetaD            dpop(1,i,j) = dpop(1,i,j)+ddela(ibin1,i)
cWTmetaD     .      *dela(1,ibin2,j)*s
cWTmetaD            dpop(2,i,j) = dpop(2,i,j)+ddela(ibin2,j)
cWTmetaD     .          *dela(1,ibin1,i)*s
cmABP            dpop(1,i,j) = dpop(1,i,j)+ddela(ibin1,i)
cmABP     .      *dela(1,ibin2,j)
cmABP            dpop(2,i,j) = dpop(2,i,j)+ddela(ibin2,j)
cmABP     .          *dela(1,ibin1,i)
         enddo
      enddo
!----------------------------------------------------------------------
98765 continue
!Add the bias forces back to system force
      fr1=0d0
      fr2=0d0
      edge=CVRESTd0 
      if(angl1.gt.edge)then
         fr1=-400d0*bolt*(angl1-edge)
      endif
      if(angl2.gt.edge)then
         fr2=-400d0*bolt*(angl2-edge)
      endif
      ff=0d0
      sp=400d0*bolt !or want 100?
      do i=1,np
         sum=0d0
         do j=1,3
            ip=i*3-(3-j)
            ff(ip) = ff(ip)+(ave1(1)+fr1)*jaco(1,ip)+
     .           (ave1(2)+fr2)*jaco(2,ip)
            sum=sum+(xx(ip)-point(1,j))*(xx(ip)-point(1,j))
! this restrains to the box exactly. i didn't like it so much.
c            if(xx(ip).ge.width)then!box restraint
c               ff(ip)=ff(ip) - 100d0*bolt*(xx(ip)-width)
c            endif
c            if(xx(ip).le.0d0)then
c               ff(ip)=ff(ip) - 100d0*bolt*(xx(ip)-0d0)
c            endif
         enddo
         sum=sqrt(sum)
c         sum=sp/sum
         if(sum.gt.radius)then!sphere restraint
            do j=1,3
               ip=i*3-(3-j)
               ff(ip) = ff(ip)-sp*(sum-radius)*(xx(ip)-point(1,j))/sum
            enddo
         endif
      enddo
!           dihedral restraint
! which we used to fix the azepine pucker
c      jaco = 0d0 !this will hold the phi-psi derivatives
c      natom1=5!5 
c      natom2=6!7
c      natom3=7!9
c      natom4=8!15
c      qa(1) = xx(3*natom1-(3-1))
c      qa(2) = xx(3*natom1-(3-2))
c      qa(3) = xx(3*natom1-(3-3))
c      qa(4) = xx(3*natom2-(3-1))
c      qa(5) = xx(3*natom2-(3-2))
c      qa(6) = xx(3*natom2-(3-3))
c      qa(7) = xx(3*natom3-(3-1))
c      qa(8) = xx(3*natom3-(3-2))
c      qa(9) = xx(3*natom3-(3-3))
c      qa(10) = xx(3*natom4-(3-1))
c      qa(11) = xx(3*natom4-(3-2))
c      qa(12) = xx(3*natom4-(3-3))
c      do i=1,3
c         do j=1,4
c            deid(i,j) = 0d0
c         enddo
c      enddo
c      angle = 0d0
c      call bmdDI(qa,deid,angle)
c      pi=acos(-1d0)
c      angle = angle*pi/180d0
c      if(istep.eq.0)dih=angle
c      do idm=1,3
c         ip = natom1*3-(3-idm)
c         jaco(1,ip)=jaco(1,ip)+deid(idm,1)
c      enddo
c      do idm=1,3
c         ip = natom2*3-(3-idm)
c         jaco(1,ip)=jaco(1,ip)+deid(idm,2)
c      enddo
c      do idm=1,3
c         ip = natom3*3-(3-idm)
c         jaco(1,ip)=jaco(1,ip)+deid(idm,3)
c      enddo
c      do idm=1,3
c         ip = natom4*3-(3-idm)
c         jaco(1,ip)=jaco(1,ip)+deid(idm,4)
c      enddo
c      do i=5,8
c         do j=1,3
c            ip=i*3-(3-j)
c            ff(ip) = ff(ip)-sp*(angle-dih)*jaco(1,ip)
c         enddo
c      enddo

!----------------------------------------------------------------------
!     Write a restart file and get convergence curve
      if(mod(istep,50000).eq.0)then
c      if(1.eq.0)then !debug
         write(81,*) np
         write(81,*)
         do i=1,np
            write(81,*) 'C ', xx(i*3-(3-1))*10d0, xx(i*3-(3-2))*10d0,
     .           xx(i*3-(3-3))*10d0
         enddo
         call flush(81)
         open(99,file='restartABP')
         do j=1,nbin
            do i=1,nbin
               write(99,111) pop(i,j), dpop(1,i,j), dpop(2,i,j), 
     .              decon(i,j)
            enddo
         enddo
         close(99)
         !Get the zero-of energy
         bang=0d0    
         do j=1,nbin
            do i=1,nbin
cWTmetaD               sf=exp(-pop(i,j)/DelT)*omega
cWTmetaD               freest = bolt*log(sf)-pop(i,j)
cmABP                  freest=-bolt*log(decon(i,j)*pop(i,j)**(bee/(1d0-bee)))               
               if(freest.lt.bang)bang=freest !Sets the zero
            enddo
         enddo
         er = 0d0
         open(89,file='freeE')
         sumd=0d0
         orml=0d0
         do j=1,nbin
            at2 = dx * (dble(j) - 0.5d0)
            do i=1,nbin           
               at1 = dx * (dble(i) - 0.5d0)
cWTmetaD               sf=exp(-pop(i,j)/DelT)*omega
cWTmetaD               freest = bolt*log(sf)-pop(i,j)
cmABP                  freest=-bolt*log(decon(i,j)*pop(i,j)**(bee/(1d0-bee)))               
               freest=freest-bang !put min to zero
               write(89,*) at1, at2, freest, decon(i,j) !write free energy to fort.89
            enddo
            write(89,*)
         enddo
         close(89)
77777    continue
!                 TimeStep, CV1, CV2, hill height
         write(88,*) istep, angl1, angl2, 
cmABP     .        cee*bee*bolt/(1d0+cee*(1d0-bee)*pop(ibin1,ibin2))  !last dist to middl
cWTmetaD     .        s  !last dist to middl
         call flush(88) !write to fort.88
      endif

 111   format(4(E12.5,1X))
      return
      end !end the ABP routine(s)-------------------------
!---------------------------------------------------------
!---------------------------------------------------------
!---------------------------------------------------------
!---------------------------------------------------------
      subroutine bmdDI(qa,deid,angle)
      implicit integer*4 (i-n)
      implicit real*8 (a-h,o-z)
      real*8 deid(3,4),qa(12),angle,radian
      
      radian = 180d0/acos(-1d0)
      do i=1,4
         do  j=1,3
            deid(j,i)=0d0
         enddo
      enddo
      
      xia = qa(1)
      yia = qa(2)
      zia = qa(3)
      xib = qa(4)
      yib = qa(5)
      zib = qa(6)
      xic = qa(7)
      yic = qa(8)
      zic = qa(9)
      xid = qa(10)
      yid = qa(11)
      zid = qa(12)
      xba = xib - xia
      yba = yib - yia
      zba = zib - zia
      xcb = xic - xib
      ycb = yic - yib
      zcb = zic - zib      
      xdc = xid - xic
      ydc = yid - yic
      zdc = zid - zic
      xt = yba*zcb - ycb*zba
      yt = zba*xcb - zcb*xba
      zt = xba*ycb - xcb*yba
      xu = ycb*zdc - ydc*zcb
      yu = zcb*xdc - zdc*xcb
      zu = xcb*ydc - xdc*ycb
      xtu = yt*zu - yu*zt
      ytu = zt*xu - zu*xt
      ztu = xt*yu - xu*yt
      rt2 = xt*xt + yt*yt + zt*zt
      ru2 = xu*xu + yu*yu + zu*zu
      if(rt2.lt. 0.000000001)then
       rt2=0.000000001
      endif
      if(ru2.lt. 0.000000001)then
       ru2 = 0.000000001
      endif
      rtru = sqrt(rt2 * ru2)
      if(rtru.ne.0d0)then
      rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb)
      cosine = (xt*xu + yt*yu + zt*zu) / rtru
      sine = (xcb*xtu + ycb*ytu + zcb*ztu) / (rcb*rtru)
      cosine = min(1.0d0,max(-1.0d0,cosine))
      angle = radian * acos(cosine)
      if (sine .lt. 0.0d0)  angle = -angle+360d0
      if(angle.ge.360d0)angle=angle-360d0
      dedphi = 1d0
      xca = xic - xia
      yca = yic - yia
      zca = zic - zia      
      xdb = xid - xib
      ydb = yid - yib
      zdb = zid - zib
      dedxt = dedphi * (yt*zcb - ycb*zt) / (rt2*rcb)
      dedyt = dedphi * (zt*xcb - zcb*xt) / (rt2*rcb)
      dedzt = dedphi * (xt*ycb - xcb*yt) / (rt2*rcb)
      dedxu = -dedphi * (yu*zcb - ycb*zu) / (ru2*rcb)
      dedyu = -dedphi * (zu*xcb - zcb*xu) / (ru2*rcb)
      dedzu = -dedphi * (xu*ycb - xcb*yu) / (ru2*rcb)
      dedxia = zcb*dedyt - ycb*dedzt
      dedyia = xcb*dedzt - zcb*dedxt
      dedzia = ycb*dedxt - xcb*dedyt
      dedxib = yca*dedzt - zca*dedyt + zdc*dedyu - ydc*dedzu
      dedyib = zca*dedxt - xca*dedzt + xdc*dedzu - zdc*dedxu
      dedzib = xca*dedyt - yca*dedxt + ydc*dedxu - xdc*dedyu
      dedxic = zba*dedyt - yba*dedzt + ydb*dedzu - zdb*dedyu
      dedyic = xba*dedzt - zba*dedxt + zdb*dedxu - xdb*dedzu
      dedzic = yba*dedxt - xba*dedyt + xdb*dedyu - ydb*dedxu
      dedxid = zcb*dedyu - ycb*dedzu
      dedyid = xcb*dedzu - zcb*dedxu
      dedzid = ycb*dedxu - xcb*dedyu
      deid(1,1) = deid(1,1) + dedxia
      deid(2,1) = deid(2,1) + dedyia 
      deid(3,1) = deid(3,1) + dedzia
      deid(1,2) = deid(1,2) + dedxib
      deid(2,2) = deid(2,2) + dedyib
      deid(3,2) = deid(3,2) + dedzib
      deid(1,3) = deid(1,3) + dedxic
      deid(2,3) = deid(2,3) + dedyic
      deid(3,3) = deid(3,3) + dedzic
      deid(1,4) = deid(1,4) + dedxid
      deid(2,4) = deid(2,4) + dedyid
      deid(3,4) = deid(3,4) + dedzid
      else
       write(*,*) 'failed dihed'
       stop
      endif
c      deid = deid !* radian
      RETURN
      END
