c      Truly locally Discontinuous Galerkin formulation for 
c      Poisson problem - Lagrangian interpolation 
c
c      
c
c     ************************************************************
c     *                                                          *
c     *            * * *   Poisson Problem      * * *            *
c     *                                                          *
c     *                                                          *
c     *     DUAL HIBRID TRULY LDG FINITE ELEMENT METHODS         *
c     *                                                          *
c     *                 LOCALLY CONSERVATIVE		             *
c     *                                                          *
c     *                                                          *
c     *                    ABIMAEL LOULA                         *
c     *                                                          *
c     *                   Setembro de 2010                       *
c     ************************************************************

c
c         Residuo no interior dos elementos:
c
c
c         (grad u , grad v ) - (f,v) +  
c
c
c         Termos de saltos nas arestas:
c
c
c         beta(u - up, v-vp)
c
c         -(grad(u).n , (v-vp))
c
c         lambda((u-up) , grad(v).n)
c
c
c.... program to set storage capacity, precision and input/output units
c
      common /bpoint/ mfirst,mlast,ilast,mtot,iprec
      common /iounit/ iin,iout,iecho,ioupp,itest1,itest2 
      character*4 ia 
      parameter (ndim=200000000) 
      common a(ndim) 
      common /dictn/ ia(10000000)
c
c        mfirst = address of first available word in blank common
c        mlast  = address of last available word in blank common
c        mtot   = total storage allocated to blank common
c        iprec  = precision flag; eq.1, single precision
c                                 eq.2, double precision
c
      iin= 8
      iout= 11
      iecho=12
      ioupp=13
	itest1=14
	itest2=15
      open(unit=iin, file= 'ldgcPkQl.dat',status='old')
      open(unit=iecho, file= 'ldgcPkQl.eco')	
      open(unit=iout, file= 'errofem-ns-c.con')      
      open(unit=ioupp, file= 'ldgherrol-ns-c.con')
      open(unit=itest1, file= 'errointer-c.con')      
      open(unit=itest2, file= 'teste2-c.con')
c
      mfirst = 1 
      ilast  = 0 
      mlast  = ndim 
      mtot   = ndim 
      iprec  = 2 
c
c        iin    = input unit number
c        iecho  = output unit of input data
c        ioupp  = output unit of post-processed displacements
c
c
         call lpgm
c
c.... system-dependent unit/file specifications
c
c
         close(iin)
         close(iecho)
	   close(iout)
	   close(ioupp)
	   close(itest1)
	   close(itest2)
c
      stop
      end
c
c**** new **********************************************************************
      subroutine lpgm
c
c.... LPGM - a linear static finite element analysis program for 
c            Petrov Galerkin methods : global driver
c
      real*8 zero,pt1667,pt25,pt5,one,two,three,four,five,six,
     & tempf
      character*4 title
c
c.... remove above card for single-precision operation
c
c.... catalog of common statements
c
      common /bpoint/ mfirst,mlast,ilast,mtot,iprec
      common /colhtc/ neqc
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
      common /genelc/ n,nel(3),incel(3),inc(3)                          
      common /genflc/ tempf(6,20),nf,numgpf,nincf(3),incf(3)
      common /info  / iexec,iprtin,irank,nsd,numnp,ndof,nlvect,
     &                numeg,nmultpd,nmultpc,nedge
      common /iounit/ iin,iout,iecho,ioupp,itest1,itest2 
      common /labels/ labeld(3),label1(16),label2(3)
      common /spoint/ mpd,mpdprj,mpx,mpidc,mpiedge,mpf,mpdiag,
     &                mpngrp,mpalhs,mpdlhs,mpbrhs,mped,index
      common /titlec/ title(20)
      character*4 ia 
      common a(1) 
      common /dictn/ ia(1) 

c
c.... input phase
c
      call echo
c
  100 continue
c
      read(iin,1000) title
      if (title(1).eq.'*end') return
      read(iin,2000) iexec,iprtin,irank,
     &           nsd,numnp,ndof,nlvect,numeg,nedge,npar
c
      write(iecho,3000) title , iexec,iprtin
      write(iecho,4000) irank ,   nsd, numnp,  ndof,
     &                 nlvect,numeg,nedge,npar
c
c.... initialization phase
c
c....    set memory pointers for static data arrays,
c        and call associated input routines 
c
      mpx     = mpoint('x       ',nsd   ,numnp  ,0     ,iprec)
      mped    = mpoint('ideg    ',ndof  ,nedge  ,0     ,1)
      mpiedge = mpoint('iedge   ',npar  ,nedge  ,0     ,1)
c
c.... input coordinate data
c
      call coord(a(mpx),nsd,numnp,iprtin)
c
c.... input boundary condition data and establish equation numbers
c
      call readbcedge(a(mped),npar,nedge,ndof,iprtin)
c
      mpngrp = mpoint('ngrp    ',numeg ,0     ,0     ,1)
c
c.... input element data
c
      call elemnt('input_el',a(mpngrp))
c
       nmultpd = nedge*npar
c
      mpd    = mpoint('d       ',ndof  ,nmultpc ,0     ,iprec)
      mpdprj = mpoint('dprj    ',ndof  ,nmultpd ,0     ,iprec)
c
      if (nlvect.eq.0) then
         mpf = 1
      else
         mpf = mpoint('f       ',ndof  ,nmultpc ,nlvect,iprec)
      endif
c
c.... input nodal force and prescribed kinematic boundary-value data
c
      if (nlvect.gt.0) call input(a(mpf),ndof,nmultpc,0,nlvect,
     &                            iprtin)
c
c.... allocate memory for id array and clear 
c
      mpidc  = mpoint('idc     ',ndof  ,nmultpc ,0     ,1)
c
c.... allocate memory for idiag array and clear 
c
      mpdiag = mpoint('idiag   ',nmultpc  ,0     ,0     ,1)
      call iclear(a(mpdiag),nmultpc)
c
c.... input side data for ipar 
c
      call elemnt('input_sd',a(mpngrp))
c
c
c
c.... determine addresses of diagonals in left-hand-side matrix
c
      call diag(a(mpdiag),neqc,nalhs)
c
c.... allocate memory for global equation system
c
      mpalhs = mpoint('alhs    ',nalhs,0,0,iprec)
      mpdlhs = mpoint('dlhs    ',nalhs,0,0,iprec)
      mpbrhs = mpoint('brhs    ',neqc ,0,0,iprec)
      meanbw = nalhs/neqc
      nwords = mtot - mlast + mfirst - 1
c
c.... write equation system data
c
      write(iecho,5000) title,neqc,nalhs,meanbw,nwords
c
c.... solution phase
c
c
      if (iexec.eq.1) call driver(neqc,nalhs)
c
c.... print memory-pointer dictionary
c
c      call prtdc
c
c
c.... print elapsed time summary
c
      go to 100
c
 1000 format(20a4)
 2000 format(16i10)
 3000 format(///,20a4///
     &' e x e c u t i o n   c o n t r o l   i n f o r m a t i o n '//5x,
     &' execution code  . . . . . . . . . . . . . . (iexec ) = ',i10//5x,
     &'    eq. 0, data check                                   ',   /5x,
     &'    eq. 1, execution                                    ',  //5x,
     &' input data print code . . . . . . . . . . . (iprtin) = ',i10//5x,
     &'    eq. 0, print nodal and element input data           ',   /5x,
     &'    eq. 1, do not print nodal and element input data    ',   /5x)
 4000 format(5x,
     &' rank check code . . . . . . . . . . . . . . (irank ) = ',i10//5x,
     &'    eq. 0, do not perform rank check                    ',   /5x,
     &'    eq. 1, print numbers of zero and nonpositive pivots ',   /5x,
     &'    eq. 2, print all pivots                             ',  //5x,
     &' number of space dimensions  . . . . . . . .(nsd   ) = ',i10//5x,
     &' number of nodal points  . . . . . . . . . .(numnp ) = ',i10//5x,
     &' number of nodal degrees-of-freedom  . . . .(ndof  ) = ',i10//5x,
     &' number of load vectors  . . . . . . . . . .(nlvect) = ',i10//5x,
     &' number of element groups  . . . . . . . . .(numeg ) = ',i10//5x,
     &' number of edge. . . .  .  . . . . . . . . .(nedges) = ',i10//5x)
4500  format(///,20a4//5x,
     &'nsd = ',i10,5x,'nmultpc= ',i10,5x,'ndof = ',i10,5x,
     &'nlvect = ',i10,5x,'numeg = ',i10)
 5000 format(///,20a4///
     &' e q u a t i o n    s y s t e m    d a t a              ',  //5x,
     &' number of equations . . . . . . . . . . . . (neqc  ) = ',i8//5x,   
     &' number of terms in left-hand-side matrix  . (nalhs ) = ',i8//5x,
     &' mean half bandwidth . . . . . . . . . . . . (meanbw) = ',i8//5x,
     &' total length of blank common required . . . (nwords) = ',i8    )
c
      end
c**** new **********************************************************************
      subroutine driver(neqc,nalhs)
c
c.... solution driver program 
c
      common /info  / iexec,iprtin,irank,nsd,numnp,ndof,nlvect,
     &                numeg,nmultpd,nmultpc,nedge

      common /iounit/ iin,iout,iecho,ioupp,itest1,itest2 
      common /spoint/ mpd,mpdprj,mpx,mpidc,mpiedge,mpf,mpdiag,
     &                mpngrp,mpalhs,mpdlhs,mpbrhs,mped,index

      character*4 ia 
      common a(1) 
      common /dictn/ ia(1)
c
c
c      clear left and right hand side
c
      call clear(a(mpalhs),nalhs)
      call clear(a(mpdlhs),nalhs)
c
      call clear(a(mpbrhs),neqc)
c
c      account the nodal forces in the r.h.s.
c
      if (nlvect.gt.0)
     &   call load(a(mpidc),a(mpf),a(mpbrhs),ndof,nmultpc,nlvect)
c
c	clear displacement array
c
      call clear(a(mpd),ndof*nmultpc)
      call clear(a(mpdprj),ndof*nmultpd)
c
      if (nlvect.gt.0)
     &   call ftod(a(mpidc),a(mpd),a(mpf),ndof,nmultpc,nlvect)
c
c      form the l.h.s and r.h.s. at element level
c
c      write(*,*) "Aqui 1"
      call elemnt('form_lrs',a(mpngrp))
c      write(*,*) "Aqui 2"
c
cc      call teste(a(mpidc),ndof,nmultpc)
c
c      factorization of the stiffnes matrix
c
      if(neqc.eq.0) go to 1111
c      write(*,*) "Aqui 5"
      call factns(a(mpalhs),a(mpdlhs),a(mpdiag),neqc)
c
c      back substitution
c
      call backns(a(mpalhs),a(mpdlhs),a(mpbrhs),a(mpdiag),neqc)
c
 1111 continue  
c
      call  btod(a(mpidc),a(mpd),a(mpbrhs),ndof,nmultpc)
c
c.... write output 
c
c
         call printd(' d i s p l a c e m e n t s                  ',
     &               a(mpd),ndof,nmultpc,iecho)
c
c
c	post-processing phase
c
      call elemnt('pos_proc',a(mpngrp))
c
c
100   continue
c
      return
      end 
c**** new **********************************************************************
      subroutine addlhs(alhs,eleffm,idiag,lm,nee,diag)
c
c.... program to add element left-hand-side matrix to
c        global left-hand-side matrix
c
c        diag = .true., add diagonal element matrix
c
c        diag = .false, add upper triangle of full element matrix
c
      implicit real*8 (a-h,o-z)
c
c.... remove above card for single-precision operation
c
      logical diag
      dimension alhs(*),eleffm(nee,*),idiag(*),lm(*)
c
      if (diag) then
c
         do 100 j=1,nee
         k = lm(j)
         if (k.gt.0) then
            l = idiag(k)
            alhs(l) = alhs(l) + eleffm(j,j)
         endif
  100    continue
c
      else
c
         do 300 j=1,nee
         k = lm(j)
         if (k.gt.0) then
c
            do 200 i=1,j
            m = lm(i)
            if (m.gt.0) then
               if (k.ge.m) then
                  l = idiag(k) - k + m
               else
                  l = idiag(m) - m + k
               endif
               alhs(l) = alhs(l) + eleffm(i,j)
            endif
  200       continue
c
         endif
  300    continue
c
      endif
c
      return
      end
c**** new **********************************************************************
      subroutine addrhs(brhs,elresf,lm,nee)
c
c.... program to add element residual-force vector to
c        global right-hand-side vector
c
      implicit real*8 (a-h,o-z)
c
c.... remove above card for single-precision operation
c
      dimension brhs(*),elresf(*),lm(*)
c
      do 100 j=1,nee
      k = lm(j)
      if (k.gt.0) brhs(k) = brhs(k) + elresf(j)
  100 continue
c
      return
      end
c**** new **********************************************************************
      subroutine back(a,b,idiag,neqc)
c
c.... program to perform forward reduction and back substitution
c
      implicit real*8 (a-h,o-z)
c
c.... remove above card for single-precision operation
c
      dimension a(*),b(*),idiag(*)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
c
c.... forward reduction
c
      jj = 0
c
      do 100 j=1,neqc
      jjlast = jj
      jj     = idiag(j)
      jcolht = jj - jjlast
      if (jcolht.gt.1) 
     &   b(j) = b(j) - coldot(a(jjlast+1),b(j-jcolht+1),jcolht-1)
  100 continue
c
c.... diagonal scaling
c
      do 200 j=1,neqc
      ajj = a(idiag(j))
      if (ajj.ne.zero) b(j) = b(j)/ajj
  200 continue
c
c.... back substitution
c
      if (neqc.eq.1) return
      jjnext = idiag(neqc)
c
      do 400 j=neqc,2,-1
      jj     = jjnext
      jjnext = idiag(j-1)
      jcolht = jj - jjnext
      if (jcolht.gt.1) then
         bj = b(j)
         istart = j - jcolht + 1
         jtemp  = jjnext - istart + 1
c
         do 300 i=istart,j-1
         b(i) = b(i) - a(jtemp+i)*bj
  300    continue
c
      endif
c
  400 continue
c
      return
      end
c**** new **********************************************************************
      subroutine readbcedge(ideg,npar,nedge,ndof,iprtin)
c
c.... program to read, generate and write boundary condition data
c        and establish equation numbers
c
      dimension ideg(ndof,*)
c
      common /iounit/ iin,iout,iecho,ioupp,itest1,itest2 
      logical pflag
c
      call iclear(ideg,ndof*nedge)
      call igen(ideg,ndof)
c
      if (iprtin.eq.0) then
         nn=0
         do 200 n=1,nedge
         pflag = .false.
c
         do 100 i=1,ndof
         if (ideg(i,n).ne.0) pflag = .true.
  100    continue
c
         if (pflag) then      
            nn = nn + 1
            if (mod(nn,50).eq.1) write(iecho,1000) (i,i=1,ndof)
            write(iecho,2000) n,(ideg(i,n),i=1,ndof)
         endif
  200    continue
      endif
c
      return
c
 1000 format(//, 10x, ' e d g e   b o u n d a r y   c o n d i t i o n  c o
     & d e s'///
     & 5x,'   node no.',3x,6(13x,'dof',i1:)//)
 1100 format(//, 10x,' n o d a l   b o u n d a r y   c o n d i t i o n  c o
     & d e s'///
     & 5x,'   node no.',3x,6(13x,'dof',i1:)//)
 2000 format(6x,i10,5x,6(5x,i10))
c
      end
c**** new **********************************************************************
c**** new **********************************************************************
      block data
c
c.... program to define output labels and numerical constants
c
      implicit real*8 (a-h,o-z)
c
c.... remove above card for single-precision operation
c
c
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
      common /labels/ labeld(3),label1(16),label2(3)
c
c        labeld(3)  = displacement, velocity and acceleration labels
c        label1(16) = output labels for element-type 1
c        label2(3)  = output labels for element-type 2
c
c.... note: add label arrays for any additional elements        
c
      data   zero,pt1667,pt25,pt5
     &      /0.0d0,0.1666666666666667d0,0.25d0,0.5d0/,
     &       one,two,three,four,five,six
     &      /1.0d0,2.0d0,3.0d0,4.0d0,5.0d0,6.0d0/
c
      data labeld/'disp','vel ','acc '/
c
      data label1/'s 11','s 22','s 12','s 33','ps 1','ps 2',
     &            'tau ','sang','e 11','e 22','g 12','e 33',
     &            'pe 1','pe 2','gam ','eang'/
c
      data label2/'strs','forc','strn'/
c
      end
c**** new **********************************************************************
      subroutine btdb(elstif,b,db,nee,nrowb,nstr)
c
c.... program to multiply b(transpose) * db taking account of symmetry
c        and accumulate into element stiffness matrix
c
      implicit real*8 (a-h,o-z)
c
c.... remove above card for single-precision operation
c
      dimension elstif(nee,*),b(nrowb,*),db(nrowb,*)
c
      do 200 j=1,nee
c
      do 100 i=1,j
      elstif(i,j) = elstif(i,j) + coldot(b(1,i),db(1,j),nstr)   
  100 continue
c
  200 continue
c
      return
      end
c**** new **********************************************************************
      subroutine teste(idc,ndof,numnp)
c
c.... program to perform transfer from r.h.s. to displacement array
c
      implicit real*8 (a-h,o-z)
c
c.... remove above card for single-precision operation
c
      dimension idc(ndof,*)
c
         do n=1,numnp
           write(88,*) n,idc(1,n)
         end do
         return
         end
c**** new **********************************************************************
      subroutine btod(idc,d,brhs,ndof,numnp)
c
c.... program to perform transfer from r.h.s. to displacement array
c
      implicit real*8 (a-h,o-z)
c
c.... remove above card for single-precision operation
c
      dimension idc(ndof,*),d(ndof,*),brhs(*)
c
         do 200 i=1,ndof
c
         do 100 j=1,numnp
         k = idc(i,j)
         if (k.gt.0) d(i,j) = brhs(k)
  100    continue
c
  200    continue
c
      return
      end
c**** new **********************************************************************
      subroutine clear(a,m)
c
c.... program to clear a floating-point array
c
      implicit real*8 (a-h,o-z)
c
c.... remove above card for single-precision operation
c
      dimension a(*)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
c
      do 100 i=1,m
      a(i) = zero
  100 continue
c
      return
      end
c**** new **********************************************************************
      function coldot(a,b,n)
c
c.... program to compute the dot product of vectors stored column-wise
c
      implicit real*8 (a-h,o-z)
c
c.... remove above card for single-precision operation
c
      dimension a(*),b(*)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
c
      coldot = zero
c
      do 100 i=1,n
      coldot = coldot + a(i)*b(i)
  100 continue
c
      return
      end
c**** new **********************************************************************
      subroutine colht(idiag,lm,ned,nen,numel,neqc)
c
c.... program to compute column heights in global left-hand-side matrix
c
      dimension idiag(*),lm(ned,nen,*)
c
      do 500 k=1,numel
      min = neqc
c
      do 200 j=1,nen
c
      do 100 i=1,ned
      num = lm(i,j,k)
      if (num.gt.0) min = min0(min,num)
  100 continue
c
  200 continue
c
      do 400 j=1,nen
c
      do 300 i=1,ned
      num = lm(i,j,k)
      if (num.gt.0) then
         m = num - min
         if (m.gt.idiag(num)) idiag(num) = m
      endif
c
  300 continue

  400 continue
c
  500 continue
c
      return
      end
c**** new **********************************************************************
      subroutine coord(x,nsd,numnp,iprtin)
c
c.... program to read, generate and write coordinate data
c
c        x(nsd,numnp) = coordinate array
c
      implicit real*8 (a-h,o-z)
c
c.... remove above card for single-precision operation
c
      dimension x(nsd,*)
      common /iounit/ iin,iout,iecho,ioupp,itest1,itest2 
c
      call genfl(x,nsd)
c
      if (iprtin.eq.1) return
c
      do 100 n=1,numnp
      if (mod(n,50).eq.1) write(iecho,1000) (i,i=1,nsd)
      write(iecho,2000) n,(x(i,n),i=1,nsd)   
  100 continue
c
      return
c
 1000 format(///,' n o d a l   c o o r d i n a t e   d a t a '///5x,
     &' node no.',3(13x,' x',i1,' ',:)//)
 2000 format(6x,i10,10x,3(1pe15.8,2x))
      end
c**** new ********************************************************************** 
      subroutine dctnry(name,ndim1,ndim2,ndim3,mpoint,ipr,mlast,ilast) 
c 
c.... program to store pointer information in dictionary 
c 
      character*4 name(2) 
      character*4 ia 
      common na(1) 
      common /dictn/ ia(1) 
      common /iounit/ iin,iout,iecho,ioupp,itest1,itest2 
c 
      mlast = mlast - 5 
      ia(ilast+1) = name(1) 
      ia(ilast+2) = name(2) 
      na(mlast+1) = mpoint 
      na(mlast+2) = ndim1 
      na(mlast+3) = ndim2 
      na(mlast+4) = ndim3 
      na(mlast+5) = ipr 
      ilast = ilast + 2 
c 
      return 
      end 
c**** new **********************************************************************
      subroutine diag(idiag,neqc,n)
c
c.... program to compute diagonal addresses of left-hand-side matrix
c
      dimension idiag(*)
c
      n = 1
      idiag(1) = 1 
      if (neqc.eq.1) return
c
      do 100 i=2,neqc
      idiag(i) = idiag(i) + idiag(i-1) + 1
  100 continue
      n = idiag(neqc)
c
      return
      end
c**** new ********************************************************************** 
      subroutine echo 
c 
c.... program to echo input data 
c 
      character*4 ia(20) 
      common /iounit/ iin,iout,iecho,ioupp,itest1,itest2 
c 
      read(iin,1000) iech 
      if (iech.eq.0) return 
c 
      write(iecho,2000) iech 
      backspace iin 
c 
      do 100 i=1,100000 
      read(iin,3000,end=200) ia 
      if (mod(i,50).eq.1) write(iecho,4000) 
      write(iecho,5000) ia 
  100 continue 
c 
  200 continue 
      rewind iin 
      read(iin,1000) iech 
c 
      return 
c 
 1000 format(16i10) 
 2000 format(///,' i n p u t   d a t a   f i l e               ',  //5x, 
     &' echo print code . . . . . . . . . . . . . . (iecho ) = ',i10//5x, 
     &'    eq. 0, no echo of input data                        ',   /5x, 
     &'    eq. 1, echo input data                              ',   ///) 
 3000 format(20a4) 
 4000 format(' ',8('123456789*'),//) 
 5000 format(' ',20a4) 
      end 
c**** new **********************************************************************
      subroutine elemnt(task,ngrp)
c
c.... program to calculate element task number
c
      character*8 task,eltask(4) 
      dimension ngrp(*)
      common /info  / iexec,iprtin,irank,nsd,numnp,ndof,nlvect,
     &                numeg,nmultpd,nmultpc,nedge

      common /iounit/ iin,iout,iecho,ioupp,itest1,itest2 
      character*4 ia 
      common na(1) 
      common /dictn/ ia(1) 
      data ntask,    eltask
     &    /    4,'input_el',
     &           'input_sd',
     &           'form_lrs',
     &           'pos_proc'/
c
      do 100 i=1,ntask
      if (task.eq.eltask(i)) itask = i
  100 continue
c
      do 200 neg=1,numeg
c
      if (itask.eq.1) then
         mpnpar = mpoint('npar    ',16   ,0,0,1)
         ngrp(neg) = mpnpar
         call elcard(na(mpnpar),neg)
      else
         mpnpar = ngrp(neg)
      endif
c
      ntype  = na(mpnpar)
      call elmlib(ntype,mpnpar,itask,neg)
  200 continue
c
      return
      end
c**** new **********************************************************************
      subroutine elcard(npar,neg)
c
c.... program to read element group control card
c
      dimension npar(*)
      common /iounit/ iin,iout,iecho,ioupp,itest1,itest2 
c
      read(iin,1000) (npar(i),i=1,16)
      write(iecho,2000) neg
c
      return
c
 1000 format(16i10)
 2000 format(//,' e l e m e n t   g r o u p   d a t a         ',  //5x,
     &' element group number  . . . . . . . . . . (neg   ) = ',i10/// )
c
      end
c***************************************************************************
c
c     Solvers
c
c**** newfrey *********************************************************
      subroutine addnsl(alhs,clhs,eleffm,idiag,lm,nee,ldiag)
c
c         program to add element left-hand-side matrix to          
c                global left-hand-side matrix                      
c                                                                  
c        ldiag = .true.,  add diagonal element matrix              
c                                                                  
c        ldiag = .false, then                                     
c        add full nonsymmetric element matrix                   
c                                                                  
c
      implicit real*8 (a-h,o-z)
c
c.... remove above card for single-precision operation
c
      logical ldiag
      dimension alhs(*),clhs(*),eleffm(nee,*),idiag(*),lm(*)
c
      if (ldiag) then
c
         do 100 j=1,nee
            k = iabs(lm(j))
            if (k.gt.0) then
               l = idiag(k)
               alhs(l) = alhs(l) + eleffm(j,j)
            endif
  100    continue
c
      else
c
         do 400 j=1,nee
            k = iabs(lm(j))
            if (k.gt.0) then
               do 200 i=1,nee
                  m = iabs(lm(i))
                  if (m.gt.0) then
                     if (k.gt.m) then
                        l = idiag(k) - k + m
                        alhs(l) = alhs(l) + eleffm(i,j)
                     else
                        l = idiag(m) - m + k
                        clhs(l) = clhs(l) + eleffm(i,j)
                     endif
                     if (k.eq.m) then
                        l = idiag(k)
                        alhs(l) = alhs(l) + eleffm(i,j)
                        clhs(l) = alhs(l)
                     endif
                  endif
  200          continue
            endif
  400    continue
c
      endif
c
      return
      end
c**** newfrey *********************************************************
      subroutine backns(a,c,b,idiag,neqc)
c
c.... program to perform forward reduction and back substitution
c
      implicit real*8 (a-h,o-z)
c
c.... deactivate above card(s) for single-precision operation
c
      dimension a(*),c(*),b(*),idiag(*)
c
c
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
c
c.... forward reduction
c
c
      jj = 0
c
      do 100 j=1,neqc
      jjlast = jj
      jj     = idiag(j)
      jcolht = jj - jjlast
      if (jcolht.gt.1) then
          b(j) = b(j) - coldot(c(jjlast+1),b(j-jcolht+1),jcolht-1)
      endif
  100 continue
c
c.... diagonal scaling
c
      do 200 j=1,neqc
      ajj = a(idiag(j))
c
c.... warning: diagonal scaling is not performed if ajj equals zero
c
      if (ajj.ne.zero) b(j) = b(j)/ajj
  200 continue
c
c.... back substitution
c
      if (neqc.eq.1) return
      jjnext = idiag(neqc)
c
      do 400 j=neqc,2,-1
      jj     = jjnext
      jjnext = idiag(j-1)
      jcolht = jj - jjnext
      if (jcolht.gt.1) then
         bj = b(j)
         istart = j - jcolht + 1
         jtemp  = jjnext - istart + 1
c
         do 300 i=istart,j-1
         b(i) = b(i) - a(jtemp+i)*bj
  300    continue
c
      endif
c
  400 continue
c
      return
      end
c**** newfrey *********************************************************
      subroutine factns(a,c,idiag,neqc)
c
c.... program to perform crout factorization: a = l * d * u
c
c        a(i):  coefficient matrix stored in compacted column form;
c               after factorization contains d and u
c
c        c(i):  non-symmetric lower triangular coefficient matrix stored in
c                compacted row form; after factorization contains l
c 
c
      implicit real*8 (a-h,o-z)
c
c.... deactivate above card(s) for single-precision operation
c
      dimension a(*),c(*),idiag(*)
c
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
c
      do i=1,neqc
        write(98,*) i, idiag(i)
      end do
      jj = 0
c
      do 300 j=1,neqc
c
      jjlast = jj
      jj     = idiag(j)
      jcolht = jj - jjlast
c
      if (jcolht.gt.2) then
c
c....... for column j and i.le.j-1, replace a(i,j) with d(i,i)*u(i,j)
c
         istart = j - jcolht + 2
         jm1    = j - 1
         ij     = jjlast + 2
         ii     = idiag(istart-1)
c
         do 100 i=istart,jm1
c
         iilast = ii
         ii     = idiag(i)
         icolht = ii - iilast
         length = min0(icolht-1,i - istart + 1)
         if (length.gt.0)  then
            a(ij) = a(ij) - coldot(a(ij-length),c(ii-length),length)
            c(ij) = c(ij) - coldot(c(ij-length),a(ii-length),length)
         endif
         ij = ij + 1
  100    continue
c
      endif
c
      if (jcolht.ge.2) then
c
c....... for column j and i.le.j-1, replace a(i,j) with u(i,j);
c           replace a(j,j) with d(j,j).
c
         jtemp = j - jj
c
         do 200 ij=jjlast+1,jj-1
c
         ii = idiag(jtemp + ij)
c
c....... warning: the following calculations are skipped 
c                 if a(ii) equals zero
c
         if (a(ii).ne.zero) then
             c(ij) = c(ij)/a(ii)
             a(jj) = a(jj) - c(ij)*a(ij)
             a(ij) = a(ij)/a(ii)
         endif
  200    continue
c
      endif
c
  300 continue
c
      return
      end
c**** new **********************************************************************
      subroutine factor(a,idiag,neqc)
c
c.... program to perform Crout factorization: a = u(transpose) * d * u
c
c        a(i):  coefficient matrix stored in compacted column form;
c               after factorization contains d and u
c
      implicit real*8 (a-h,o-z)
c
c.... remove above card for single-precision operation
c
      dimension a(*),idiag(*)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
c
      jj = 0
c
      do 300 j=1,neqc
c
      jjlast = jj
      jj     = idiag(j)
      jcolht = jj - jjlast
c
      if (jcolht.gt.2) then
c
c....... for column j and i.le.j-1, replace a(i,j) with d(i,i)*u(i,j)
c
         istart = j - jcolht + 2
         jm1    = j - 1
         ij     = jjlast + 2
         ii     = idiag(istart-1)
c
         do 100 i=istart,jm1
c
         iilast = ii
         ii     = idiag(i)
         icolht = ii - iilast
         jlngth = i - istart + 1
         length = min0(icolht-1,jlngth)
         if (length.gt.0) 
     &      a(ij) = a(ij) - coldot(a(ii-length),a(ij-length),length)
         ij = ij + 1
  100    continue
c
      endif
c
      if (jcolht.ge.2) then
c
c....... for column j and i.le.j-1, replace a(i,j) with u(i,j);
c           replace a(j,j) with d(j,j).
c
         jtemp = j - jj
c
         do 200 ij=jjlast+1,jj-1
c
         ii = idiag(jtemp + ij)
         if (a(ii).ne.zero) then
            temp  = a(ij)
            a(ij) = temp/a(ii)
            a(jj) = a(jj) - temp*a(ij)
         endif
  200    continue
c
      endif
c
  300 continue
c
      return
      end
c**** new **********************************************************************
      subroutine formlm (idc,ien,lm,ndof,ned,nen,numel)
c
c.... program to form lm array
c
      dimension idc(ndof,*),ien(nen,*),lm(ned,nen,*)
c
      do 300 k=1,numel
c
      do 200 j=1,nen
      node=ien(j,k)
c
      do 100 i=1,ndof
      lm(i,j,k) = idc(i,node)
  100 continue
c
  200 continue
c
  300 continue
c
      return
      end
c**** new **********************************************************************
      subroutine ftod(idc,d,f,ndof,numnp,nlvect)
c
c.... program to compute displacement boundary conditions
c
      implicit real*8 (a-h,o-z)
c
c.... remove above card for single-precision operation
c
      dimension idc(ndof,*),d(ndof,*),f(ndof,numnp,*)
c
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
c
      do 300 i=1,ndof
c
            do 200 j=1,numnp
c
            k = idc(i,j)
            if (k.gt.0) go to 200
            val = zero
                  do 100 lv=1,nlvect
                  val = val + f(i,j,lv)
100               continue
c
            d(i,j) = val
c
  200       continue
c
 300  continue
      return
      end
c******************************************************************************
      subroutine genelad(lado,nside)                                     
c                                                                       
c.... program to read and generate element node and material numbers    
c                                                                       
c         lado(nside,numel) = element node numbers                         
c         mat(numel) e,    = element material numbers                     
c         nen            = number of element nodes (le.27)              
c         n              = element number                               
c         ng             = generation parameter                         
c         nel(i)         = number of elements in direction i            
c         incel(i)       = element number increment for direction i     
c         inc(i)         = node number increment for direction i        
c                                                                       
      dimension lado(nside,*),itemp(27)                             
      common /iounit/ iin,iout,iecho,ioupp,itest1,itest2 
      common /genelc/ n,nel(3),incel(3),inc(3)                          
c
  100 continue
      read(iin,1000) n,(itemp(i),i=1,nside),ng  
      if (n.eq.0) return  
c
	do no =1,nside
	  lado(no,n) = itemp(no) 
	end do
c
      if (ng.ne.0) then
c                                                                       
c....... generate data                                                     
c                                                                       
         read(iin,1000) (nel(i),incel(i),inc(i),i=1,3)  
         call genelad1(lado,nside)                                          
      endif
c
      go to 100                                                         
c
c                                                                       
 1000 format(20i10)                                             
c                                                                       
      end                                                               
c******************************************************************************
c******************************************************************************
      subroutine genelpar(ipar,ien,lado,idside,iedge,indedg,idlsd,
     &              nen,nside,nodsp,numel,npars)      
c                                                                       
c.... gera numeração dos parametros dos multiplicadores poe elemento   
c                                                                       
c         ipar(npars*nside,numel) = element parameter numbers 
c         nside          = number of element sides                        
c         nodsp          = number of element parameters (le.27)              
c         n              = element number                               
c         ng             = generation parameter                         
c         nel(i)         = number of elements in direction i            
c         incel(i)       = element number increment for direction i     
c         inc(i)         = node number increment for direction i        
c                                                                       
      dimension lado(nside,*),ipar(nodsp,*),iedge(npars,*),indedg(*) 
      dimension ien(nen,*),idside(nside,*),idlsd(*)   
      common /iounit/ iin,iout,iecho,ioupp,itest1,itest2 
c
c   	
      do 2002 nel=1,numel

c
c     nós dos vertices   
c
      do 1230 ns=1,nside
      ns1=idside(ns,1)
      ns2=idside(ns,2)
      nl1=ien(ns1,nel)
      nl2=ien(ns2,nel)
c      
      if(nl2.gt.nl1) then 
        idlsd(1) = 1
        idlsd(2) = 2
	  id3 = npars-2
	  if(id3.gt.0) then
          do  nn=3,npars
            idlsd(nn) = nn
	    end do
	  end if 
	else
	  idlsd(1) = 2
	  idlsd(2) = 1
	  id3 = npars-2
	  if(id3.gt.0) then
	    do il=3,npars
	      idlsd(il) = npars + 3 - il
	    end do 
	  end if
      end if 
c

	  la = lado(ns,nel)
	  kedg = indedg(la)
c
 	     neds = idlsd(1)  
           lada = iedge(neds,la)
           nsp = idside(ns,1)
	     ipar(nsp,nel) = lada 
c
 1230  continue


c
c    nós sobre os lados      
c 
      if(npars.gt.2) then         
	do 2001 ns=1,nside
c
      ns1=idside(ns,1)
      ns2=idside(ns,2)
      nl1=ien(ns1,nel)
      nl2=ien(ns2,nel)
c      
      if(nl2.gt.nl1) then 
        idlsd(1) = 1
        idlsd(2) = 2
	  id3 = npars-2
	  if(id3.gt.0) then
          do  nn=3,npars
            idlsd(nn) = nn
	    end do
	  end if 
	else
	  idlsd(1) = 2
	  idlsd(2) = 1
	  id3 = npars-2
	  if(id3.gt.0) then
	    do il=3,npars
	      idlsd(il) = npars + 3 - il
	    end do 
	  end if
      end if 
c
	  la = lado(ns,nel)
	  kedg = indedg(la)
c
       do np=3,npars
 	     neds = idlsd(np)  
           lada = iedge(neds,la)
           nsp = idside(ns,np)
	     ipar(nsp,nel) = lada 
	 end do
c
 2001  continue
       end if
 2002  continue 
c
       return
c                                                                      
c                                                                       
      end                                                               
c******************************************************************************
c******************************************************************************
      subroutine genelad1(ipar,nodsp)                                    
c                                                                       
c.... program to generate element node and material numbers             
c                                                                       
      dimension ipar(nodsp,*)                                      
      common /genelc/ n,nel(3),incel(3),inc(3)                          
c
c                                                                       
c.... set defaults                                                      
c                                                                       
      call geneld                                                       
c                                                                       
c.... generation algorithm                                              
c                                                                       
      ie = n                                                            
      je = n                                                            
      ke = n                                                            
c                                                                      
      ii = nel(1)                                                       
      jj = nel(2)                                                       
      kk = nel(3)                                                       
c   
c
      do 300 k=1,kk                                                     
c
      do 200 j=1,jj                                                     
c
      do 100 i=1,ii                                                     
c                                                                       
      if (i.ne.ii) then
         le = ie                                                           
         ie = le + incel(1)                                                
       call genelip(ipar(1,ie),ipar(1,le),inc(1),nodsp)                       
      endif
  100 continue                                                          
c                                                                       
      if (j.ne.jj) then
         le = je                                                           
         je = le + incel(2)                                                
       call genelip(ipar(1,je),ipar(1,le),inc(2),nodsp)                       
         ie = je                                                           
      endif
  200 continue                                                          
c                                                                       
      if (k.ne.kk) then
         le = ke                                                           
         ke = le + incel(3)                                                
       call genelip(ipar(1,ke),ipar(1,le),inc(3),nodsp)                       
         ie = ke                                                           
      endif
  300 continue                                                          
c                                                                       
      return                                                            
      end                                                               
c******************************************************************************
      subroutine genelip(ien2,ien1,inc,nodsp)                              
c                                                                       
c.... program to increment element node numbers                         
c                                                                       
      dimension ien1(*),ien2(*)                                         
c
c                                                                       
      do 100 i=1,nodsp                                                    
      if (ien1(i).eq.0) then
         ien2(i) = 0
      else
         ien2(i) = ien1(i) + inc                         
      endif
  100 continue                                                          
c                                                                       
      return                                                            
      end                                                               
c******************************************************************************
      subroutine genel(ien,mat,nen)                                     
c                                                                       
c.... program to read and generate element node and material numbers    
c                                                                       
c         ien(nen,numel) = element node numbers                         
c         mat(numel)     = element material numbers                     
c         nen            = number of element nodes (le.27)              
c         n              = element number                               
c         ng             = generation parameter                         
c         nel(i)         = number of elements in direction i            
c         incel(i)       = element number increment for direction i     
c         inc(i)         = node number increment for direction i        
c                                                                       
      dimension ien(nen,*),mat(*),itemp(27)                             
      common /genelc/ n,nel(3),incel(3),inc(3)                          
      common /iounit/ iin,iout,iecho,ioupp,itest1,itest2 
c
  100 continue
      read(iin,1000) n,m,(itemp(i),i=1,nen),ng  
      if (n.eq.0) return
      call imove(ien(1,n),itemp,nen)                                    
      mat(n)=m                                                          
      if (ng.ne.0) then
c                                                                       
c....... generate data                                                     
c                                                                       
       read(iin,1000) (nel(i),incel(i),inc(i),i=1,3)                     
       call genel1(ien,mat,nen) 
      endif
      go to 100                                                         
c                                                                       
 1000 format(20i10)                                             
c                                                                       
      end                                                               
c******************************************************************************
      subroutine genel1(ien,mat,nen)                                    
c                                                                       
c.... program to generate element node and material numbers             
c                                                                       
      dimension ien(nen,*),mat(*)                                       
      common /genelc/ n,nel(3),incel(3),inc(3)                          
c
 1000 format(16i10)                                             
c                                                                       
c.... set defaults                                                      
c                                                                       
      call geneld                                                       
c                                                                       
c.... generation algorithm                                              
c
      ie = n                                                            
      je = n                                                            
      ke = n                                                            
c                                                                      
      ii = nel(1)                                                       
      jj = nel(2)                                                       
      kk = nel(3)                                                       
c
c                                                                       
      do 300 k=1,kk                                                     
c
      do 200 j=1,jj                                                     
c
      do 100 i=1,ii                                                     
c
c                                                                       
      if (i.ne.ii) then
         le = ie            	                                                 
         ie = le + incel(1)                                                
         call geneli(ien(1,ie),ien(1,le),inc(1),nen)  
         mat(ie) = mat(le)                                                 
      endif
  100 continue                                                          
c                                                                       
      if (j.ne.jj) then
         le = je                                                           
         je = le + incel(2)                                                
         call geneli(ien(1,je),ien(1,le),inc(2),nen)                       
         mat(je) = mat(le)                                                 
         ie = je                                                           
      endif
  200 continue                                                          
c                                                                             
      if (k.ne.kk) then
         le = ke                                                           
         ke = le + incel(3)                                                
         call geneli(ien(1,ke),ien(1,le),inc(3),nen)                       
         mat(ke) = mat(le)                                                 
         ie = ke                                                           
      endif
  300 continue                                                          
c                                                                       
      return                                                            
      end                                                               
c******************************************************************************
      subroutine geneld                                                 
c                                                                       
c.... program to set defaults for element node       
c        and material number generation                              
c                                                                       
      common /genelc/ n,nel(3),incel(3),inc(3)                          
c                                                                       
      if (nel(1).eq.0) nel(1) = 1                                       
      if (nel(2).eq.0) nel(2) = 1                                       
      if (nel(3).eq.0) nel(3) = 1                                       
c                                                                       
      if (incel(1).eq.0) incel(1) = 1                                   
      if (incel(2).eq.0) incel(2) = nel(1)                              
      if (incel(3).eq.0) incel(3) = nel(1)*nel(2)                       
c                                                                       
      if (inc(1).eq.0) inc(1) = 1                                       
      if (inc(2).eq.0) inc(2) = (1+nel(1))*inc(1)                       
      if (inc(3).eq.0) inc(3) = (1+nel(2))*inc(2)                       
c                                                                       
      return                                                            
      end                                                               
c******************************************************************************
      subroutine geneli(ien2,ien1,inc,nen)                              
c                                                                       
c.... program to increment element node numbers                         
c                                                                       
      dimension ien1(*),ien2(*)                                         
c                                                                       
      do 100 i=1,nen                                                    
      if (ien1(i).eq.0) then
         ien2(i) = 0
      else
         ien2(i) = ien1(i) + inc                         
      endif
  100 continue                                                          
c                                                                       
      return                                                            
      end                                                               
c******************************************************************************
      subroutine genfl(a,nra)                                           
c                                                                       
c.... program to read and generate floating-point nodal data            
c                                                                       
c         a       = input array                                         
c         nra     = number of rows in a (le.6)                          
c         n       = node number                                         
c         numgp   = number of generation points                         
c         ninc(i) = number of increments for direction i                
c         inc(i)  = increment for direction i                           
c                                                                       
      implicit real*8(a-h,o-z)                                          
c                                                                       
c.... remove above card for single-precision operation               
c                                                                       
      dimension a(nra,*)                                                
      common /iounit/ iin,iout,iecho,ioupp,itest1,itest2 
      common /genflc/ temp(6,20),n,numgp,ninc(3),inc(3)                 
c                                                                       
  100 continue                                                          
      read(iin,1000) n,numgp,(temp(i,1),i=1,nra)                        
      if (n.eq.0) return                                                
      call move(a(1,n),temp,nra)                                       
      if (numgp.ne.0) then
         do 200 j=2,numgp                                                  
c                                                                       
         read(iin,1000) m,mgen,(temp(i,j),i=1,nra)                         
         if (mgen.ne.0) call move(temp(1,j),a(1,m),nra) 
c                                                                       
  200    continue                                                          
         read(iin,2000) (ninc(i),inc(i),i=1,3)                             
         call genfl1(a,nra)                                                
      endif
      go to 100                                                         
c                                                                       
 1000 format(2i10,6f10.0)                                                
 2000 format(16i10)                                                      
c                                                                       
      end                                                               
c******************************************************************************
      subroutine genfl1(a,nra)                                          
c                                                                       
c.... program to generate floating-point nodal data 
c        via isoparametric interpolation         
c                                                                       
c         iopt = 1, generation along a line                             
c              = 2, generation over a surface                           
c              = 3, generation within a volume                            
c                                                                       
      implicit real*8(a-h,o-z)                                          
c                                                                       
c.... remove above card for single-precision operation                  
c                                                                       
      dimension a(nra,*),sh(20)                                         
      common /genflc/ temp(6,20),n,numgp,ninc(3),inc(3)                 
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
c
      iopt = 3                                                          
      if (ninc(3).eq.0) iopt = 2                                        
      if (ninc(2).eq.0) iopt = 1                                        
c                                                                       
      dr = zero                                         
      ds = zero                                                        
      dt = zero                                   
c                                                                       
      if (ninc(1).ne.0) dr = two/ninc(1)                                
      if (ninc(2).ne.0) ds = two/ninc(2)                                
      if (ninc(3).ne.0) dt = two/ninc(3)                                
c                                                                       
      ii = ninc(1)+1                                                    
      jj = ninc(2)+1                                                    
      kk = ninc(3)+1                                                    
c                                                                       
      ni = n                                                            
      nj = n                                                            
      nk = n                                                            
c                                                                       
      t = -one                                                          
      do 300 k=1,kk                                                     
c
      s = -one                                                          
      do 200 j=1,jj                                                     
c
      r = -one                                                          
      do 100 i=1,ii                                                     
c                                                                       
      call gensh(r,s,t,sh,numgp,iopt)                                   
      call multab(temp,sh,a(1,ni),6,20,nra,numgp,nra,1,1)               
      ni = ni + inc(1)                                                      
      r = r + dr                                                            
  100 continue                                                          
c                                                                       
      nj = nj + inc(2)                                                      
      ni = nj                                                             
      s = s + ds                                                            
  200 continue                                                          
c                                                                       
      nk = nk + inc(3)                                                      
      ni = nk                                                             
      t = t + dt                                                            
  300 continue                                                          
c                                                                       
      return                                                            
      end                                                               
c****************************************************************************** 
    
      subroutine gensh(r,s,t,sh,numgp,iopt)                             
c                                                                       
c.... program to call shape function routines         
c        for isoparametric generation         
c                                                                       
      implicit real*8(a-h,o-z)                                      
c                                                                       
c.... modify above card for single-precision operation               
c                                                                       
      dimension sh(*)                                                   
c                                                                       
      go to (100,200,300),iopt                                                
c                                                                       
  100 call gensh1(r,sh,numgp)                                           
      return                                                            
c                                                                       
  200 call gensh2(r,s,sh,numgp)                                         
      return                                                            
c                                                                       
  300 call gensh3(r,s,t,sh,numgp)                                       
      return                                                            
c                                                                       
      end                                                               
c******************************************************************************
      subroutine gensh1(r,sh,n)                                         
c                                                                       
c.... program to compute 1d shape functions           
c        for isoparametric generation                     
c                                                                       
      implicit real*8(a-h,o-z)                                          
c                                                                       
c.... modify above card(s) for single-precision operation               
c                                                                       
      dimension sh(*)                                                   
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
c                                                                       
      sh(2) = pt5*r                                                       
      sh(1) = pt5 - sh(2)                                                   
      sh(2) = pt5 + sh(2)                                                   
      if (n.eq.3) then
         sh(3) = one - r*r                                                     
         sh(1) = sh(1) - pt5*sh(3)                                             
         sh(2) = sh(2) - pt5*sh(3)                                             
      endif
c                                                                       
      return                                                            
      end                                                               
c******************************************************************************
      subroutine gensh2(r,s,sh,n)                                       
c                                                                       
c.... program to compute 2d shape functions 
c        for isoparametric generation    
c                                                                       
      implicit real*8(a-h,o-z)                                          
c                                                                       
c.... modify above card for single-precision operation               
c                                                                       
      dimension sh(*)                                      
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
c
      r2 = pt5*r                                                          
      r1 = pt5 - r2                                                         
      r2 = pt5 + r2                                                         
      s2 = pt5*s                                                          
      s1 = pt5 - s2                                                         
      s2 = pt5 + s2                                                         
      sh(1) = r1*s1                                                       
      sh(2) = r2*s1                                                       
      sh(3) = r2*s2                                                       
      sh(4) = r1*s2                                                       
      if (n.eq.4) return                                                
c                                                                       
      r3 = one - r*r                                                        
      s3 = one - s*s                                                        
      sh(5) = r3*s1                                                       
      sh(6) = s3*r2                                                       
      sh(7) = r3*s2                                                       
      sh(8) = s3*r1                                                       
      sh(1) = sh(1) - pt5*(sh(5) + sh(8))
      sh(2) = sh(2) - pt5*(sh(6) + sh(5))
      sh(3) = sh(3) - pt5*(sh(7) + sh(6))
      sh(4) = sh(4) - pt5*(sh(8) + sh(7))
c                                                                       
      return                                                            
      end                                                               
c******************************************************************************
      subroutine gensh3(r,s,t,sh,n)                                     
c                                                                       
c.... program to compute 3d shape functions            
c        for isoparametric generation   
c                                                                       
      implicit real*8(a-h,o-z)                                          
c                                                                       
c.... modify above card for single-precision operation               
c                                                                       
      dimension sh(*)                                                   
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
c                                                                       
      r2 = pt5*r
      r1 = pt5 - r2                                                         
      r2 = pt5 + r2                                                         
      s2 = pt5*s
      s1 = pt5 - s2                                                         
      s2 = pt5 + s2                                                         
      t2 = pt5*t
      t1 = pt5 - t2                                                         
      t2 = pt5 + t2                                                         
c                                                                       
      rs1 = r1*s1                                                         
      rs2 = r2*s1                                                         
      rs3 = r2*s2                                                         
      rs4 = r1*s2                                                         
      sh(1) = rs1*t1                                                      
      sh(2) = rs2*t1                                                      
      sh(3) = rs3*t1                                                      
      sh(4) = rs4*t1                                                      
      sh(5) = rs1*t2                                                      
      sh(6) = rs2*t2                                                      
      sh(7) = rs3*t2                                                      
      sh(8) = rs4*t2                                                      
      if (n.eq.8) return                                                 
c                                                                       
      r3 = one - r*r                                                        
      s3 = one - s*s                                                        
      t3 = one - t*t                                                        
      sh(17) = t3*rs1                                                     
      sh(18) = t3*rs2                                                     
      sh(19) = t3*rs3                                                     
      sh(20) = t3*rs4                                                     
      rs1 = r3*s1                                                         
      rs2 = s3*r2                                                         
      rs3 = r3*s2                                                         
      rs4 = s3*r1                                                         
      sh( 9) = rs1*t1                                                     
      sh(10) = rs2*t1                                                     
      sh(11) = rs3*t1                                                     
      sh(12) = rs4*t1                                                     
      sh(13) = rs1*t2                                                     
      sh(14) = rs2*t2                                                     
      sh(15) = rs3*t2                                                     
      sh(16) = rs4*t2                                                     
c                                                                       
      sh(1) = sh(1) - pt5*(sh( 9) + sh(12) + sh(17))
      sh(2) = sh(2) - pt5*(sh( 9) + sh(10) + sh(18))
      sh(3) = sh(3) - pt5*(sh(10) + sh(11) + sh(19))
      sh(4) = sh(4) - pt5*(sh(11) + sh(12) + sh(20))
      sh(5) = sh(5) - pt5*(sh(13) + sh(16) + sh(17))
      sh(6) = sh(6) - pt5*(sh(13) + sh(14) + sh(18))
      sh(7) = sh(7) - pt5*(sh(14) + sh(15) + sh(19))
      sh(8) = sh(8) - pt5*(sh(15) + sh(16) + sh(20))
c                                                                       
      return                                                            
      end                                                               
c**** new **********************************************************************
      subroutine iclear(ia,m)
c
c.... program to clear an integer array
c
      dimension ia(*)
c
      do 100 i=1,m
      ia(i) = 0
  100 continue
c
      return
      end
c**** new **********************************************************************
      subroutine igen(ia,m)
c
c.... program to read and generate integer nodal data
c
c        ia = input array
c         m = number of rows in ia
c         n = node number
c        ne = end node in generation sequence
c        ng = generation increment
c
      dimension ia(m,*),ib(13)
      common /iounit/ iin,iout,iecho,ioupp,itest1,itest2 
c
  100 continue
      read(iin,1000) n,ne,ng,(ib(i),i=1,m)
      if (n.eq.0) return
      if (ng.eq.0) then
         ne = n
         ng = 1
      else
         ne = ne - mod(ne-n,ng)
      endif
c
      do 200 i=n,ne,ng
      call imove(ia(1,i),ib,m)
  200 continue
c
      go to 100
c
 1000 format(16i10)
      end
c**** new *********************************************************************
      subroutine imove(ia,ib,n)
c
c.... program to move an integer array 
c
      dimension ia(*),ib(*)
c
      do 100 i=1,n
      ia(i)=ib(i)
  100 continue

      return
      end
c**** new **********************************************************************
      subroutine input(f,ndof,numnp,j,nlvect,iprtin)
c
c.... program to read, generate and write nodal input data
c
c        f(ndof,numnp,nlvect) = prescribed forces/kinematic data (j=0)
c                             = nodal body forces(j=1)
c
      implicit real*8 (a-h,o-z)
c
c.... remove above card for single-precision operation
c
      logical lzero
      dimension f(ndof,numnp,*)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
      common /iounit/ iin,iout,iecho,ioupp,itest1,itest2 
c
      call clear(f,nlvect*numnp*ndof)
c
      do 100 nlv=1,nlvect
      call genfl(f(1,1,nlv),ndof)
      call ztest(f(1,1,nlv),ndof*numnp,lzero)
c
      if (iprtin.eq.0) then
c
         if (lzero) then
            if (j.eq.0) write(iecho,1000) nlv
            if (j.eq.1) write(iecho,2000)
         else
            if (j.eq.0) call printf(f,ndof,numnp,nlv)
c
            if (j.eq.1) 
     &      call printd(' n o d a l  b o d y  f o r c e s  ',
     &                  f,ndof,numnp,iecho)
c
         endif
      endif
c
  100 continue
c
      return
 1000 format(/////,' there are no nonzero prescribed forces and ',
     &    'kinematic boundary conditions for load vector number ',i10)
 2000 format(/////,' there are no nonzero nodal body forces')
      end
c**** new **********************************************************************
      subroutine interp(x,y,xx,yy,n)
c
c.... program to perform linear interpolation
c
c        x(i) = abscissas
c        y(i) = ordinates
c          xx = input abscissa
c          yy = output ordinate
c           n = total number of data points (1.le.i.le.n)
c
      implicit real*8 (a-h,o-z)
c
c.... remove above card for single-precision operation
c
      dimension x(*),y(*)
c
      if (xx.le.x(1)) then
         yy = y(1)
         return
      endif
c
      if (xx.ge.x(n)) then
         yy = y(n)
         return
      endif
c
      do 100 i=1,n
      if (x(i).ge.xx) then
         yy = y(i-1) + (xx - x(i-1))*(y(i) - y(i-1))/(x(i) - x(i-1))
         return
      endif
  100 continue
c
      end
c**** new **********************************************************************
      subroutine kdbc(eleffm,elresf,dl,nee)
c
c.... program to adjust load vector for prescribed displacement
c     boundary condition
c
      implicit real*8 (a-h,o-z)
c
c.... remove above card for single-precision operation
c
      dimension eleffm(nee,*),elresf(*),dl(*)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
c
      do 200 j=1,nee
c
      val=dl(j)
      if(val.eq.zero) go to 200
c
      do 100 i=1,nee
      elresf(i)=elresf(i)-eleffm(i,j)*val
100   continue
c
200   continue
c
      return
      end
c**** new **********************************************************************
      subroutine load(idc,f,brhs,ndof,numnp,nlvect)
c
c.... program to accumulate nodal forces and transfer into
c        right-hand-side vector
c
      implicit real*8 (a-h,o-z)
c
c.... remove above card for single-precision operation
c
      dimension idc(ndof,*),f(ndof,numnp,*),brhs(*)
c
      do 300 i=1,ndof
c
      do 200 j=1,numnp
      k = idc(i,j)
      if (k.gt.0) then
c
         do 100 nlv=1,nlvect
         brhs(k) = brhs(k) + f(i,j,nlv)
  100    continue
c
      endif
c
  200 continue
c
  300 continue
c
      return
      end
c**** new **********************************************************************
      subroutine local(ien,x,xl,nen,nrowx,nrowxl)
c
c.... program to localize a global array
c
c        note: it is assumed nrowxl.le.nrowx
c
      implicit real*8 (a-h,o-z)
c
c.... remove above card for single-precision operation
c
      dimension ien(*),x(nrowx,*),xl(nrowxl,*)
c
      do 200 j=1,nen
      node = ien(j)
c
      do 100 i=1,nrowxl
      xl(i,j)= x(i,node)
  100 continue
c
  200 continue
c
      return
      end
c**** new **********************************************************************
      function lout(i,j)
c
c.... program to determine logical switch
c
      logical lout
c
      lout = .false.
      if (j.eq.0) return
      if (mod(i,j).eq.0) lout = .true.
c
      return
      end
c**** new **********************************************************************
      subroutine matadd(a,b,c,ma,mb,mc,m,n,iopt)
c
c.... program to add rectangular matrices
c
      implicit real*8 (a-h,o-z)
c
c.... remove above card for single-precision operation
c
      dimension a(ma,*),b(mb,*),c(mc,*)
c
      go to (1000,2000,3000),iopt
c
c.... iopt = 1, add entire matrices
c
 1000 do 1200 j=1,n
c
      do 1100 i=1,m 
      c(i,j) = a(i,j) + b(i,j)
 1100 continue
c
 1200 continue
      return
c
c.... iopt = 2, add lower triangular and diagonal elements
c
 2000 do 2200 j=1,n
c
      do 2100 i=j,m 
      c(i,j) = a(i,j) + b(i,j)
 2100 continue
c
 2200 continue
      return
c
c.... iopt = 3, add upper triangular and diagonal elements
c
 3000 do 3200 j=1,n
c
      do 3100 i=1,j 
      c(i,j) = a(i,j) + b(i,j)
 3100 continue
c
 3200 continue
      return
c
      end
c**** new **********************************************************************
      subroutine minmax(x,xmax,xmin,l,m,n)
c
c.... program to compute the min and max in the row of a matrix
c
c        x = matrix
c        l = number of rows in x
c        m = number of columns in x
c        n = row number
c
      dimension x(l,*)
c
      xmax = x(n,1)
      xmin = x(n,1)
c
      do 100 i = 2,m
        if (x(n,i).gt.xmax) xmax = x(n,i)
        if (x(n,i).lt.xmin) xmin = x(n,i)
  100 continue
c
      return
      end
c**** new **********************************************************************
      subroutine move(a,b,n)
c
c.... program to move a floating-point array
c
      implicit real*8(a-h,o-z)
c
c.... remove above card for single-precision operation
c
      dimension a(*),b(*)
c
      do 100 i=1,n
      a(i) = b(i)
  100 continue
c
      return
      end
c**** new ********************************************************************** 
      function mpoint(name,ndim1,ndim2,ndim3,ipr) 
c 
c.... program to calculate storage pointer 
c 
      character*4 name(2) 
      common /bpoint/ mfirst,mlast,ilast,mtot,iprec
c 
      mpoint = mfirst 
      if ( iprec.eq.2 .and. mod(mpoint,2).eq.0 ) mpoint = mpoint + 1 
      call dctnry(name,ndim1,ndim2,ndim3,mpoint,ipr,mlast,ilast) 
      mfirst = mpoint + ndim1*max0(1,ndim2)*max0(1,ndim3)*ipr 
      if (mfirst.ge.mlast) call serror(name,mfirst-mlast) 
c 
      return 
      end 
c**** new **********************************************************************
      subroutine multab(a,b,c,ma,mb,mc,l,m,n,iopt)
c
c.... program to multiply two matrices
c
c        l = range of dot-product index
c        m = number of active rows in c
c        n = number of active columns in c
c
      implicit real*8 (a-h,o-z)
c                                                                       
c.... remove above card for single-precision operation               
c                                                                       
      dimension a(ma,*),b(mb,*),c(mc,*)
c
      go to (1000,2000,3000,4000),iopt
c
c.... iopt = 1, c(i,j) = a(i,k)*b(k,j) , (c = a * b)
c
 1000 do 1200 i=1,m
c
      do 1100 j=1,n
      c(i,j) = rcdot(a(i,1),b(1,j),ma,l)
 1100 continue
c
 1200 continue
      return
c                                            t
c.... iopt = 2, c(i,j) = a(k,i)*b(k,j) (c = a  * b)
c
 2000 do 2200 i=1,m
c
      do 2100 j=1,n
      c(i,j) = coldot(a(1,i),b(1,j),l)
 2100 continue
c
 2200 continue
      return
c                                                t
c.... iopt = 3, c(i,j) = a(i,k)*b(j,k) (c = a * b )
c
 3000 do 3200 i=1,m
c
      do 3100 j=1,n
      c(i,j) = rowdot(a(i,1),b(j,1),ma,mb,l)
 3100 continue
c
 3200 continue
      return
c                                            t    t
c.... iopt = 4, c(i,j) = a(k,i)*b(j,k) (c = a  * b )
c
 4000 do 4200 i=1,m
c
      do 4100 j=1,n
      c(i,j) = rcdot(b(j,1),a(1,i),mb,l)
 4100 continue
c
 4200 continue
c
      return
      end
c**** new **********************************************************************
      subroutine pivots(a,idiag,neqc,nsq,*)
c
c.... program to determine the number of zero and negative terms in
c        array d of factorization a = u(transpose) * d * u
c
      implicit real*8 (a-h,o-z)
c                                                                       
c.... remove above card for single-precision operation               
c                                                                       
      dimension a(*),idiag(*)
      common /iounit/ iin,iout,iecho,ioupp,itest1,itest2 
c
      iz = 0
      in = 0
c
      do 100 n=1,neqc
      i = idiag(n)
      if (a(i).eq.0.) iz = iz + 1
      if (a(i).lt.0.) in = in + 1
  100 continue
c
      write(iecho,1000) nsq,iz,in
c
      return 1
c
 1000 format(' ',
     &' zero and/or negative pivots encountered               ', ///5x,
     &' time sequence number   . . . . . . . . . . (nsq  ) = ',i10//5x,
     &' number of zeroes . . . . . . . . . . . . . . . . . . = ',i10//5x,
     &' number of negatives  . . . . . . . . . . . . . . . = ',i10//5x)
c
      end

c**** new **********************************************************************
      subroutine princ(n,s,p)
c
c.... program to compute principal values of symmetric 2nd-rank tensor
c     
c        s = symmetric second-rank tensor stored as a vector
c        n = number of dimensions (2 or 3)
c        p = vector of principal values 
c
c.... the components of s must be stored in the following orders
c
c        2-d problems: s11,s22,s12
c        3-d problems: s11,s22,s33,s12,s23,s31
c
      implicit real*8 (a-h,o-z)
c
c.... remove above card for single precision operation
c
      dimension s(*),p(*)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
      data rt2/1.41421356237309d0/,pi23/2.09439510239321d0/
c
      if (n.eq.2) then
c
c....... 2-d problem
c
         a = 22.5d0/datan(one)
         x = pt5*(s(1) + s(2))
         y = pt5*(s(1) - s(2))
         r = dsqrt(y*y + s(3)*s(3))
         p(1) = x + r
         p(2) = x - r
         p(3) = r
         p(4) = 45.0d0
         if (y.ne.zero.or.s(3).ne.zero) p(4) = a*atan2(s(3),y)
      endif
c
      if (n.eq.3) then
c
c....... 3-d problem
c
  100    r = zero
         x = (s(1) + s(2) + s(3))/three
         y = s(1)*(s(2) + s(3)) + s(2)*s(3)
     &       - s(4)*s(4) - s(6)*s(6) - s(5)*s(5) 
         z = s(1)*s(2)*s(3) - two*s(4)*s(6)*s(5) - s(1)*s(5)*s(5)
     &       - s(2)*s(6)*s(6) - s(3)*s(4)*s(4)
         t = three*x*x - y
         u = zero
         if (t.ne.zero) then
            u = dsqrt(two*t/three)
            a = (z + (t - x*x)*x)*rt2/u**3
            r = dsqrt(dabs(one - a*a))
            r = datan2(r,a)/three
         endif
         p(1) = x + u*rt2*cos(r)
         p(2) = x + u*rt2*cos(r - pi23)
         p(3) = x + u*rt2*cos(r + pi23)
      endif
c
      return
      end
c**** new **********************************************************************
      subroutine printd(name,dva,ndof,numnp,icode)
c
c
c.... program to print kinematic data
c
      implicit real*8(a-h,o-z)
c
c.... remove above card for single precision operation
c
      logical lzero
	character *4 name(11)
      dimension dva(ndof,*)
c
      nn = 0
c
      do 100 n=1,numnp
      call ztest(dva(1,n),ndof,lzero)
      if (.not.lzero) then
         nn = nn + 1
         if (mod(nn,50).eq.1) 
     &      write(icode,1000) name,(i,i=1,ndof)
         write(icode,2000) n,(dva(i,n),i=1,ndof)
      endif
  100 continue
c
      return
c
 1000 format(///,11a4//6xx,'node',6(11x,'dof',i1)/)
 2000 format(1x,i10,2x,6(1pe13.6,2x))
      end
c**** new **********************************************************************
      subroutine printf(f,ndof,numnp,nlv)
c
c.... program to print prescribed force and boundary condition data
c
      implicit real*8(a-h,o-z)
c
c.... remove above card for single precision operation
c
      logical lzero
      dimension f(ndof,numnp,*)
      common /iounit/ iin,iout,iecho,ioupp,itest1,itest2 
c
      nn = 0
c
      do 100 n=1,numnp
      call ztest(f(1,n,nlv),ndof,lzero)
      if (.not.lzero) then
         nn = nn + 1
         if (mod(nn,50).eq.1) 
     &      write(iecho,1000) nlv,(i,i=1,ndof)
         write(iecho,2000) n,(f(i,n,nlv),i=1,ndof)
      endif
  100 continue
c
      return
c
 1000 format(///,
     &' p r e s c r i b e d   f o r c e s   a n d   k i n e m a t i c ',
     &'  b o u n d a r y   c o n d i t i o n s'//5x,
     &' load vector number = ',i10///5x,
     &' node no.',6(13x,'dof',i1,:)/)
 2000 format(6x,i10,10x,6(1pe15.8,2x))
      end
c**** new **********************************************************************
      subroutine printp(a,idiag,neqc,nsq,*)
c
c.... program to print array d after Crout factorization 
c        a = u(transpose) * d * u
c
      implicit real*8(a-h,o-z)
c
c.... remove above card for single precision operation
c
      dimension a(*),idiag(*)
      common /iounit/ iin,iout,iecho,ioupp,itest1,itest2 
c
      do 100 n=1,neqc
      if (mod(n,50).eq.1) write(iecho,1000) nsq
      write(iecho,1000)
      i = idiag(n)
      write(iecho,2000) n,a(i)
  100 continue
c
      return 1
c
 1000 format(///,' array d of factorization',/
     &' a = u(transpose) * d * u ',                               //5x,
     &' time sequence number   . . . . . . . . . . . (nsq) = ',i10//5x)
 2000 format(1x,i10,4x,1pe20.8)
      end
c**** new **********************************************************************
      subroutine prntels(mat,ien,nen,numel)
c
c.... program to print data for element with "nen" nodes
c
c        note: presently the label formats are limited to
c              elements with one to nine nodes
c
      dimension mat(*),ien(nen,*)
      common /iounit/ iin,iout,iecho,ioupp,itest1,itest2 
c
      do 100 n=1,numel
      if (mod(n,50).eq.1) write(iecho,1000) (i,i=1,nen)
      write(iecho,2000) n,mat(n),(ien(i,n),i=1,nen)
  100 continue
c
      return
c
 1000 format(///,
     &' d a t a   f o r    e l e m e n t   s i d e s ',//1x,
     &20x, ' g e o m e t r i c    d  a t a ',//1x,
     &' element      material',16('      side',i2))
 2000 format(1x,i10,20(2x,i10))
      end
c**** new **********************************************************************
      subroutine rnprntels(ien,nen,numel)
c
c.... program to print data for element with "nen" nodes
c
c        note: presently the label formats are limited to
c              elements with one to nine nodes
c
      dimension ien(nen,*)
      common /iounit/ iin,iout,iecho,ioupp,itest1,itest2 
c
      do 100 n=1,numel
      if (mod(n,50).eq.1) write(iecho,1000) (i,i=1,nen)
      write(iecho,2000) n,(ien(i,n),i=1,nen)
  100 continue
c
      return
c
 1000 format(///,
     &' d a t a   f o r    e l e m e n t   s i d e s ',//1x,
     7  20x, 'r e n u m b e r e d ',//1x,
     &' element    ',16('      side',i2))
 2000 format(1x,i10,20(2x,i10))
      end
c**** new **********************************************************************
      subroutine prntelp(mat,ien,nen,numel)
c
c.... program to print data for element with "nen" nodes
c
c        note: presently the label formats are limited to
c              elements with one to nine nodes
c
      dimension mat(*),ien(nen,*)
      common /iounit/ iin,iout,iecho,ioupp,itest1,itest2 
c
      do 100 n=1,numel
      if (mod(n,50).eq.1) then
	  if(nen.le.16) write(iecho,1000) (i,i=1,nen)
	  if(nen.eq.20) write(iecho,1020) (i,i=1,nen)
	  if(nen.eq.24) write(iecho,1024) (i,i=1,nen)
	  if(nen.eq.28) write(iecho,1028) (i,i=1,nen)
	  if(nen.eq.32) write(iecho,1032) (i,i=1,nen)
	end if
      write(iecho,2000) n,mat(n),(ien(i,n),i=1,nen)
  100 continue
c
      return
c
 1000 format(///,
     &' d a t a   f o r   e l e m e n t   p a r a m e t e r s ',//1x,
     &' element      material',16('      node',i2))
 1020 format(///,
     &' d a t a   f o r   e l e m e n t   p a r a m e t e r s ',//1x,
     &' element      material',20('      node',i2))
 1024 format(///,
     &' d a t a   f o r   e l e m e n t   p a r a m e t e r s ',//1x,
     &' element      material',24('      node',i2))
 1028 format(///,
     &' d a t a   f o r   e l e m e n t   p a r a m e t e r s ',//1x,
     &' element      material',28('      node',i2))
 1032 format(///,
     &' d a t a   f o r   e l e m e n t   p a r a m e t e r s ',//1x,
     &' element      material',32('      node',i2))
2000  format(1x,i10,64(2x,i10))
      end
c**** new ********************************************************************** 
      subroutine prntel(mat,ien,nen,numel)
c
c.... program to print data for element with "nen" nodes
c
c        note: presently the label formats are limited to
c              elements with one to nine nodes
c
      dimension mat(*),ien(nen,*)
      common /iounit/ iin,iout,iecho,ioupp,itest1,itest2 
c
      do 100 n=1,numel
      if (mod(n,50).eq.1) write(iecho,1000) (i,i=1,nen)
      write(iecho,2000) n,mat(n),(ien(i,n),i=1,nen)
  100 continue
c
      return
c
 1000 format(///,
     &' d a t a   f o r   e l e m e n t   n o d e s ',//1x,
     &' element      material',16('      node',i2))
 2000 format(1x,i10,20(2x,i10))
      end
c**** new ********************************************************************** 
      subroutine prtdc 
c 
c.... program to print memory-pointer dictionary 
c 
      common /bpoint/ mfirst,mlast,ilast,mtot,iprec 
      common /iounit/ iin,iout,iecho,ioupp,itest1,itest2 
      character *4 ia 
      common na(1) 
      common /dictn/ ia(1) 
c 
      n = (mtot-mlast)/5 
      j = mtot + 1 
c 
      k = 1 
      do 100 i=1,n 
      if (mod(i,50).eq.1) write(iout,1000) 
      j = j - 5 
      call prtdc1(i,ia(k),na(j),na(j+1),na(j+2),na(j+3),na(j+4)) 
      k = k + 2 
 100  continue 
c 
      return 
c 
 1000 format(///, 
     &' d y n a m i c   s t o r a g e    a l l o c a t i o n', 
     &'   i n f o r m a t i o n '// 
     &  12x,'array no.',5x,'array',8x,'address',6x,'dim1',6x,'dim2', 
     &  6x, 'dim3',6x,'prec.'/) 
c 
      end 
c**** new ********************************************************************** 
      subroutine prtdc1(i,iname,iadd,ndim1,ndim2,ndim3,ipr) 
c 
c.... program to print memory-pointer information for an array 
c 
      character *4 iname(2) 
      common /iounit/ iin,iout,iecho,ioupp,itest1,itest2 
c 
      if (i.eq.1) neg = 1 
      if (iname(1).eq.'npar') then 
        write (iout,1000) neg 
        neg = neg + 1 
      endif 
      write(iout,2000) i,iname,iadd,ndim1,ndim2,ndim3,ipr 
c 
      return 
c 
 1000 format(/14x,'*****',7x,'begin element group number',i10/' ') 
 2000 format(14x,i10,7x,2a4,1x,6i10) 
      end 
c**** new **********************************************************************
      function rcdot(a,b,ma,n)
c
c.... program to compute the dot product of a vector stored row-wise
c        with a vector stored column-wise
c
      implicit real*8 (a-h,o-z)
c                                                                       
c.... remove above card for single-precision operation               
c                                                                       
      dimension a(ma,*),b(*)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
c
      rcdot = zero
c
      do 100 i=1,n
      rcdot = rcdot + a(1,i)*b(i)
  100 continue
c
      return
      end
c**** new **********************************************************************
      function rowdot(a,b,ma,mb,n)
c
c.... program to compute the dot product of vectors stored row-wise
c
      implicit real*8 (a-h,o-z)
c
c.... remove above card for single precision operation
c
      dimension a(ma,*),b(mb,*)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
c
      rowdot = zero
c
      do 100 i=1,n
      rowdot = rowdot + a(1,i)*b(1,i)
  100 continue
c
      return
      end    
      subroutine serror(name,i) 
c 
c.... program to print error message if available storage is exceeded 
c 
      character*4 name(2) 
      common /iounit/ iin,iout,iecho,ioupp,itest1,itest2 
c 
      call prtdc 
      write(*,1000) i,name
      pause 
      stop 
c 
 1000 format(1x,5('*'),'storage exceeded by ',i10, 
     &' words in attempting to store array ',2a4) 
      end 

c***** end *********************************************************************
      subroutine setupd(c,dmat,const,nstr,nrowb)
c
c.... program to calculate the d matrix
c
      implicit real*8 (a-h,o-z)
c
c.... remove above card for single precision operation
c
      dimension c(nrowb,*),dmat(nrowb,*)
c
      do 200 j=1,nstr
c
      do 100 i=1,j
      dmat(i,j) = const*c(i,j)
      dmat(j,i) = dmat(i,j)
  100 continue
c
  200 continue
c
      return
      end
c**** new **********************************************************************
      subroutine shgq4(xl,det,shl,shg,nint,nel,neg,quad)
c
c.... program to calculate global derivatives of shape functions and
c        jacobian determinants for a 4-node quadrilateral element
c
c        xl(j,i)    = global coordinates
c        det(l)     = jacobian determinant
c        shl(1,i,l) = local ("xi") derivative of shape function
c        shl(2,i,l) = local ("eta") derivative of shape function
c        shl(3,i,l) = local  shape function
c        shg(1,i,l) = x-derivative of shape function
c        shg(2,i,l) = y-derivative of shape function
c        shg(3,i,l) = shl(3,i,l)
c        xs(i,j)    = jacobian matrix
c                 i = local node number or global coordinate number
c                 j = global coordinate number
c                 l = integration-point number
c              nint = number of integration points, eq. 1 or 4
c
      implicit real*8 (a-h,o-z)
c
c.... remove above card for single precision operation
c
      logical quad
      dimension xl(2,*),det(*),shl(3,4,*),shg(3,4,*),xs(2,2)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
      common /iounit/ iin,iout,iecho,ioupp,itest1,itest2 
c
      call move(shg,shl,12*nint)
c
      do 700 l=1,nint
c
      if (.not.quad) then
         do 100 i=1,3
         shg(i,3,l) = shl(i,3,l) + shl(i,4,l)
         shg(i,4,l) = zero
  100    continue
      endif
c
      do 300 j=1,2
      do 200 i=1,2
      xs(i,j) = rowdot(shg(i,1,l),xl(j,1),3,2,4)
  200 continue
  300 continue
c
      det(l) = xs(1,1)*xs(2,2)-xs(1,2)*xs(2,1)
      if (det(l).le.zero) then
         write(iecho,1000) nel,neg
         stop
      endif
c
      do 500 j=1,2
      do 400 i=1,2
      xs(i,j) = xs(i,j)/det(l) 
  400 continue
  500 continue
c
      do 600 i=1,4
        temp = xs(2,2)*shg(1,i,l) - xs(1,2)*shg(2,i,l)
        shg(2,i,l) = - xs(2,1)*shg(1,i,l) + xs(1,1)*shg(2,i,l)
        shg(1,i,l) = temp
  600 continue
c
  700 continue
c
      return
c
 1000 format(///,'non-positive determinant in element number  ',i10,
     &          ' in element group  ',i10)
      end
c**** new **********************************************************************
      subroutine shlq4(shl,w,nint)
c
c.... program to calculate integration-rule weights, shape functions
c        and local derivatives for a four-node quadrilateral element
c
c               s,t = local element coordinates ("xi", "eta", resp.)       
c        shl(1,i,l) = local ("xi") derivative of shape function
c        shl(2,i,l) = local ("eta") derivative of shape function
c        shl(3,i,l) = local  shape function
c              w(l) = integration-rule weight
c                 i = local node number
c                 l = integration point number
c              nint = number of integration points, eq. 1 or 4
c
      implicit real*8 (a-h,o-z)
c
c.... remove above card for single precision operation
c
      dimension shl(3,4,*),w(*),ra(4),sa(4)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
      data ra/-0.5d0,0.5d0,0.5d0,-0.5d0/,sa/-0.5d0,-0.5d0,0.5d0,0.5d0/
c
      g = zero
      w(1) = four
      if (nint.eq.4) then
         g = two/dsqrt(three)
         w(1) = one
         w(2) = one
         w(3) = one
         w(4) = one
      endif
c
      do 200 l=1,nint
      r = g*ra(l)
      s = g*sa(l)
c
      do 100 i=1,4
      tempr = pt5 + ra(i)*r
      temps = pt5 + sa(i)*s
      shl(1,i,l) = ra(i)*temps
      shl(2,i,l) = tempr*sa(i)
      shl(3,i,l) = tempr*temps
  100 continue
c
  200 continue
c
      return
      end
c**** new **********************************************************************
      subroutine smult(a,b,c,mb,mc,m,n,iopt)
c
c.... program to perform scalar multiplication of a matrix
c
c        c(i,j) = a*b(i,j)
c
      implicit real*8 (a-h,o-z)
c
c.... remove above card for single precision operation
c
      dimension b(mb,*),c(mc,*)
c
      go to (1000,2000,3000),iopt
c
c.... iopt = 1, multiply entire matrix
c
 1000 do 1200 j=1,n
c
      do 1100 i=1,m
      c(i,j) = a*b(i,j)
 1100 continue
c
 1200 continue
      return
c
c.... iopt = 2, multiply lower triangular and diagonal elements
c
 2000 do 2200 j=1,n
c
      do 2100 i=j,m
      c(i,j) = a*b(i,j)
 2100 continue
c
 2200 continue
      return
c
c.... iopt = 3, multiply upper triangular and diagonal elements
c
 3000 do 3200 j=1,n
c
      do 3100 i=1,j
      c(i,j) = a*b(i,j)
 3100 continue
c
 3200 continue
      return
c
      end
c**** new **********************************************************************
      subroutine ztest(a,n,lzero)
c
c.... program to determine if an array contains only zero entries
c
      implicit real*8 (a-h,o-z)
c
c.... remove above card for single precision operation
c
      dimension a(*)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
      logical lzero
c
      lzero = .true.
c
      do 100 i=1,n
      if (a(i).ne.zero) then
         lzero = .false.
         return
      endif
  100 continue
c
      end     
c
c********************************************************************
      subroutine shap2m(s,xl,det,sh,nen,ien,nesd)
c********************************************************************
      implicit real*8 (a-h,o-z)
      dimension xl(nesd,*),sh(2,*),ien(*)
c
c     shape function
c
      sh(2,1)=(1.d00-s)/2.d00
      sh(2,2)=(1.d00+s)/2.d00
c
c
      sh(1,1)=-.5d00
      sh(1,2)=.5d00
c
c     3 node correction
c
      if(nen.eq.3) then
      corr=1.d00-s*s
      corrh=corr/2
      sh(2,1)=sh(2,1)-corrh
      sh(2,2)=sh(2,2)-corrh
      sh(2,3)=corr
c
      corr=-2.d00*s
      corrh=-s
      sh(1,1)=sh(1,1)-corrh
      sh(1,2)=sh(1,2)-corrh
      sh(1,3)=corr
c
      end if
c
      det=0.d00
      do 100 l=1,nen
      det=det+xl(1,l)*sh(1,l)
100   continue
c
c     global derivatives
c
      do 200 l=1,nen
      sh(1,l)=sh(1,l)/det
200   continue
      return
      end
c************************************************************************
       subroutine shap20(s,t,x,det,sh,sh2,nen,inc,ien,quad)
c************************************************************************
c      program to compute shape functions for quadrilateral		*
c									*
c        s,t        = natural coordinates				*
c        sh(nsd,i)  = first derivatives of shape functions		*
c        sh(3,i)    = shape functions					*
c        sh2(3,i)   = second derivatives of shape functions		*
c        xs(nsd,nsd)= jacobian matrix					*
c        det        = jacobian determinant				*
c        x(nsd,nen) = global coordinates				*
c************************************************************************
      implicit real*8(a-h,o-z)
      logical quad
      dimension sa(4),ta(4),sh(3,*),sh2(3,*),x(2,*),xs(2,2)
      dimension ien(1)
      data sa/-0.5d0,0.5d0,0.5d0,-0.5d0/,ta/-0.5d0,-0.5d0,0.5d0,0.5d0/
c
      do 10 i=1,4
      sh(3,i)=(0.5d0+sa(i)*s)*(0.5d0+ta(i)*t)
      sh(1,i)=sa(i)*(0.5d0+ta(i)*t)
      sh(2,i)=ta(i)*(0.5d0+sa(i)*s)
10    continue
c
      if(quad) goto 30
      do 20 i=1,3
      sh(i,3)=sh(i,3)+sh(i,4)
      sh(i,4)=0.d0
20    continue
c
30    if(nen.eq.8.or.nen.eq.9) call shap21(s,t,sh,nen,ien)
c
      do 40 i=1,2
      do 40 j=1,2
      xs(i,j)=0.d0
      do 40 k=1,nen
      xs(i,j)=xs(i,j)+x(i,k)*sh(j,k)
40    continue
      det=xs(1,1)*xs(2,2)-xs(1,2)*xs(2,1)
c
      do 50 i=1,2
      do 50 j=1,2
      xs(i,j)=xs(i,j)/det
50    continue

c      call shap22(s,t,xs,x,sh,sh2,ien,nen)

      do 60 i=1,nen
      temp=xs(2,2)*sh(1,i)-xs(2,1)*sh(2,i)
      sh(2,i)=-xs(1,2)*sh(1,i)+xs(1,1)*sh(2,i)
      sh(1,i)=temp
60    continue

      if(inc.eq.0) return
c
c.... incompatible modes
c
      if(quad) goto 80
      do 70 i=1,3
      do 70 j=5,6
   70 sh(i,j)=0.d0
      return
   80 sh(1,5)=-s-s
      sh(2,5)=0.d0
      sh(3,5)=1.d0-s*s
      sh(1,6)=0.d0
      sh(2,6)=-t-t
      sh(3,6)=1.d0-t*t
      xs(1,1)=0.25d0*(-x(1,1)+x(1,2)+x(1,3)-x(1,4))
      xs(1,2)=0.25d0*(-x(1,1)-x(1,2)+x(1,3)+x(1,4))
      xs(2,1)=0.25d0*(-x(2,1)+x(2,2)+x(2,3)-x(2,4))
      xs(2,2)=0.25d0*(-x(2,1)-x(2,2)+x(2,3)+x(2,4))
      do 90 i=5,6
      temp=(xs(2,2)*sh(1,i)-xs(2,1)*sh(2,i))/det
      sh(2,i)=(-xs(1,2)*sh(1,i)+xs(1,1)*sh(2,i))/det
   90 sh(1,i)=temp
      return
      end



  
      subroutine shap21 (s,t,sh,nen,ien)
c************************************************************************
c   program to compute shape functions and local derivatives 		*
c   for quadrilateral -  shape functions  5  to  9			*
c************************************************************************
      implicit real*8(a-h,o-z)
      dimension sh(3,*),ien(*)
      data zero/0.d0/
c
      ss=(1.d0-s*s)/2.d0
      tt=(1.d0-t*t)/2.d0
      do 10 i=5,8
      do 10 j=1,3
   10 sh(j,i)=zero
      s19=zero
      s29=zero
      s39=zero
c
      if(nen.ne.9) goto 15
      if(ien(9).eq.0) goto 15
      s19=-2.d0*tt*s
      s29=-2.d0*ss*t
      s39=2.d0*tt*ss
      sh(1,9)=2.d0*s19
      sh(2,9)=2.d0*s29
      sh(3,9)=2.d0*s39
      do 20 i=1,4
      do 20 j=1,3
   20 sh(j,i)=sh(j,i)-sh(j,9)/4.d0
c
   15 if(ien(5).eq.0) goto 30
      sh(1,5)=-s*(1.d0-t)-s19
      sh(2,5)=-ss-s29
      sh(3,5)=ss*(1.d0-t)-s39
c
   30 if(nen.lt.6) goto 60
      if(ien(6).eq.0) goto 40
      sh(1,6)=tt-s19
      sh(2,6)=-t*(1.d0+s)-s29
      sh(3,6)=tt*(1.d0+s)-s39
c
   40 if(nen.lt.7) goto 60
      if(ien(7).eq.0) goto 50
      sh(1,7)=-s*(1.d0+t)-s19
      sh(2,7)=ss-s29
      sh(3,7)=ss*(1.d0+t)-s39
c
   50 if(nen.lt.8) goto 60
      if(ien(8).eq.0) goto 60
      sh(1,8)=-tt-s19
      sh(2,8)=-t*(1.d0-s)-s29
      sh(3,8)=tt*(1.d0-s)-s39
c
   60 k=8
      do 70 i=1,4
      l=i+4
      do 80 j=1,3
   80 sh(j,i)=sh(j,i)-0.5d0*(sh(j,k)+sh(j,l))
   70 k=l
c
      return
      end

c************************************************************************
      subroutine shap22(s,t,xs,x,sh,sh2,ien,nen)
c************************************************************************
c   compute second derivatives of shape functions 1 to 9 for		*
c   quadrilateral 							*
c          eh(1,i)  =  d2(Ni)/de2					*
c          eh(2,i)  =  d2(Ni)/dn2					*
c          eh(3,i)  =  d2(Ni)/dedn					*
c************************************************************************
      implicit real*8(a-h,o-z)
      dimension sh(3,*),sh2(3,*),ien(*),xs(2,2),x(2,*),eh(3,9)
      call clear(eh,27)
c
c.....local second derivatives of shape functions
c
      eh(3,1)= 0.25d0
      eh(3,2)=-0.25d0
      eh(3,3)= 0.25d0
      eh(3,4)=-0.25d0
c
      if(nen.le.4) goto 110
c
      if(nen.ne.9) goto 30
      if(ien(9).eq.0) goto 30
      eh(1,9)=-2.d0*(1.d0-t*t)
      eh(2,9)=-2.d0*(1.d0-s*s)
      eh(3,9)=4.d0*s*t
      do 40 i=1,4
      do 40 j=1,3
   40 eh(j,i)=eh(j,i)-eh(j,9)/4.d0
   30 e19=0.5d0*eh(1,9)
      e29=0.5d0*eh(2,9)
      e39=0.5d0*eh(3,9)
c
      if(ien(5).eq.0) goto 50
      eh(1,5)=-1.d0+t-e19
      eh(2,5)=-e29
      eh(3,5)=s-e39
c
   50 if(nen.lt.6) goto 80
      if(ien(6).eq.0) goto 60
      eh(1,6)=-e19
      eh(2,6)=-1.d0-s-e29
      eh(3,6)=-t-e39
c
   60 if(nen.lt.7) goto 80
      if (ien(7).eq.0) goto 70
      eh(1,7)=-1.d0-t-e19
      eh(2,7)=-e29
      eh(3,7)=-s-e39
c
   70 if(nen.lt.8) goto 80
      if(ien(8).eq.0) goto 80
      eh(1,8)=-e19
      eh(2,8)=-1.d0+s-e29
      eh(3,8)=t-e39
c
   80 k=8
      do 90 i=1,4
      l=i+4
      do 100 j=1,3
  100 eh(j,i)=eh(j,i)-0.5d0*(eh(j,k)+eh(j,l))
   90 k=l
c
c.....global second derivatives
c
110   call shap23(xs,x,eh,sh,sh2,nen)
      return
      end



      subroutine shap23(xs,x,eh,sh,sh2,nen)
c************************************************************************
c   transform second derivatives from natural coordinates to		*
c   global coordinates 							*
c         sh2(1,i)  =  d2(Ni)/dx2					*
c         sh2(2,i)  =  d2(Ni)/dy2					*
c         sh2(3,i)  =  d2(Ni)/dxdy					*
c************************************************************************
      implicit real*8(a-h,o-z)
      dimension xs(2,2),x(2,*),eh(3,*),sh(3,*),sh2(3,*)
      dimension t2(3,3),c1(3,2),t1(3,2),xj(2,2)
c
      call clear(c1,6)
      call clear(t1,6)
      call clear(sh2,3*nen)
c
c.... form j inverse of jacobian matrix
c
      xj(1,1)=xs(2,2)
      xj(2,2)=xs(1,1)
      xj(1,2)=-xs(2,1)
      xj(2,1)=-xs(1,2)
c
c.... form t2
c
      do 10 i=1,2
      t2(i,3)=2.d0*xj(i,1)*xj(i,2)
      t2(3,i)=xj(1,i)*xj(2,i)
      do 10 j=1,2
      t2(i,j)=xj(i,j)**2
10    continue
      t2(3,3)=xj(1,1)*xj(2,2)+xj(1,2)*xj(2,1)
c
c.... form c1
c
      do 20 n=1,3
      do 20 i=1,2
      do 20 j=1,nen
      c1(n,i)=c1(n,i)+eh(n,j)*x(i,j)
20    continue
c
c.... form t1
c
      do 30 i=1,3
      do 30 j=1,2
      do 30 k=1,3
      t1(i,j)=t1(i,j)-t2(i,k)*(c1(k,1)*xj(1,j)+c1(k,2)*xj(2,j))
30    continue
c
c.... transformation from natural coor. to global coor.
c
      do 50 n=1,nen
      do 50 l=1,3
         do 60 i=1,2
         sh2(l,n)=sh2(l,n)+t1(l,i)*sh(i,n)
60       continue
         do 70 i=1,3
         sh2(l,n)=sh2(l,n)+t2(l,i)*eh(i,n)
70       continue
50    continue
c
      return
      end

c**** new *************************************************************
      subroutine oneshl(shl,w,nint,nen)
c
c.... program to calculate integration-rule weights, shape functions 
c        and local derivatives for a two, three or four  node, 
c        one-dimensional element
c
c                 r = local element coordinate ("xi")
c        shl(1,i,l) = local ("xi") derivative of shape function
c        shl(2,i,l) = shape function
c              w(l) = integration-rule weight
c                 i = local node number 
c                 l = integration-point number
c              nint = number of integration points, eq. 1, 2, 3, 4, 6, 7 ,8
c
      implicit real*8 (a-h,o-z)
c
c.... remove above card for single precision operation
c
      dimension shl(3,nen,*),w(*),ra(10),xa(10)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
      data  five9/0.5555555555555555d0/,eight9/0.8888888888888888d0/
c
c     
c
      if (nint.eq.1) then
         w(1)  = two
         ra(1) = zero
      endif
c
      if (nen.eq.1) xa(1) = zero
c
c
c
      if (nint.eq.2) then
         w(1) = one
         w(2) = one
         ra(1)=-.577350269189626
         ra(2)= .577350269189625
      endif
c
      if (nen.eq.2) then
         xa(1) = -one
         xa(2) =  one
      endif
c
c
      if (nint.eq.3) then
         w(1) = five9
         w(2) = eight9
         w(3) = five9
         ra(1)=-.774596669241483
         ra(2)= zero
         ra(3)= .774596669241483
      endif
c
c 
      if(nen.eq.3) then
         xa(1)= -one
         xa(2)= one
         xa(3)= zero
      endif
c
c
      if (nint.eq.4) then
         w(1) = .347854845137454
         w(2) = .652145154862546
         w(3) = .652145154862546
         w(4) = .347854845137454
         ra(1)=-.861136311594053
         ra(2)=-.339981043584856
         ra(3)= .339981043584856
         ra(4)= .861136311594053
      endif
c
      if (nen.eq.4) then
         xa(1) = -one
         xa(2) = one
         xa(3) = -.333333333333333
         xa(4) =  .333333333333333
         endif
c
c
c
       if(nint.eq.5) then
        w(1) = .236926885056189
        w(2) = .478628670499366
        w(3) = .568888888888888
        w(4) = .478628670499366
        w(5) = .236926885056189
        ra(1)=-.906179845938664
        ra(2)=-.538469310105683
        ra(3)= zero
        ra(4)= .538469310105683
        ra(5)= .906179845938664
       endif
c
       if(nen.eq.5) then
         xa(1)= -one 
         xa(2)=  one 
         xa(3)= -pt5
         xa(4)= zero
         xa(5)= pt5
       endif
c
       if(nint.eq.6) then
         w(1) = .171324492397170
         w(2) = .360761573048139
         w(3) = .467913934572691
         w(4) = .467913934572691
         w(5) = .360761573048139
         w(6) = .171324492397170
c         
         ra(1)=-.932469514203152
         ra(2)=-.661209386466265
         ra(3)=-.238619186083197
         ra(4)= .238619186083197
         ra(5)= .661209386466365
         ra(6)= .932469514203152
        endif
c
        if(nen.eq.6) then
         xa(1) = -one
         xa(2) =  one
         xa(3) = -.600000000000000
         xa(4) = -.200000000000000
         xa(5) =  .200000000000000
         xa(6) =  .600000000000000
        endif
c
       if(nint.eq.7) then
         w(1) = .129484966168870
         w(2) = .279705391489277 
         w(3) = .381830050505119
         w(4) = .417959183673469
         w(5) = .381830050505119
         w(6) = .279705391489277
         w(7) = .129484966168870
c         
         ra(1)=-.949107912342759
         ra(2)=-.741531185599394
         ra(3)=-.405845151377397
         ra(4)= zero
         ra(5)= .405845151377397
         ra(6)= .741531185599394
         ra(7)= .949107912342759
        endif
c
        if(nen.eq.7) then
         xa(1) = -one
         xa(2) =  one
         xa(3) = -.666666666666666
         xa(4) = -.333333333333333
         xa(5) = zero
         xa(6) =  .333333333333333
         xa(7) =  .666666666666666
        endif
c
       if(nint.eq.8) then
         w(1) = .101228536290376
         w(2) = .222381034453374
         w(3) = .313706645877887
         w(4) = .362683783378362
         w(5) = .362683783378362
         w(6) = .313706645877887
         w(7) = .222381034453374
         w(8) = .101228536290376
c
         ra(1)=-.960289856497536
         ra(2)=-.796666477413627
         ra(3)=-.525532409916329
         ra(4)=-.183434642495650
         ra(5)= .183434642495650
         ra(6)= .525532409916329
         ra(7)= .796666477413627
         ra(8)= .960289856497536
        endif
        if(nen.eq.8) then
         xa(1) = -one
         xa(2) =  one
         xa(3) = -0.71428571428571
         xa(4) = -0.42857142857143
         xa(5) = -0.14285714285714
         xa(6) = 0.14285714285714
         xa(7) = 0.42857142857143
         xa(8) = 0.71428571428571
        endif
c
      do 100 l = 1, nint
         r = ra(l)
c
        if(nen.eq.1) then
        shl(1,1,l) = zero
        shl(2,1,l) = one
        go to 100
        endif
c
        do 50 i = 1, nen
         aa = one
         bb = one
         aax = zero
         do 40 j =1, nen
          daj = one
          if (i .ne. j)then
          aa = aa * ( r - xa(j))
          bb = bb * ( xa(i) - xa(j))
          do 30 k = 1, nen
           if(k .ne. i .and. k .ne. j) daj = daj * ( r - xa(k))
   30     continue
          aax =aax + daj
          endif
   40    continue
        shl(2,i,l) = aa/bb
        shl(1,i,l) = aax/bb
   50  continue
c
  100  continue      
       return
       end
c***********************************************************************
      subroutine oneshg(xl,det,shl,shg,nen,nint,nesd,ns,nel,neg)
c
c.... program to calculate global derivatives of shape functions 
c        and jacobian determinants for the bi-dimensional,
c        elastic beam element
c
c           xl(j,l) = global coordinates of integration points
c        shl(1,i,l) = local ("xi") derivative of shape function
c        shl(2,i,l) = shape function
c        shg(1,i,l) = global ("arc-length") derivative of shape ftn
c        shg(2,i,l) = shl(2,i,l)
c            det(l) = euclidean length 
c                 i = local node number 
c                 j = global coordinate number
c                 l = integration-point number
c              nint = number of integration points
c
      implicit real*8 (a-h,o-z)
c
c.... remove above card for single precision operation
c
      dimension xl(nesd,*),det(*),shl(3,nen,*),shg(3,nen,*)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
      common /iounit/ iin,iout,iecho,ioupp,itest1,itest2 
c
      do 400 l=1,nint
c
      det(l)=zero
      x1=0.d0
      x2=0.d0
      do 100 j=1,nen
      x1=x1+shl(1,j,l)*xl(1,j)
      x2=x2+shl(1,j,l)*xl(2,j)
100   continue
      det(l)=dsqrt(x1*x1+x2*x2)
c
      if (det(l).le.zero) then
         write(iecho,1000) ns,nel,neg
         stop
      endif
c
      do 300 i=1,nen
      shg(1,i,l)=shl(1,i,l)/det(l)
      shg(2,i,l)=shl(2,i,l)
300   continue
c
  400 continue
 1000 format(///,' oneshg - non-positive determinant in side ',i10,/,
     &          ' in element ',i10,5x,' in element group  ',i10)
c
      return
      end
c***********************************************************************
      subroutine oneshgp(xl,det,shl,shlp,shgp,
     &                   nen,npars,nint,nesd,ns,nel,neg)
c
c.... program to calculate global derivatives of shape functions 
c        and jacobian determinants for the bi-dimensional,
c        elastic beam element
c
c           xl(j,l) = global coordinates of integration points
c        shl(1,i,l) = local ("xi") derivative of shape function
c        shl(2,i,l) = shape function
c        shg(1,i,l) = global ("arc-length") derivative of shape ftn
c        shg(2,i,l) = shl(2,i,l)
c            det(l) = euclidean length 
c                 i = local node number 
c                 j = global coordinate number
c                 l = integration-point number
c              nint = number of integration points
c
      implicit real*8 (a-h,o-z)
c
c.... remove above card for single precision operation
c
      dimension xl(nesd,*),det(*),shl(3,nen,*)
	dimension shlp(3,npars,*),shgp(3,npars,*)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
      common /iounit/ iin,iout,iecho,ioupp,itest1,itest2 
c
c
      do 400 l=1,nint
c
      det(l)=zero
      x1=0.d0
      x2=0.d0
      do 100 j=1,nen
      x1=x1+shl(1,j,l)*xl(1,j)
      x2=x2+shl(1,j,l)*xl(2,j)
100   continue
      det(l)=dsqrt(x1*x1+x2*x2)
c
      if (det(l).le.zero) then
         write(iecho,1000) ns,nel,neg
         stop
      endif
c
      do 300 i=1,npars
      shgp(1,i,l)=shlp(1,i,l)/det(l)
      shgp(2,i,l)=shlp(2,i,l)
300   continue
c
  400 continue
 1000 format(///,' oneshgp - non-positive determinant in side ',i10,/,
     &          ' in element ',i10,5x,'in oneshgp group  ',i10)
c
      return
      end
c**** new **********************************************************************
      subroutine shgq(xl,det,shl,shg,nint,nel,neg,quad,nen)
c
c.... program to calculate global derivatives of shape functions and
c        jacobian determinants for a  quadrilateral element
c
c        xl(j,i)    = global coordinates
c        det(l)     = jacobian determinant
c        shl(1,i,l) = local ("xi") derivative of shape function
c        shl(2,i,l) = local ("eta") derivative of shape function
c        shl(3,i,l) = local  shape function
c        shg(1,i,l) = x-derivative of shape function
c        shg(2,i,l) = y-derivative of shape function
c        shg(3,i,l) = shl(3,i,l)
c        xs(i,j)    = jacobian matrix
c                 i = local node number or global coordinate number
c                 j = global coordinate number
c                 l = integration-point number
c              nint = number of integration points, eq. 1 or 4
c
      implicit real*8 (a-h,o-z)
c
c.... remove above card for single precision operation
c
      logical quad
      dimension xl(2,*),det(*),shl(3,nen,*),shg(3,nen,*),xs(2,2)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
      common /iounit/ iin,iout,iecho,ioupp,itest1,itest2 
c
      call move(shg,shl,3*nen*nint)
c
      do 700 l=1,nint
c
      if (.not.quad) then
         do 100 i=1,3
         shg(i,3,l) = shl(i,3,l) + shl(i,4,l)
         shg(i,4,l) = zero
  100    continue
      endif
c
      do 300 j=1,2
      do 200 i=1,2
      xs(i,j) = rowdot(shg(i,1,l),xl(j,1),3,2,nen)
  200 continue
  300 continue
c
      det(l) = xs(1,1)*xs(2,2)-xs(1,2)*xs(2,1)
      if (det(l).le.zero) then
         write(iecho,1000) nel,neg
         stop
      endif
c
      do 500 j=1,2
      do 400 i=1,2
      xs(i,j) = xs(i,j)/det(l) 
  400 continue
  500 continue
c
      do 600 i=1,nen
        temp = xs(2,2)*shg(1,i,l) - xs(1,2)*shg(2,i,l)
        shg(2,i,l) = - xs(2,1)*shg(1,i,l) + xs(1,1)*shg(2,i,l)
        shg(1,i,l) = temp
  600 continue
c
  700 continue
c
      return
c
 1000 format(////,' shgq - non-positive determinant - element ',i10,
     &          ' in element group  ',i10)
      end
c******* nova subrotina ***************************************************************
c                                              
      subroutine shlt(shl,w,nint,nen) 
c
c     Objetivo: calcular pesos, funcoes de interpolacao e derivadas locais
c               para elementos triangulares
c
c----------------------------------------------------------------------
c
c 
c.... program to calculate integration-rule weights, shape functions 
c        and local derivatives for a triangular element 
c 
c        c1, c2, c3 = local element coordinates ("l1", "l2", "l3".) 
c        shl(j,i,l) = local ("j") derivative of shape function 
c        shl(3,i,l) = local  shape function 
c              w(l) = integration-rule weight 
c                 i = local node number 
c                 l = integration point number 
c              nint = number of integration points, eq. 1 or 4 
c 
      implicit real*8 (a-h,o-z) 
c 
c.... remove above card for single precision operation 
c 
      dimension shl(3,nen,*),w(*),cl1(16),cl2(16),cl3(16)
      data   zero,pt1667,pt25,pt5 
     &      /0.0d0,0.1666666666666667d0,0.25d0,0.5d0/, 
     &       one,two,three,four,five,six 
     &      /1.0d0,2.0d0,3.0d0,4.0d0,5.0d0,6.0d0/ 
      data r1/0.33333333333333333333d00/,w1/1.d00/, 
     &     r2/0.5d00                   /,w2/0.3333333333333333333d00/, 
     &     r3a/0.3333333333333333333d00/,w3a/-0.5625d00/, 
     &     r3b1/0.6d00                 /,w3b/0.520833333333333333d00/, 
     &     r3b2/0.2d00                 /, 
     &     r7a/0.3333333333333333333d00/,w7a/0.225d00/, 
     &     r7b/0.0597158717d00         /,w7b/0.1323941527d00/, 
     &     r7c/0.4701420641d00         /, 
     &     r7d/0.7974269853d00         /,w7d/0.1259391805d00/, 
     &     r7e/0.1012865073d00         /, 
     &     r13a/0.3333333333333333333d00/,w13a/-0.1495700444676d0/, 
     &     r13b/0.479308067841923d0         /,w13b/0.175615257433204d0/, 
     &     r13c/0.260345966079038d0         /,w13c/0.053347235608839d0/ 
     &     r13d/0.869739794195568d0         /,w13d/0.077113760890257d0/, 
     &     r13e/0.065130102902216d0        /, 
     &     r13f/0.638444188568809d0         /, 
     &     r13g/0.312865496004875d0         / 
     &     r13h/0.048690315425316d0         / 
c 
      if (nint.eq.1) then 
            w(1)=w1/two 
            cl1(1)=r1 
            cl2(1)=r1 
            cl3(1)=one-r1-r1 
      end if 
c 
      if(nint.eq.3) then 
            do 10 i=1,3 
                  w(i)=w2/two 
   10       continue 
            cl1(1)=r2 
            cl2(1)=r2 
            cl3(1)=zero 
            cl1(2)=zero 
            cl2(2)=r2 
            cl3(2)=r2 
            cl1(3)=r2 
            cl2(3)=zero 
            cl3(3)=r2 
      end if 
c 
      if(nint.eq.4) then 
      w(1)= w3a/two 
            do 20 i=2,4 
                  w(i)=w3b/two 
   20       continue 
            cl1(1)=r3a 
            cl2(1)=r3a 
            cl3(1)=one - r3a - r3a 
            cl1(2)=r3b1 
            cl2(2)=r3b2 
            cl3(2)=r3b2 
            cl1(3)=r3b2 
            cl2(3)=r3b1 
            cl3(3)=r3b2 
            cl1(4)=r3b2 
            cl2(4)=r3b2 
            cl3(4)=r3b1 
      end if 
c 
      if(nint.eq.7) then 
      w(1)= w7a/two 
            do 30 i=2,4 
                  w(i)=w7b/two 
   30       continue 
            do 40 i=5,7 
                  w(i)=w7d/two 
   40       continue 
            cl1(1)=r7a 
            cl2(1)=r7a 
            cl3(1)=r7a 
      do 50 i=2,4 
                  cl1(i)=r7c 
                  cl2(i)=r7c 
                  cl3(i)=r7c 
   50 continue 
            cl1(2)=r7b 
            cl2(3)=r7b 
            cl3(4)=r7b 
      do 60 i=5,7 
                  cl1(i)=r7e 
                  cl2(i)=r7e 
                  cl3(i)=r7e 
   60 continue 
            cl1(5)=r7d 
            cl2(6)=r7d 
            cl3(7)=r7d 
      end if 
      if (nint.eq.13) then 
        
      w(1)= w13a/two 
            do 31 i=2,4 
                  w(i)=w13b/two 
   31       continue 
            do 41 i=5,7 
                  w(i)=w13c/two 
   41       continue 
            do 51 i=8,13 
                  w(i)=w13d/two 
   51       continue 
            cl1(1)=r13a 
            cl2(1)=r13a 
            cl3(1)=r13a 
      do 52 i=2,4 
                  cl1(i)=r13c 
                  cl2(i)=r13c 
                  cl3(i)=r13c 
   52 continue 
 
            cl1(2)=r13b 
            cl2(3)=r13b 
            cl3(4)=r13b 
      do 53 i=5,7 
                  cl1(i)=r13e 
                  cl2(i)=r13e 
                  cl3(i)=r13e 
   53 continue 
 
            cl1(5)=r13d 
            cl2(6)=r13d 
            cl3(7)=r13d 
 
            cl1(8)=r13f 
            cl2(8)=r13g 
            cl3(8)=r13h 
            cl1(9)=r13h 
            cl2(9)=r13f 
            cl3(9)=r13g 
            cl1(10)=r13g 
            cl2(10)=r13h 
            cl3(10)=r13f 
 
 
            cl1(11)=r13f 
            cl2(11)=r13h 
            cl3(11)=r13g 
            cl1(12)=r13g 
            cl2(12)=r13f 
            cl3(12)=r13h 
            cl1(13)=r13h 
            cl2(13)=r13g 
            cl3(13)=r13f 
       endif 
 
c 
      do 200 l=1,nint 
c 
            c1 = cl1(l) 
            c2 = cl2(l) 
            c3 = cl3(l) 
            shl(1,1,l)= one 
            shl(2,1,l)= zero 
            shl(3,1,l)= c1 
            shl(1,2,l)= zero 
            shl(2,2,l)= one 
            shl(3,2,l)= c2 
            shl(1,3,l)=-one 
            shl(2,3,l)=-one 
            shl(3,3,l)= c3 
            if(nen.eq.6) then 
                  shl(1,4,l)= four * c2 
                  shl(2,4,l)= four * c1 
                  shl(3,4,l)= four * c1 * c2 
                  shl(1,5,l)=-four * c2 
                  shl(2,5,l)= four * (c3 - c2) 
                  shl(3,5,l)= four * c2 * c3 
                  shl(1,6,l)= four * (c3 - c1) 
                  shl(2,6,l)=-four * c1 
                  shl(3,6,l)= four * c3 * c1 
c 
                  do 70 i=1,3 
                        shl(i,1,l)=shl(i,1,l) 
     &                  -pt5*(shl(i,4,l)+shl(i,6,l)) 
                        shl(i,2,l)=shl(i,2,l) 
     &                  -pt5*(shl(i,4,l)+shl(i,5,l)) 
                        shl(i,3,l)=shl(i,3,l) 
     &                  -pt5*(shl(i,5,l)+shl(i,6,l)) 
   70             continue 
            end if 
 
            if (nen.eq.10) then 
 
                  shl(1,1,l)= pt5*((three*c1-one)*(three*c1-two) 
     &            +c1*three*(three*c1-two)+c1*three*(three*c1-one)) 
                  shl(2,1,l)= zero 
                  shl(3,1,l)= pt5*c1*(three*c1-one)*(three*c1-two) 
 
                  shl(1,2,l)=zero 
                  shl(2,2,l)=  pt5*((three*c2-one)*(three*c2-two) 
     &            +c2*three*(three*c2-two)+c2*three*(three*c2-one)) 
                  shl(3,2,l)=  pt5*c2*(three*c2-one)*(three*c2-two) 
 
                  shl(1,3,l)= -pt5*((three*c3-one)*(three*c3-two) 
     &            +c3*three*(three*c3-two)+c3*three*(three*c3-one)) 
                  shl(2,3,l)=-pt5*((three*c3-one)*(three*c3-two) 
     &            +c3*three*(three*c3-two)+c3*three*(three*c3-one)) 
                  shl(3,3,l)=  pt5*c3*(three*c3-one)*(three*c3-two) 
 
 
 
 
                  shl(1,4,l) = (9.0d0/two)*(c2*(three*c1-one) 
     &            +c1*c2*three) 
                  shl(2,4,l) = (9.0d0/two)*c1*(three*c1-one) 
                  shl(3,4,l) = (9.0d0/two)*c1*c2*(three*c1-one) 
 
                  shl(1,5,l) = (9.0d0/two)*c2*(three*c2-one) 
                  shl(2,5,l) = (9.0d0/two)*(c1*(three*c2-one) 
     &            +c1*c2*three) 
                  shl(3,5,l) = (9.0d0/two)*c1*c2*(three*c2-one) 
 
 
 
                 shl(1,6,l) =- (9.0d0/two)*c2*(three*c2-one) 
                 shl(2,6,l) =  (9.0d0/two)*(c3*(three*c2-one) 
     &           +c2*c3*three-c2*(three*c2-one)) 
                 shl(3,6,l) = (9.0d0/two)*c3*c2*(three*c2-one) 
 
 
                  shl(1,7,l) = (9.0d0/two)*(-c2*(three*c3-one) 
     &            -c2*c3*three) 
                  shl(2,7,l) =  (9.0d0/two)*(c3*(three*c3-one) 
     &            -c2*c3*three -c2*(three*c3-one)) 
                  shl(3,7,l) = (9.0d0/two)*c2*c3*(three*c3-one) 
 
                  shl(1,9,l)= (9.0d0/two)*(-c1*(three*c1-one) 
     &            +c1*c3*three +c3*(three*c1-one)) 
                  shl(2,9,l) =  (9.0d0/two)*(-c1*(three*c1-one)) 
                  shl(3,9,l) = (9.0d0/two)*c3*c1*(three*c1-one) 
 
 
                  shl(1,8,l) = (9.0d0/two)*(-c1*(three*c3-one) 
     &            -c1*c3*three +c3*(three*c3-one)) 
                  shl(2,8,l) =  (9.0d0/two)*( 
     &            -c1*c3*three -c1*(three*c3-one)) 
                  shl(3,8,l) = (9.0d0/two)*c1*c3*(three*c3-one) 
 
                 shl(1,10,l) = 27.0d0*(c2*c3-c1*c2) 
                 shl(2,10,l) = 27.0d0*(c1*c3-c1*c2) 
                 shl(3,10,l) = 27.0d0*c1*c2*c3 
          endif 
  200      continue 
c 
           return 
           end
C 
c		**************************************************************************
c******* nova subrotina ***************************************************************
      subroutine shg2q(xl,shl,shl2,shg2,nint,nel,neg,quad,nen)
c
c.... program to calculate global second derivatives of shape functions
c        for a  quadrilateral element
c
c        xl(j,i)    = global coordinates
c        shl(1,i,l) = local ("xi") derivative of shape function
c        shl(2,i,l) = local ("eta") derivative of shape function
c        shl(3,i,l) = local  shape function
c        shl2(1,i,l) = local second ("xi") derivative of shape function
c        shl2(2,i,l) = local second ("eta") derivative of shape function
c        shl2(3,i,l) = local second ("xi*eta") derivative of shape function
c        shg2(1,i,l) = second x-derivative of shape function
c        shg2(2,i,l) = second y-derivative of shape function
c        shg2(3,i,l) = second xy-derivative of shape function
c        xs(i,j)    = jacobian matrix
c                 i = local node number or global coordinate number
c                 j = global coordinate number
c                 l = integration-point number
c              nint = number of integration points, eq. 1 or 4
c
      implicit real*8 (a-h,o-z)
c
c.... remove above card for single precision operation
c
      logical quad
      dimension xl(2,*),shl(3,nen,*),shl2(3,nen,*),shg2(3,nen,*),
     &          xs(2,2),t1(3,2),t2(3,3),c1(3,2)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
      common /iounit/ iin,iout,iecho,ioupp,itest1,itest2 
c
      do 960 l=1,nint
c
      if (.not.quad) then
         do 100 i=1,3
         shl(i,3,l) = shl(i,3,l) + shl(i,4,l)
         shl(i,4,l) = zero
         shl2(i,3,l) = shl2(i,3,l) + shl2(i,4,l)
         shl2(i,4,l) = zero
  100    continue
      endif
c
      do 300 j=1,2
      do 200 i=1,2
      xs(i,j) = rowdot(shl(i,1,l),xl(j,1),3,2,nen)
  200 continue
  300 continue
c
      det = xs(1,1)*xs(2,2)-xs(1,2)*xs(2,1)
      if (det.le.zero) then
         write(iecho,1000) nel,neg
         stop
      endif
c
      do 500 j=1,2
      do 400 i=1,2
      xs(i,j) = xs(i,j)/det
  400 continue
  500 continue
c
c      jacobian elements
c
      dxirx=xs(2,2)
      dxiry=xs(1,2)
      detarx=xs(2,1)
      detary=xs(1,1)
      xs(1,1)=dxirx
      xs(1,2)=detarx
      xs(2,1)=dxiry
      xs(2,2)=detary
c
c      calculation of global second order derivatives
c      
c      {d2global}=[t1]{dlocal}+[t2]{d2local}
c      [t1]=-[t2][c1][xs]
c
c.... form t2
c
      t2(1,1)=xs(1,1)**2
      t2(1,2)=xs(1,2)**2
      t2(1,3)=two*xs(1,1)*xs(1,2)
      t2(2,1)=xs(2,1)**2
      t2(2,2)=xs(2,2)**2
      t2(2,3)=two*xs(2,1)*xs(2,2)
      t2(3,1)=xs(1,1)*xs(2,1)
      t2(3,2)=xs(1,2)*xs(2,2)
      t2(3,3)=xs(1,1)*xs(2,2)+xs(1,2)*xs(2,1)
c
c.... form c1
c
      do 600 i=1,3
      do 600 j=1,2
      c1(i,j) = rowdot(shl2(i,1,l),xl(j,1),3,2,nen)
600   continue
c
c.... form t1
c
      do 700 i=1,3
      do 700 j=1,2
c--------------------      
      t1(i,j)=0.d00
c--------------------
      do 700 k=1,3
      t1(i,j)=t1(i,j)-t2(i,k)*(c1(k,1)*xs(1,j)+c1(k,2)*xs(2,j))
700    continue
c
c.... transformation from natural coor. to global coor.
c
      do 950 j=1,nen
      do 950 i=1,3
c------------------------
      shg2(i,j,l)=0.d00
c------------------------
         do 800 k=1,2
         shg2(i,j,l)=shg2(i,j,l)+t1(i,k)*shl(k,j,l)
800       continue
         do 900 k=1,3
         shg2(i,j,l)=shg2(i,j,l)+t2(i,k)*shl2(k,j,l)
900       continue
950    continue
c
  960 continue
c
      return
c
 1000 format(///,'non-positive determinant in element number  ',i10,
     &          ' in element group  ',i10)
      end
c**** new **********************************************************************
      subroutine shl2q(shl2,nint,nen)
c
c.... program to calculate local second  derivatives 
c     for a four-node quadrilateral element
c
c               s,t = local element coordinates ("xi", "eta", resp.)       
c        shl2(1,i,l) = local second ("xi") derivative of shape function
c        shl2(2,i,l) = local second ("eta") derivative of shape function
c        shl2(3,i,l) = local second ("xi*eta") derivative of shape function
c                 i = local node number
c                 l = integration point number
c              nint = number of integration points, eq. 1 or 4
c
      implicit real*8 (a-h,o-z)
c
c.... remove above card for single precision operation
c
      dimension shl2(3,nen,*),ra(16),sa(16)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
      data r1/0.d00/,
     &     r2/0.577350269189626d00/,
     &     r3a/0.774596669241483d00/,
     &     r3b/0.d00/,
     &     r4a/0.861136311594053d00/,
     &     r4b/0.339981043584856d00/
c
      if (nint.eq.1) then
            ra(1)=r1
            sa(1)=r1
      end if
c
      if(nint.eq.4) then
            ra(1)=r2
            sa(1)=r2
            ra(2)=-r2
            sa(2)=r2
            ra(3)=-r2
            sa(3)=-r2
            ra(4)=r2
            sa(4)=-r2
      end if
c
      if(nint.eq.9) then
            ra(1)=r3a
            sa(1)=r3a
            ra(2)=-r3a
            sa(2)=r3a
            ra(3)=-r3a
            sa(3)=-r3a
            ra(4)=r3a
            sa(4)=-r3a
c
            ra(5)=r3a
            sa(5)=r3b
            ra(6)=r3b
            sa(6)=r3a
            ra(7)=-r3a
            sa(7)=r3b
            ra(8)=r3b
            sa(8)=-r3a
c
            ra(9)=r3b
            sa(9)=r3b
      end if
c
      if(nint.eq.16) then
            ra(1)=r4a
            sa(1)=r4a
            ra(2)=-r4a
            sa(2)=r4a
            ra(3)=-r4a
            sa(3)=-r4a
            ra(4)=r4a
            sa(4)=-r4a
c
            ra(5)=r4b
            sa(5)=r4b
            ra(6)=-r4b
            sa(6)=r4b
            ra(7)=-r4b
            sa(7)=-r4b
            ra(8)=r4b
            sa(8)=-r4b
c
            ra(9)=r4b
            sa(9)=r4a
            ra(10)=-r4b
            sa(10)=r4a
            ra(11)=-r4a
            sa(11)=r4b
            ra(12)=-r4a
            sa(12)=-r4b
            ra(13)=-r4b
            sa(13)=-r4a
            ra(14)=r4b
            sa(14)=-r4a
            ra(15)=r4a
            sa(15)=-r4b
            ra(16)=r4a
            sa(16)=r4b
      end if
c
      do 200 l=1,nint
c
            r=ra(l)
            s=sa(l)
            shl2(1,1,l)=zero
            shl2(2,1,l)=zero
            shl2(3,1,l)=pt25
            shl2(1,2,l)=zero
            shl2(2,2,l)=zero
            shl2(3,2,l)=-pt25
            shl2(1,3,l)=zero
            shl2(2,3,l)=zero
            shl2(3,3,l)=pt25
            shl2(1,4,l)=zero
            shl2(2,4,l)=zero
            shl2(3,4,l)=-pt25
            if(nen.eq.9) then
                  onepr=one+r
                  onemr=one-r
                  oneps=one+s
                  onems=one-s
		  onemrs=one-r*r
		  onemss=one-s*s
                  shl2(1,5,l)=-oneps
                  shl2(2,5,l)=zero
                  shl2(3,5,l)=-r
                  shl2(1,6,l)=zero
                  shl2(2,6,l)=-onemr
                  shl2(3,6,l)=s
                  shl2(1,7,l)=-onems
                  shl2(2,7,l)=zero
                  shl2(3,7,l)=r
                  shl2(1,8,l)=zero
                  shl2(2,8,l)=-onepr
                  shl2(3,8,l)=-s
                  shl2(1,9,l)=-two*onemss
                  shl2(2,9,l)=-two*onemrs
                  shl2(3,9,l)=four*r*s
c
                  do 1111 k=5,8
                        do 2222 i=1,3
                              shl2(i,k,l)=shl2(i,k,l)-pt5*shl2(i,9,l)
2222                    continue
1111              continue
c
                  do 3333 i=1,3
                        shl2(i,1,l)=shl2(i,1,l)
     &                  -pt5*(shl2(i,5,l)+shl2(i,8,l))-pt25*shl2(i,9,l)
                        shl2(i,2,l)=shl2(i,2,l)
     &                  -pt5*(shl2(i,6,l)+shl2(i,5,l))-pt25*shl2(i,9,l)
                        shl2(i,3,l)=shl2(i,3,l)
     &                  -pt5*(shl2(i,7,l)+shl2(i,6,l))-pt25*shl2(i,9,l)
                        shl2(i,4,l)=shl2(i,4,l)
     &                  -pt5*(shl2(i,8,l)+shl2(i,7,l))-pt25*shl2(i,9,l)
3333              continue
            end if
            if (nen.eq.16) then
                  onemrsq=one-r*r
                  onemssq=one-s*s
                  onep3r=one+three*r
                  onem3r=one-three*r
                  onep3s=one+three*s
                  onem3s=one-three*s
c
c	not inplemented
c
c
            end if
  200 continue
c
      return
      end
c*******************************************************************
      subroutine shgqs(xl,det,shl,shg,nint,nel,neg,quad,nen,
     &                 shlnode,nenode)
c
c.... program to calculate global derivatives of shape functions and
c        jacobian determinants for a  quadrilateral element
c
c        xl(j,i)    = global coordinates
c        det(l)     = jacobian determinant
c        shl(1,i,l) = local ("xi") derivative of shape function
c        shl(2,i,l) = local ("eta") derivative of shape function
c        shl(3,i,l) = local  shape function
c        shg(1,i,l) = x-derivative of shape function
c        shg(2,i,l) = y-derivative of shape function
c        shg(3,i,l) = shl(3,i,l)
c        xs(i,j)    = jacobian matrix
c                 i = local node number or global coordinate number
c                 j = global coordinate number
c                 l = integration-point number
c              nint = number of integration points, eq. 1 or 4
c
      implicit real*8 (a-h,o-z)
c
c.... remove above card for single precision operation
c
      logical quad
      dimension xl(2,*),det(*),shl(3,nen,*),shg(3,nen,*),xs(2,2),
     &          shlnode(3,nenode,*)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
      common /iounit/ iin,iout,iecho,ioupp,itest1,itest2 
c
      call move(shg,shl,3*nen*nint)
c
      do 700 l=1,nint
c
      if (.not.quad) then
         do 100 i=1,3
         shg(i,3,l) = shl(i,3,l) + shl(i,4,l)
         shg(i,4,l) = zero
  100    continue
      endif
c
      do 300 j=1,2
      do 200 i=1,2
      xs(i,j) = rowdot(shlnode(i,1,l),xl(j,1),3,2,nenode)
  200 continue
  300 continue
c
      det(l) = xs(1,1)*xs(2,2)-xs(1,2)*xs(2,1)
      if (det(l).le.zero) then
         write(iecho,1000) nel,neg
         stop
      endif
c
      do 500 j=1,2
      do 400 i=1,2
      xs(i,j) = xs(i,j)/det(l) 
  400 continue
  500 continue
c
      do 600 i=1,nen
        temp = xs(2,2)*shg(1,i,l) - xs(1,2)*shg(2,i,l)
        shg(2,i,l) = - xs(2,1)*shg(1,i,l) + xs(1,1)*shg(2,i,l)
        shg(1,i,l) = temp
  600 continue
c
  700 continue
c
      return
c
 1000 format(///,'shgqs - non-positive determinant - element ',i10,
     &          ' in element group  ',i10,'  shgqs')
      end
c**** new *********************************************************************
c**** new *********************************************************************
      subroutine shgtqsd(ien,xl,shl,shg,
     &                 nside,nint,nints,nel,neg,nen,
     &                 shlnode,nenode)
c
c.... program to calculate global derivatives of shape functions and
c        jacobian determinants for a  quadrilateral element
c
c        xl(j,i)    = global coordinates
c        det(l)     = jacobian determinant
c        shl(1,i,l) = local ("xi") derivative of shape function
c        shl(2,i,l) = local ("eta") derivative of shape function
c        shl(3,i,l) = local  shape function
c        shg(1,i,l) = x-derivative of shape function
c        shg(2,i,l) = y-derivative of shape function
c        shg(3,i,l) = shl(3,i,l)
c        xs(i,j)    = jacobian matrix
c                 i = local node number or global coordinate number
c                 j = global coordinate number
c                 l = integration-point number
c              nint = number of integration points, eq. 1 or 4
c
      implicit real*8 (a-h,o-z)
c
c.... remove above card for single precision operation
c
      dimension xl(2,*),shl(3,nen,*),shg(3,nen,*),xs(2,2),
     &          shlnode(3,nenode,*),ien(*)
      dimension igas(100)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
      common /iounit/ iin,iout,iecho,ioupp,itest1,itest2 
c
cc      call move(shg,shl,3*nen*nint)
c
c   renumeracao
c
      ns = 1
      nl1=ien(1)
      nl2=ien(2)
c      
      if(nl2.gt.nl1) then 
        do nn=1,nints
           ng = (ns-1)*nints + nn
           igas(ng) = ng
	  end do
	else
        do  nn=1,nints
            ng = (ns-1)*nints + nn
            ngs = (ns-1)*nints
	      igas(ng) = ngs + nints + 1 - nn
	    end do 
      end if 
c
c
      ns = 2
      nl1=ien(2)
      nl2=ien(3)
c      
      if(nl2.gt.nl1) then 
        do  nn=1,nints
           ng = (ns-1)*nints + nn
           igas(ng) = ng
	  end do
	else
        do  nn=1,nints
            ng = (ns-1)*nints + nn
            ngs = (ns-1)*nints
	      igas(ng) = ngs + nints + 1 - nn
	    end do 
      end if 
c
c
      if(nside.eq.3) then 
        ns = 3
        nl1=ien(3)
        nl2=ien(1)
c      
        if(nl2.gt.nl1) then 
          do  nn=1,nints
             ng = (ns-1)*nints + nn
             igas(ng) = ng
	    end do
	   else
           do  nn=1,nints
             ng = (ns-1)*nints + nn
             ngs = (ns-1)*nints
	       igas(ng) = ngs + nints + 1 - nn
	     end do 
         end if 
c
      end if 
c
c
      if(nside.eq.4) then 
        ns = 3
        nl1=ien(3)
        nl2=ien(4)
c      
        if(nl2.gt.nl1) then 
          do  nn=1,nints
             ng = (ns-1)*nints + nn
             igas(ng) = ng
	    end do
	  else
          do  nn=1,nints
              ng = (ns-1)*nints + nn
              ngs = (ns-1)*nints
	        igas(ng) = ngs + nints + 1 - nn
	    end do 
         end if 
c 
        ns = 4
        nl1=ien(4)
        nl2=ien(1)
c      
        if(nl2.gt.nl1) then 
          do  nn=1,nints
             ng = (ns-1)*nints + nn
             igas(ng) = ng
	    end do
	  else
          do  nn=1,nints
              ng = (ns-1)*nints + nn
              ngs = (ns-1)*nints
	        igas(ng) = ngs + nints + 1 - nn
	    end do 
         end if 
c 
      end if
c 
c
      do 700 l=1,nint
c
c
      do 300 j=1,2
      do 200 i=1,2
      xs(i,j) = rowdot(shlnode(i,1,l),xl(j,1),3,2,nenode)
  200 continue
  300 continue
c
      deta = xs(1,1)*xs(2,2)-xs(1,2)*xs(2,1)
      if (deta.le.zero) then
         write(iecho,1000) nel,neg
         stop
      endif
c
      do 500 j=1,2
      do 400 i=1,2
      xs(i,j) = xs(i,j)/deta 
  400 continue
  500 continue
c
      do 600 i=1,nen
c
        lga = igas(l)
        shg(3,i,lga) = shl(3,i,l)        
        temp = xs(2,2)*shl(1,i,l) - xs(1,2)*shl(2,i,l)
        shg(2,i,lga) = - xs(2,1)*shl(1,i,l) + xs(1,1)*shl(2,i,l)
        shg(1,i,lga) = temp
  600 continue
c
  700 continue
c      
      return
c
 1000 format(///,'shgtqsd - non-positive determinant - element ',i10,
     &          ' in element group  ',i10,'  shgtqsd')
      end
c**** new *********************************************************************
c**** new **********************************************************************
      subroutine elmlib(ntype,mpnpar,itask,neg)
c
c.... program to call element routines
c
      common a(1)
c
      go to (10,20) ntype
c
c    POISSON PROBLEM - KINEMATIC FORMULATION
c
  10  continue
      call pflux(itask,a(mpnpar),a(mpnpar+16),neg)
  20  continue
      return
      end
c**** new **********************************************************************
      subroutine pflux(itask,npar,mp,neg)
c___________________________________________________________
c
c..... program to set storage and call tasks for the 
c             primal mixed Poisson  problem
c     with continuous temperature and discontinuous flux
c___________________________________________________________
      dimension npar(*),mp(*)
      common /bpoint/ mfirst,mlast,ilast,mtot,iprec
      common /iounit/ iin,iout,iecho,ioupp,itest1,itest2 
      common /info  / iexec,iprtin,irank,nsd,numnp,ndof,nlvect,
     &                numeg,nmultpd,nmultpc,nedge

      common /spoint/ mpd,mpdprj,mpx,mpidc,mpiedge,mpf,mpdiag,
     &                mpngrp,mpalhs,mpdlhs,mpbrhs,mped,index

      character*4 ia
      common a(1)
      common /dictn/ ia(1)
c
c
      mw     = 1
      mdet   = 2
      mshl   = 3
      mshg   = 4
c
      mc     = 5
      mgrav  = 6
      mien   = 7
      mmat   = 8
      mlm    = 9
c
      mxl    = 10
      mdl    = 11
c
c	pointers for condensation
c
      mipar  = 12
      mlado  = 13
c
      melefd = 14
      melred = 15
	mdlp   = 16
c
c    inetrais nas arestas
c
      mdetpn = 17
      mshlpn = 18
      mshgpn = 19
      mwpn   = 20
c
c.....pointers for boundary itegrals
c
      mdside = 21
      mxls   = 22
      midlsd = 23
c
c    geoemetria das arstas
c      
      mdetn  = 24
      mshln  = 25
      mshgn  = 26
      mwn    = 27
c
c    valores na fronteira
c
      mdetb  = 28
      mshlb  = 29
      mshgb  = 30
	mddis  = 31
c
      mdetp  = 32
      mshlp  = 33
      mshgp  = 34
      mwp    = 35
c
c    matrizes e vetores da formulação hibrida
c
	melmbb  = 36
	melmcb  = 37
c
      melfbb  = 38
      melfcb  = 39
c
c
      mshlpsd = 40
      mshgpsd = 41
c
      melmdb  = 42
c
      mshlsd = 43
      mshgsd = 44
c
      maelm = 45 
      mfelm = 46
      mxvla = 47
      mxvlb = 48
      mxvle = 49
c
	melmbc  = 50
c
	mddisa = 51
c
      mindno = 52
      mindedg= 53
      miedgto= 54
c
      mdlprj = 55
c
      ntype  = npar( 1)
      numel  = npar( 2)
      numat  = npar( 3)
      nint   = npar( 4)
      nen    = npar( 5)
	nenp   = npar( 7)
	npars  = npar( 8)
      nints  = npar( 9)
c
      if(nen.eq.0) nen=4
c
      if(nen.eq.4)  nenlad = 2
      if(nen.eq.8)  nenlad = 3
      if(nen.eq.9)  nenlad = 3
      if(nen.eq.16) nenlad = 4
c
      if(nen.eq.3)  nenlad = 2
      if(nen.eq.6)  nenlad = 3
      if(nen.eq.10) nenlad = 4
c
      if(nenp.eq.4) then
          nnods = 2
          nside = 4
       end if
c
      if(nenp.eq.8.or.nenp.eq.9) then
          nnods = 3
          nside = 4
      end if
c
      if(nenp.eq.16) then
          nnods = 4
          nside = 4
      end if
      if(nenp.eq.25) then
          nnods = 5
          nside = 4
      end if
c
      if(nenp.eq.36) then
          nnods = 6
          nside = 4
      end if
      if(nenp.eq.49) then
          nnods = 7
          nside = 4
      end if
      if(nenp.eq.64) then
          nnods = 8
          nside = 4
      end if
c
c
      if(nenp.eq.3) then
	    nnods  = 2
	    nside  = 3
      end if
c
      if(nenp.eq.6) then
	    nnods  = 3
	    nside  = 3
      end if
c
      if(nenp.eq.10) then
	    nnods  = 4
	    nside  = 3
      end if
c
      nintb=nside*nints
c
c
c      parameters for post processing
c
	 nodsp=(npars-1)*nside ! para o multiplicador continuo
	 nodspd = npars*nside
c
c.... set element parameters
c
      ndimc  = 8
      ned    = 1
      ncon   = 2
      nee    = nodsp*ned
	neep   = nenp*ned
cc      necon  = nencon*ncon
      neesq  = nee*nee
      nesd   = 2
      nrowsh = 3
      ngrav  = 3

c
       if (itask.eq.1) then
c
c....... set memory pointers
c
c
c        note:  the mp array is stored directly after the npar array,
c               beginning at location mpnpar + 16 of blank common.
c               the variable "junk" is not used subsequently.
c
         junk       = mpoint('mp      ',80     ,0     ,0     ,1)
c
         mp(mw    ) = mpoint('w       ',nint   ,0     ,0     ,iprec)
         mp(mdet  ) = mpoint('det     ',nint   ,0     ,0     ,iprec)
         mp(mshl  ) = mpoint('shl     ',nrowsh ,nen   ,nint  ,iprec)
         mp(mshg  ) = mpoint('shg     ',nrowsh ,nen   ,nint  ,iprec)
         mp(mc    ) = mpoint('c       ',ndimc  ,numat ,0     ,iprec)
         mp(mgrav ) = mpoint('grav    ',ngrav  ,0     ,0     ,iprec)
c
         mp(mien  ) = mpoint('ien     ',nen    ,numel ,0     ,1)
         mp(mmat  ) = mpoint('mat     ',numel  ,0     ,0     ,1)
c
         mp(mlm   ) = mpoint('lm      ',ned    ,nodsp ,numel ,1)
         mp(mxl   ) = mpoint('xl      ',nesd   ,nen   ,0     ,iprec)
         mp(mdl   ) = mpoint('dl      ',ned    ,nodsp ,0     ,iprec)
c
         mp(mipar ) = mpoint('ipar    ',nodsp  ,numel ,0,     1)
         mp(mlado)  = mpoint('lado    ',nside  ,numel ,0     ,1)
c
c	
         mp(melefd) = mpoint('eleffd  ',nee    ,nee   ,0     ,iprec)
         mp(melred) = mpoint('elresd  ',nee    ,0     ,0     ,iprec)
         mp(mdlp  ) = mpoint('dlp     ',ned    ,nenp  ,0     ,iprec)
c
         mp(mdetpn) = mpoint('detpn   ',nints  ,0     ,0     ,iprec)
         mp(mshlpn) = mpoint('shlpn   ',nrowsh ,npars ,nints ,iprec)
         mp(mshgpn) = mpoint('shgpn   ',nrowsh ,npars ,nints ,iprec)
         mp(mwpn  ) = mpoint('wpn     ',nints  ,0     ,0     ,iprec)
c
         mp(mdside) = mpoint('idside  ',nside  ,npars ,0     ,1    )
         mp(mxls  ) = mpoint('xls     ',nesd   ,nenlad,0     ,iprec)
         mp(midlsd) = mpoint('idlsd   ',npars  ,0     ,0     ,1    )
         mp(mdetn ) = mpoint('detn    ',nints  ,0     ,0     ,iprec)
         mp(mshln ) = mpoint('shln    ',3      ,nnods ,nints ,iprec)
         mp(mshgn ) = mpoint('shgn    ',3      ,nnods ,nints ,iprec)
         mp(mwn   ) = mpoint('wn      ',nints  ,0     ,0     ,iprec)
c
         mp(mdetb ) = mpoint('detb    ',nints  ,0     ,0     ,iprec)
         mp(mshlb ) = mpoint('shlb    ',3      ,nenlad,nints ,iprec)
         mp(mshgb ) = mpoint('shgb    ',3      ,nenlad,nints ,iprec)
c
         mp(mddis ) = mpoint('ddis    ',ned    ,nenp  ,numel ,iprec )
         mp(mwp   ) = mpoint('wp      ',nint   ,0     ,0     ,iprec)
         mp(mdetp ) = mpoint('detp    ',nint   ,0     ,0     ,iprec)
         mp(mshlp ) = mpoint('shlp    ',nrowsh ,nenp  ,nint  ,iprec)
         mp(mshgp ) = mpoint('shgp    ',nrowsh ,nenp  ,nint  ,iprec)
c
c
         mp(melmbb) = mpoint('elmbb   ',neep   ,neep  ,0     ,iprec)
         mp(melmcb) = mpoint('elmcb   ',neep   ,nee   ,0     ,iprec)
         mp(melfbb) = mpoint('elfbb   ',neep   ,0     ,0     ,iprec)
         mp(melfcb) = mpoint('elfcb   ',nee    ,0     ,0     ,iprec)
c
c
         mp(mshlpsd) = mpoint('shlpsd  ',nrowsh ,nenp  ,nintb ,iprec)
         mp(mshgpsd) = mpoint('shgpsd  ',nrowsh ,nenp  ,nintb ,iprec)
c
         mp(melmdb) = mpoint('elmdb   ',nee    ,neep   ,0     ,iprec)
c
         mp(mshlsd) = mpoint('shlsd   ',nrowsh ,nen   ,nintb ,iprec)
         mp(mshgsd) = mpoint('shgsd   ',nrowsh ,nen   ,nintb ,iprec)
c
         mp(maelm) = mpoint('aelm    ',neep   ,neep  ,numel ,iprec)
         mp(mfelm) = mpoint('felm    ',neep   ,numel ,0     ,iprec)
         mp(mxvla) = mpoint('xvla    ',nints  ,nedge ,0     ,iprec)
         mp(mxvlb) = mpoint('xvlb    ',nints  ,nedge ,0     ,iprec)
         mp(mxvle) = mpoint('xvle    ',nints  ,nedge ,0     ,iprec)
         mp(melmbc)= mpoint('elmbc   ',nee    ,neep  ,0     ,iprec)
         mp(mddisa)= mpoint('ddisa   ',ned    ,nenp  ,numel ,iprec )
c
         mp(mindno) =  mpoint('indno   ',numnp  ,0     ,0     ,1    )
         mp(mindedg) = mpoint('indedg  ',nedge  ,0     ,0     ,1    )
         mp(miedgto) = mpoint('iedgto  ',nedge  ,0     ,0     ,1    )
c
         mp(mdlprj ) = mpoint('dlprj   ',ned    ,nodspd,0     ,iprec)
c
      endif
c
c.... task calls
c
      if (itask.gt.4) return
      go to (100,200,300,400),itask
c
  100 continue
c
c.... input element data ('input___')
c
      call flux1el(a(mp(mshl  )),a(mp(mw    )),
     &           a(mp(mc    )),a(mp(mgrav )),
     &           a(mp(mien  )),a(mp(mmat  )),
     &           a(mpidc     ),a(mp(mlm   )),
     &           a(mpdiag    ),a(mp(mipar )),
     &           a(mpx       ),a(mp(mlado)),
     &           a(mp(mshlpn)),a(mp(mwpn  )),
     &           a(mp(mshln )),a(mp(mwn   )),
     &           a(mp(mshlb )),a(mp(mshlp )),
     &           a(mp(mwp   )),a(mp(mdside)),
c     
     &           a(mp(mshlpsd)),
     &           a(mp(mshlsd)),
c
     &           a(mp(mindno)),
     &           a(mp(mindedg)),
     &           a(mpiedge)  ,
     &           a(mp(miedgto)),
     &           a(mp(midlsd)),
c     
     &           ntype ,numel ,numat ,
     &           nint  ,nrowsh,nesd  ,
     &           nen   ,ndof  ,ned   ,
     &           iprtin,numnp ,ncon  ,
     &           nencon,necon ,nints ,
     &           nnods ,nenlad,npars ,
     &           nenp  ,nside ,nodsp ,
     &           nedge ,nmultpd,nmultpc)
c
c
      return
c
  200 continue
c
      call flux1sd(a(mp(mien  )),a(mp(mmat  )),
     &           a(mp(mlm   )),
     &           a(mpidc     ),a(mped      ),
     &           a(mpdiag    ),a(mp(mipar )),
c
     &           a(mp(mlado)) , 
c
     &           a(mp(mdside)),
     &           a(mp(mindno)),
     &           a(mp(mindedg)),
     &           a(mpiedge)  ,
     &           a(mp(miedgto)),
     &           a(mp(midlsd)),
c     
     &           ntype ,numel ,numat ,
     &           nint  ,nrowsh,nesd  ,
     &           nen   ,ndof  ,ned   ,
     &           iprtin,numnp ,ncon  ,
     &           nencon,necon ,nints ,
     &           nnods ,nenlad,npars ,
     &           nenp  ,nside ,nodsp ,
     &           nedge ,nmultpd,nmultpc)
c
      return
c
  300 continue
c
c
c.... form element stiffnes matrix and force vector and
c     assemble them into global left and right hand-side
c     ('form_lrs')
c
      call flux2(a(mp(mien  )),a(mpx       ),a(mp(mxl   )),
     &           a(mpd       ),a(mp(mdl   )),a(mp(mmat  )),
     &           a(mp(mdet  )),a(mp(mshl  )),a(mp(mshg  )),
     &           a(mp(mw    )),a(mp(mc    )),
     &           a(mpalhs    ),a(mpdlhs    ),
     &           a(mpbrhs    ),a(mpdiag    ),a(mp(mlm   )),
     &           a(mp(mgrav )),a(mp(mipar )),a(mp(mlado )),
     &           a(mp(melefd)),a(mp(melred)),a(mp(mshln )),
     &           a(mp(mshgn )),a(mp(mwn   )),a(mp(mdetn )),
     &           a(mp(mdetb )),a(mp(mshlb )),a(mp(mshgb )),
     &           a(mp(mdetpn)),a(mp(mshlpn)),a(mp(mshgpn)),
     &           a(mp(mdside)),a(mp(mxls  )),a(mp(midlsd)),
     &           a(mp(mddis )),a(mp(mdetp )),
     &           a(mp(mshlp )),a(mp(mshgp )),
c     
     &           a(mp(melmbb)),
     &           a(mp(melmcb)),
     &           a(mp(melmbc)),
     &           a(mp(melfbb)),a(mp(melfcb)),
     &           a(mp(melmdb)),
c
     $           a(mped      ),
c
     &           a(mp(mshlpsd)),
     &           a(mp(mshgpsd)),
c
     &           a(mp(mshlsd)),
     &           a(mp(mshgsd)),
c
c
     &           a(mp(maelm)),
     &           a(mp(mfelm)),
     &           a(mp(mxvla)),
     &           a(mp(mxvlb)),
c     
     &           a(mpdprj  ),
c
     &           a(mpiedge)  ,
     &           a(mpidc     ),
c     
     &                 numel ,neesq ,nen   ,nsd   ,
     &                 nesd  ,nint  ,neg   ,nrowsh,
     &                 ned   ,nee   ,numnp ,ndof  ,
     &                 ncon  ,nencon,necon ,neep  ,
     &                 nints ,nnods ,nenlad,npars ,
     &                 nside ,nenp  ,nedge ,nodsp ,
     &                 index ,nmultpd,nmultpc)

c
      return
c
  400 continue
      call flux3(a(mp(mien  )),a(mpx       ),a(mp(mxl   )),
     &           a(mpd       ),a(mp(mdl   )),a(mp(mmat  )),
     &           a(mp(mdet  )),a(mp(mshl  )),a(mp(mshg  )),
     &           a(mp(mw    )),a(mp(mc    )),
     &           a(mp(mgrav )),a(mp(mipar )),a(mp(mlado )),
     &           a(mp(melefd)),a(mp(melred)),a(mp(mshln )),
     &           a(mp(mshgn )),a(mp(mwn   )),a(mp(mdetn )),
     &           a(mp(mdetb )),a(mp(mshlb )),a(mp(mshgb )),
     &           a(mp(mdetpn)),a(mp(mshlpn)),a(mp(mshgpn)),
     &           a(mp(mdside)),a(mp(mxls  )),a(mp(midlsd)),
     &           a(mp(mddis )),a(mp(mdetp )),
     &           a(mp(mshlp )),a(mp(mshgp )),
c     
     &           a(mp(melmbb)),
     &           a(mp(melmcb)),
     &           a(mp(melfbb)),a(mp(melfcb)),
     &           a(mp(melmdb)),
c
c
     $             a(mped      ),
c
c
     &           a(mp(mshlpsd)),
     &           a(mp(mshgpsd)),
c
     &           a(mp(mshlsd)),
     &           a(mp(mshgsd)),
c
     &           a(mp(maelm)),
     &           a(mp(mfelm)),
     &           a(mp(mxvla)),
     &           a(mp(mxvlb)),
c     
     &           a(mpiedge)  ,
c     
     &           numel ,neesq ,nen   ,
     &           nsd   ,nesd  ,nint  ,
     &           neg   ,nrowsh,ned   ,
     &           nee   ,numnp ,ndof  ,
     &           ncon  ,nencon,necon ,
     &           neep  ,nints ,nnods ,
     &           nenlad,npars ,nside ,
     &           nenp  ,nodsp ,index ,
     &                  nedge)
c
      read(iin,*) index,iwrite
c
      if(index.gt.0) then
c      write(*,*) "Aqui 3"
      call flnormc(a(mpdp(mien )),a(mpx       ),a(mp(mxl   )),
     &            a(mpd      ),a(mp(mdl   )),a(mp(mmat  )),
     &            a(mp(mc   )),a(mp(mipar )),
     &            a(mp(mdlp )),
     &            a(mp(mdet  )),
     &            a(mp(mshl )),a(mp(mshg  )),a(mp(mw    )),
     &            a(mp(mddis )),a(mp(mdetp)),a(mp(mshlp )),
     &            a(mp(mshgp )),
c
     &            a(mp(mshln )),a(mp(mshgn )),
     &            a(mp(mdetb )),a(mp(mshlb )),a(mp(mshgb )),
     &            a(mp(mdetpn)),a(mp(mshlpn)),a(mp(mshgpn)),
     &            a(mp(mdside)),a(mp(mxls  )),a(mp(midlsd)),
     &            a(mp(mgrav )),a(mp(mwn   )),
c
     &            numel ,neesq ,nen   ,nsd   ,
     &            nesd  ,nint  ,neg   ,nrowsh,
     &            ned   ,nee   ,numnp ,ndof  ,
     &            ncon  ,nencon,necon ,index ,
     &            nints ,iwrite,iout  ,nenp  ,
     &            nside ,nnods ,nenlad,npars ,
     &            nmultpd,nmultpc,nodsp )
c
            go to 12345
c
      call flninter(a(mp(mien )),a(mpx       ),a(mp(mxl   )),
     &            a(mpd      ),a(mp(mdl   )),a(mp(mmat  )),
     &            a(mp(mc   )),a(mp(mipar )),
     &            a(mp(mdlp )),
     &            a(mp(mdet  )),
     &            a(mp(mshl )),a(mp(mshg  )),a(mp(mw    )),
     &            a(mp(mddis )),a(mp(mdetp)),a(mp(mshlp )),
     &            a(mp(mshgp )),
c
     &            a(mp(mshln )),a(mp(mshgn )),
     &            a(mp(mdetb )),a(mp(mshlb )),a(mp(mshgb )),
     &            a(mp(mdetpn)),a(mp(mshlpn)),a(mp(mshgpn)),
     &            a(mp(mdside)),a(mp(mxls  )),a(mp(midlsd)),
     &            a(mp(mgrav )),a(mp(mwn   )),
c
     &            numel ,neesq ,nen   ,nsd   ,
     &            nesd  ,nint  ,neg   ,nrowsh,
     &            ned   ,nee   ,numnp ,ndof  ,
     &            ncon  ,nencon,necon ,index ,
     &            nints ,iwrite,itest1,nenp  ,
     &            nside ,nnods ,nenlad,npars ,
     &            nmultpd,nmultpc,nodsp )
       end if
c
12345 continue
c
c      projecalo local
c
c      write(*,*) "Aqui 4"
      call flux30(a(mp(mien  )),a(mpx       ),a(mp(mxl   )),
     &           a(mpdprj    ),a(mp(mdlprj )),a(mp(mmat  )),
     &           a(mp(mdet  )),a(mp(mshl  )),a(mp(mshg  )),
     &           a(mp(mw    )),a(mp(mc    )),
     &           a(mp(mgrav )),a(mp(mipar )),a(mp(mlado )),
     &           a(mp(melefd)),a(mp(melred)),a(mp(mshln )),
     &           a(mp(mshgn )),a(mp(mwn   )),a(mp(mdetn )),
     &           a(mp(mdetb )),a(mp(mshlb )),a(mp(mshgb )),
     &           a(mp(mdetpn)),a(mp(mshlpn)),a(mp(mshgpn)),
     &           a(mp(mdside)),a(mp(mxls  )),a(mp(midlsd)),
     &           a(mp(mddis )),a(mp(mdetp )),
     &           a(mp(mshlp )),a(mp(mshgp )),
c     
     &           a(mp(melmbb)),
     &           a(mp(melmcb)),
     &           a(mp(melfbb)),a(mp(melfcb)),
     &           a(mp(melmdb)),
c
c
     $             a(mped      ),
c
c
     &           a(mp(mshlpsd)),
     &           a(mp(mshgpsd)),
c
     &           a(mp(mshlsd)),
     &           a(mp(mshgsd)),
c
     &           a(mp(maelm)),
     &           a(mp(mfelm)),
     &           a(mp(mxvla)),
     &           a(mp(mxvlb)),
c     
     &           numel ,neesq ,nen   ,
     &           nsd   ,nesd  ,nint  ,
     &           neg   ,nrowsh,ned   ,
     &           nee   ,numnp ,ndof  ,
     &           ncon  ,nencon,necon ,
     &           neep  ,nints ,nnods ,
     &           nenlad,npars ,nside ,
     &           nenp  ,nodsp ,index ,
     &                  nedge)
c
      if(index.gt.0) then
      call flnorm(a(mp(mien )),a(mpx       ),a(mp(mxl   )),
     &            a(mpdprj  ),a(mp(mdlprj )),a(mp(mmat  )),
     &            a(mp(mc   )),a(mp(mipar )),
     &            a(mp(mdlp )),
     &            a(mp(mdet  )),
     &            a(mp(mshl )),a(mp(mshg  )),a(mp(mw    )),
     &            a(mp(mddis )),a(mp(mdetp)),a(mp(mshlp )),
     &            a(mp(mshgp )),
c
     &            a(mp(mshln )),a(mp(mshgn )),
     &            a(mp(mdetb )),a(mp(mshlb )),a(mp(mshgb )),
     &            a(mp(mdetpn)),a(mp(mshlpn)),a(mp(mshgpn)),
     &            a(mp(mdside)),a(mp(mxls  )),a(mp(midlsd)),
     &            a(mp(mgrav )),a(mp(mwn   )),
c
     &            a(mp(mlado )),
c
     &            numel ,neesq ,nen   ,nsd   ,
     &            nesd  ,nint  ,neg   ,nrowsh,
     &            ned   ,nee   ,numnp ,ndof  ,
     &            ncon  ,nencon,necon ,index ,
     &            nints ,iwrite,ioupp  ,nenp  ,
     &            nside ,nnods ,nenlad,npars ,
     &            nmultpd,nmultpc,nodsp,nodspd )
      end if
c
      return
      end
c
c**** new **********************************************************************
      subroutine flux1el(shl   ,w     ,       
     &                 c     ,grav  ,
     &                 ien   ,mat   ,
     &                 idc   ,lm    ,
     &                 idiag ,ipar  ,
     &                 x     ,lado  ,
     &                 shlpn ,wpn   ,
     &                 shln  ,wn    ,
     &                 shlb  ,shlp  ,
     &                 wp    ,idside,
c     
     &                 shlpsd,     
     &                 shlsd,
c
     &                 indno ,
     &                 indedg,
     &                 iedge ,
     &                 iedgto,
     &                 idlsd ,
c     
c          
     &           ntype ,numel ,numat ,
     &           nint  ,nrowsh,nesd  ,
     &           nen   ,ndof  ,ned   ,
     &           iprtin,numnp ,ncon  ,
     &           nencon,necon ,nints ,
     &           nnods ,nenlad,npars ,
     &           nenp  ,nside ,nodsp ,
     &           nedge ,nmultpd,nmultpc)
c
c.... program to read, generate and write data for the
c        four-node quadrilateral, elastic continuum element
c
      implicit real*8 (a-h,o-z)
c                                                                       
c.... remove above card for single-precision operation               
c                                                                       
      dimension shl(nrowsh,nen,*),w(*),
     &          c(8,*),grav(*),ien(nen,*),mat(*),
     &          idc(ndof,*),lm(ned,nodsp,*),idiag(*),
     &          ipar(nodsp,*),idside(nside,*),
     &          x(nesd,*),lado(nside,*),
     $          shlpn(nrowsh,npars,*),
     &          wpn(*),shlp(nrowsh,nenp,*),wp(*),
     &          shln(nrowsh,nnods,*),wn(*),shlb(nrowsh,nenlad,*)
      dimension shlpsd(nrowsh,nenp,*)
      dimension shlsd(nrowsh,nen,*)
      dimension indno(*),indedg(*),iedgto(*),iedge(npars,*),idlsd(*)
c
      common /iounit/ iin,iout,iecho,ioupp,itest1,itest2 
      common /colhtc / neqc
c
      write(iecho,1000) ntype,numel,numat
      write(iecho,2000) nint
c
c     geomeria 1D / Lagrange
c
c      generation of local shape functions and weight values
c
      call oneshl(shlb,wp,nints,nenlad)
c
      call oneshl(shln,wn,nints,nnods)
      call oneshl(shlpn,wpn,nints,npars)
c
	if(nen.eq.3.or.nen.eq.6.or.nen.eq.10) then
        call shlt(shl,w,nint,nen)
        call shlt(shlp,wp,nint,nenp)
        call shltpbk(shlpsd,nenp,nside,nnods,nints)
        call shltpbk(shlsd,nen,nside,nenlad,nints)
	else
        call shlq(shl,w,nint,nen)
        call shlqpk(shlp,wp,nint,nenp)
        call shlqpbk(shlpsd,nenp,nside,nnods,nints)
        call shlqpbk(shlsd,nen,nside,nenlad,nints)
	end if
c
      nintb=nside*nints
c
c      read material properties
c
      call fluxmx(c,numat)
c
c	constant body forces
c
      read (iin,5000) (grav(i),i=1,3)
      write (iecho,6000) (grav(i),i=1,3)
c
c    generation of conectivities
c
      call genel(ien,mat,nen)
c
      if (iprtin.eq.0) call prntel(mat,ien,nen,numel)
c
c   generation of conectivety for element multipliers
c
      nenlb = nenp
      call genside(idside,nside,nenlb)
      call genelad(lado,nside)
      if (iprtin.eq.0) call prntels(mat,lado,nside,numel)
c
c   geraçao das conectividades para o multiplicador continuo
c
      do i=1,numnp
        indno(i) = 0
      end do
c
      do n=1,nedge
        indedg(n) = 0
      end do
c
      keq = 0
      kedg = 0
c
c
      do 1111 nel=1,numel
c
      do 1100 ns=1,nside
        neledg  = lado(ns,nel)
c
c
        ns1 = idside(ns,1)
        ns2 = idside(ns,2)
        nl1 = ien(ns1,nel)
        nl2 = ien(ns2,nel)
c      
      if(nl2.gt.nl1) then
        idlsd(1) = 1
        idlsd(2) = 2 
 	  id3 = npars-2
	  if(id3.gt.0) then
          do  nn=3,npars
           idlsd(nn) = nn
	    end do
	   end if 
	else
	  idlsd(1) = 2
	  idlsd(2) = 1
	  id3 = npars-2
	  if(id3.gt.0) then
	    do il=3,npars
	      idlsd(il) = npars + 3 - il
	    end do 
	  end if
      end if 
c
      i1 = idlsd(1)
      i2 = idlsd(2)
c
      if(indedg(neledg).eq.0) then
c
      kedg = kedg + 1
      indedg(neledg) = kedg
      iedgto(kedg) = neledg
c
         if(indno(nl1).eq.0) then
            keq = keq + 1
            indno(nl1) = keq
            iedge(i1,neledg) = keq
         else
c
            iedge(i1,neledg) = indno(nl1)
c
         end if
c
         if(indno(nl2).eq.0) then
              keq = keq + 1
              indno(nl2) = keq 
              iedge(i2,neledg) = keq
         else
              iedge(i2,neledg) = indno(nl2)
         end if
c
         if(npars.gt.2) then
             do i=3,npars
                keq = keq + 1 
                isv = idlsd(i)
                iedge(isv,neledg) = keq
             end do
          end if
c
      end if
 1100 continue
c
 1111 continue
c
      write(44,*) numel,nedge,npars
      do nd = 1,nedge
       write(44,4444) nd,indedg(nd),iedgto(nd),(iedge(i,nd),i=1,npars)
4444   format(10 i8)
      end do 
c
      do nel=1,numel 
c
      do ns=1,nside
        nd  = lado(ns,nel)
        write(54,*) nel,ns,(iedge(i,nd),i=1,npars)
      end do
      end do
c
      nmultpc = keq      
c
      return
c
 1000 format(///,
     &' d u a l   h i b r i d   m i x e d   f o r m u l a t i o n ',///
     &//5x,' element type number . . . . . . . . . . .(ntype ) = ',i10
     &//5x,' number of elements  . . . . . . . . . . .(numel ) = ',i10
     &//5x,' number of element material sets . . . . .(numat ) = ',i10)
2000  format(
     &//5x,' numerical integration points  . . . . . .(nint  ) = ',i10)
 5000 format(8f10.0)
 6000 format(////' ',
     &' g r a v i t y   v e c t o r   c o m p o n e n t s      ',//5x,
     &' W-1 direction (sinxsiny)  . . . . . . . = ',      1pe15.8//5x,
     &' W-2 direction (cosxcosy)  . . . . . . . = ',      1pe15.8//5x,
     &' W-3 direction (Polinomial). . . . . . . = ',      1pe15.8//5x)
7000  format(i10)
c
      end
c**** new **********************************************************************
      subroutine flux1sd(ien   ,mat   ,
     &                 lm    ,
     &                 idc   ,ideg  ,
     &                 idiag ,ipar  ,
c
     &                 lado  ,
c
     &                 idside,
     &                 indno ,
     &                 indedg,
     &                 iedge ,
     &                 iedgto,
     &                 idlsd ,
c     
     &                 ntype ,numel ,numat ,
     &                 nint  ,nrowsh,nesd  ,
     &                 nen   ,ndof  ,ned   ,
     &                 iprtin,numnp ,ncon  ,
     &                 nencon,necon ,nints ,
     &                 nnods,nenlad ,npars ,
     &                 nenp ,nside  ,nodsp ,
     &                 nedge,nmultpd,nmultpc)
c
c.... program to read, generate and write data for the
c        four-node quadrilateral, elastic continuum element
c
      implicit real*8 (a-h,o-z)
c                                                                       
c.... remove above card for single-precision operation               
c                                                                       
      dimension ien(nen,*),mat(*),
     &          lm(ned,nodsp,*),idiag(*),
     &          ipar(nodsp,*),idside(nside,*)
c
      dimension lado(nside,*),iedge(npars,*)
      dimension indno(*),indedg(*),iedgto(*),idlsd(*)
      dimension idc(ndof,*),ideg(ndof,*)
c
      common /iounit/ iin,iout,iecho,ioupp,itest1,itest2 
      common /colhtc / neqc
c
c condições de contorno para o multiplicador continuo
c
      call setbcedgec(ideg,idc,iedge,indno,indedg,iedgto,
     &               npars,nedge,ndof,nmultpc,neqc,iprtin)
c
c   generation of multiplier equations
c
      call genelpar(ipar,ien,lado,idside,iedge,indedg,idlsd,
     &              nen,nside,nodsp,numel,npars)
c
      if (iprtin.eq.0) call prntelp(mat,ipar,nodsp,numel)
c
c     generation of lm array
c
      call formlm(idc,ipar,lm,ndof,ned,nodsp,numel)
c
c     modification of idiag array
c
      call colht(idiag,lm,ned,nodsp,numel,neqc)
c
c     Neuman B.c. not implemented YET
c
c
      return
c
 1000 format(///,
     &' d u a l   h i b r i d   m i x e d   f o r m u l a t i o n ',///
     &//5x,' element type number . . . . . . . . . . .(ntype ) = ',i10
     &//5x,' number of elements  . . . . . . . . . . .(numel ) = ',i10
     &//5x,' number of element material sets . . . . .(numat ) = ',i10)
2000  format(
     &//5x,' numerical integration points  . . . . . .(nint  ) = ',i10)
 5000 format(8f10.0)
 6000 format(////' ',
     &' g r a v i t y   v e c t o r   c o m p o n e n t s      ',//5x,
     &' W-1 direction (sinxsiny)  . . . . . . . = ',      1pe15.8//5x,
     &' W-2 direction (cosxcosy)  . . . . . . . = ',      1pe15.8//5x,
     &' W-3 direction (Polinomial). . . . . . . = ',      1pe15.8//5x)
7000  format(i10)
c
      end
c**** new **********************************************************************
      subroutine setbcedgec(ideg,idc,iedge,indno,indedg,iedgto,
     &                      npar,nedge,ndof,nmultpc,neqc,iprtin)
c
c.... program to read, generate and write boundary condition data
c        and establish equation numbers
c
      dimension idc(ndof,*),ideg(ndof,*),iedge(npar,*)
      dimension indno(*),indedg(*),iedgto(*)
c
      common /iounit/ iin,iout,iecho,ioupp,itest1,itest2 
      logical pflag
c
      call iclear(idc,ndof*nmultpc)
c
c    idc - prescribed dof
c
        do n=1,nedge 
           kedg = iedgto(n)
           do j=1,ndof
           if(ideg(j,kedg).ne.0) then
	       do is=1,npar
	         kk = iedge(is,kedg)
	         idc(j,kk) = ideg(j,kedg)
	       end do
	     end if
	     end do
	  end do
c
      if (iprtin.eq.0) then
         nn=0
         do 220 n=1,nmultpc
         pflag = .false.
c
         do 110 i=1,ndof
         if (idc(i,n).ne.0) pflag = .true.
  110    continue
c
         if (pflag) then      
            nn = nn + 1
            if (mod(nn,50).eq.1) write(55,1100) (i,i=1,ndof)
            write(55,2000) n,(idc(i,n),i=1,ndof)
         endif
  220    continue
      endif
c
c.... establish equation numbers
c
      neqc = 0
c
      do 400 n=1,nmultpc
c
      do 300 i=1,ndof
      if (idc(i,n).eq.0) then
         neqc = neqc + 1
         idc(i,n) = neqc
      else
         idc(i,n) = 1 - idc(i,n)
      endif
c
  300 continue
c
  400 continue
c
      return
c
 1000 format(//,10x,' e d g e   b o u n d a r y   c o n d i t i o n '//
     & 5x,'   node no.',3x,6(13x,'dof',i1:)//)
 1100 format(//,10x,' n o d a l   b o u n d a r y   c o n d i t i o n'//
     & 5x,'   node no.',3x,6(13x,'dof',i1:)//)
 2000 format(6x,i10,5x,6(5x,i10))
c
      end
c**** new **********************************************************************
c******************************************************************************
      subroutine genside(idside,nside,nen)                            
c                                                                       
c.... program to read and generate element node and material numbers    
c                                                                       
c         idside(nside,nnods) = element sides                        
c                                                                       
      dimension idside(nside,*)                            
c
c.....define idside
c
           if(nen.eq.4) then
              idside(1,1) = 1
              idside(1,2) = 2
c
              idside(2,1) = 2
              idside(2,2) = 3
c
              idside(3,1) = 3
              idside(3,2) = 4
c
              idside(4,1) = 4
              idside(4,2) = 1
           end if
c
            if(nen.eq.9) then
              idside(1,1) = 1
              idside(1,2) = 2
              idside(1,3) = 5
c
              idside(2,1) = 2
              idside(2,2) = 3
              idside(2,3) = 6
c
              idside(3,1) = 3
              idside(3,2) = 4
              idside(3,3) = 7
c
              idside(4,1) = 4
              idside(4,2) = 1
              idside(4,3) = 8
             end if
c
           if(nen.eq.16) then
c
c  lado 1
c
              idside(1,1) = 1
              idside(1,2) = 2
              idside(1,3) = 5
              idside(1,4) = 6
c
c   lado 2
c
              idside(2,1) = 2
              idside(2,2) = 3
              idside(2,3) = 7
              idside(2,4) = 8
c
c   lado 3
c             
              idside(3,1) = 3
              idside(3,2) = 4
              idside(3,3) = 9
              idside(3,4) = 10
c
c   lado 4
c              
              idside(4,1) = 4
              idside(4,2) = 1
              idside(4,3) = 11
              idside(4,4) = 12
c
c
           end if
c
           if(nen.eq.25) then
c
c     lado 1
c
              idside(1,1) = 1
              idside(1,2) = 2
              idside(1,3) = 5
              idside(1,4) = 6
              idside(1,5) = 7
c
c     lado 2
c
              idside(2,1) = 2
              idside(2,2) = 3
              idside(2,3) = 8
              idside(2,4) = 9
              idside(2,5) = 10
c
c     lado 3
c
              idside(3,1) = 3
              idside(3,2) = 4
              idside(3,3) = 11
              idside(3,4) = 12
              idside(3,5) = 13

c
c     lado 4
c
              idside(4,1) = 4
              idside(4,2) = 1
              idside(4,3) = 14
              idside(4,4) = 15
              idside(4,5) = 16

c
           end if
c
c
           if(nen.eq.36) then
c
c     lado 1
c
              idside(1,1) = 1
              idside(1,2) = 2
              idside(1,3) = 5
              idside(1,4) = 6
              idside(1,5) = 7
              idside(1,6) = 8
c
c     lado 2
c
              idside(2,1) = 2
              idside(2,2) = 3
              idside(2,3) = 9
              idside(2,4) = 10
              idside(2,5) = 11
              idside(2,6) = 12
c
c     lado 3
c
              idside(3,1) = 3
              idside(3,2) = 4
              idside(3,3) = 13
              idside(3,4) = 14
              idside(3,5) = 15
              idside(3,6) = 16

c
c     lado 4
c
              idside(4,1) = 4
              idside(4,2) = 1
              idside(4,3) = 17
              idside(4,4) = 18
              idside(4,5) = 19
              idside(4,6) = 20
c
           end if
c
c
           if(nen.eq.49) then
c
c     lado 1
c
              idside(1,1) = 1
              idside(1,2) = 2
              idside(1,3) = 5
              idside(1,4) = 6
              idside(1,5) = 7
              idside(1,6) = 8
              idside(1,7) = 9
c
c     lado 2
c
              idside(2,1) = 2
              idside(2,2) = 3
              idside(2,3) = 10
              idside(2,4) = 11
              idside(2,5) = 12
              idside(2,6) = 13
              idside(2,7) = 14
c
c     lado 3
c
              idside(3,1) = 3
              idside(3,2) = 4
              idside(3,3) = 15
              idside(3,4) = 16
              idside(3,5) = 17
              idside(3,6) = 18
              idside(3,7) = 19

c
c     lado 4
c
              idside(4,1) = 4
              idside(4,2) = 1
              idside(4,3) = 20
              idside(4,4) = 21
              idside(4,5) = 22
              idside(4,6) = 23
              idside(4,7) = 24
c
           end if
c
c
           if(nen.eq.64) then
c
c     lado 1
c
              idside(1,1) = 1
              idside(1,2) = 2
              idside(1,3) = 5
              idside(1,4) = 6
              idside(1,5) = 7
              idside(1,6) = 8
              idside(1,7) = 9
              idside(1,8) = 10
c
c     lado 2
c
              idside(2,1) = 2
              idside(2,2) = 3
              idside(2,3) = 11
              idside(2,4) = 12
              idside(2,5) = 13
              idside(2,6) = 14
              idside(2,7) = 15
              idside(2,8) = 16
c
c     lado 3
c
              idside(3,1) = 3
              idside(3,2) = 4
              idside(3,3) = 17
              idside(3,4) = 18
              idside(3,5) = 19
              idside(3,6) = 20
              idside(3,7) = 21
              idside(3,8) = 22

c
c     lado 4
c
              idside(4,1) = 4
              idside(4,2) = 1
              idside(4,3) = 23
              idside(4,4) = 24
              idside(4,5) = 25
              idside(4,6) = 26
              idside(4,7) = 27
              idside(4,8) = 28
c
           end if
c
c
          if(nen.eq.3) then
             idside(1,1) = 1
             idside(1,2) = 2
             
             idside(2,1) = 2
             idside(2,2) = 3
             
             idside(3,1) = 3
             idside(3,2) = 1
          end if
c
          if(nen.eq.6) then
             idside(1,1) = 1
             idside(1,2) = 2
             idside(1,3) = 4
            
             idside(2,1) = 2
             idside(2,2) = 3
             idside(2,3) = 5
            
             idside(3,1) = 3
             idside(3,2) = 1
             idside(3,3) = 6
          end if
c
          if(nen.eq.10) then
c         
             idside(1,1) = 1
             idside(1,2) = 2
             idside(1,3) = 4
             idside(1,4) = 5
c          
             idside(2,1) = 2
             idside(2,2) = 3
             idside(2,3) = 6
             idside(2,4) = 7
c             
             idside(3,1) = 3
             idside(3,2) = 1
             idside(3,3) = 8
             idside(3,4) = 9
c
c
          end if
c
       return 
c                                                                       
      end                                                               
c**** new **********************************************************************
      subroutine flux2(ien   ,x     ,xl    ,
     &                 d     ,dl    ,mat   ,
     &                 det   ,shl   ,shg   ,
     &                 w     ,c     ,
     &                 alhs  ,dlhs  ,
     &                 brhs  ,idiag ,lm    ,
     &                 grav  ,ipar  ,lado  ,
     &                 eleffd,elresd,shln  ,
     &                 shgn  ,wn    ,detn  ,
     &                 detb  ,shlb  ,shgb  ,
     &                 detpn ,shlpn ,shgpn ,
     &                 idside,xls   ,idlsd ,
     &                 ddis  ,detp  ,
     &                 shlp  ,shgp  , !
     &                 elmbb ,
     &                 elmcb ,
     &                 elmbc ,
     &                 elfbb ,elfcb ,
     &                 elmdb ,
c     
     &                 ideg  ,
c     
     &                 shlpsd,
     &                 shgpsd,
c     
     &                 shlsd,
     &                 shgsd,
c
c
     &                 aelm,
     &                 felm,
     &                 xvla,
     &                 xvlb,
c 
     &                 dproj,    
c
     &                 iedge,
     &                 idc  ,
c     
     &                 numel ,neesq ,nen   ,nsd   ,
     &                 nesd  ,nint  ,neg   ,nrowsh,
     &                 ned   ,nee   ,numnp ,ndof  ,
     &                 ncon  ,nencon,necon ,neep  ,
     &                 nints ,nnods ,nenlad,npars ,
     &                 nside ,nenp  ,nedge ,nodsp ,
     &                 index ,nmultpd,nmultpc)     
c
c
c.... program to calculate stifness matrix and force array for the
c        plane elasticity element and
c        assemble into the global left-hand-side matrix
c        and right-hand side vector
c
      implicit real*8 (a-h,o-z)
c                                                                       
c.... remove above card for single-precision operation               
c                                                                       
      logical diag,quad,zerodl
      dimension elmbb(neep,*),elmcb(neep,*),elmbc(nee,*)
      dimension elmdb(nee,*)
      dimension elfbb(*),elfcb(*)
      dimension ien(nen,*),x(nsd,*),xl(nesd,*),d(ndof,*),dl(ned,*),
     &          mat(*),det(*),shl(nrowsh,nen,*),shg(nrowsh,nen,*),
     &          w(*),c(8,*),alhs(*),dlhs(*),
     &          brhs(*),idiag(*),lm(ned,nodsp,*),
     &          grav(*),ipar(nodsp,*),lado(nside,*),
     $          eleffd(nee,*),elresd(*)
      dimension shln(3,nnods,*),wn(*),detn(*),shgn(3,nnods,*),
     &          detb(*),shlb(3,nenlad,*),shgb(3,nenlad,*),
     &          detpn(*),shlpn(3,npars,*),shgpn(3,npars,*),
     &          idside(nside,*),xls(nesd,*),idlsd(*),
     &          ddis(ned,nenp,*) 
	dimension dls(12),detp(*),shlp(3,nenp,*),shgp(3,nenp,*)
c
      dimension ideg(ndof,*),iedge(npars,*),idc(ndof,*)
c
      dimension dproj(ndof,*)
c
      dimension shlpsd(nrowsh,nenp,*),
     &          shgpsd(nrowsh,nenp,*)
c
      dimension shlsd(nrowsh,nen,*),
     &          shgsd(nrowsh,nen,*)
c
      dimension aelm(neep,neep,*),felm(neep,*),
     &          xvla(nints,*),xvlb(nints,*)
c
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
      common /iounit/ iin,iout,iecho,ioupp,itest1,itest2  
c      write(*,*) "Aqui 1.1"
c
c      consistent matrix
c      
      diag = .false.
      pi=4.d00*datan(1.d00)
      gf1=grav(1)
      gf2=grav(2)
      gf3=grav(3)      
c
      do 500 nel=1,numel
c
c      clear stiffness matrix and force array
c
      call clear(eleffd,nee*nee)   ! (lambda,mu)
      call clear(elmbb,neep*neep) ! A(U,V)
      call clear(elmcb,neep*nee) 
      call clear(elmbc,nee*neep)
c
      call clear(elresd,nee)
      call clear(elfbb,neep)
c
c      localize coordinates and Dirichlet b.c.
c
      call local(ien(1,nel),x,xl,nen,nsd,nesd)
      call local(ipar(1,nel),d,dl,nodsp,ndof,ned)
c
      m = mat(nel)
      quad = .true.
      if (nen.eq.4.and.ien(3,nel).eq.ien(4,nel)) quad = .false.
c
c
      call shgqs(xl,detp,shlp,shgp,nint,nel,neg,quad,nenp,shl,nen)
c
      nintb=nside*nints
      call shgtqsd(ien(1,nel),xl,shlpsd,shgpsd,
     &     nside,nintb,nints,nel,neg,nenp,shlsd,nen)
c
c....... form stiffness matrix
c
c      calculo de h - caracteristico para cada elemento
c
      if(nen.eq.3.or.nen.eq.6) then
c
      h2=xl(1,2)*xl(2,3)+xl(1,1)*xl(2,2)+xl(1,3)*xl(2,1)
     &  -xl(1,1)*xl(2,3)-xl(1,2)*xl(2,1)-xl(1,3)*xl(2,2)
          h2=h2*pt5
      h=dsqrt(h2)
c
      else
c
      h2=0
      do 100 l=1,2
      h2=h2+(xl(1,l)-xl(1,l+2))**2+(xl(2,l)-xl(2,l+2))**2
 100  continue
      h=dsqrt(h2)/2.d00
      h2=h*h
c
      end if 
c
c      set up material properties
c
 	delta1=c(1,m)
 	delta2=c(2,m)
      betah = delta1*h**delta2
	gama= c(3,m)
	teta=c(5,m)
	gamas=gama
	if(index.eq.3) gamas=-gama
	fi = pi/6.d00
c
c
c.....loop on integration points
c
      do 400 l=1,nint
      c1=detp(l)*w(l)
c      
c.... f=gf2*cos(¹x)*cos(¹y) + gf1*sin(¹x)*sin(¹y)
c
c
      xx=0.d00
      yy=0.d00
c
      do 2323 i=1,nen
      xx=xx+shl(3,i,l)*xl(1,i)
      yy=yy+shl(3,i,l)*xl(2,i)
2323  continue
c
      pix=pi*xx
      piy=pi*yy
      sx=dsin(pix)
      sy=dsin(piy)
      cx=dcos(pix)
      cy=dcos(piy)
c
      pi2=pi*pi
      co=1.d00/(pi2*2.d00+gama)
c
      ax=xx-xx*xx
      by=yy-yy*yy
      pss=gf2*cx*cy+gf1*sx*sy
c
c      loop computer volume integrals
c
      do 3000 j=1,nenp
      djx=shgp(1,j,l)*c1
      djy=shgp(2,j,l)*c1
      djn=shgp(3,j,l)*c1
c     
c.... source terms      
c
      nbj=ned*(j-1)
      nbj1=nbj+1
c
c     (grad p , grad q) + gamas (p,q) = (f,q)
c
      elfbb(nbj1)=elfbb(nbj1) + djn*pss 
c
c.... element stiffness
c      
      do 3003 i=1,nenp
      nbi=ned*(i-1)
      nbi1=nbi+1
c
	dix=shgp(1,i,l)
	diy=shgp(2,i,l)
      din=shgp(3,i,l)
c
c    (grad p, grad q) + gamas (p,q) = (f,q)
c
      elmbb(nbi1,nbj1) = elmbb(nbi1,nbj1) 
     &                 + (diy*djy+dix*djx)
     &                 + gamas*din*djn
 3003 continue
c
 3000 continue
  400 continue
c
c     armazena o vetor de carga por elementos
c
       do i=1,neep
         felm(i,nel) = elfbb(i)
       end do
c
c......boundary terms - due to symmetrization
c      typical of discontinuous Galerkin methods 
c      with SIP 
c
      do 4000 ns=1,nside
c
c-----localiza os nos do lado ns
c
c
      ns1=idside(ns,1)
      ns2=idside(ns,2)
      nl1=ien(ns1,nel)
      nl2=ien(ns2,nel)
c      
      if(nl2.gt.nl1) then 
	  sign = 1.d00
        do  nn=1,npars
           idlsd(nn)=idside(ns,nn)
	  end do
	else
	  sign = -1.d00
	  idlsd(1) = idside(ns,2)
	  idlsd(2) = idside(ns,1)
	  id3 = npars-2
	  if(id3.gt.0) then
	    do il=3,npars
	      idlsd(il) = idside(ns,npars+3-il)
	    end do 
	  end if
      end if 
c
c
      do 2000 nn=1,nenlad
      nl=idlsd(nn)
      xls(1,nn)=xl(1,nl)
      xls(2,nn)=xl(2,nl)
 2000 continue
c
      call oneshgp(xls,detn,shlb,shln,shgn,
     &             nenlad,nnods,nints,nesd,ns,nel,neg)      
      call oneshgp(xls,detpn,shlb,shlpn,shgpn,
     &             nenlad,npars,nints,nesd,ns,nel,neg)      
c
c     Dirichlet boundary donditios and local projection
c
      ndgs = lado(ns,nel)
c
c   projecao local
c
      call drchbc(shgpn,shlb,detpn,wn,
     &            gf1,gf2,gf3,xls,dls,gama,pi,fi,
     &            nints,nenlad,npars)
c
c    projeção descontinua
c
	ngsd = npars*(ndgs-1)
c
      do i=1,npars
	  dproj(1,ngsd + i) = dls(i)
	end do
c
c     weak Dirichlet bc 
c
	if(ideg(1,ndgs).eq.1) then
      call drchbc(shgpn,shlb,detpn,wn,
     &            gf1,gf2,gf3,xls,dls,gama,pi,fi,
     &            nints,nenlad,npars)
c
c     multiplicador continuo
c
      do i=1,npars
        ngs = iedge(i,ndgs)
	  d(1,ngs) = dls(i)
c
	  nls = idlsd(i)
	  dl(1,nls ) = dls(i)
	end do
	end if
c
c.....compute boundary integral
c
      do 1000 ls=1,nints
c
      lb = (ns-1)*nints + ls
c
c     geometria
c
        x1 =0.d00
        x2 =0.d00
        dx1=0.d00
        dx2=0.d00
c
        do i=1,nenlad
          x1 =x1 +xls(1,i)*shlb(2,i,ls)
          x2 =x2 +xls(2,i)*shlb(2,i,ls)
          dx1=dx1+xls(1,i)*shlb(1,i,ls)
          dx2=dx2+xls(2,i)*shlb(1,i,ls)
        end do
           dxx=dsqrt(dx1*dx1+dx2*dx2)
           xn1= sign*dx2/dxx
           xn2=-sign*dx1/dxx
c
c      jump terms
c
       do 1110 j=1,nenp
          nbj = ned*(j-1)
          nbj1=nbj+1
          djn=shgpsd(3,j,lb)*detn(ls)*wn(ls)
          djx=shgpsd(1,j,lb)*detn(ls)*wn(ls)
          djy=shgpsd(2,j,lb)*detn(ls)*wn(ls)
          gjn=djx*xn1 + djy*xn2
c
       do i=1,nenp
          nbi = ned*(i-1)
          nbi1= nbi+1
          din=shgpsd(3,i,lb)
          dix=shgpsd(1,i,lb)
          diy=shgpsd(2,i,lb)
          gin=dix*xn1 + diy*xn2
          elmbb(nbi1,nbj1) = elmbb(nbi1,nbj1) + teta*gin*djn-din*gjn
       end do
c       
 1110 continue
c
c
        do 1100 j=1,npars
c
cc          ncj1 = (ns-1)*npars + j
c
          ncj1 = idlsd(j)
c
          djn=shgpn(2,j,ls)*detpn(ls)*wn(ls)
c
c      no source term
c
          do i=1,nenp
             nbi=ned*(i-1)
             nbi1=nbi+1
c
             din=shgpsd(3,i,lb)
             dix=shgpsd(1,i,lb)
             diy=shgpsd(2,i,lb)
             gin=dix*xn1+diy*xn2
c
             elmcb(nbi1,ncj1)=elmcb(nbi1,ncj1) - teta*gin*djn
             elmbc(ncj1,nbi1)=elmbc(ncj1,nbi1) + gin*djn
c             
          end do
c
 1100 continue
 1000 continue
 4000 continue
c
c
c
c......boundary terms - penalti
c
      do 8000 ns=1,nside
c
c-----localiza os nos do lado ns
c
c
      ns1=idside(ns,1)
      ns2=idside(ns,2)
      nl1=ien(ns1,nel)
      nl2=ien(ns2,nel)
c      
      if(nl2.gt.nl1) then 
	  sign = 1.d00
        do  nn=1,npars
           idlsd(nn)=idside(ns,nn)
	  end do
	else
	  sign = -1.d00
	  idlsd(1) = idside(ns,2)
	  idlsd(2) = idside(ns,1)
	  id3 = npars-2
	  if(id3.gt.0) then
	    do il=3,npars
	      idlsd(il) = idside(ns,npars+3-il)
	    end do 
	  end if
      end if 
c
c
      do 8200 nn=1,nenlad
      nl=idlsd(nn)
      xls(1,nn)=xl(1,nl)
      xls(2,nn)=xl(2,nl)
 8200 continue
c
      call oneshgp(xls,detn,shlb,shln,shgn,
     &             nenlad,nnods,nints,nesd,ns,nel,neg)      
      call oneshgp(xls,detpn,shlb,shlpn,shgpn,
     &             nenlad,npars,nints,nesd,ns,nel,neg)      
c
c     Dirichlet boundary donditiod
c
      ndgs = lado(ns,nel)
	if(ideg(1,ndgs).eq.1) then
      call drchbc(shgpn,shlb,detpn,wn,
     &            gf1,gf2,gf3,xls,dls,gama,pi,fi,
     &            nints,nenlad,npars)
c
c
c     multiplicador continuo
c
      do i=1,npars
        ngs = iedge(i,ndgs)
	  d(1,ngs) = dls(i)
c
	  nls = idlsd(i)
	  dl(1,nls ) = dls(i)
	end do
c
	end if
c
c.....compute boundary integral
c
      do 8010 ls=1,nints
c
      lb = (ns-1)*nints + ls
c
c     geometria
c
        x1 =0.d00
        x2 =0.d00
        dx1=0.d00
        dx2=0.d00
c
        do i=1,nenlad
          x1 =x1 +xls(1,i)*shlb(2,i,ls)
          x2 =x2 +xls(2,i)*shlb(2,i,ls)
          dx1=dx1+xls(1,i)*shlb(1,i,ls)
          dx2=dx2+xls(2,i)*shlb(1,i,ls)
        end do
           dxx=dsqrt(dx1*dx1+dx2*dx2)
           xn1= sign*dx2/dxx
           xn2=-sign*dx1/dxx
c
c      jump terms
c
       do 8110 j=1,nenp
          nbj = ned*(j-1)
          nbj1=nbj+1
          djn=shgpsd(3,j,lb)*detn(ls)*wn(ls)
          djx=shgpsd(1,j,lb)*detn(ls)*wn(ls)
          djy=shgpsd(2,j,lb)*detn(ls)*wn(ls)
          gjn=djx*xn1 + djy*xn2
c
       do i=1,nenp
          nbi = ned*(i-1)
          nbi1= nbi+1
          din=shgpsd(3,i,lb)
          dix=shgpsd(1,i,lb)
          diy=shgpsd(2,i,lb)
          gin=dix*xn1 + diy*xn2
          elmbb(nbi1,nbj1) = elmbb(nbi1,nbj1) + betah*din*djn
       end do
c       
 8110 continue
c
c
        do 8100 j=1,npars
c
          ncj1 = idlsd(j)
c
          djn=shgpn(2,j,ls)*detpn(ls)*wn(ls)
c
c      no source term
c
          do i=1,nenp
             nbi=ned*(i-1)
             nbi1=nbi+1
c
             din=shgpsd(3,i,lb)
             dix=shgpsd(1,i,lb)
             diy=shgpsd(2,i,lb)
             gin=dix*xn1+diy*xn2
c             
             elmcb(nbi1,ncj1)=elmcb(nbi1,ncj1) - betah*din*djn
             elmbc(ncj1,nbi1)=elmbc(ncj1,nbi1) - betah*din*djn
          end do
c
c
          do i=1,npars
            nci1 = idlsd(i)
            din = shgpn(2,i,ls) 
            eleffd(nci1,ncj1) = eleffd(nci1,ncj1) + betah*din*djn
          end do
c
 8100 continue
 8010 continue
 8000 continue
c
c   Condensação
c
       call condtdg(elmbb,elmcb,elmbc,elmdb,
     &            elfbb,eleffd,elresd,
     &            neep,nee)
c
c     Armazena a inversa da matriz do elemento
c
      do i=1,neep
      do j=1,neep
        aelm(i,j,nel) = elmbb(i,j)
      end do
      end do
c
      if(nel.le.2) then  
	 write(23,*) nel,nee 
       do i=1,nee
         write(23,91) (eleffd(i,j),j=1,nee),elresd(i)
       end do 
  91   format(30e15.5) 
c
       end if 
c
c
c      computation of Dirichlet b.c. contribution
c
       call ztest(dl,nee,zerodl)
c
c
      if(.not.zerodl)
     & call kdbc(eleffd,elresd,dl,nee)
c
c
c.... assemble element stifness matrix and force array into global
c        left-hand-side matrix and right-hand side vector
c
      call addnsl(alhs,dlhs,eleffd,idiag,lm(1,1,nel),nee,diag)
c
      call addrhs(brhs,elresd,lm(1,1,nel),nee)
c
c
c
  500 continue
c
      return
      end
c**** new ********************************************************************** 
       subroutine condtdg(elmbb,elmcb,elmbc,elmdb,
     &            elfbb,eleffd,elresd,
     &            neep,nee)
c
      implicit real*8 (a-h,o-z)
c                                                                       
      dimension elmbb(neep,*),elmcb(neep,*),elmdb(nee,*),
     &          eleffd(nee,*),elmbc(nee,*)
      dimension elfbb(*),elresd(*)
      common /iounit/ iin,iout,iecho,ioupp,itest1,itest2 
c
c
c   condensação do sitema
c
c
c   Bb Xb   + Cb Xc    = Fbb
c
c   Bc Xb  + Cc Xc   = Fc
c
c   com eliminação das inciggnitas  Xb
c
      call invmb(elmbb,neep,neep)            
c
c    \br(D) = Bc * B^{1}
c
      do i=1,nee
      do j=1,neep
        elmdb(i,j) = 0.d00
        do k=1,neep
          elmdb(i,j) = elmdb(i,j) + elmbc(i,k)*elmbb(k,j)
        end do
      end do
      end do
c
c  eleffd = eleffd - Bc B^{1}C
c
      do i=1,nee
      do j=1,nee
        do k=1,neep
          eleffd(i,j) = eleffd(i,j) - elmdb(i,k)*elmcb(k,j)
        end do
      end do
      end do
c
c  elresd = elresd - Bc B^{1}Fb
c
      do i=1,nee
        do k=1,neep
          elresd(i) = elresd(i) - elmdb(i,k)*elfbb(k)
        end do
      end do

      return
c
      end
c**** new **********************************************************************
      subroutine flux3(ien   ,x     ,xl    ,    
     &                 d     ,dl    ,mat   ,
     &                 det   ,shl   ,shg   ,
     &                 w     ,c     ,
     &                 grav  ,ipar  ,lado  ,
     &                 eleffd,elresd,shln  ,
     &                 shgn  ,wn    ,detn  ,
     &                 detb  ,shlb  ,shgb  ,
     &                 detpn ,shlpn ,shgpn ,
     &                 idside,xls   ,idlsd ,
     &                 ddis  ,detp  ,
     &                 shlp  ,shgp  ,
     &                 elmbb ,
     &                 elmcb ,
     &                 elfbb ,elfcb ,
     &                 elmdb ,
c     
     &                 ideg  ,
c     
c     
     &                 shlpsd,
     &                 shgpsd,
c     
     &                 shlsd,
     &                 shgsd,
c
     &                 aelm,
     &                 felm,
     &                 xvla,
     &                 xvlb,
c
     &                 iedge,
c     
     &                 numel ,neesq ,nen   ,
     &                 nsd   ,nesd  ,nint  ,
     &                 neg   ,nrowsh,ned   ,
     &                 nee   ,numnp ,ndof  ,
     &                 ncon  ,nencon,necon ,
     &                 neep  ,nints ,nnods ,
     &                 nenlad,npars ,nside ,
     &                 nenp  ,nodsp ,index ,
     &                 nedge)
c
c
c     
c
c
c.... program to calculate stifness matrix and force array for the
c        plane elasticity element and
c        assemble into the global left-hand-side matrix
c        and right-hand side vector
c
      implicit real*8 (a-h,o-z)
c                                                                       
c.... remove above card for single-precision operation               
c                                                                       
      logical diag,quad,zerodl
      dimension elmbb(neep,*),elmcb(neep,*)
      dimension elmdb(nee,*)
      dimension elfbb(*),elfcb(*)
      dimension ideg(ndof,*),iedge(npars,*)
      dimension ien(nen,*),
     &          x(nsd,*),xl(nesd,*),d(ndof,*),dl(ned,*),
     &          mat(*),det(*),shl(nrowsh,nen,*),shg(nrowsh,nen,*),
     &          w(*),c(8,*),
     &          grav(*),ipar(nodsp,*),
     &          lado(nside,*),
     &          eleffd(nee,*),elresd(*)
      dimension shln(3,nnods,*),wn(*),detn(*),shgn(3,nnods,*),
     &          detb(*),shlb(3,nenlad,*),shgb(3,nenlad,*),
     &          detpn(*),shlpn(3,npars,*),shgpn(3,npars,*),
     &          idside(nside,*),xls(nesd,*),idlsd(*),
     &          ddis(ned,nenp,*) 
	dimension dls(12),detp(*),shlp(3,nenp,*),shgp(3,nenp,*)
c
      dimension shlpsd(nrowsh,nenp,*),
     &          shgpsd(nrowsh,nenp,*)
c
      dimension shlsd(nrowsh,nen,*),
     &          shgsd(nrowsh,nen,*)
c
      dimension aelm(neep,neep,*),felm(neep,*),
     &          xvla(nints,*),xvlb(nints,*),fab(200)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
      common /iounit/ iin,iout,iecho,ioupp,itest1,itest2 
c      write(*,*) "Aqui 1.2"
c

c      consistent matrix
c
      diag = .false.
      pi=4.d00*datan(1.d00)
      gf1=grav(1)
      gf2=grav(2)
      gf3=grav(3)      
c
c
      do 500 nel=1,numel
c
c      clear stiffness matrix and force array
c
c
      call clear(elfbb,neep)
c
c      localize coordinates and Dirichlet b.c.
c
      call local(ien(1,nel),x,xl,nen,nsd,nesd)
      call local(ipar(1,nel),d,dl,nodsp,ndof,ned)
c      
c
      m = mat(nel)
      quad = .true.
      if (nen.eq.4.and.ien(3,nel).eq.ien(4,nel)) quad = .false.
c
      call shgqs(xl,detp,shlp,shgp,nint,nel,neg,quad,nenp,shl,nen)
c
      nintb=nside*nints
      call shgtqsd(ien(1,nel),xl,shlpsd,shgpsd,
     &     nside,nintb,nints,nel,neg,nenp,shlsd,nen)
c
c....... form stiffness matrix
c
c      calculo de h - caracteristico para cada elemento
c
      if(nen.eq.3.or.nen.eq.6) then
c
      h2=xl(1,2)*xl(2,3)+xl(1,1)*xl(2,2)+xl(1,3)*xl(2,1)
     &  -xl(1,1)*xl(2,3)-xl(1,2)*xl(2,1)-xl(1,3)*xl(2,2)
          h2=h2*pt5
      h=dsqrt(h2)
c
      else
c
      h2=0
      do 100 l=1,2
      h2=h2+(xl(1,l)-xl(1,l+2))**2+(xl(2,l)-xl(2,l+2))**2
 100  continue
      h=dsqrt(h2)/2.d00
      h2=h*h
c
      end if
c
c      set up material properties
c
 	delta1=c(1,m)
 	delta2=c(2,m)
      betah = delta1*h**delta2
	gama= c(3,m)
	teta=c(5,m)
	gamas=gama
	if(index.eq.3) gamas=-gama
	fi = pi/6.d00
c
c     recupera o vetor de carga por elementos
c
       do i=1,neep
         elfbb(i) = felm(i,nel)
       end do
c
c......boundary terms - symmetrization
c
      do 4000 ns=1,nside
c
c-----localiza os nos do lado ns
c
c
      ns1=idside(ns,1)
      ns2=idside(ns,2)
      nl1=ien(ns1,nel)
      nl2=ien(ns2,nel)
c      
      if(nl2.gt.nl1) then 
	  sign = 1.d00
        do  nn=1,npars
           idlsd(nn)=idside(ns,nn)
	  end do
	else
	  sign = -1.d00
	  idlsd(1) = idside(ns,2)
	  idlsd(2) = idside(ns,1)
	  id3 = npars-2
	  if(id3.gt.0) then
	    do il=3,npars
	      idlsd(il) = idside(ns,npars+3-il)
	    end do 
	  end if
      end if 
c
c
c-----localiza os parametros do lado ns
c
      do nn=1,npars
	  nld = idlsd(nn)
	  dls(nn) = dl(1,nld)
      end do
c
      do 2000 nn=1,nenlad
      nl=idlsd(nn)
      xls(1,nn)=xl(1,nl)
      xls(2,nn)=xl(2,nl)
 2000 continue
c
c
      call oneshgp(xls,detn,shlb,shln,shgn,
     &             nenlad,nnods,nints,nesd,ns,nel,neg)      
      call oneshgp(xls,detpn,shlb,shlpn,shgpn,
     &             nenlad,npars,nints,nesd,ns,nel,neg)      
c
c.....compute boundary integral
c
      do 1000 ls=1,nints
       lb=(ns-1)*nints+ls
c
c
c    valores dos parametros do multiplicador
c
	  dhs=0.d00
        do i=1,npars
	    dhs = dhs + dls(i)*shlpn(2,i,ls)
        end do
c
        nla = lado(ns,nel)
        xvla(ls,nla) = dhs
c
c     geometria
c
        x1 =0.d00
        x2 =0.d00
        dx1=0.d00
        dx2=0.d00
c
        do i=1,nenlad
          x1 =x1 +xls(1,i)*shlb(2,i,ls)
          x2 =x2 +xls(2,i)*shlb(2,i,ls)
          dx1=dx1+xls(1,i)*shlb(1,i,ls)
          dx2=dx2+xls(2,i)*shlb(1,i,ls)
        end do
           dxx=dsqrt(dx1*dx1+dx2*dx2)
           xn1= sign*dx2/dxx
           xn2=-sign*dx1/dxx
c 
c    valores exatos dos multiplicadores
c
      pix=pi*x1
      piy=pi*x2
      sx=dsin(pix)
      sy=dsin(piy)
      cx=dcos(pix)
      cy=dcos(piy)
c
      pi2=pi*pi
      co=1.d00/(pi2*2.d00+gamas)
c
      sigma = dsqrt(gama)
      xc = x1*sigma*dcos(fi)
      yc = x2*sigma*dsin(fi)
c
c    valor exato do multiplicador
c
      dhse=gf2*co*cx*cy + gf1*co*sx*sy 
       
cc       write(47,*) ns,dhs,dhse
c
cc      dhs = dhse
c
       do 1110 j=1,nenp
          nbj = ned*(j-1)
          nbj1=nbj+1
          djn=shgpsd(3,j,lb)*detn(ls)*wn(ls)
          djx=shgpsd(1,j,lb)*detn(ls)*wn(ls)
          djy=shgpsd(2,j,lb)*detn(ls)*wn(ls)
          gjn=djx*xn1 + djy*xn2
c
c
          elfbb(nbj1) = elfbb(nbj1) + teta*dhs*gjn
c       
 1110 continue
c
 1000 continue
 4000 continue
c
c......boundary terms - penalti
c
      do 8000 ns=1,nside
c
c-----localiza os nos do lado ns
c
c
      ns1=idside(ns,1)
      ns2=idside(ns,2)
      nl1=ien(ns1,nel)
      nl2=ien(ns2,nel)
c      
      if(nl2.gt.nl1) then 
	  sign = 1.d00
        do  nn=1,npars
           idlsd(nn)=idside(ns,nn)
	  end do
	else
	  sign = -1.d00
	  idlsd(1) = idside(ns,2)
	  idlsd(2) = idside(ns,1)
	  id3 = npars-2
	  if(id3.gt.0) then
	    do il=3,npars
	      idlsd(il) = idside(ns,npars+3-il)
	    end do 
	  end if
      end if 
c
c
      do 8200 nn=1,nenlad
      nl=idlsd(nn)
      xls(1,nn)=xl(1,nl)
      xls(2,nn)=xl(2,nl)
c
 8200 continue
c
c
c-----localiza os parametros do lado ns
c
      do nn=1,npars
	  nld = idlsd(nn)
	  dls(nn) = dl(1,nld)
      end do
c
      call oneshgp(xls,detn,shlb,shln,shgn,
     &             nenlad,nnods,nints,nesd,ns,nel,neg)      
      call oneshgp(xls,detpn,shlb,shlpn,shgpn,
     &             nenlad,npars,nints,nesd,ns,nel,neg)      
c
c.....compute boundary integral
c
      do 8010 ls=1,nints
       lb=(ns-1)*nints+ls
c
c
c    valores dos parametros do multiplicador
c
	  dhs=0.d00
        do i=1,npars
	    dhs = dhs + dls(i)*shlpn(2,i,ls)
        end do
c
c     geometria
c
        x1 =0.d00
        x2 =0.d00
        dx1=0.d00
        dx2=0.d00
c
        do i=1,nenlad
          x1 =x1 +xls(1,i)*shlb(2,i,ls)
          x2 =x2 +xls(2,i)*shlb(2,i,ls)
          dx1=dx1+xls(1,i)*shlb(1,i,ls)
          dx2=dx2+xls(2,i)*shlb(1,i,ls)
        end do
           dxx=dsqrt(dx1*dx1+dx2*dx2)
           xn1= sign*dx2/dxx
           xn2=-sign*dx1/dxx
c
c
c 
c    valores exatos dos multiplicadores
c
      pix=pi*x1
      piy=pi*x2
      sx=dsin(pix)
      sy=dsin(piy)
      cx=dcos(pix)
      cy=dcos(piy)
c
      pi2=pi*pi
      co=1.d00/(pi2*2.d00+gamas)
c
c
c    valor exato do multiplicador
c
      dhse=gf2*co*cx*cy + gf1*co*sx*sy 
c
cc      dhs = dhse
c
       do 8110 j=1,nenp
          nbj = ned*(j-1)
          nbj1=nbj+1
          djn=shgpsd(3,j,lb)*detn(ls)*wn(ls)
          djx=shgpsd(1,j,lb)*detn(ls)*wn(ls)
          djy=shgpsd(2,j,lb)*detn(ls)*wn(ls)
          gjn=djx*xn1 + djy*xn2
c
          elfbb(nbj1) = elfbb(nbj1) + betah*dhs*djn
c       
 8110 continue
c
 8010 continue
 8000 continue
c
c	Local TLDG-solution (potencial) 
c
c     recupera a inversa da matriz do elemento
c
      do i=1,neep
        fab(i) = elfbb(i)
      do j=1,neep
        elmbb(i,j) = aelm(i,j,nel)
      end do
      end do
c
c     calcula o vetor elfbb
c
      do i=1,neep
       elfbb(i) = 0.d00
      do j=1,neep
       elfbb(i) = elfbb(i) + elmbb(i,j)*fab(j)
      end do
      end do 
c          
c      valores nodais descontinuos
c
      do j=1,nenp
	  ddis(1,j,nel) = elfbb(j)
      end do
c
c
  500 continue
c
c
      return
      end
c**** new **********************************************************************
      subroutine flux30(ien   ,x     ,xl    ,    
     &                 d     ,dl    ,mat   ,
     &                 det   ,shl   ,shg   ,
     &                 w     ,c     ,
     &                 grav  ,ipar  ,lado  ,
     &                 eleffd,elresd,shln  ,
     &                 shgn  ,wn    ,detn  ,
     &                 detb  ,shlb  ,shgb  ,
     &                 detpn ,shlpn ,shgpn ,
     &                 idside,xls   ,idlsd ,
     &                 ddis  ,detp  ,
     &                 shlp  ,shgp  ,
     &                 elmbb ,
     &                 elmcb ,
     &                 elfbb ,elfcb ,
     &                 elmdb ,
c     
     &                 ideg  ,
c     
c     
     &                 shlpsd,
     &                 shgpsd,
c     
     &                 shlsd,
     &                 shgsd,
c
     &                 aelm,
     &                 felm,
     &                 xvla,
     &                 xvlb,
c     
     &                 numel ,neesq ,nen   ,
     &                 nsd   ,nesd  ,nint  ,
     &                 neg   ,nrowsh,ned   ,
     &                 nee   ,numnp ,ndof  ,
     &                 ncon  ,nencon,necon ,
     &                 neep  ,nints ,nnods ,
     &                 nenlad,npars ,nside ,
     &                 nenp  ,nodsp ,index ,
     &                 nedge)
c
c
c     
c
c
c.... program to calculate stifness matrix and force array for the
c        plane elasticity element and
c        assemble into the global left-hand-side matrix
c        and right-hand side vector
c
      implicit real*8 (a-h,o-z)
c                                                                       
c.... remove above card for single-precision operation               
c                                                                       
      logical diag,quad,zerodl
      dimension elmbb(neep,*),elmcb(neep,*)
      dimension elmdb(nee,*)
      dimension elfbb(*),elfcb(*)
      dimension ideg(ndof,*)
      dimension ien(nen,*),
     &          x(nsd,*),xl(nesd,*),d(ndof,*),dl(ned,*),
     &          mat(*),det(*),shl(nrowsh,nen,*),shg(nrowsh,nen,*),
     &          w(*),c(8,*),
     &          grav(*),ipar(nodsp,*),
     &          lado(nside,*),
     &          eleffd(nee,*),elresd(*)
      dimension shln(3,nnods,*),wn(*),detn(*),shgn(3,nnods,*),
     &          detb(*),shlb(3,nenlad,*),shgb(3,nenlad,*),
     &          detpn(*),shlpn(3,npars,*),shgpn(3,npars,*),
     &          idside(nside,*),xls(nesd,*),idlsd(*),
     &          ddis(ned,nenp,*) 
	dimension dls(12),detp(*),shlp(3,nenp,*),shgp(3,nenp,*)
c
      dimension shlpsd(nrowsh,nenp,*),
     &          shgpsd(nrowsh,nenp,*)
c
      dimension shlsd(nrowsh,nen,*),
     &          shgsd(nrowsh,nen,*)
c
      dimension aelm(neep,neep,*),felm(neep,*),
     &          xvla(nints,*),xvlb(nints,*),fab(200)
c
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
      common /iounit/ iin,iout,iecho,ioupp,itest1,itest2 
c
c      consistent matrix
c
      diag = .false.
      pi=4.d00*datan(1.d00)
      gf1=grav(1)
      gf2=grav(2)
      gf3=grav(3)      
c
c
      do 500 nel=1,numel
c
c      clear stiffness matrix and force array
c
c
      call clear(elfbb,neep)
c
c      localize coordinates and Dirichlet b.c.
c
      call local(ien(1,nel),x,xl,nen,nsd,nesd)
cc      call local(ipar(1,nel),d,dl,nodsp,ndof,ned)
c
      do ns=1,nside 
         ndgs = lado(ns,nel)
         nls = (ns-1)*npars
	   ngsd = npars*(ndgs-1)
          do i=1,npars
	     dl(1,nls+i) = d(1,ngsd + i)
	    end do
	end do
c
      m = mat(nel)
      quad = .true.
      if (nen.eq.4.and.ien(3,nel).eq.ien(4,nel)) quad = .false.
c
      call shgqs(xl,detp,shlp,shgp,nint,nel,neg,quad,nenp,shl,nen)
c
      nintb=nside*nints
      call shgtqsd(ien(1,nel),xl,shlpsd,shgpsd,
     &     nside,nintb,nints,nel,neg,nenp,shlsd,nen)
c
c....... form stiffness matrix
c
c      calculo de h - caracteristico para cada elemento
c
      if(nen.eq.3.or.nen.eq.6) then
c
      h2=xl(1,2)*xl(2,3)+xl(1,1)*xl(2,2)+xl(1,3)*xl(2,1)
     &  -xl(1,1)*xl(2,3)-xl(1,2)*xl(2,1)-xl(1,3)*xl(2,2)
          h2=h2*pt5
      h=dsqrt(h2)
c
      else
c
      h2=0
      do 100 l=1,2
      h2=h2+(xl(1,l)-xl(1,l+2))**2+(xl(2,l)-xl(2,l+2))**2
 100  continue
      h=dsqrt(h2)/2.d00
      h2=h*h
c
      end if
c
c      set up material properties
c
 	delta1=c(1,m)
 	delta2=c(2,m)
      betah = delta1*h**delta2
	gama= c(3,m)
	teta=c(5,m)
	gamas=gama
	if(index.eq.3) gamas=-gama
	fi = pi/6.d00
c
c     recupera o vetor de carga por elementos
c
       do i=1,neep
         elfbb(i) = felm(i,nel)
       end do
c
c......boundary terms - symmetrization
c
      do 4000 ns=1,nside
c
c-----localiza os parametros do lado ns
c
      do nn=1,npars
	  nld = (ns-1)*npars + nn
	  dls(nn) = dl(1,nld)
      end do
c
c-----localiza os nos do lado ns
c
c
      ns1=idside(ns,1)
      ns2=idside(ns,2)
      nl1=ien(ns1,nel)
      nl2=ien(ns2,nel)
c      
      if(nl2.gt.nl1) then 
	  sign = 1.d00
        do  nn=1,npars
           idlsd(nn)=idside(ns,nn)
	  end do
	else
	  sign = -1.d00
	  idlsd(1) = idside(ns,2)
	  idlsd(2) = idside(ns,1)
	  id3 = npars-2
	  if(id3.gt.0) then
	    do il=3,npars
	      idlsd(il) = idside(ns,npars+3-il)
	    end do 
	  end if
      end if 
c
c
      do 2000 nn=1,nenlad
      nl=idlsd(nn)
      xls(1,nn)=xl(1,nl)
      xls(2,nn)=xl(2,nl)
c
 2000 continue
c
c
      call oneshgp(xls,detn,shlb,shln,shgn,
     &             nenlad,nnods,nints,nesd,ns,nel,neg)      
      call oneshgp(xls,detpn,shlb,shlpn,shgpn,
     &             nenlad,npars,nints,nesd,ns,nel,neg)      
c
c.....compute boundary integral
c
      do 1000 ls=1,nints
       lb=(ns-1)*nints+ls
c
        nla = lado(ns,nel)
        xvla(ls,nla) = dhs
c
c     geometria
c
        x1 =0.d00
        x2 =0.d00
        dx1=0.d00
        dx2=0.d00
c
        do i=1,nenlad
          x1 =x1 +xls(1,i)*shlb(2,i,ls)
          x2 =x2 +xls(2,i)*shlb(2,i,ls)
          dx1=dx1+xls(1,i)*shlb(1,i,ls)
          dx2=dx2+xls(2,i)*shlb(1,i,ls)
        end do
           dxx=dsqrt(dx1*dx1+dx2*dx2)
           xn1= sign*dx2/dxx
           xn2=-sign*dx1/dxx
c 
c    valores exatos dos multiplicadores
c
      pix=pi*x1
      piy=pi*x2
      sx=dsin(pix)
      sy=dsin(piy)
      cx=dcos(pix)
      cy=dcos(piy)
c
      pi2=pi*pi
      co=1.d00/(pi2*2.d00+gamas)
c
      sigma = dsqrt(gama)
      xc = x1*sigma*dcos(fi)
      yc = x2*sigma*dsin(fi)
c
c    valor exato do multiplicador
c
      dhs=gf2*co*cx*cy + gf1*co*sx*sy 
       
cc       write(47,*) ns,dhs,dhse
c
cc      dhs = dhse
c
       do 1110 j=1,nenp
          nbj = ned*(j-1)
          nbj1=nbj+1
          djn=shgpsd(3,j,lb)*detn(ls)*wn(ls)
          djx=shgpsd(1,j,lb)*detn(ls)*wn(ls)
          djy=shgpsd(2,j,lb)*detn(ls)*wn(ls)
          gjn=djx*xn1 + djy*xn2
c
c
          elfbb(nbj1) = elfbb(nbj1) + teta*dhs*gjn
c       
 1110 continue
c
 1000 continue
 4000 continue
c
c......boundary terms - penalti
c
      do 8000 ns=1,nside
c
c-----localiza os parametros do lado ns
c
      do nn=1,npars
	  nld = (ns-1)*npars + nn
	  dls(nn) = dl(1,nld)
      end do
c
c-----localiza os nos do lado ns
c
c
      ns1=idside(ns,1)
      ns2=idside(ns,2)
      nl1=ien(ns1,nel)
      nl2=ien(ns2,nel)
c      
      if(nl2.gt.nl1) then 
	  sign = 1.d00
        do  nn=1,npars
           idlsd(nn)=idside(ns,nn)
	  end do
	else
	  sign = -1.d00
	  idlsd(1) = idside(ns,2)
	  idlsd(2) = idside(ns,1)
	  id3 = npars-2
	  if(id3.gt.0) then
	    do il=3,npars
	      idlsd(il) = idside(ns,npars+3-il)
	    end do 
	  end if
      end if 
c
c
      do 8200 nn=1,nenlad
      nl=idlsd(nn)
      xls(1,nn)=xl(1,nl)
      xls(2,nn)=xl(2,nl)
c
 8200 continue
c
c
      call oneshgp(xls,detn,shlb,shln,shgn,
     &             nenlad,nnods,nints,nesd,ns,nel,neg)      
      call oneshgp(xls,detpn,shlb,shlpn,shgpn,
     &             nenlad,npars,nints,nesd,ns,nel,neg)      
c
c.....compute boundary integral
c
      do 8010 ls=1,nints
       lb=(ns-1)*nints+ls
c
c     geometria
c
        x1 =0.d00
        x2 =0.d00
        dx1=0.d00
        dx2=0.d00
c
        do i=1,nenlad
          x1 =x1 +xls(1,i)*shlb(2,i,ls)
          x2 =x2 +xls(2,i)*shlb(2,i,ls)
          dx1=dx1+xls(1,i)*shlb(1,i,ls)
          dx2=dx2+xls(2,i)*shlb(1,i,ls)
        end do
           dxx=dsqrt(dx1*dx1+dx2*dx2)
           xn1= sign*dx2/dxx
           xn2=-sign*dx1/dxx
c
c
c 
c    valores exatos dos multiplicadores
c
      pix=pi*x1
      piy=pi*x2
      sx=dsin(pix)
      sy=dsin(piy)
      cx=dcos(pix)
      cy=dcos(piy)
c
      pi2=pi*pi
      co=1.d00/(pi2*2.d00+gamas)
c
c
c    valor exato do multiplicador
c
      dhs=gf2*co*cx*cy + gf1*co*sx*sy 
c
cc      dhs = dhse
c
       do 8110 j=1,nenp
          nbj = ned*(j-1)
          nbj1=nbj+1
          djn=shgpsd(3,j,lb)*detn(ls)*wn(ls)
          djx=shgpsd(1,j,lb)*detn(ls)*wn(ls)
          djy=shgpsd(2,j,lb)*detn(ls)*wn(ls)
          gjn=djx*xn1 + djy*xn2
c
          elfbb(nbj1) = elfbb(nbj1) + betah*dhs*djn
c       
 8110 continue
c
 8010 continue
 8000 continue
c
c	Local TLDG-solution (potencial) 
c
c     recupera a inversa da matriz do elemento
c
      do i=1,neep
        fab(i) = elfbb(i)
      do j=1,neep
        elmbb(i,j) = aelm(i,j,nel)
      end do
      end do
c
c     calcula o vetor elfbb
c
      do i=1,neep
       elfbb(i) = 0.d00
      do j=1,neep
       elfbb(i) = elfbb(i) + elmbb(i,j)*fab(j)
      end do
      end do 
c          
c      valores nodais descontinuos
c
      do j=1,nenp
	  ddis(1,j,nel) = elfbb(j)
      end do
c
      do j=1,nenp
	  ddis(1,j,nel) = elfbb(j)
      end do
c
  500 continue
c
c
      return
      end
c**** new **********************************************************************
      subroutine drchbc(shgpn,shlb,detpn,wn,
     &            gf1,gf2,gf3,xls,dls,gama,pi,fi,
     &            nints,nenlad,npars)                     
c                                                                        
c     Dirichlet boundary conditions                                             
c                                                                        
      implicit real*8(a-h,o-z)                                  
      dimension shlb(3,nenlad,*),shgpn(3,npars,*),xls(2,*)
	dimension detpn(*),wn(*),dls(*)
	dimension aa(10,10),bb(10)
c
      do i=1,npars
	 bb(i) = 0.d00
	do j=1,npars
	 aa(i,j) = 0.d00
      end do
	end do
c
      do 1000 ls=1,nints
c
c     geometria
c
        x1 =0.d00
        x2 =0.d00
c
        do i=1,nenlad
          x1 =x1 +xls(1,i)*shlb(2,i,ls)
          x2 =x2 +xls(2,i)*shlb(2,i,ls)
        end do
c 
      pix=pi*x1
      piy=pi*x2
      sx=dsin(pix)
      sy=dsin(piy)
      cx=dcos(pix)
      cy=dcos(piy)
c
      pi2=pi*pi
      co=1.d00/(pi2*2.d00+gama)
      sigma=dsqrt(gama)
      xc=x1*sigma*dcos(fi)
      yc=x2*sigma*dsin(fi)
c
c    valor exato do multiplicador
c
      dhse = gf2*co*cx*cy + gf1*co*sx*sy + gf3*dcos(xc + yc)
c
        do 1100 j=1,npars
c
c
         djn=shgpn(2,j,ls)*detpn(ls)*wn(ls)
c
c    source term
c
         bb(j) = bb(j) + dhse*djn
c
        do i=1,npars
c
          din=shgpn(2,i,ls)
c
c     L2-projection at the boundary (side)
c             
          aa(i,j) = aa(i,j) + din*djn
c
        end do
c
 1100 continue
 1000 continue
c
      if(npars.gt.1) call invmb(aa,10,npars)
      if(npars.eq.1) aa(1,1) = 1.d00/aa(1,1)
c
	do i=1,npars
	  dls(i) = 0.d00
	do j=1,npars
	  dls(i) = dls(i) + aa(i,j)*bb(j)
	end do
	end do
      return                                                             
      end                                                                
c**** new **********************************************************************
      subroutine solvetdg(elmbb,elfbb,neep)
c
      implicit real*8 (a-h,o-z)
c                                                                       
      dimension elmbb(neep,*),elfbb(*)
      dimension fab(200)
      common /iounit/ iin,iout,iecho,ioupp,itest1,itest2 
c
c   resolve o sistema
c
c   B Xb = Fb
c
      call invmb(elmbb,neep,neep)            
c
c    Aa = A^{-1} Fa
c 
      do j=1,neep
        fab(j)=0.d00
        do k=1,neep
          fab(j) = fab(j) + elmbb(j,k)*elfbb(k)
        end do
      end do
c
c  xb=fab
c
      do i=1,neep
      elfbb(i) = fab(i)
      end do      
c
      return
c
      end
c**** new **********************************************************************
      subroutine invmb(am,ndim,m)                      
c                                                                        
c     subrotina de inversao                                              
c                                                                        
      implicit real*8(a-h,o-z)                                  
      dimension ipi(200),ind(200,2),piv(200),dis(200,1),am(ndim,*)
      ncoln=0                                                            
      det=1.0                                                            
      do 20 j=1,m                                                        
   20 ipi(j)=0                                                           
      do 550 i=1,m                                                       
      amax=0.0                                                           
      do 105 j=1,m                                                       
      if(ipi(j)-1)60,105,60                                              
   60 do 100 k=1,m                                                       
      if(ipi(k)-1) 80,100,740                                            
   80 if(dabs(amax)-dabs(am(j,k)))85,100,100                             
   85 irow=j                                                             
      ico=k                                                              
      amax=am(j,k)                                                       
  100 continue                                                           
  105 continue                                                           
      ipi(ico)=ipi(ico)+1                                                
      if(irow-ico)140,260,140                                            
  140 det=-det                                                           
      do 200 l=1,m                                                       
      swap=am(irow,l)                                                    
      am(irow,l)=am(ico,l)                                               
  200 am(ico,l)=swap                                                     
      if(ncoln) 260,260,210                                              
  210 do 250 l=1,ncoln                                                   
      swap=dis(irow,l)                                                   
      dis(irow,l)=dis(ico,l)                                             
  250 dis(ico,l)=swap                                                    
  260 ind(i,1)=irow                                                      
      piv(i)=am(ico,ico)                                                 
      ind(i,2)=ico                                                       
      det=det*piv(i)                                                     
      am(ico,ico)=1.0                                                    
      do 350 l=1,m                                                       
  350 am(ico,l)=am(ico,l)/piv(i)                                         
      if(ncoln) 380,380,360                                              
  360 do 370 l=1,ncoln                                                   
  370 dis(ico,l)=dis(ico,l)/piv(i)                                       
  380 do 550 lz=1,m                                                      
      if(lz-ico)400,550,400                                              
  400 t=am(lz,ico)                                                       
      am(lz,ico)=0.0                                                     
      do 450 l=1,m                                                       
  450 am(lz,l)=am(lz,l)-am(ico,l)*t                                      
      if(ncoln)550,550,460                                               
  460 do 500 l=1,ncoln                                                   
  500 dis(lz,l)=dis(lz,l)-dis(ico,l)*t                                   
  550 continue                                                           
      do 710 i=1,m                                                       
      l=m+1-i                                                            
      if(ind(l,1)-ind(l,2))630,710,630                                   
  630 jrow=ind(l,1)                                                      
      jco=ind(l,2)                                                       
      do 705 k=1,m                                                       
      swap=am(k,jrow)                                                    
      am(k,jrow)=am(k,jco)                                               
      am(k,jco)=swap                                                     
  705 continue                                                           
  710 continue                                                           
  740 continue                                                           
      return                                                             
      end                                                                
c**** new **********************************************************************
      subroutine fluxmx(c,numat)
c
c.... program to read, write and store properties
c     for plane stres mixed elements
c
c
      implicit real*8 (a-h,o-z)
c
c.... remove above card for single precision operation
c
      dimension c(8,*)
      common /iounit/ iin,iout,iecho,ioupp,itest1,itest2 
c
      do 100 n=1,numat
      if (mod(n,50).eq.1) write(iecho,1000) numat
c
      read(iin,2000) m,del1,del2,del3,del4,del5
      c(1,m)=del1
      c(2,m)=del2
      c(3,m)=del3
      c(4,m)=del4
      c(5,m)=del5
      write(iecho,3000) m,del1,del2,del3,del4,del5
c
     
  100 continue
c
      return
c
 1000 format(///,
     &' m a t e r i a l   s e t   d a t a       '   //5x,
     &' number of material sets . . . . . (numat ) = ',i10//,
     & 7x,'set',7x,'del1',7x,'del2',
     & 7x,'del3',7x,'del4',7x,'del5',/)
 2000 format(i10,8f10.0)
 3000 format(i10,2x,8(1x,1pe10.3))
      end
c**** new **********************************************************************
c-----------------------------------------------------
      subroutine uexafx(x,y,ue,duex,duey,gama,fi,index)
c-----------------------------------------------------
      implicit real*8(a-h,o-z)
      dimension ue(*),duex(*),duey(*)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
c
c
c
c
       go to (100,200,300) index

c
c        PROBLEM 1 - T=sin¹x*sin¹y 
c
 100   continue
      pi=4.d00*datan(1.d00)
      pi2=pi*pi
      px=pi*x
      py=pi*y
      sx=dsin(px)
      sy=dsin(py)
      cx=dcos(px)
      cy=dcos(py)
      co=1.d00/(pi2*2.d00 + gama)
c
      ue(1) = -pi*co*cx*sy
      ue(2) = -pi*co*sx*cy
      ue(3) = co*sx*sy
c      
      duex(1) =  pi2*co*sx*sy
      duex(2) = -pi2*co*cx*cy
      duex(3) =  pi*co*cx*sy
c
      duey(1) = -pi2*co*cx*cy
      duey(2) =  pi2*co*sx*sy
      duey(3) =  pi*co*sx*cy
      go to 1001
c
c     2  PROBLEM 2 - T=cos¹xcos¹y
c
 200   continue
      pi=4.d00*datan(1.d00)
      pi2=pi*pi
      px=pi*x
      py=pi*y
      sx=dsin(px)
      sy=dsin(py)
      cx=dcos(px)
      cy=dcos(py)
      co=1.d00/(pi2*2.d00 + gama)
c      
      ue(1) = co*pi*sx*cy
      ue(2) = co*pi*cx*sy
      ue(3) = co*cx*cy
c      
      duex(1)= co*pi2*cx*cy
      duex(2)=-co*pi2*sx*sy
      duex(3)=-co*pi*sx*cy
c      
      duey(1)=-co*pi2*sx*sy
      duey(2)= co*pi2*cx*cy
      duey(3)=-co*pi*cx*sy
      go to 1001
c     
c     Polinomial
c
300   continue
c
      sigma = dsqrt(gama)
      xc = x*sigma*dcos(fi)
      yc = y*sigma*dsin(fi)
c         
      ue(1) = sigma*dcos(fi)*dsin(xc + yc)
      ue(2) = sigma*dsin(fi)*dsin(xc + yc)
      ue(3) = dcos(xc + yc) 
c      
      duex(1)= sigma**2*dcos(fi)*dcos(fi)*dcos(xc + yc)         
      duex(2)= sigma**2*dsin(fi)*dcos(fi)*dcos(xc + yc)
      duex(3)=-sigma*dcos(fi)*dsin(xc + yc) 
c      
      duey(1)= sigma**2*dsin(fi)*dcos(fi)*dcos(xc + yc)
      duey(2)= sigma**2*dsin(fi)*dsin(fi)*dcos(xc + yc) 
      duey(3)=-sigma*dsin(fi)*dsin(xc + yc)
c      
1001  continue
      return
      end
c-----------------------------------------------------
      subroutine elemdlp(xl,dlp,gama,fi,ned,nenp,index)
c-----------------------------------------------------
      implicit real*8(a-h,o-z)
      dimension xl(2,*),dlp(ned,*)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
c
      pi=4.d00*datan(1.d00)
      pi2=pi*pi
c
      do 1001 n=1,nenp 
       x = xl(1,n)
       y = xl(2,n)
       go to (100,200,300) index
c
c        PROBLEM 1 - T=sin¹x*sin¹y 
c
 100   continue
      px=pi*x
      py=pi*y
      sx=dsin(px)
      sy=dsin(py)
      cx=dcos(px)
      cy=dcos(py)
      co=1.d00/(pi2*2.d00 + gama)
c
      dlp(1,n) = co*sx*sy
c
      go to 1001
c
c     2  PROBLEM 2 - T=cos¹xcos¹y
c
 200   continue
      px=pi*x
      py=pi*y
      sx=dsin(px)
      sy=dsin(py)
      cx=dcos(px)
      cy=dcos(py)
      co=1.d00/(pi2*2.d00 + gama)
c      
      dlp(1,n) = co*cx*cy
c      
      go to 1001
c     
c     Plane wave
c
300   continue
c
      sigma = dsqrt(gama)
      xc = x*sigma*dcos(fi)
      yc = y*sigma*dsin(fi)
c         
      dlp(1,n) = dcos(xc + yc) 
c      
1001  continue
      return
      end
c**** new ********************************************************************** 
      subroutine shlq(shl,w,nint,nen)
c
c.... program to calculate integration-rule weights, shape functions
c        and local derivatives for a four-node quadrilateral element
c
c               s,t = local element coordinates ("xi", "eta", resp.)       
c        shl(1,i,l) = local ("xi") derivative of shape function
c        shl(2,i,l) = local ("eta") derivative of shape function
c        shl(3,i,l) = local  shape function
c              w(l) = integration-rule weight
c                 i = local node number
c                 l = integration point number
c              nint = number of integration points, eq. 1 or 4
c
      implicit real*8 (a-h,o-z)
c
c.... remove above card for single precision operation
c
      dimension wone(8),raone(8)
      dimension shl(3,nen,*),w(*),ra(64),sa(64)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
      data  five9/0.5555555555555555d0/,eight9/0.8888888888888888d0/
      data r1/0.d00/,w1/2.d00/,
     &     r2/0.577350269189626d00/,w2/1.d00/,
     &     r3a/0.774596669241483d00/,w3a/0.555555555555556d00/,
     &     r3b/0.d00/,w3b/0.888888888888889d00/,
     &     r4a/0.861136311594053d00/,w4a/0.347854845137454d00/,
     &     r4b/0.339981043584856d00/,w4b/0.652145154862546d00/
c
      if (nint.eq.1) then
         wone(1)  = two
         raone(1) = zero
	nintx=1
	ninty=1
      endif
c
c
      if (nint.eq.4) then
         wone(1) = one
         wone(2) = one
         raone(1)=-.577350269189626
         raone(2)= .577350269189625
	nintx=2
	ninty=2
      endif
c
c
      if (nint.eq.9) then
         wone(1) = five9
         wone(2) = five9
         wone(3) = eight9
         raone(1)=-.774596669241483
         raone(2)= .774596669241483
         raone(3)= zero
	nintx=3
	ninty=3
      endif
c
      if (nint.eq.16) then
         wone(1) = .347854845137454
         wone(2) = .347854845137454
         wone(3) = .652145154862546
         wone(4) = .652145154862546
         raone(1)=-.861136311594053
         raone(2)= .861136311594053
         raone(3)=-.339981043584856
         raone(4)= .339981043584856
	nintx=4
	ninty=4
      endif
c
c
       if(nint.eq.25) then
        wone(1) = .236926885056189
        wone(2) = .236926885056189
        wone(3) = .478628670499366
        wone(4) = .478628670499366
        wone(5) = .568888888888888
        raone(1)=-.906179845938664
        raone(2)= .906179845938664
        raone(3)=-.538469310105683
        raone(4)= .538469310105683
        raone(5)= zero
	nintx=5
	ninty=5
       endif
c
       if(nint.eq.36) then
         wone(1) = .171324492397170
         wone(2) = .171324492397170
         wone(3) = .360761573048139
         wone(4) = .360761573048139
         wone(5) = .467913934572691
         wone(6) = .467913934572691
         raone(1)=-.932469514203152
         raone(2)= .932469514203152
         raone(3)=-.661209386466265
         raone(4)= .661209386466365
         raone(5)=-.238619186083197
         raone(6)= .238619186083197
	nintx=6
	ninty=6
        endif
c
c
       if(nint.eq.49) then
         wone(1) = .129484966168870
         wone(2) = .129484966168870
         wone(3) = .279705391489277 
         wone(4) = .279705391489277
         wone(5) = .381830050505119
         wone(6) = .381830050505119
         wone(7) = .417959183673469
         raone(1)=-.949107912342759
         raone(2)= .949107912342759
         raone(3)=-.741531185599394
         raone(4)= .741531185599394
         raone(5)=-.405845151377397
         raone(6)= .405845151377397
         raone(7)= zero
	nintx=7
	ninty=7
        endif
c
c
       if(nint.eq.64) then
         wone(1) = .101228536290376
         wone(2) = .101228536290376
         wone(3) = .222381034453374
         wone(4) = .222381034453374
         wone(5) = .313706645877887
         wone(6) = .313706645877887
         wone(7) = .362683783378362
         wone(8) = .362683783378362
         raone(1)=-.960289856497536
         raone(2)= .960289856497536
         raone(3)=-.796666477413627
         raone(4)= .796666477413627
         raone(5)=-.525532409916329
         raone(6)= .525532409916329
         raone(7)=-.183434642495650
         raone(8)= .183434642495650
	nintx=8
	ninty=8
        endif
c
c
         l=0
	   do ly=1,ninty
	   do lx=1,nintx
	   l = l+1
	      w(l) = wone(lx)*wone(ly)
	      ra(l) = raone(lx)
	      sa(l) = raone(ly)
         end do
	   end do
c
      do 200 l=1,nint
c
            r=ra(l)
            s=sa(l)
c
	if(nen.eq.4) then
            f1 = pt5*(one-r)
            f2 = pt5*(one+r)
	      fx1=-pt5
	      fx2= pt5
            g1 = pt5*(one-s)
            g2 = pt5*(one+s)
	      gx1=-pt5
	      gx2= pt5
            shl(1,1,l)=fx1*g1
            shl(2,1,l)=f1*gx1
            shl(3,1,l)=f1*g1
            shl(1,2,l)=fx2*g1
            shl(2,2,l)=f2*gx1
            shl(3,2,l)=f2*g1
            shl(1,3,l)=fx2*g2
            shl(2,3,l)=f2*gx2
            shl(3,3,l)=f2*g2
            shl(1,4,l)=fx1*g2
            shl(2,4,l)=f1*gx2
            shl(3,4,l)=f1*g2
	   end if
c
c
         if(nen.eq.9) then
            f1 = -pt5*(one-r)*r
            f2 =  pt5*(one+r)*r
	      f3 = (one+r)*(one-r)
c
c	      
	      f1x = pt5*r - pt5*(one-r)
	      f2x = pt5*r + pt5*(one+r)
	      f3x =-two*r

            g1 =-pt5*(one-s)*s
            g2 = pt5*(one+s)*s
	      g3 = (one+s)*(one-s)
c
	      g1x = pt5*s - pt5*(one-s)
	      g2x = pt5*s + pt5*(one+s)
	      g3x =-two*s
c
c
            shl(3,1,l)=f1*g1
            shl(3,2,l)=f2*g1
            shl(3,3,l)=f2*g2
            shl(3,4,l)=f1*g2
c
                  shl(3,5,l)=f3*g1
                  shl(3,6,l)=f2*g3
                  shl(3,7,l)=f3*g2
                  shl(3,8,l)=f1*g3
                  shl(3,9,l)=f3*g3
c
            shl(1,1,l)=f1x*g1
            shl(1,2,l)=f2x*g1
            shl(1,3,l)=f2x*g2
            shl(1,4,l)=f1x*g2
c
                  shl(1,5,l)=f3x*g1
                  shl(1,6,l)=f2x*g3
                  shl(1,7,l)=f3x*g2
                  shl(1,8,l)=f1x*g3
                  shl(1,9,l)=f3x*g3
c
            shl(2,1,l)=f1*g1x
            shl(2,2,l)=f2*g1x
            shl(2,3,l)=f2*g2x
            shl(2,4,l)=f1*g2x
c
                  shl(2,5,l)=f3*g1x
                  shl(2,6,l)=f2*g3x
                  shl(2,7,l)=f3*g2x
                  shl(2,8,l)=f1*g3x
                  shl(2,9,l)=f3*g3x
c
            end if
c
c
            if (nen.eq.16) then
                  onemrsq=one-r*r
                  onemssq=one-s*s
                  onep3r=one+three*r
                  onem3r=one-three*r
                  onep3s=one+three*s
                  onem3s=one-three*s
		  f1=-1.d00/16.d00*(9.d00*r*r-1.d00)*(r-1.d00)
		  f2=9.d00/16.d00*(1.d00-r*r)*onem3r
		  f3=9.d00/16.d00*(1.d00-r*r)*onep3r
		  f4=1.d00/16.d00*(9.d00*r*r-1.d00)*(r+1.d00)
c
		  f1x=-1.d00/16.d00*(18.d00*r)*(r-1.d00)
     &               -1.d00/16.d00*(9.d00*r*r-1.d00)
		  f2x=9.d00/16.d00*(-2.d00*r)*onem3r
     &               -9.d00/16.d00*(1.d00-r*r)*3.d00
		  f3x=9.d00/16.d00*(-2.d00*r)*onep3r
     &               -9.d00/16.d00*(r*r-1.d00)*3.d00
      	  f4x=1.d00/16.d00*(18.d00*r)*(r+1.d00)
     &		       +1.d00/16.d00*(9.d00*r*r-1.d00)
c
            g1=-1.d00/16.d00*(9.d00*s*s-1.d00)*(s-1.d00)
		  g2=9.d00/16.d00*(1.d00-s*s)*onem3s
		  g3=9.d00/16.d00*(1.d00-s*s)*onep3s
		  g4=1.d00/16.d00*(9.d00*s*s-1.d00)*(s+1.d00)
c
		  g1x=-1.d00/16.d00*(18.d00*s)*(s-1.d00)
     &               -1.d00/16.d00*(9.d00*s*s-1.d00)
		  g2x=9.d00/16.d00*(-2.d00*s)*onem3s
     &               -9.d00/16.d00*(1.d00-s*s)*3.d00
		  g3x=9.d00/16.d00*(-2.d00*s)*onep3s
     &               -9.d00/16.d00*(s*s-1.d00)*3.d00
            g4x=1.d00/16.d00*(18.d00*s)*(s+1.d00)
     &		       +1.d00/16.d00*(9.d00*s*s-1.d00)
c
           shl(3,1,l)=f1*g1
	     shl(3,2,l)=f4*g1
	     shl(3,3,l)=f4*g4
	     shl(3,4,l)=f1*g4
c	     
           shl(3,5,l)=f2*g1
	     shl(3,6,l)=f3*g1
	     shl(3,7,l)=f4*g2
	     shl(3,8,l)=f4*g3
	     shl(3,9,l)=f3*g4
	     shl(3,10,l)=f2*g4
	     shl(3,11,l)=f1*g3
	     shl(3,12,l)=f1*g2
	     shl(3,13,l)=f2*g2
	     shl(3,14,l)=f3*g2
	     shl(3,15,l)=f3*g3
	     shl(3,16,l)=f2*g3
c
c
           shl(1,1,l)=f1x*g1
	     shl(1,2,l)=f4x*g1
	     shl(1,3,l)=f4x*g4
	     shl(1,4,l)=f1x*g4
c	     
           shl(1,5,l)=f2x*g1
	     shl(1,6,l)=f3x*g1
	     shl(1,7,l)=f4x*g2
	     shl(1,8,l)=f4x*g3
	     shl(1,9,l)=f3x*g4
	     shl(1,10,l)=f2x*g4
	     shl(1,11,l)=f1x*g3
	     shl(1,12,l)=f1x*g2
	     shl(1,13,l)=f2x*g2
	     shl(1,14,l)=f3x*g2
	     shl(1,15,l)=f3x*g3
	     shl(1,16,l)=f2x*g3
c
c
           shl(2,1,l)=f1*g1x
	     shl(2,2,l)=f4*g1x
	     shl(2,3,l)=f4*g4x
	     shl(2,4,l)=f1*g4x
c	     
           shl(2,5,l)=f2*g1x
	     shl(2,6,l)=f3*g1x
	     shl(2,7,l)=f4*g2x
	     shl(2,8,l)=f4*g3x
	     shl(2,9,l)=f3*g4x
	     shl(2,10,l)=f2*g4x
	     shl(2,11,l)=f1*g3x
	     shl(2,12,l)=f1*g2x
	     shl(2,13,l)=f2*g2x
	     shl(2,14,l)=f3*g2x
	     shl(2,15,l)=f3*g3x
	     shl(2,16,l)=f2*g3x
c
            end if
  200 continue
c
      return
      end

c**** new ********************************************************************** 
      subroutine legdre1d(shlone,wone,nint,nen)
c
c.... program to calculate integration-rule weights, shape functions
c        and local derivatives for a four-node quadrilateral element
c
c               s,t = local element coordinates ("xi", "eta", resp.)       
c        shl(1,i,l) = local ("xi") derivative of shape function
c        shl(2,i,l) = local ("eta") derivative of shape function
c        shl(3,i,l) = local  shape function
c              w(l) = integration-rule weight
c                 i = local node number
c                 l = integration point number
c              nint = number of integration points, eq. 1 or 4
c
      implicit real*8 (a-h,o-z)
c
c.... remove above card for single precision operation
c
      dimension wone(*),raone(8)
	dimension shlone(3,nen,*)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
      data  five9/0.5555555555555555d0/,eight9/0.8888888888888888d0/
      data r1/0.d00/,w1/2.d00/,
     &     r2/0.577350269189626d00/,w2/1.d00/,
     &     r3a/0.774596669241483d00/,w3a/0.555555555555556d00/,
     &     r3b/0.d00/,w3b/0.888888888888889d00/,
     &     r4a/0.861136311594053d00/,w4a/0.347854845137454d00/,
     &     r4b/0.339981043584856d00/,w4b/0.652145154862546d00/
c
      if (nint.eq.1) then
         wone(1)  = two
         raone(1) = zero
	nintx=1
	ninty=1
      endif
c
c
      if (nint.eq.2) then
         wone(1) = one
         wone(2) = one
         raone(1)=-.577350269189626
         raone(2)= .577350269189625
      endif
c
c

      if (nint.eq.3) then
         wone(1) = five9
         wone(2) = five9
         wone(3) = eight9
         raone(1)=-.774596669241483
         raone(2)= .774596669241483
         raone(3)= zero
      endif
c 
c
      if (nint.eq.4) then
         wone(1) = .347854845137454
         wone(2) = .347854845137454
         wone(3) = .652145154862546
         wone(4) = .652145154862546
         raone(1)=-.861136311594053
         raone(2)= .861136311594053
         raone(3)=-.339981043584856
         raone(4)= .339981043584856
      endif
c
c
       if(nint.eq.5) then
        wone(1) = .236926885056189
        wone(2) = .236926885056189
        wone(3) = .478628670499366
        wone(4) = .478628670499366
        wone(5) = .568888888888888
        raone(1)=-.906179845938664
        raone(2)= .906179845938664
        raone(3)=-.538469310105683
        raone(4)= .538469310105683
        raone(5)= zero
       endif
c
c
       if(nint.eq.6) then
         wone(1) = .171324492397170
         wone(2) = .171324492397170
         wone(3) = .360761573048139
         wone(4) = .360761573048139
         wone(5) = .467913934572691
         wone(6) = .467913934572691
         raone(1)=-.932469514203152
         raone(2)= .932469514203152
         raone(3)=-.661209386466265
         raone(4)= .661209386466365
         raone(5)=-.238619186083197
         raone(6)= .238619186083197
        endif
c
c
       if(nint.eq.7) then
         wone(1) = .129484966168870
         wone(2) = .129484966168870
         wone(3) = .279705391489277 
         wone(4) = .279705391489277
         wone(5) = .381830050505119
         wone(6) = .381830050505119
         wone(7) = .417959183673469
         raone(1)=-.949107912342759
         raone(2)= .949107912342759
         raone(3)=-.741531185599394
         raone(4)= .741531185599394
         raone(5)=-.405845151377397
         raone(6)= .405845151377397
         raone(7)= zero
        endif
c
c
       if(nint.eq.8) then
         wone(1) = .101228536290376
         wone(2) = .101228536290376
         wone(3) = .222381034453374
         wone(4) = .222381034453374
         wone(5) = .313706645877887
         wone(6) = .313706645877887
         wone(7) = .362683783378362
         wone(8) = .362683783378362
         raone(1)=-.960289856497536
         raone(2)= .960289856497536
         raone(3)=-.796666477413627
         raone(4)= .796666477413627
         raone(5)=-.525532409916329
         raone(6)= .525532409916329
         raone(7)=-.183434642495650
         raone(8)= .183434642495650
        endif
c
c    polinomios de Legendre
c
      do 100 l = 1, nint
         r = raone(l)
c
      shlone(1,1,l) = zero
      shlone(2,1,l) = one
        if(nen.eq.1) go to 100
      shlone(1,2,l) = one
      shlone(2,2,l) = r
        if(nen.eq.2) go to 100
      shlone(1,3,l) = 3.D0 * r
      shlone(2,3,l) = 0.3D1 / 0.2D1 * r ** 2 
     #              - 0.1D1 / 0.2D1
        if(nen.eq.3) go to 100
      shlone(1,4,l) = 0.15D2 / 0.2D1 * r ** 2 
     #              - 0.3D1 / 0.2D1
      shlone(2,4,l) = 0.5D1 / 0.2D1 * r ** 3 
     #              - 0.3D1 / 0.2D1 * r
        if(nen.eq.4) go to 100
      shlone(1,5,l) = 0.35D2 / 0.2D1 * r ** 3 
     #              - 0.15D2 / 0.2D1 * r
      shlone(2,5,l) = 0.35D2 / 0.8D1 * r ** 4 
     #              - 0.15D2 / 0.4D1 * r ** 2 
     #              + 0.3D1 / 0.8D1
        if(nen.eq.5) go to 100
      shlone(1,6,l) = 0.315D3 / 0.8D1 * r ** 4 
     #              - 0.105D3 / 0.4D1 * r ** 2 
     #              + 0.15D2 / 0.8D1
      shlone(2,6,l) = 0.63D2 / 0.8D1 * r ** 5 
     #              - 0.35D2 / 0.4D1 * r ** 3 
     #              + 0.15D2 / 0.8D1 * r
        if(nen.eq.6) go to 100
      shlone(1,7,l) = 0.693D3 / 0.8D1 * r ** 5 
     #              - 0.315D3 / 0.4D1 * r ** 3 
     #              + 0.105D3 / 0.8D1 * r
      shlone(2,7,l) = 0.231D3 / 0.16D2 * r ** 6 
     #              - 0.315D3 / 0.16D2 * r ** 4 
     #              + 0.105D3 / 0.16D2 * r ** 2 
     #              - 0.5D1 / 0.16D2
        if(nen.eq.7) go to 100
      shlone(1,8,l) = 0.3003D4 / 0.16D2 * r ** 6 
     #              - 0.3465D4 / 0.16D2 * r ** 4 
     #              + 0.945D3 / 0.16D2 * r ** 2 
     #              - 0.35D2 / 0.16D2
      shlone(2,8,l) = 0.429D3 / 0.16D2 * r ** 7 
     #              - 0.693D3 / 0.16D2 * r ** 5 
     #              + 0.315D3 / 0.16D2 * r ** 3 
     #              - 0.35D2 / 0.16D2 * r
c
  100  continue      
c
c
      return
      end
c**** new ********************************************************************** 
c**** new ********************************************************************** 
      subroutine shlegdre(shl,w,nint,nen)
c
c.... program to calculate integration-rule weights, shape functions
c        and local derivatives for a four-node quadrilateral element
c
c               s,t = local element coordinates ("xi", "eta", resp.)       
c        shl(1,i,l) = local ("xi") derivative of shape function
c        shl(2,i,l) = local ("eta") derivative of shape function
c        shl(3,i,l) = local  shape function
c              w(l) = integration-rule weight
c                 i = local node number
c                 l = integration point number
c              nint = number of integration points, eq. 1 or 4
c
      implicit real*8 (a-h,o-z)
c
c.... remove above card for single precision operation
c
      dimension wone(8),raone(8)
	dimension shlone(2,8,8)
      dimension shl(3,nen,*),w(*)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
      data  five9/0.5555555555555555d0/,eight9/0.8888888888888888d0/
      data r1/0.d00/,w1/2.d00/,
     &     r2/0.577350269189626d00/,w2/1.d00/,
     &     r3a/0.774596669241483d00/,w3a/0.555555555555556d00/,
     &     r3b/0.d00/,w3b/0.888888888888889d00/,
     &     r4a/0.861136311594053d00/,w4a/0.347854845137454d00/,
     &     r4b/0.339981043584856d00/,w4b/0.652145154862546d00/
c
      if (nint.eq.1) then
         wone(1)  = two
         raone(1) = zero
	nintx=1
	ninty=1
      endif
c
      if (nen.eq.1) then
	nenx=1
	neny=1
      endif
c
c
      if (nint.eq.4) then
         wone(1) = one
         wone(2) = one
         raone(1)=-.577350269189626
         raone(2)= .577350269189625
	nintx=2
	ninty=2
      endif
c
      if (nen.eq.4) then
	nenx=2
	neny=2
      endif
c
c
c
      if (nint.eq.9) then
         wone(1) = five9
         wone(2) = five9
         wone(3) = eight9
         raone(1)=-.774596669241483
         raone(2)= .774596669241483
         raone(3)= zero
	nintx=3
	ninty=3
      endif
c 
      if(nen.eq.9) then
	nenx=3
	neny=3
      endif
c
      if (nint.eq.16) then
         wone(1) = .347854845137454
         wone(2) = .347854845137454
         wone(3) = .652145154862546
         wone(4) = .652145154862546
         raone(1)=-.861136311594053
         raone(2)= .861136311594053
         raone(3)=-.339981043584856
         raone(4)= .339981043584856
	nintx=4
	ninty=4
      endif
c
      if (nen.eq.16) then
	nenx=4
	neny=4
         endif
c
c
       if(nint.eq.25) then
        wone(1) = .236926885056189
        wone(2) = .236926885056189
        wone(3) = .478628670499366
        wone(4) = .478628670499366
        wone(5) = .568888888888888
        raone(1)=-.906179845938664
        raone(2)= .906179845938664
        raone(3)=-.538469310105683
        raone(4)= .538469310105683
        raone(5)= zero
	nintx=5
	ninty=5
       endif
c
       if(nen.eq.25) then
	nenx=5
	neny=5
       endif
c
       if(nint.eq.36) then
         wone(1) = .171324492397170
         wone(2) = .171324492397170
         wone(3) = .360761573048139
         wone(4) = .360761573048139
         wone(5) = .467913934572691
         wone(6) = .467913934572691
         raone(1)=-.932469514203152
         raone(2)= .932469514203152
         raone(3)=-.661209386466265
         raone(4)= .661209386466365
         raone(5)=-.238619186083197
         raone(6)= .238619186083197
	nintx=6
	ninty=6
        endif
c
        if(nen.eq.36) then
	nenx=6
	neny=6
        endif
c
       if(nint.eq.49) then
         wone(1) = .129484966168870
         wone(2) = .129484966168870
         wone(3) = .279705391489277 
         wone(4) = .279705391489277
         wone(5) = .381830050505119
         wone(6) = .381830050505119
         wone(7) = .417959183673469
         raone(1)=-.949107912342759
         raone(2)= .949107912342759
         raone(3)=-.741531185599394
         raone(4)= .741531185599394
         raone(5)=-.405845151377397
         raone(6)= .405845151377397
         raone(7)= zero
	nintx=7
	ninty=7
        endif
c
        if(nen.eq.49) then
	nenx=7
	neny=7
        endif
c
       if(nint.eq.64) then
         wone(1) = .101228536290376
         wone(2) = .101228536290376
         wone(3) = .222381034453374
         wone(4) = .222381034453374
         wone(5) = .313706645877887
         wone(6) = .313706645877887
         wone(7) = .362683783378362
         wone(8) = .362683783378362
         raone(1)=-.960289856497536
         raone(2)= .960289856497536
         raone(3)=-.796666477413627
         raone(4)= .796666477413627
         raone(5)=-.525532409916329
         raone(6)= .525532409916329
         raone(7)=-.183434642495650
         raone(8)= .183434642495650
	nintx=8
	ninty=8
        endif
        if(nen.eq.64) then
	nenx=8
	neny=8
        endif
c
c    polinomios de Legendre
c
      do 100 l = 1, nintx
         r = raone(l)
c
      shlone(1,1,l) = zero
      shlone(2,1,l) = one
        if(nenx.eq.1) go to 100
      shlone(1,2,l) = one
      shlone(2,2,l) = r
        if(nenx.eq.2) go to 100
      shlone(1,3,l) = 3.D0 * r
      shlone(2,3,l) = 0.3D1 / 0.2D1 * r ** 2 - 0.1D1 / 0.2D1
        if(nenx.eq.3) go to 100
      shlone(1,4,l) = 0.15D2 / 0.2D1 * r ** 2 - 0.3D1 / 0.2D1
      shlone(2,4,l) = 0.5D1 / 0.2D1 * r ** 3 - 0.3D1 / 0.2D1 * r
        if(nenx.eq.4) go to 100
      shlone(1,5,l) = 0.35D2 / 0.2D1 * r ** 3 - 0.15D2 / 0.2D1 * r
      shlone(2,5,l) = 0.35D2 / 0.8D1 * r ** 4 - 0.15D2 / 0.4D1 * r ** 2 
     #              + 0.3D1 / 0.8D1
        if(nenx.eq.5) go to 100
      shlone(1,6,l) = 0.315D3 / 0.8D1 * r ** 4 
     #              - 0.105D3 / 0.4D1 * r ** 2 + 0.15D2 / 0.8D1
      shlone(2,6,l) = 0.63D2 / 0.8D1 * r ** 5 - 0.35D2 / 0.4D1 * r ** 3 
     #              + 0.15D2 / 0.8D1 * r
        if(nenx.eq.6) go to 100
      shlone(1,7,l) = 0.693D3 / 0.8D1 * r ** 5 
     #              - 0.315D3 / 0.4D1 * r ** 3 + 0.105D3 / 0.8D1 * r
      shlone(2,7,l) = 0.231D3 / 0.16D2 * r ** 6 
     #              - 0.315D3 / 0.16D2 * r ** 4 
     #              + 0.105D3 / 0.16D2 * r ** 2 - 0.5D1 / 0.16D2
        if(nenx.eq.7) go to 100
      shlone(1,8,l) = 0.3003D4 / 0.16D2 * r ** 6 
     #              - 0.3465D4 / 0.16D2 * r ** 4 
     #              + 0.945D3 / 0.16D2 * r ** 2 - 0.35D2 / 0.16D2
      shlone(2,8,l) = 0.429D3 / 0.16D2 * r ** 7 
     #              - 0.693D3 / 0.16D2 * r ** 5 
     #              + 0.315D3 / 0.16D2 * r ** 3 - 0.35D2 / 0.16D2 * r
c
  100  continue      
c
         l=0
	   do ly=1,ninty
	   do lx=1,nintx
	   l = l+1
	      w(l) = wone(lx)*wone(ly)
         end do
	   end do
c
      l=0
	do ly=1,ninty
      do lx=1,nintx
	  l = l+1
	  j=0
	 do iy=1,neny
	 do ix=1,nenx
         j = j + 1
	   shl(1,j,l) = shlone(1,ix,lx)*shlone(2,iy,ly)
	   shl(2,j,l) = shlone(2,ix,lx)*shlone(1,iy,ly)
	   shl(3,j,l) = shlone(2,ix,lx)*shlone(2,iy,ly)
	 end do
	 end do
	end do
	end do
c
      return
      end
c**** new ********************************************************************** 
      subroutine shlqpbk(shl,nen,nside,nnods,nints)
c
c.... program to calculate integration-rule weights, shape functions
c        and local derivatives for a four-node quadrilateral element
c
c               s,t = local element coordinates ("xi", "eta", resp.)       
c        shl(1,i,l) = local ("xi") derivative of shape function
c        shl(2,i,l) = local ("eta") derivative of shape function
c        shl(3,i,l) = local  shape function
c              w(l) = integration-rule weight
c                 i = local node number
c                 l = integration point number
c              nint = number of integration points, eq. 1 or 4
c
      implicit real*8 (a-h,o-z)
c
c.... remove above card for single precision operation
c
      dimension raone(8),xaone(8),paone(8)
	dimension shlx(2,8),shly(2,8),inod(8,8)
      dimension shl(3,nen,*)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
      data  five9/0.5555555555555555d0/,eight9/0.8888888888888888d0/
      data r1/0.d00/,w1/2.d00/,
     &     r2/0.577350269189626d00/,w2/1.d00/,
     &     r3a/0.774596669241483d00/,w3a/0.555555555555556d00/,
     &     r3b/0.d00/,w3b/0.888888888888889d00/,
     &     r4a/0.861136311594053d00/,w4a/0.347854845137454d00/,
     &     r4b/0.339981043584856d00/,w4b/0.652145154862546d00/
c

      if (nints.eq.1) then
         raone(1) = zero
      endif
c
      if(nnods.eq.1) then
        xaone(1) = zero
        paone(1) = zero
      end if
c
      if (nints.eq.2) then
         raone(1)=-.577350269189626
         raone(2)= .577350269189625
c
         paone(1)= .577350269189626
         paone(2)=-.577350269189625
c
      endif
c
      if (nnods.eq.2) then
         xaone(1) = -one
         xaone(2) =  one
      endif
c
c
      if (nints.eq.3) then
         raone(1)=-.774596669241483
         raone(2)= zero
         raone(3)= .774596669241483
c
         paone(1)= .774596669241483
         paone(2)= zero
         paone(3)=-.774596669241483
c
      endif
c 
      if(nnods.eq.3) then
         xaone(1)= -one
         xaone(2)= one
         xaone(3)= zero
      endif
c
      if (nints.eq.4) then
         raone(1)=-.861136311594053
         raone(2)=-.339981043584856
         raone(3)= .339981043584856
         raone(4)= .861136311594053
c
         paone(1)= .861136311594053
         paone(2)= .339981043584856
         paone(3)=-.339981043584856
         paone(4)=-.861136311594053
c
      endif
c
      if (nnods.eq.4) then
         xaone(1) = -one
         xaone(2) = one
         xaone(3) = -.333333333333333
         xaone(4) =  .333333333333333
      endif
c
c
       if(nints.eq.5) then
        raone(1)=-.906179845938664
        raone(2)=-.538469310105683
        raone(3)= zero
        raone(4)= .538469310105683
        raone(5)= .906179845938664
c
        paone(1)= .906179845938664
        paone(2)= .538469310105683
        paone(3)= zero
        paone(4)=-.538469310105683
        paone(5)=-.906179845938664
c
       endif
c
       if(nnods.eq.5) then
         xaone(1)= -one 
         xaone(2)=  one 
         xaone(3)= -pt5
         xaone(4)= zero
         xaone(5)= pt5
       endif
c
       if(nints.eq.6) then
         raone(1)=-.932469514203152
         raone(2)=-.661209386466265
         raone(3)=-.238619186083197
         raone(4)= .238619186083197
         raone(5)= .661209386466365
         raone(6)= .932469514203152
c
         paone(1)= .932469514203152
         paone(2)= .661209386466265
         paone(3)= .238619186083197
         paone(4)=-.238619186083197
         paone(5)=-.661209386466365
         paone(6)=-.932469514203152
c
        endif
c
        if(nnods.eq.6) then
         xaone(1) = -one
         xaone(2) =  one
         xaone(3) = -.600000000000000
         xaone(4) = -.200000000000000
         xaone(5) =  .200000000000000
         xaone(6) =  .600000000000000
        endif
c
       if(nints.eq.7) then
         raone(1)=-.949107912342759
         raone(2)=-.741531185599394
         raone(3)=-.405845151377397
         raone(4)= zero
         raone(5)= .405845151377397
         raone(6)= .741531185599394
         raone(7)= .949107912342759
c
         paone(1)= .949107912342759
         paone(2)= .741531185599394
         paone(3)= .405845151377397
         paone(4)= zero
         paone(5)=-.405845151377397
         paone(6)=-0.741531185599394
         paone(7)=-.949107912342759
c
        endif
c
        if(nnods.eq.7) then
         xaone(1) = -one
         xaone(2) =  one
         xaone(3) = -.666666666666666
         xaone(4) = -.333333333333333
         xaone(5) = zero
         xaone(6) =  .333333333333333
         xaone(7) =  .666666666666666
        endif
c
       if(nints.eq.8) then
         raone(1)=-.960289856497536
         raone(2)=-.796666477413627
         raone(3)=-.525532409916329
         raone(4)=-.183434642495650
         raone(5)= .183434642495650
         raone(6)= .525532409916329
         raone(7)= .796666477413627
         raone(8)= .960289856497536
c
         paone(1)= .960289856497536
         paone(2)= .796666477413627
         paone(3)= .525532409916329
         paone(4)= .183434642495650
         paone(5)=-.183434642495650
         paone(6)=-.525532409916329
         paone(7)=-.796666477413627
         paone(8)=-.960289856497536
c
        endif
        if(nnods.eq.8) then
         xaone(1) = -one
         xaone(2) =  one
         xaone(3) = -0.71428571428571
         xaone(4) = -0.42857142857143
         xaone(5) = -0.14285714285714
         xaone(6) = 0.14285714285714
         xaone(7) = 0.42857142857143
         xaone(8) = 0.71428571428571
        endif
c
c
	if(nen.eq.1) inod(1,1) = 1
	if(nen.eq.4) then
	  inod(1,1) = 1
	  inod(2,1) = 2
	  inod(1,2) = 4
	  inod(2,2) = 3
	end if
c
	if(nen.eq.9) then
	  inod(1,1) = 1
	  inod(2,1) = 2
	  inod(3,1) = 5
c
	  inod(1,2) = 4
	  inod(2,2) = 3
	  inod(3,2) = 7
c
	  inod(1,3) = 8
	  inod(2,3) = 6
	  inod(3,3) = 9
	end if
c
	if(nen.eq.16) then
	  inod(1,1) = 1
	  inod(2,1) = 2
	  inod(3,1) = 5
	  inod(4,1) = 6
c
	  inod(1,2) = 4
	  inod(2,2) = 3
	  inod(3,2) = 10
	  inod(4,2) = 9
c
	  inod(1,3) = 12
	  inod(2,3) = 7
	  inod(3,3) = 13
	  inod(4,3) = 14
c
	  inod(1,4) = 11
	  inod(2,4) = 8
	  inod(3,4) = 16
	  inod(4,4) = 15
	end if
c
c
	if(nen.eq.25) then
	  inod(1,1) = 1
	  inod(2,1) = 2
	  inod(3,1) = 5
	  inod(4,1) = 6
	  inod(5,1) = 7
c
	  inod(1,2) = 4
	  inod(2,2) = 3
	  inod(3,2) = 13
	  inod(4,2) = 12
	  inod(5,2) = 11
c
	  inod(1,3) = 16
	  inod(2,3) = 8
	  inod(3,3) = 17
	  inod(4,3) = 18
	  inod(5,3) = 19
c
	  inod(1,4) = 15
	  inod(2,4) = 9
	  inod(3,4) = 24
	  inod(4,4) = 25
	  inod(5,4) = 20
c	  
	  inod(1,5) = 14
	  inod(2,5) = 10
	  inod(3,5) = 23
	  inod(4,5) = 22
	  inod(5,5) = 21
	end if
c
	if(nen.eq.36) then
	  inod(1,1) = 1
	  inod(2,1) = 2
	  inod(3,1) = 5
	  inod(4,1) = 6
	  inod(5,1) = 7
	  inod(6,1) = 8
c
	  inod(1,2) = 4
	  inod(2,2) = 3
	  inod(3,2) = 16
	  inod(4,2) = 15
	  inod(5,2) = 14
	  inod(6,2) = 13
c
	  inod(1,3) = 20
	  inod(2,3) = 9
	  inod(3,3) = 21
	  inod(4,3) = 22
	  inod(5,3) = 23
	  inod(6,3) = 24
c
	  inod(1,4) = 19
	  inod(2,4) = 10
	  inod(3,4) = 32
	  inod(4,4) = 33
	  inod(5,4) = 34
	  inod(6,4) = 25
c	  
	  inod(1,5) = 18
	  inod(2,5) = 11
	  inod(3,5) = 31
	  inod(4,5) = 36
	  inod(5,5) = 35
	  inod(6,5) = 26
c
	  inod(1,6) = 17
	  inod(2,6) = 12
	  inod(3,6) = 30
	  inod(4,6) = 29
	  inod(5,6) = 28
	  inod(6,6) = 27
c	  
	end if
c
	if(nen.eq.49) then
	  inod(1,1) = 1
	  inod(2,1) = 2
	  inod(3,1) = 5
	  inod(4,1) = 6
	  inod(5,1) = 7
	  inod(6,1) = 8
	  inod(7,1) = 9
c
	  inod(1,2) = 4
	  inod(2,2) = 3
	  inod(3,2) = 19
	  inod(4,2) = 18
	  inod(5,2) = 17
	  inod(6,2) = 16
	  inod(7,2) = 15
c
	  inod(1,3) = 24
	  inod(2,3) = 10
	  inod(3,3) = 25
	  inod(4,3) = 26
	  inod(5,3) = 27
	  inod(6,3) = 28
	  inod(7,3) = 29
c
	  inod(1,4) = 23
	  inod(2,4) = 11
	  inod(3,4) = 40
	  inod(4,4) = 41
	  inod(5,4) = 42
	  inod(6,4) = 43
	  inod(7,4) = 30
c	  
	  inod(1,5) = 22
	  inod(2,5) = 12
	  inod(3,5) = 39
	  inod(4,5) = 48
	  inod(5,5) = 49
	  inod(6,5) = 44
	  inod(7,5) = 31
c
	  inod(1,6) = 21
	  inod(2,6) = 13
	  inod(3,6) = 38
	  inod(4,6) = 47
	  inod(5,6) = 46
	  inod(6,6) = 45
	  inod(7,6) = 32
c
	  inod(1,7) = 20
	  inod(2,7) = 14
	  inod(3,7) = 37
	  inod(4,7) = 36
	  inod(5,7) = 35
	  inod(6,7) = 34
	  inod(7,7) = 33
c	  
	end if
c
c
	if(nen.eq.64) then
	  inod(1,1) = 1
	  inod(2,1) = 2
	  inod(3,1) = 5
	  inod(4,1) = 6
	  inod(5,1) = 7
	  inod(6,1) = 8
	  inod(7,1) = 9
	  inod(8,1) = 10
c
	  inod(1,2) = 4
	  inod(2,2) = 3
	  inod(3,2) = 22
	  inod(4,2) = 21
	  inod(5,2) = 20
	  inod(6,2) = 19
	  inod(7,2) = 18
	  inod(8,2) = 17
c
	  inod(1,3) = 28
	  inod(2,3) = 11
	  inod(3,3) = 29
	  inod(4,3) = 30
	  inod(5,3) = 31
	  inod(6,3) = 32
	  inod(7,3) = 33
	  inod(8,3) = 34
c
	  inod(1,4) = 27
	  inod(2,4) = 12
	  inod(3,4) = 48
	  inod(4,4) = 49
	  inod(5,4) = 50
	  inod(6,4) = 51
	  inod(7,4) = 52
	  inod(8,4) = 35
c	  
	  inod(1,5) = 26
	  inod(2,5) = 13
	  inod(3,5) = 47
	  inod(4,5) = 60
	  inod(5,5) = 61
	  inod(6,5) = 62
	  inod(7,5) = 53
	  inod(8,5) = 36
c
	  inod(1,6) = 25
	  inod(2,6) = 14
	  inod(3,6) = 46
	  inod(4,6) = 59
	  inod(5,6) = 64
	  inod(6,6) = 63
	  inod(7,6) = 54
	  inod(8,6) = 37
c
	  inod(1,7) = 24
	  inod(2,7) = 15
	  inod(3,7) = 45
	  inod(4,7) = 58
	  inod(5,7) = 57
	  inod(6,7) = 56
	  inod(7,7) = 55
	  inod(8,7) = 38
c
	  inod(1,8) = 23
	  inod(2,8) = 16
	  inod(3,8) = 44
	  inod(4,8) = 43
	  inod(5,8) = 42
	  inod(6,8) = 41
	  inod(7,8) = 40
	  inod(8,8) = 39
c	  
	end if
c
      lb=0
      do 200 ns=1,nside    
      do 200 l = 1, nints
      lb = lb + 1
      if(ns.eq.1) then
         r = raone(l)
         s = -one
      end if 
      if(ns.eq.2) then
         r = one
         s = raone(l)
      end if 
c
      if(ns.eq.3) then
         r = paone(l)
         s = one
      end if 
c
      if(ns.eq.4) then
         r = -one
         s = paone(l)
      end if 
c
        shlx(1,1) = zero
        shly(1,1) = zero
        shlx(2,1) = one
        shly(2,1) = one
      if(nnods.eq.1) go to 100
c
        do 50 i = 1, nnods
         aa = one
         bb = one
         cc = one
         aax= zero
         aay= zero
         do 40 j =1, nnods
          daj = one
          caj = one
          if (i .ne. j)then
          aa = aa * ( r - xaone(j))
          cc = cc * ( s - xaone(j))
          bb = bb * ( xaone(i) - xaone(j))
          do 30 k = 1, nnods
           if(k .ne. i .and. k .ne. j) daj = daj * ( r - xaone(k))
           if(k .ne. i .and. k .ne. j) caj = caj * ( s - xaone(k))
   30     continue
          aax = aax + daj
          aay = aay + caj
          endif
   40    continue
        shlx(1,i) = aax/bb
        shly(1,i) = aay/bb
        shlx(2,i) = aa/bb
        shly(2,i) = cc/bb
   50  continue
c
c
  100  continue 
c
	 do iy=1,nnods
	 do ix=1,nnods
         j = inod(ix,iy)
	   shl(1,j,lb) = shlx(1,ix)*shly(2,iy)
	   shl(2,j,lb) = shlx(2,ix)*shly(1,iy)
	   shl(3,j,lb) = shlx(2,ix)*shly(2,iy)
	 end do
	 end do
 200  continue
c
      return
      end
c**** new **********************************************************************
c**** new ********************************************************************** 
      subroutine shlqpk(shl,w,nint,nen)
c
c.... program to calculate integration-rule weights, shape functions
c        and local derivatives for a four-node quadrilateral element
c
c               s,t = local element coordinates ("xi", "eta", resp.)       
c        shl(1,i,l) = local ("xi") derivative of shape function
c        shl(2,i,l) = local ("eta") derivative of shape function
c        shl(3,i,l) = local  shape function
c              w(l) = integration-rule weight
c                 i = local node number
c                 l = integration point number
c              nint = number of integration points, eq. 1 or 4
c
      implicit real*8 (a-h,o-z)
c
c.... remove above card for single precision operation
c
      dimension wone(8),raone(8),xaone(8)
	dimension shlone(2,8,8),inod(8,8)
      dimension shl(3,nen,*),w(*),ra(64),sa(64)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
      data  five9/0.5555555555555555d0/,eight9/0.8888888888888888d0/
      data r1/0.d00/,w1/2.d00/,
     &     r2/0.577350269189626d00/,w2/1.d00/,
     &     r3a/0.774596669241483d00/,w3a/0.555555555555556d00/,
     &     r3b/0.d00/,w3b/0.888888888888889d00/,
     &     r4a/0.861136311594053d00/,w4a/0.347854845137454d00/,
     &     r4b/0.339981043584856d00/,w4b/0.652145154862546d00/
c
      if (nint.eq.1) then
         wone(1)  = two
         raone(1) = zero
	nintx=1
	ninty=1
      endif
c
      if (nen.eq.1) xaone(1) = zero
c
      if (nint.eq.4) then
         wone(1) = one
         wone(2) = one
         raone(1)=-.577350269189626
         raone(2)= .577350269189625
	nintx=2
	ninty=2
      endif
c
      if (nen.eq.4) then
         xaone(1) = -one
         xaone(2) =  one
	nenx=2
	neny=2
      endif
c
c
c
      if (nint.eq.9) then
         wone(1) = five9
         wone(2) = five9
         wone(3) = eight9
         raone(1)=-.774596669241483
         raone(2)= .774596669241483
         raone(3)= zero
	nintx=3
	ninty=3
      endif
c 
      if(nen.eq.9) then
         xaone(1)= -one
         xaone(2)= one
         xaone(3)= zero
	nenx=3
	neny=3
      endif
c
      if (nint.eq.16) then
         wone(1) = .347854845137454
         wone(2) = .347854845137454
         wone(3) = .652145154862546
         wone(4) = .652145154862546
         raone(1)=-.861136311594053
         raone(2)= .861136311594053
         raone(3)=-.339981043584856
         raone(4)= .339981043584856
	nintx=4
	ninty=4
      endif
c
      if (nen.eq.16) then
         xaone(1) = -one
         xaone(2) = one
         xaone(3) = -.333333333333333
         xaone(4) =  .333333333333333
	nenx=4
	neny=4
         endif
c
c
c
       if(nint.eq.25) then
        wone(1) = .236926885056189
        wone(2) = .236926885056189
        wone(3) = .478628670499366
        wone(4) = .478628670499366
        wone(5) = .568888888888888
        raone(1)=-.906179845938664
        raone(2)= .906179845938664
        raone(3)=-.538469310105683
        raone(4)= .538469310105683
        raone(5)= zero
	nintx=5
	ninty=5
       endif
c
       if(nen.eq.25) then
         xaone(1)= -one 
         xaone(2)=  one 
         xaone(3)= -pt5
         xaone(4)= zero
         xaone(5)= pt5
	nenx=5
	neny=5
       endif
c
       if(nint.eq.36) then
         wone(1) = .171324492397170
         wone(2) = .171324492397170
         wone(3) = .360761573048139
         wone(4) = .360761573048139
         wone(5) = .467913934572691
         wone(6) = .467913934572691
         raone(1)=-.932469514203152
         raone(2)= .932469514203152
         raone(3)=-.661209386466265
         raone(4)= .661209386466365
         raone(5)=-.238619186083197
         raone(6)= .238619186083197
	nintx=6
	ninty=6
        endif
c
        if(nen.eq.36) then
         xaone(1) = -one
         xaone(2) =  one
         xaone(3) = -.600000000000000
         xaone(4) = -.200000000000000
         xaone(5) =  .200000000000000
         xaone(6) =  .600000000000000
	nenx=6
	neny=6
        endif
c
       if(nint.eq.49) then
         wone(1) = .129484966168870
         wone(2) = .129484966168870
         wone(3) = .279705391489277 
         wone(4) = .279705391489277
         wone(5) = .381830050505119
         wone(6) = .381830050505119
         wone(7) = .417959183673469
         raone(1)=-.949107912342759
         raone(2)= .949107912342759
         raone(3)=-.741531185599394
         raone(4)= .741531185599394
         raone(5)=-.405845151377397
         raone(6)= .405845151377397
         raone(7)= zero
	nintx=7
	ninty=7
        endif
c
        if(nen.eq.49) then
         xaone(1) = -one
         xaone(2) =  one
         xaone(3) = -.666666666666666
         xaone(4) = -.333333333333333
         xaone(5) = zero
         xaone(6) =  .333333333333333
         xaone(7) =  .666666666666666
	nenx=7
	neny=7
        endif
c
       if(nint.eq.64) then
         wone(1) = .101228536290376
         wone(2) = .101228536290376
         wone(3) = .222381034453374
         wone(4) = .222381034453374
         wone(5) = .313706645877887
         wone(6) = .313706645877887
         wone(7) = .362683783378362
         wone(8) = .362683783378362
         raone(1)=-.960289856497536
         raone(2)= .960289856497536
         raone(3)=-.796666477413627
         raone(4)= .796666477413627
         raone(5)=-.525532409916329
         raone(6)= .525532409916329
         raone(7)=-.183434642495650
         raone(8)= .183434642495650
	nintx=8
	ninty=8
        endif
        if(nen.eq.64) then
         xaone(1) = -one
         xaone(2) =  one
         xaone(3) = -0.71428571428571
         xaone(4) = -0.42857142857143
         xaone(5) = -0.14285714285714
         xaone(6) = 0.14285714285714
         xaone(7) = 0.42857142857143
         xaone(8) = 0.71428571428571
	nenx=8
	neny=8
        endif
c
c
	if(nen.eq.1) inod(1,1) = 1
	if(nen.eq.4) then
	  inod(1,1) = 1
	  inod(2,1) = 2
	  inod(1,2) = 4
	  inod(2,2) = 3
	end if
c
	if(nen.eq.9) then
	  inod(1,1) = 1
	  inod(2,1) = 2
	  inod(3,1) = 5
c
	  inod(1,2) = 4
	  inod(2,2) = 3
	  inod(3,2) = 7
c
	  inod(1,3) = 8
	  inod(2,3) = 6
	  inod(3,3) = 9
	end if
c
c
c
	if(nen.eq.16) then
	  inod(1,1) = 1
	  inod(2,1) = 2
	  inod(3,1) = 5
	  inod(4,1) = 6
c
	  inod(1,2) = 4
	  inod(2,2) = 3
	  inod(3,2) = 10
	  inod(4,2) = 9
c
	  inod(1,3) = 12
	  inod(2,3) = 7
	  inod(3,3) = 13
	  inod(4,3) = 14
c
	  inod(1,4) = 11
	  inod(2,4) = 8
	  inod(3,4) = 16
	  inod(4,4) = 15
	end if
c
c
	if(nen.eq.25) then
	  inod(1,1) = 1
	  inod(2,1) = 2
	  inod(3,1) = 5
	  inod(4,1) = 6
	  inod(5,1) = 7
c
	  inod(1,2) = 4
	  inod(2,2) = 3
	  inod(3,2) = 13
	  inod(4,2) = 12
	  inod(5,2) = 11
c
	  inod(1,3) = 16
	  inod(2,3) = 8
	  inod(3,3) = 17
	  inod(4,3) = 18
	  inod(5,3) = 19
c
	  inod(1,4) = 15
	  inod(2,4) = 9
	  inod(3,4) = 24
	  inod(4,4) = 25
	  inod(5,4) = 20
c	  
	  inod(1,5) = 14
	  inod(2,5) = 10
	  inod(3,5) = 23
	  inod(4,5) = 22
	  inod(5,5) = 21
	end if
c
c
	if(nen.eq.36) then
	  inod(1,1) = 1
	  inod(2,1) = 2
	  inod(3,1) = 5
	  inod(4,1) = 6
	  inod(5,1) = 7
	  inod(6,1) = 8
c
	  inod(1,2) = 4
	  inod(2,2) = 3
	  inod(3,2) = 16
	  inod(4,2) = 15
	  inod(5,2) = 14
	  inod(6,2) = 13
c
	  inod(1,3) = 20
	  inod(2,3) = 9
	  inod(3,3) = 21
	  inod(4,3) = 22
	  inod(5,3) = 23
	  inod(6,3) = 24
c
	  inod(1,4) = 19
	  inod(2,4) = 10
	  inod(3,4) = 32
	  inod(4,4) = 33
	  inod(5,4) = 34
	  inod(6,4) = 25
c	  
	  inod(1,5) = 18
	  inod(2,5) = 11
	  inod(3,5) = 31
	  inod(4,5) = 36
	  inod(5,5) = 35
	  inod(6,5) = 26
c
	  inod(1,6) = 17
	  inod(2,6) = 12
	  inod(3,6) = 30
	  inod(4,6) = 29
	  inod(5,6) = 28
	  inod(6,6) = 27
c	  
	end if
c
c
	if(nen.eq.49) then
	  inod(1,1) = 1
	  inod(2,1) = 2
	  inod(3,1) = 5
	  inod(4,1) = 6
	  inod(5,1) = 7
	  inod(6,1) = 8
	  inod(7,1) = 9
c
	  inod(1,2) = 4
	  inod(2,2) = 3
	  inod(3,2) = 19
	  inod(4,2) = 18
	  inod(5,2) = 17
	  inod(6,2) = 16
	  inod(7,2) = 15
c
	  inod(1,3) = 24
	  inod(2,3) = 10
	  inod(3,3) = 25
	  inod(4,3) = 26
	  inod(5,3) = 27
	  inod(6,3) = 28
	  inod(7,3) = 29
c
	  inod(1,4) = 23
	  inod(2,4) = 11
	  inod(3,4) = 40
	  inod(4,4) = 41
	  inod(5,4) = 42
	  inod(6,4) = 43
	  inod(7,4) = 30
c	  
	  inod(1,5) = 22
	  inod(2,5) = 12
	  inod(3,5) = 39
	  inod(4,5) = 48
	  inod(5,5) = 49
	  inod(6,5) = 44
	  inod(7,5) = 31
c
	  inod(1,6) = 21
	  inod(2,6) = 13
	  inod(3,6) = 38
	  inod(4,6) = 47
	  inod(5,6) = 46
	  inod(6,6) = 45
	  inod(7,6) = 32
c
	  inod(1,7) = 20
	  inod(2,7) = 14
	  inod(3,7) = 37
	  inod(4,7) = 36
	  inod(5,7) = 35
	  inod(6,7) = 34
	  inod(7,7) = 33
c	  
	end if
c
c
	if(nen.eq.64) then
	  inod(1,1) = 1
	  inod(2,1) = 2
	  inod(3,1) = 5
	  inod(4,1) = 6
	  inod(5,1) = 7
	  inod(6,1) = 8
	  inod(7,1) = 9
	  inod(8,1) = 10
c
	  inod(1,2) = 4
	  inod(2,2) = 3
	  inod(3,2) = 22
	  inod(4,2) = 21
	  inod(5,2) = 20
	  inod(6,2) = 19
	  inod(7,2) = 18
	  inod(8,2) = 17
c
	  inod(1,3) = 28
	  inod(2,3) = 11
	  inod(3,3) = 29
	  inod(4,3) = 30
	  inod(5,3) = 31
	  inod(6,3) = 32
	  inod(7,3) = 33
	  inod(8,3) = 34
c
	  inod(1,4) = 27
	  inod(2,4) = 12
	  inod(3,4) = 48
	  inod(4,4) = 49
	  inod(5,4) = 50
	  inod(6,4) = 51
	  inod(7,4) = 52
	  inod(8,4) = 35
c	  
	  inod(1,5) = 26
	  inod(2,5) = 13
	  inod(3,5) = 47
	  inod(4,5) = 60
	  inod(5,5) = 61
	  inod(6,5) = 62
	  inod(7,5) = 53
	  inod(8,5) = 36
c
	  inod(1,6) = 25
	  inod(2,6) = 14
	  inod(3,6) = 46
	  inod(4,6) = 59
	  inod(5,6) = 64
	  inod(6,6) = 63
	  inod(7,6) = 54
	  inod(8,6) = 37
c
	  inod(1,7) = 24
	  inod(2,7) = 15
	  inod(3,7) = 45
	  inod(4,7) = 58
	  inod(5,7) = 57
	  inod(6,7) = 56
	  inod(7,7) = 55
	  inod(8,7) = 38
c
	  inod(1,8) = 23
	  inod(2,8) = 16
	  inod(3,8) = 44
	  inod(4,8) = 43
	  inod(5,8) = 42
	  inod(6,8) = 41
	  inod(7,8) = 40
	  inod(8,8) = 39
c	  
	end if
c
c
      do 100 l = 1, nintx
         r = raone(l)
c
        if(nenx.eq.1) then
        shlone(1,1,l) = zero
        shlone(2,1,l) = one
        go to 100
        endif
c
        do 50 i = 1, nenx
         aa = one
         bb = one
         aax = zero
         do 40 j =1, nenx
          daj = one
          if (i .ne. j)then
          aa = aa * ( r - xaone(j))
          bb = bb * ( xaone(i) - xaone(j))
          do 30 k = 1, nenx
           if(k .ne. i .and. k .ne. j) daj = daj * ( r - xaone(k))
   30     continue
          aax =aax + daj
          endif
   40    continue
        shlone(2,i,l) = aa/bb
        shlone(1,i,l) = aax/bb
   50  continue
c
  100  continue      
c
         l=0
	   do ly=1,ninty
	   do lx=1,nintx
	   l = l+1
	      w(l) = wone(lx)*wone(ly)
         end do
	   end do
c
      l=0
	do ly=1,ninty
      do lx=1,nintx
	  l = l+1
	  r=ra(l)
	  s=sa(l)
	 do iy=1,neny
	 do ix=1,nenx
         j = inod(ix,iy)
	   shl(1,j,l) = shlone(1,ix,lx)*shlone(2,iy,ly)
	   shl(2,j,l) = shlone(2,ix,lx)*shlone(1,iy,ly)
	   shl(3,j,l) = shlone(2,ix,lx)*shlone(2,iy,ly)
	 end do
	 end do
	end do
	end do
c
      return
      end
c**** new **********************************************************************
      subroutine flnormc(ien   ,x     ,xl   ,    
     &                 d     ,dl    ,mat   ,
     &                 c     ,ipar  ,
     &                 dlp   ,
     &                 det   ,
     &                 shl   , shg  ,wt    ,
     &                 ddis  ,detp  ,shlp  , 
     &                 shgp ,  
c
     &                 shln  ,shgn  ,
     &                 detn  ,shlb  ,shgb  ,
     &                 detpn ,shlpn ,shgpn ,
     &                 idside,xls   ,idlsd ,
     &                 grav  ,wn    ,
c                  
     &                 numel ,neesq ,nen   ,nsd   ,
     &                 nesd  ,nint  ,neg   ,nrowsh,
     &                 ned   ,nee   ,numnp ,ndof  ,
     &                 ncon  ,nencon,necon ,index ,
     &                 nints,iwrite ,iplt ,nenp   ,
     &                 nside,nnods ,nenlad,npars  , 
     &                 nmultpd, nmultpc,nodsp)
c------------------------------------------------------------------
c     Some variables used in this subroutine
c            j: degree of freedom (1,...ncon)
c         u(j): finite element solution
c        du(j): derivative of finite element solution
c        ue(j): exact solution
c       due(j): derivative of exact solution
c       el2(j): error in L2
c      epri(j): error in the seminorm of H1 (L2 of derivatives)
c     el2el(j): error in L2 in the element domain
c    epriel(j): error in the seminorm of H1 in the element domain
c----------------------------------------------------------------
c
c
c    Program to calculate and print the L2 and H1 seminorm of the error
c    for each degree of freedom in the finite element solution .
c    The trapezoidal rule is used for integration in each element .
c    The number of integration points is given by nints .
c
c    This version is only applicable to 2D.
c
      implicit real*8 (a-h,o-z)
c
c.... remove above card for single-precision operation
c
      character *1 tab
      dimension ien(nen,*),x(nsd,*),xl(nesd,*),d(ndof,*),dl(ned,*),
     &          mat(*),c(8,*),ipar(nodsp,*),
     &          ddis(ned,nenp,*)
      dimension u(6),dux(6),duy(6),ue(6),duex(6),duey(6)
      dimension shl(3,nen,*),shg(3,nen,*),det(*),wt(*)
      dimension el2(6),epri(6),eprix(6),epriy(6),el2el(6),
     &          eprxel(6),epryel(6)
	dimension detp(*),shlp(3,nenp,*),shgp(3,nenp,*),dlp(ned,*)
c
      dimension detn(*),shlb(3,nenlad,*),shgb(3,nenlad,*),
     &          detpn(*),shlpn(3,npars,*),shgpn(3,npars,*),
     &          idside(nside,*),xls(nesd,*),idlsd(*)
	dimension dls(12),grav(*),wn(*)
      dimension shln(3,nnods,*),shgn(3,nnods,*)
c
      common /colhtc/ neqc
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
      common /iounit/ iin,iout,iecho,ioupp,itest1,itest2 
c
      pi=4.*datan(1.d00)
      dpi=2.d00*pi
c
      gf1=grav(1)
      gf2=grav(2)
      gf3=grav(3) 
      beta0=c(1,1) 
	gama=c(3,1)  
	gamas=gama
	if(index.eq.3) gamas=-gama
	fi = pi/6.d00
         tab=char(9)
         xint=0.d00 
         yint=0.d00
c
cc      call clear ( el2, ncon )
cc      call clear (eprix, ncon )
cc      call clear (epriy, ncon )
c
      divu2=0.d00
	edp = 0.d00
	edpx= 0.d00
	edpy= 0.d00
	xmlt= 0.d00
c
c.....  loop on elements
c 
      do 50 n=1,numel
c
      call local(ien(1,n),x,xl,nen,nsd,nesd)
      call local(ipar(1,n),d,dl,nodsp,ndof,ned)
c
c
      do i=1,nenp
	 dlp(1,i) = ddis(1,i,n)
	end do
c
c
c
cc      call clear ( el2el, ncon )
cc      call clear (eprxel, ncon )
cc      call clear (epryel, ncon )
      divue=0.d00
	pmede=0.d00
	dpe = 0.d00
	dpex= 0.d00
	dpey= 0.d00
c
c.....  loop on integration points
c
c
c.....Triangles or quadrilaterals
c
c
c
      call shgqs(xl,detp,shlp,shgp,nint,n,neg,.true.,nenp,shl,nen)
c
c
        do 4040 l=1,nint
        ct=detp(l)*wt(l)
        call clear ( u, ncon )
        call clear (dux, ncon )
        call clear (duy, nncon )
        xint=0.d0
        yint=0.d0
         do 3030 i=1,nen
            xint=xint+shl(3,i,l)*xl(1,i)
            yint=yint+shl(3,i,l)*xl(2,i)
3030     continue
c
c
         call uexafx(xint,yint,ue,duex,duey,gama,fi,index)
c
c    pressao descontinua
c
       pe = 0.d00
	pex = 0.d00
	pey = 0.d00
	do i=1,nenp
	  pe = pe + shgp(3,i,l)*dlp(1,i)
	  pex= pex+ shgp(1,i,l)*dlp(1,i)
	  pey= pey+ shgp(2,i,l)*dlp(1,i)
	end do
	  dpe = dpe + ct*(pe - ue(3))**2
	  dpex= dpex + ct*(pex-duex(3))**2
	  dpey= dpey + ct*(pey-duey(3))**2
c
      if(iwrite.ne.0) then
         write(43,2424) n,l,pe,ue(3)
2424   format(2i10,6e15.5)
222    format(8e20.9)
       end if
4040  continue
c
c
	edp = edp + dpe
	edpx= edpx+ dpex
	edpy= edpy+ dpey
c
c
c......boundary terms - multiplier
c
	xmlte = 0.d00
      do 4000 ns=1,nside
c
c-----localiza os nos do lado ns
c
c
      ns1=idside(ns,1)
      ns2=idside(ns,2)
      nl1=ien(ns1,n)
      nl2=ien(ns2,n)
c      
      if(nl2.gt.nl1) then 
	  sign = 1.d00
        do  nn=1,npars
           idlsd(nn)=idside(ns,nn)
	  end do
	else
	  sign = -1.d00
	  idlsd(1) = idside(ns,2)
	  idlsd(2) = idside(ns,1)
	  id3 = npars-2
	  if(id3.gt.0) then
	    do il=3,npars
	      idlsd(il) = idside(ns,npars+3-il)
	    end do 
	  end if
      end if 
c
c
c-----localiza os parametros do lado ns
c
      do nn=1,npars
	  nld = idlsd(nn)
	  dls(nn) = dl(1,nld)
      end do
c
      do 2000 nn=1,nenlad
      nl=idlsd(nn)
      xls(1,nn)=xl(1,nl)
      xls(2,nn)=xl(2,nl)
c
 2000 continue
c
c
      call oneshgp(xls,detn,shlb,shln,shgn,
     &             nenlad,nnods,nints,nesd,ns,n,neg)      
c
      call oneshgp(xls,detpn,shlb,shlpn,shgpn,
     &             nenlad,npars,nints,nesd,ns,n,neg)      
c  
c
c.....compute boundary integral
c

      do 1000 ls=1,nints
c
      cwn = wn(ls)*detn(ls)
c
c    valores dos parametros do multiplicador
c
	  dhs=0.d00
        do i=1,npars
	    dhs = dhs + dls(i)*shgpn(2,i,ls)
        end do
c
c     geometria
c
        x1 =0.d00
        x2 =0.d00
        dx1=0.d00
        dx2=0.d00
c
        do i=1,nenlad
          x1 =x1 +xls(1,i)*shlb(2,i,ls)
          x2 =x2 +xls(2,i)*shlb(2,i,ls)
          dx1=dx1+xls(1,i)*shlb(1,i,ls)
          dx2=dx2+xls(2,i)*shlb(1,i,ls)
        end do
           dxx=dsqrt(dx1*dx1+dx2*dx2)
           xn1= sign*dx2/dxx
           xn2=-sign*dx1/dxx
c
c    valores exatos dos multiplicadores
c
c
      pix=pi*x1
      piy=pi*x2
      sx=dsin(pix)
      sy=dsin(piy)
      cx=dcos(pix)
      cy=dcos(piy)
c
      pi2=pi*pi
      co=1.d00/(pi2*2.d00+gama)
c
c
      sigma = dsqrt(gama)
      xc = x1*sigma*dcos(fi)
      yc = x2*sigma*dsin(fi)
c
c    valor exato do multiplicador
c
      dhse=gf2*co*cx*cy + gf1*co*sx*sy 
     &   +gf3*dcos(xc + yc)
c
c
      xmlte = xmlte + cwn*(dhse-dhs)**2
c
 1000 continue
 4000 continue
c
      xmlt = xmlt + xmlte
c
  50   continue
c
      divu = dsqrt(divu2)
c
c...... error in all domain
c

	el2p = dsqrt(edp)/co
	eh1p = dsqrt(edpx+edpy)/co
c
	edp0 = dlog10(el2p)
	gradp0= dlog10(eh1p)
c
c      
	edp = dlog10(dsqrt(edp))
	gradp= dlog10(dsqrt(edpx+edpy))
	xmlt= dlog10(dsqrt(xmlt))
c

       xel = dfloat(numel)/dfloat(nnods-1)**2
       xnp = neqc
       xnp = dlog10(xnp)
       xel = dlog10(xel)/2.d00

      write(iplt,2101) xel,xnp,edp,gradp,xmlt
c
      iplt0 = iplt+11
      write(iplt0,2101) xel,xnp,edp0,gradp0
c
      ipl = iplt+20	
      write(ipl,2102) beta0,el2p,eh1p
       return
c
2100  format(5x,i3,4x,e15.7,5x,e15.7)
2101  format(9(e15.5))
2102  format(4(f13.5))
       return
       end
c**** new **********************************************************************
      subroutine flnorm(ien   ,x     ,xl   ,    
     &                 d     ,dl    ,mat   ,
     &                 c     ,ipar  ,
     &                 dlp   ,
     &                 det   ,
     &                 shl   , shg  ,wt    ,
     &                 ddis  ,detp  ,shlp  , 
     &                 shgp ,  
c
     &                 shln  ,shgn  ,
     &                 detn  ,shlb  ,shgb  ,
     &                 detpn ,shlpn ,shgpn ,
     &                 idside,xls   ,idlsd ,
     &                 grav  ,wn    ,
c
     &                 lado  ,
c                  
     &                 numel ,neesq ,nen   ,nsd   ,
     &                 nesd  ,nint  ,neg   ,nrowsh,
     &                 ned   ,nee   ,numnp ,ndof  ,
     &                 ncon  ,nencon,necon ,index ,
     &                 nints,iwrite ,iplt ,nenp   ,
     &                 nside,nnods ,nenlad,npars  , 
     &                 nmultpd, nmultpc,nodsp,nospd)
c------------------------------------------------------------------
c     Some variables used in this subroutine
c            j: degree of freedom (1,...ncon)
c         u(j): finite element solution
c        du(j): derivative of finite element solution
c        ue(j): exact solution
c       due(j): derivative of exact solution
c       el2(j): error in L2
c      epri(j): error in the seminorm of H1 (L2 of derivatives)
c     el2el(j): error in L2 in the element domain
c    epriel(j): error in the seminorm of H1 in the element domain
c----------------------------------------------------------------
c
c
c    Program to calculate and print the L2 and H1 seminorm of the error
c    for each degree of freedom in the finite element solution .
c    The trapezoidal rule is used for integration in each element .
c    The number of integration points is given by nints .
c
c    This version is only applicable to 2D.
c
      implicit real*8 (a-h,o-z)
c
c.... remove above card for single-precision operation
c
      character *1 tab
      dimension ien(nen,*),x(nsd,*),xl(nesd,*),d(ndof,*),dl(ned,*),
     &          mat(*),c(8,*),ipar(nodsp,*),
     &          ddis(ned,nenp,*)
      dimension u(6),dux(6),duy(6),ue(6),duex(6),duey(6)
      dimension shl(3,nen,*),shg(3,nen,*),det(*),wt(*)
      dimension el2(6),epri(6),eprix(6),epriy(6),el2el(6),
     &          eprxel(6),epryel(6)
	dimension detp(*),shlp(3,nenp,*),shgp(3,nenp,*),dlp(ned,*)
c
      dimension detn(*),shlb(3,nenlad,*),shgb(3,nenlad,*),
     &          detpn(*),shlpn(3,npars,*),shgpn(3,npars,*),
     &          idside(nside,*),xls(nesd,*),idlsd(*)
	dimension dls(12),grav(*),wn(*)
      dimension shln(3,nnods,*),shgn(3,nnods,*)
c
      dimension lado(nside,*)
c
      common /colhtc/ neqc
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
      common /iounit/ iin,iout,iecho,ioupp,itest1,itest2 
c
      pi=4.*datan(1.d00)
      dpi=2.d00*pi
c
      gf1=grav(1)
      gf2=grav(2)
      gf3=grav(3) 
      beta0=c(1,1) 
	gama=c(3,1)  
	gamas=gama
	if(index.eq.3) gamas=-gama
	fi = pi/6.d00
         tab=char(9)
         xint=0.d00 
         yint=0.d00
c
cc      call clear ( el2, ncon )
cc      call clear (eprix, ncon )
cc      call clear (epriy, ncon )
c
      divu2=0.d00
	edp = 0.d00
	edpx= 0.d00
	edpy= 0.d00
	xmlt= 0.d00
c
c.....  loop on elements
c      
      do 50 n=1,numel
c
      call local(ien(1,n),x,xl,nen,nsd,nesd)
cc      call local(ipar(1,n),d,dl,nodsp,ndof,ned)
      do ns=1,nside
         ndgs = lado(ns,n)
         nls = (ns-1)*npars
	   ngsd = npars*(ndgs-1)
          do i=1,npars
	     dl(1,nls+i) = d(1,ngsd + i)
	    end do
	end do
c
      do i=1,nenp
	 dlp(1,i) = ddis(1,i,n)
	end do
c
c
c
cc      call clear ( el2el, ncon )
cc      call clear (eprxel, ncon )
cc      call clear (epryel, ncon )
      divue=0.d00
	pmede=0.d00
	dpe = 0.d00
	dpex= 0.d00
	dpey= 0.d00
c
c.....  loop on integration points
c
c
c.....Triangles or quadrilaterals
c
c
c
      call shgqs(xl,detp,shlp,shgp,nint,n,neg,.true.,nenp,shl,nen)
c
c
        do 4040 l=1,nint
        ct=detp(l)*wt(l)
        call clear ( u, ncon )
        call clear (dux, ncon )
        call clear (duy, nncon )
        xint=0.d0
        yint=0.d0
         do 3030 i=1,nen
            xint=xint+shl(3,i,l)*xl(1,i)
            yint=yint+shl(3,i,l)*xl(2,i)
3030     continue
c
c
         call uexafx(xint,yint,ue,duex,duey,gama,fi,index)
c
c    pressao descontinua
c
       pe = 0.d00
	pex = 0.d00
	pey = 0.d00
	do i=1,nenp
	  pe = pe + shgp(3,i,l)*dlp(1,i)
	  pex= pex+ shgp(1,i,l)*dlp(1,i)
	  pey= pey+ shgp(2,i,l)*dlp(1,i)
	end do
	  dpe = dpe + ct*(pe - ue(3))**2
	  dpex= dpex + ct*(pex-duex(3))**2
	  dpey= dpey + ct*(pey-duey(3))**2
c
      if(iwrite.ne.0) then
         write(43,2424) n,l,pe,ue(3)
2424   format(2i10,6e15.5)
222    format(8e20.9)
       end if
4040  continue
c
c
	edp = edp + dpe
	edpx= edpx+ dpex
	edpy= edpy+ dpey
c
c
c......boundary terms - multiplier
c
	xmlte = 0.d00
      do 4000 ns=1,nside
c
c-----localiza os parametros do lado ns
c
      do nn=1,npars
	  nld = (ns-1)*npars + nn
	  dls(nn) = dl(1,nld)
      end do
c
c-----localiza os nos do lado ns
c
c
      ns1=idside(ns,1)
      ns2=idside(ns,2)
      nl1=ien(ns1,n)
      nl2=ien(ns2,n)
c      
      if(nl2.gt.nl1) then 
	  sign = 1.d00
        do  nn=1,npars
           idlsd(nn)=idside(ns,nn)
	  end do
	else
	  sign = -1.d00
	  idlsd(1) = idside(ns,2)
	  idlsd(2) = idside(ns,1)
	  id3 = npars-2
	  if(id3.gt.0) then
	    do il=3,npars
	      idlsd(il) = idside(ns,npars+3-il)
	    end do 
	  end if
      end if 
c
c
      do 2000 nn=1,nenlad
      nl=idlsd(nn)
      xls(1,nn)=xl(1,nl)
      xls(2,nn)=xl(2,nl)
 2000 continue
c
c
      call oneshgp(xls,detn,shlb,shln,shgn,
     &             nenlad,nnods,nints,nesd,ns,n,neg)      
c
      call oneshgp(xls,detpn,shlb,shlpn,shgpn,
     &             nenlad,npars,nints,nesd,ns,n,neg)      
c  
c
c.....compute boundary integral
c

      do 1000 ls=1,nints
c
      cwn = wn(ls)*detn(ls)
c
c    valores dos parametros do multiplicador
c
	  dhs=0.d00
        do i=1,npars
	    dhs = dhs + dls(i)*shgpn(2,i,ls)
        end do
c
c     geometria
c
        x1 =0.d00
        x2 =0.d00
        dx1=0.d00
        dx2=0.d00
c
        do i=1,nenlad
          x1 =x1 +xls(1,i)*shlb(2,i,ls)
          x2 =x2 +xls(2,i)*shlb(2,i,ls)
          dx1=dx1+xls(1,i)*shlb(1,i,ls)
          dx2=dx2+xls(2,i)*shlb(1,i,ls)
        end do
           dxx=dsqrt(dx1*dx1+dx2*dx2)
           xn1= sign*dx2/dxx
           xn2=-sign*dx1/dxx
c
c    valores exatos dos multiplicadores
c
c
      pix=pi*x1
      piy=pi*x2
      sx=dsin(pix)
      sy=dsin(piy)
      cx=dcos(pix)
      cy=dcos(piy)
c
      pi2=pi*pi
      co=1.d00/(pi2*2.d00+gama)
c
c
      sigma = dsqrt(gama)
      xc = x1*sigma*dcos(fi)
      yc = x2*sigma*dsin(fi)
c
c    valor exato do multiplicador
c
      dhse=gf2*co*cx*cy + gf1*co*sx*sy 
     &   +gf3*dcos(xc + yc)
c
c
      xmlte = xmlte + cwn*(dhse-dhs)**2
c
 1000 continue
 4000 continue
c
      xmlt = xmlt + xmlte
c
  50   continue
c
      divu = dsqrt(divu2)
c
c...... error in all domain
c

	el2p = dsqrt(edp)/co
	eh1p = dsqrt(edpx+edpy)/co
c
	edp0 = dlog10(el2p)
	gradp0= dlog10(eh1p)
c
c      
	edp = dlog10(dsqrt(edp))
	gradp= dlog10(dsqrt(edpx+edpy))
	xmlt = dlog10(dsqrt(xmlt))
c

       xel = dfloat(numel)/dfloat(nnods-1)**2
       xnp = neqc
       xnp = dlog10(xnp)
       xel = dlog10(xel)/2.d00

      write(iplt,2101) xel,xnp,edp,gradp,xmlt
c
      iplt0 = iplt+11
      write(iplt0,2101) xel,xnp,edp0,gradp0
c
      ipl = iplt+20	
      write(ipl,2102) beta0,el2p,eh1p
       return
c
2100  format(5x,i3,4x,e15.7,5x,e15.7)
2101  format(9(e15.5))
2102  format(4(f13.5))
       return
       end
c**** new **********************************************************************
      subroutine flninter(ien   ,x     ,xl   ,    
     &                 d     ,dl    ,mat   ,
     &                 c     ,ipar  ,
     &                 dlp   ,
     &                 det   ,
     &                 shl   , shg  ,wt    ,
     &                 ddis  ,detp  ,shlp  , 
     &                 shgp ,  
c
     &                 shln  ,shgn  ,
     &                 detn  ,shlb  ,shgb  ,
     &                 detpn ,shlpn ,shgpn ,
     &                 idside,xls   ,idlsd ,
     &                 grav  ,wn    ,
c                  
     &                 numel ,neesq ,nen   ,nsd   ,
     &                 nesd  ,nint  ,neg   ,nrowsh,
     &                 ned   ,nee   ,numnp ,ndof  ,
     &                 ncon  ,nencon,necon ,index ,
     &                 nints,iwrite ,iplt ,nenp   ,
     &                 nside,nnods ,nenlad,npars  , 
     &                 nmultpd, nmultpc,nodsp)
c------------------------------------------------------------------
c     Some variables used in this subroutine
c            j: degree of freedom (1,...ncon)
c         u(j): finite element solution
c        du(j): derivative of finite element solution
c        ue(j): exact solution
c       due(j): derivative of exact solution
c       el2(j): error in L2
c      epri(j): error in the seminorm of H1 (L2 of derivatives)
c     el2el(j): error in L2 in the element domain
c    epriel(j): error in the seminorm of H1 in the element domain
c----------------------------------------------------------------
c
c
c    Program to calculate and print the L2 and H1 seminorm of the error
c    for each degree of freedom in the finite element solution .
c    The trapezoidal rule is used for integration in each element .
c    The number of integration points is given by nints .
c
c    This version is only applicable to 2D.
c
      implicit real*8 (a-h,o-z)
c
c.... remove above card for single-precision operation
c
      character *1 tab
      dimension ien(nen,*),x(nsd,*),xl(nesd,*),d(ndof,*),dl(ned,*),
     &          mat(*),c(8,*),ipar(nodsp,*),
     &         ddis(ned,nenp,*)
      dimension u(6),dux(6),duy(6),ue(6),duex(6),duey(6)
      dimension shl(3,nen,*),shg(3,nen,*),det(*),wt(*)
      dimension el2(6),epri(6),eprix(6),epriy(6),el2el(6),
     &          eprxel(6),epryel(6)
	dimension detp(*),shlp(3,nenp,*),shgp(3,nenp,*),dlp(ned,*)
c
      dimension detn(*),shlb(3,nenlad,*),shgb(3,nenlad,*),
     &          detpn(*),shlpn(3,npars,*),shgpn(3,npars,*),
     &          idside(nside,*),xls(nesd,*),idlsd(*)
	dimension dls(12),grav(*),wn(*)
      dimension shln(3,nnods,*),shgn(3,nnods,*)
      dimension sxlhp(64,64),sxlhv(64,64),sxlhs(8),xlpn(2,64)
      dimension shlsd(8,8),xlps(2,8)
c
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
      common /iounit/ iin,iout,iecho,ioupp,itest1,itest2 
c
      pi=4.*datan(1.d00)
      dpi=2.d00*pi
c
      gf1=grav(1)
      gf2=grav(2)
      gf3=grav(3)  
	gama=c(3,1) 
	gamas=gama
	if(index.eq.3) gamas=-gama
	fi = pi/6.d00
         tab=char(9)
         xint=0.d00 
         yint=0.d00
c
      divu2=0.d00
	edp = 0.d00
	edpx= 0.d00
	edpy= 0.d00
	xmlt= 0.d00
c
c
       call shapesd(shlsd,nenlad,npars)
       call shapen(sxlhp,nen,nenp)
c
c.....  loop on elements
c      
      do 50 n=1,numel
c
       call local(ien(1,n),x,xl,nen,nsd,nesd)
       call local(ipar(1,n),d,dl,nodsp,ndof,ned)
c
      do i=1,nenp
      xlpn(1,i) = 0.d00
      xlpn(2,i) = 0.d00
      do j=1,nen
        xlpn(1,i) = xlpn(1,i) + sxlhp(j,i)*xl(1,j)
        xlpn(2,i) = xlpn(2,i) + sxlhp(j,i)*xl(2,j)
       end do
       end do
c       
       call elemdlp(xlpn,dlp,gama,fi,ned,nenp,index)	
c
c
	pmede=0.d00
	dpe = 0.d00
	dpex= 0.d00
	dpey= 0.d00
c
c.....  loop on integration points
c
c
c.....Triangles or quadrilaterals
c
c
      call shgqs(xl,detp,shlp,shgp,nint,n,neg,.true.,nenp,shl,nen)
c
c
        do 4040 l=1,nint
        ct=detp(l)*wt(l)
        call clear ( u, ncon )
        call clear (dux, ncon )
        call clear (duy, ncon )
        xint=0.d0
        yint=0.d0
         do 3030 i=1,nen
            xint=xint+shl(3,i,l)*xl(1,i)
            yint=yint+shl(3,i,l)*xl(2,i)
3030     continue
c
c
         call uexafx(xint,yint,ue,duex,duey,gama,fi,index)
c
c    pressao descontinua
c
       pe = 0.d00
	pex = 0.d00
	pey = 0.d00
	do i=1,nenp
	  pe = pe + shgp(3,i,l)*dlp(1,i)
	  pex= pex+ shgp(1,i,l)*dlp(1,i)
	  pey= pey+ shgp(2,i,l)*dlp(1,i)
	end do
	  dpe = dpe + ct*(pe - ue(3))**2
	  dpex= dpex + ct*(pex-duex(3))**2
	  dpey= dpey + ct*(pey-duey(3))**2
c
      if(iwrite.ne.0) then
         write(43,2424) n,l,pe,ue(3)
2424   format(2i10,6e15.5)
222    format(8e20.9)
       end if
c
4040  continue
c
	edp = edp + dpe
	edpx= edpx+ dpex
	edpy= edpy+ dpey
c
c......boundary terms - multiplier
c
	xmlte = 0.d00
      do 4000 ns=1,nside
c
c-----localiza os parametros do lado n
c
      do nn=1,npars
	  nld = (ns-1)*npars + nn
	  dls(nn) = dl(1,nld)
      end do
c
c-----localiza os nos do lado ns
c
c
      ns1=idside(ns,1)
      ns2=idside(ns,2)
      nl1=ien(ns1,n)
      nl2=ien(ns2,n)
c      
      if(nl2.gt.nl1) then 
	  sign = 1.d00
        do  nn=1,npars
           idlsd(nn)=idside(ns,nn)
	  end do
	else
	  sign = -1.d00
	  idlsd(1) = idside(ns,2)
	  idlsd(2) = idside(ns,1)
	  id3 = npars-2
	  if(id3.gt.0) then
	    do il=3,npars
	      idlsd(il) = idside(ns,npars+3-il)
	    end do 
	  end if
      end if 
c
c
      do 2000 nn=1,nenlad
      nl=idlsd(nn)
      xls(1,nn)=xl(1,nl)
      xls(2,nn)=xl(2,nl)
 2000 continue
c
      do nn=1,npars
      xlps(1,nn) = 0.d00
      xlps(2,nn) = 0.d00
      do jj=1,nenlad
        xlps(1,nn) = xlps(1,nn) + shlsd(jj,nn)*xls(1,jj)
        xlps(2,nn) = xlps(2,nn) + shlsd(jj,nn)*xls(2,jj)
      end do
      end do
c
       call elemdlp(xlps,dls,gama,fi,ned,npars,index)	
c
      call oneshgp(xls,detn,shlb,shln,shgn,
     &             nenlad,nnods,nints,nesd,ns,n,neg)      
c
      call oneshgp(xls,detpn,shlb,shlpn,shgpn,
     &             nenlad,npars,nints,nesd,ns,n,neg)      
c  
c
c.....compute boundary integral
c

      do 1000 ls=1,nints
c
      cwn = wn(ls)*detn(ls)
c
c    valores dos parametros do multiplicador
c
	  dhs=0.d00
        do i=1,npars
	    dhs = dhs + dls(i)*shgpn(2,i,ls)
        end do
c
c     geometria
c
        x1 =0.d00
        x2 =0.d00
        dx1=0.d00
        dx2=0.d00
c
        do i=1,nenlad
          x1 =x1 +xls(1,i)*shlb(2,i,ls)
          x2 =x2 +xls(2,i)*shlb(2,i,ls)
          dx1=dx1+xls(1,i)*shlb(1,i,ls)
          dx2=dx2+xls(2,i)*shlb(1,i,ls)
        end do
           dxx=dsqrt(dx1*dx1+dx2*dx2)
           xn1= sign*dx2/dxx
           xn2=-sign*dx1/dxx
c
c    valores exatos dos multiplicadores
c
c
      pix=pi*x1
      piy=pi*x2
      sx=dsin(pix)
      sy=dsin(piy)
      cx=dcos(pix)
      cy=dcos(piy)
c
      pi2=pi*pi
      co=1.d00/(pi2*2.d00+gama)
c
c
      sigma = dsqrt(gama)
      xc = x1*sigma*dcos(fi)
      yc = x2*sigma*dsin(fi)
c
c    valor exato do multiplicador
c
      dhse=gf2*co*cx*cy + gf1*co*sx*sy 
     &   +gf3*dcos(xc + yc)
c
c
      xmlte = xmlte + cwn*(dhse-dhs)**2
c
 1000 continue
 4000 continue
c
      xmlt = xmlt + xmlte
c
  50   continue
c
      divu = dsqrt(divu2)
c
c...... error in all domain
c
	edp = dlog10(dsqrt(edp))
	gradp= dlog10(dsqrt(edpx+edpy))

	 xmlt = dlog10(dsqrt(xmlt))
       xel = dfloat(numel)/dfloat(nnods-1)**2
       xnp = nmultpc
       xel = dlog10(xel)/2.d00
       xnp = dlog10(xnp)/2.d00

      write(iplt,2101) xel,xnp,edp,gradp,xmlt
       return
c
2100  format(5x,i3,4x,e15.7,5x,e15.7)
2101  format(9(e15.4))
       end
c**** new **********************************************************************
      subroutine flnorma(ien   ,x     ,xl   ,    
     &                 d     ,dl    ,mat   ,
     &                 c     ,ipar  ,
     &                 dlp   ,
     &                 det   ,
     &                 shl   , shg  ,wt    ,
     &                 ddis  ,detp  ,shlp  , 
     &                 shgp  ,
c                  
     &                 numel ,neesq ,nen   ,nsd   ,
     &                 nesd  ,nint  ,neg   ,nrowsh,
     &                 ned   ,nee   ,numnp ,ndof  ,
     &                 ncon  ,nencon,necon ,index ,
     &                 nints,iwrite ,iplt ,nenp   ,
     &                 nside,nnods ,nenlad,npars  ,
     &                nmultpd,nmultpc,nodsp )
c------------------------------------------------------------------
c     Some variables used in this subroutine
c            j: degree of freedom (1,...ncon)
c         u(j): finite element solution
c        du(j): derivative of finite element solution
c        ue(j): exact solution
c       due(j): derivative of exact solution
c       el2(j): error in L2
c      epri(j): error in the seminorm of H1 (L2 of derivatives)
c     el2el(j): error in L2 in the element domain
c    epriel(j): error in the seminorm of H1 in the element domain
c----------------------------------------------------------------
c
c
c    Program to calculate and print the L2 and H1 seminorm of the error
c    for each degree of freedom in the finite element solution .
c    The trapezoidal rule is used for integration in each element .
c    The number of integration points is given by nints .
c
c    This version is only applicable to 2D.
c
      implicit real*8 (a-h,o-z)
c
c.... remove above card for single-precision operation
c
     
      dimension ien(nen,*),x(nsd,*),xl(nesd,*),d(ndof,*),dl(ned,*),
     &          mat(*),c(8,*),ipar(nodsp,*),
     &          ddis(ned,nenp,*)
      dimension u(6),dux(6),duy(6),ue(6),duex(6),duey(6)
      dimension shl(3,nen,*),shg(3,nen,*),det(*),wt(*)
      dimension el2(6),epri(6),eprix(6),epriy(6),el2el(6),
     &          eprxel(6),epryel(6)
	dimension detp(*),shlp(3,nenp,*),shgp(3,nenp,*),dlp(ned,*)
c
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
      common /iounit/ iin,iout,iecho,ioupp,itest1,itest2 
c
         pi=4.*datan(1.d00)
         dpi=2.d00*pi
         gama=c(3,1)
         gamas=gama
         if(index.eq.3) gamas=-gama
	   fi = pi/6.d00
         xint=0.d00 
         yint=0.d00
c
	edp = 0.d00
	edpx= 0.d00
	edpy= 0.d00
	xmlt= 0.d00
c
c.....  loop on elements
c
      do 50 n=1,numel
c
      call local(ien(1,n),x,xl,nen,nsd,nesd)
      call local(ipar(1,n),d,dl,nodsp,ndof,ned)
c
      do i=1,nenp
	 dlp(1,i) = ddis(1,i,n)
	end do
c
c
	dpe = 0.d00
	dpex= 0.d00
	dpey= 0.d00
c
c.....  loop on integration points
c
c
c.....Triangles or quadrilaterals
c
c
      call shgqs(xl,detp,shlp,shgp,nint,n,neg,.true.,nenp,shl,nen)
c
c
        do 4040 l=1,nint
        ct=detp(l)*wt(l)
        call clear ( u, ncon )
        call clear (dux, ncon )
        call clear (duy, ncon )
        xint=0.d0
        yint=0.d0
         do 3030 i=1,nen
            xint=xint+shl(3,i,l)*xl(1,i)
            yint=yint+shl(3,i,l)*xl(2,i)
3030     continue
c
         call uexafx(xint,yint,ue,duex,duey,gama,fi,index)
c
c    pressao descontinua
c
       pe = 0.d00
	pex = 0.d00
	pey = 0.d00
	do i=1,nenp
	  pe = pe + shgp(3,i,l)*dlp(1,i)
	  pex= pex+ shgp(1,i,l)*dlp(1,i)
	  pey= pey+ shgp(2,i,l)*dlp(1,i)
	end do
	  dpe = dpe + ct*(pe - ue(3))**2
	  dpex= dpex + ct*(pex-duex(3))**2
	  dpey= dpey + ct*(pey-duey(3))**2
c
      if(iwrite.ne.0) then
         write(43,2424) n,l,pe,ue(3)
2424   format(2i10,6e15.5)
222    format(8e20.9)
       end if
c
4040  continue
c
c
	edp = edp + dpe
	edpx= edpx+ dpex
	edpy= edpy+ dpey
c
  50   continue
c
      divu = dsqrt(divu2)

c
c...... error in all domain
c
	 edp  = dlog10(dsqrt(edp))
	 gradp= dlog10(dsqrt(edpx+edpy))

       xel = dfloat(numel)/dfloat(nnods-1)**2
       xnp = nmultpc
       xel = dlog10(xel)/2.d00
       xnp = dlog10(xnp)/2.d00

      write(iplt,2101) xel,xnp,edp,gradp
       return
c
2100  format(5x,i3,4x,e15.7,5x,e15.7)
2101  format(9(e15.4))
       end
c**** new ********************************************************************** 
      subroutine shapen(shl,nen,nenc)
c
c
      implicit real*8 (a-h,o-z)
c
c.... remove above card for single precision operation
c
      dimension raone(8),xaone(8)
      dimension shlone(8,8),inod(8,8),inodc(8,8)
      dimension shl(64,64)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
c
      if (nen.eq.1) xaone(1) = zero
      if (nenc.eq.1) raone(1) = zero
c
c
      if (nen.eq.4) then
         xaone(1) = -one
         xaone(2) =  one
	nenx=2
	neny=2
      endif
c
      if (nenc.eq.4) then
         raone(1) = -one
         raone(2) =  one
	nenrx=2
	nenry=2
      endif
c
c 
      if(nen.eq.9) then
         xaone(1)= -one
         xaone(2)= one
         xaone(3)= zero
	nenx=3
	neny=3
      endif
c 
      if(nenc.eq.9) then
         raone(1)= -one
         raone(2)= one
         raone(3)= zero
	nenrx=3
	nenry=3
      endif
c
c
      if (nen.eq.16) then
         xaone(1) = -one
         xaone(2) = one
         xaone(3) = -.333333333333333
         xaone(4) =  .333333333333333
	nenx=4
	neny=4
         endif
c
      if (nenc.eq.16) then
         raone(1) = -one
         raone(2) = one
         raone(3) = -.333333333333333
         raone(4) =  .333333333333333
	nenrx=4
	nenry=4
         endif
c
c
       if(nen.eq.25) then
         xaone(1)= -one 
         xaone(2)=  one 
         xaone(3)= -pt5
         xaone(4)= zero
         xaone(5)= pt5
	nenx=5
	neny=5
       endif
c
       if(nenc.eq.25) then
         raone(1)= -one 
         raone(2)=  one 
         raone(3)= -pt5
         raone(4)= zero
         raone(5)= pt5
	nenrx=5
	nenry=5
       endif
c
c
        if(nen.eq.36) then
         xaone(1) = -one
         xaone(2) =  one
         xaone(3) = -.600000000000000
         xaone(4) = -.200000000000000
         xaone(5) =  .200000000000000
         xaone(6) =  .600000000000000
	nenx=6
	neny=6
        endif
c
        if(nenc.eq.36) then
         raone(1) = -one
         raone(2) =  one
         raone(3) = -.600000000000000
         raone(4) = -.200000000000000
         raone(5) =  .200000000000000
         raone(6) =  .600000000000000
	nenrx=6
	nenry=6
        endif
c
c
        if(nen.eq.49) then
         xaone(1) = -one
         xaone(2) =  one
         xaone(3) = -.666666666666666
         xaone(4) = -.333333333333333
         xaone(5) = zero
         xaone(6) =  .333333333333333
         xaone(7) =  .666666666666666
	nenx=7
	neny=7
        endif
c
        if(nenc.eq.49) then
         raone(1) = -one
         raone(2) =  one
         raone(3) = -.666666666666666
         raone(4) = -.333333333333333
         raone(5) = zero
         raone(6) =  .333333333333333
         raone(7) =  .666666666666666
	nenrx=7
	nenry=7
        endif
c
c
        if(nen.eq.64) then
         xaone(1) = -one
         xaone(2) =  one
         xaone(3) = -0.71428571428571
         xaone(4) = -0.42857142857143
         xaone(5) = -0.14285714285714
         xaone(6) = 0.14285714285714
         xaone(7) = 0.42857142857143
         xaone(8) = 0.71428571428571
	nenx=8
	neny=8
        endif
c
        if(nenc.eq.64) then
         raone(1) = -one
         raone(2) =  one
         raone(3) = -0.71428571428571
         raone(4) = -0.42857142857143
         raone(5) = -0.14285714285714
         raone(6) = 0.14285714285714
         raone(7) = 0.42857142857143
         raone(8) = 0.71428571428571
	nenrx=8
	nenry=8
        endif
c
c
	if(nen.eq.1) inod(1,1) = 1
	if(nen.eq.4) then
	  inod(1,1) = 1
	  inod(2,1) = 2
	  inod(1,2) = 4
	  inod(2,2) = 3
	end if
c
	if(nen.eq.9) then
	  inod(1,1) = 1
	  inod(2,1) = 2
	  inod(3,1) = 5
c
	  inod(1,2) = 4
	  inod(2,2) = 3
	  inod(3,2) = 7
c
	  inod(1,3) = 8
	  inod(2,3) = 6
	  inod(3,3) = 9
	end if
c
c
c
	if(nen.eq.16) then
	  inod(1,1) = 1
	  inod(2,1) = 2
	  inod(3,1) = 5
	  inod(4,1) = 6
c
	  inod(1,2) = 4
	  inod(2,2) = 3
	  inod(3,2) = 10
	  inod(4,2) = 9
c
	  inod(1,3) = 12
	  inod(2,3) = 7
	  inod(3,3) = 13
	  inod(4,3) = 14
c
	  inod(1,4) = 11
	  inod(2,4) = 8
	  inod(3,4) = 16
	  inod(4,4) = 15
	end if
c
c
	if(nen.eq.25) then
	  inod(1,1) = 1
	  inod(2,1) = 2
	  inod(3,1) = 5
	  inod(4,1) = 6
	  inod(5,1) = 7
c
	  inod(1,2) = 4
	  inod(2,2) = 3
	  inod(3,2) = 13
	  inod(4,2) = 12
	  inod(5,2) = 11
c
	  inod(1,3) = 16
	  inod(2,3) = 8
	  inod(3,3) = 17
	  inod(4,3) = 18
	  inod(5,3) = 19
c
	  inod(1,4) = 15
	  inod(2,4) = 9
	  inod(3,4) = 24
	  inod(4,4) = 25
	  inod(5,4) = 20
c	  
	  inod(1,5) = 14
	  inod(2,5) = 10
	  inod(3,5) = 23
	  inod(4,5) = 22
	  inod(5,5) = 21
	end if
c
c
	if(nen.eq.36) then
	  inod(1,1) = 1
	  inod(2,1) = 2
	  inod(3,1) = 5
	  inod(4,1) = 6
	  inod(5,1) = 7
	  inod(6,1) = 8
c
	  inod(1,2) = 4
	  inod(2,2) = 3
	  inod(3,2) = 16
	  inod(4,2) = 15
	  inod(5,2) = 14
	  inod(6,2) = 13
c
	  inod(1,3) = 20
	  inod(2,3) = 9
	  inod(3,3) = 21
	  inod(4,3) = 22
	  inod(5,3) = 23
	  inod(6,3) = 24
c
	  inod(1,4) = 19
	  inod(2,4) = 10
	  inod(3,4) = 32
	  inod(4,4) = 33
	  inod(5,4) = 34
	  inod(6,4) = 25
c	  
	  inod(1,5) = 18
	  inod(2,5) = 11
	  inod(3,5) = 31
	  inod(4,5) = 36
	  inod(5,5) = 35
	  inod(6,5) = 26
c
	  inod(1,6) = 17
	  inod(2,6) = 12
	  inod(3,6) = 30
	  inod(4,6) = 29
	  inod(5,6) = 28
	  inod(6,6) = 27
c	  
	end if
c
c
	if(nen.eq.49) then
	  inod(1,1) = 1
	  inod(2,1) = 2
	  inod(3,1) = 5
	  inod(4,1) = 6
	  inod(5,1) = 7
	  inod(6,1) = 8
	  inod(7,1) = 9
c
	  inod(1,2) = 4
	  inod(2,2) = 3
	  inod(3,2) = 19
	  inod(4,2) = 18
	  inod(5,2) = 17
	  inod(6,2) = 16
	  inod(7,2) = 15
c
	  inod(1,3) = 24
	  inod(2,3) = 10
	  inod(3,3) = 25
	  inod(4,3) = 26
	  inod(5,3) = 27
	  inod(6,3) = 28
	  inod(7,3) = 29
c
	  inod(1,4) = 23
	  inod(2,4) = 11
	  inod(3,4) = 40
	  inod(4,4) = 41
	  inod(5,4) = 42
	  inod(6,4) = 43
	  inod(7,4) = 30
c	  
	  inod(1,5) = 22
	  inod(2,5) = 12
	  inod(3,5) = 39
	  inod(4,5) = 48
	  inod(5,5) = 49
	  inod(6,5) = 44
	  inod(7,5) = 31
c
	  inod(1,6) = 21
	  inod(2,6) = 13
	  inod(3,6) = 38
	  inod(4,6) = 47
	  inod(5,6) = 46
	  inod(6,6) = 45
	  inod(7,6) = 32
c
	  inod(1,7) = 20
	  inod(2,7) = 14
	  inod(3,7) = 37
	  inod(4,7) = 36
	  inod(5,7) = 35
	  inod(6,7) = 34
	  inod(7,7) = 33
c	  
	end if
c
c
	if(nen.eq.64) then
	  inod(1,1) = 1
	  inod(2,1) = 2
	  inod(3,1) = 5
	  inod(4,1) = 6
	  inod(5,1) = 7
	  inod(6,1) = 8
	  inod(7,1) = 9
	  inod(8,1) = 10
c
	  inod(1,2) = 4
	  inod(2,2) = 3
	  inod(3,2) = 22
	  inod(4,2) = 21
	  inod(5,2) = 20
	  inod(6,2) = 19
	  inod(7,2) = 18
	  inod(8,2) = 17
c
	  inod(1,3) = 28
	  inod(2,3) = 11
	  inod(3,3) = 29
	  inod(4,3) = 30
	  inod(5,3) = 31
	  inod(6,3) = 32
	  inod(7,3) = 33
	  inod(8,3) = 34
c
	  inod(1,4) = 27
	  inod(2,4) = 12
	  inod(3,4) = 48
	  inod(4,4) = 49
	  inod(5,4) = 50
	  inod(6,4) = 51
	  inod(7,4) = 52
	  inod(8,4) = 35
c	  
	  inod(1,5) = 26
	  inod(2,5) = 13
	  inod(3,5) = 47
	  inod(4,5) = 60
	  inod(5,5) = 61
	  inod(6,5) = 62
	  inod(7,5) = 53
	  inod(8,5) = 36
c
	  inod(1,6) = 25
	  inod(2,6) = 14
	  inod(3,6) = 46
	  inod(4,6) = 59
	  inod(5,6) = 64
	  inod(6,6) = 63
	  inod(7,6) = 54
	  inod(8,6) = 37
c
	  inod(1,7) = 24
	  inod(2,7) = 15
	  inod(3,7) = 45
	  inod(4,7) = 58
	  inod(5,7) = 57
	  inod(6,7) = 56
	  inod(7,7) = 55
	  inod(8,7) = 38
c
	  inod(1,8) = 23
	  inod(2,8) = 16
	  inod(3,8) = 44
	  inod(4,8) = 43
	  inod(5,8) = 42
	  inod(6,8) = 41
	  inod(7,8) = 40
	  inod(8,8) = 39
c	  
	end if
c
c
c
	if(nenc.eq.1) inodc(1,1) = 1
	if(nenc.eq.4) then
	  inodc(1,1) = 1
	  inodc(2,1) = 2
	  inodc(1,2) = 4
	  inodc(2,2) = 3
	end if
c
	if(nenc.eq.9) then
	  inodc(1,1) = 1
	  inodc(2,1) = 2
	  inodc(3,1) = 5
c
	  inodc(1,2) = 4
	  inodc(2,2) = 3
	  inodc(3,2) = 7
c
	  inodc(1,3) = 8
	  inodc(2,3) = 6
	  inodc(3,3) = 9
	end if
c
c
c
	if(nenc.eq.16) then
	  inodc(1,1) = 1
	  inodc(2,1) = 2
	  inodc(3,1) = 5
	  inodc(4,1) = 6
c
	  inodc(1,2) = 4
	  inodc(2,2) = 3
	  inodc(3,2) = 10
	  inodc(4,2) = 9
c
	  inodc(1,3) = 12
	  inodc(2,3) = 7
	  inodc(3,3) = 13
	  inodc(4,3) = 14
c
	  inodc(1,4) = 11
	  inodc(2,4) = 8
	  inodc(3,4) = 16
	  inodc(4,4) = 15
	end if
c
c
	if(nenc.eq.25) then
	  inodc(1,1) = 1
	  inodc(2,1) = 2
	  inodc(3,1) = 5
	  inodc(4,1) = 6
	  inodc(5,1) = 7
c
	  inodc(1,2) = 4
	  inodc(2,2) = 3
	  inodc(3,2) = 13
	  inodc(4,2) = 12
	  inodc(5,2) = 11
c
	  inodc(1,3) = 16
	  inodc(2,3) = 8
	  inodc(3,3) = 17
	  inodc(4,3) = 18
	  inodc(5,3) = 19
c
	  inodc(1,4) = 15
	  inodc(2,4) = 9
	  inodc(3,4) = 24
	  inodc(4,4) = 25
	  inodc(5,4) = 20
c	  
	  inodc(1,5) = 14
	  inodc(2,5) = 10
	  inodc(3,5) = 23
	  inodc(4,5) = 22
	  inodc(5,5) = 21
	end if
c
c
	if(nenc.eq.36) then
	  inodc(1,1) = 1
	  inodc(2,1) = 2
	  inodc(3,1) = 5
	  inodc(4,1) = 6
	  inodc(5,1) = 7
	  inodc(6,1) = 8
c
	  inodc(1,2) = 4
	  inodc(2,2) = 3
	  inodc(3,2) = 16
	  inodc(4,2) = 15
	  inodc(5,2) = 14
	  inodc(6,2) = 13
c
	  inodc(1,3) = 20
	  inodc(2,3) = 9
	  inodc(3,3) = 21
	  inodc(4,3) = 22
	  inodc(5,3) = 23
	  inodc(6,3) = 24
c
	  inodc(1,4) = 19
	  inodc(2,4) = 10
	  inodc(3,4) = 32
	  inodc(4,4) = 33
	  inodc(5,4) = 34
	  inodc(6,4) = 25
c	  
	  inodc(1,5) = 18
	  inodc(2,5) = 11
	  inodc(3,5) = 31
	  inodc(4,5) = 36
	  inodc(5,5) = 35
	  inodc(6,5) = 26
c
	  inodc(1,6) = 17
	  inodc(2,6) = 12
	  inodc(3,6) = 30
	  inodc(4,6) = 29
	  inodc(5,6) = 28
	  inodc(6,6) = 27
c	  
	end if
c
c
	if(nenc.eq.49) then
	  inodc(1,1) = 1
	  inodc(2,1) = 2
	  inodc(3,1) = 5
	  inodc(4,1) = 6
	  inodc(5,1) = 7
	  inodc(6,1) = 8
	  inodc(7,1) = 9
c
	  inodc(1,2) = 4
	  inodc(2,2) = 3
	  inodc(3,2) = 19
	  inodc(4,2) = 18
	  inodc(5,2) = 17
	  inodc(6,2) = 16
	  inodc(7,2) = 15
c
	  inodc(1,3) = 24
	  inodc(2,3) = 10
	  inodc(3,3) = 25
	  inodc(4,3) = 26
	  inodc(5,3) = 27
	  inodc(6,3) = 28
	  inodc(7,3) = 29
c
	  inodc(1,4) = 23
	  inodc(2,4) = 11
	  inodc(3,4) = 40
	  inodc(4,4) = 41
	  inodc(5,4) = 42
	  inodc(6,4) = 43
	  inodc(7,4) = 30
c	  
	  inodc(1,5) = 22
	  inodc(2,5) = 12
	  inodc(3,5) = 39
	  inodc(4,5) = 48
	  inodc(5,5) = 49
	  inodc(6,5) = 44
	  inodc(7,5) = 31
c
	  inodc(1,6) = 21
	  inodc(2,6) = 13
	  inodc(3,6) = 38
	  inodc(4,6) = 47
	  inodc(5,6) = 46
	  inodc(6,6) = 45
	  inodc(7,6) = 32
c
	  inodc(1,7) = 20
	  inodc(2,7) = 14
	  inodc(3,7) = 37
	  inodc(4,7) = 36
	  inodc(5,7) = 35
	  inodc(6,7) = 34
	  inodc(7,7) = 33
c	  
	end if
c
c
	if(nenc.eq.64) then
	  inodc(1,1) = 1
	  inodc(2,1) = 2
	  inodc(3,1) = 5
	  inodc(4,1) = 6
	  inodc(5,1) = 7
	  inodc(6,1) = 8
	  inodc(7,1) = 9
	  inodc(8,1) = 10
c
	  inodc(1,2) = 4
	  inodc(2,2) = 3
	  inodc(3,2) = 22
	  inodc(4,2) = 21
	  inodc(5,2) = 20
	  inodc(6,2) = 19
	  inodc(7,2) = 18
	  inodc(8,2) = 17
c
	  inodc(1,3) = 28
	  inodc(2,3) = 11
	  inodc(3,3) = 29
	  inodc(4,3) = 30
	  inodc(5,3) = 31
	  inodc(6,3) = 32
	  inodc(7,3) = 33
	  inodc(8,3) = 34
c
	  inodc(1,4) = 27
	  inodc(2,4) = 12
	  inodc(3,4) = 48
	  inodc(4,4) = 49
	  inodc(5,4) = 50
	  inodc(6,4) = 51
	  inodc(7,4) = 52
	  inodc(8,4) = 35
c	  
	  inodc(1,5) = 26
	  inodc(2,5) = 13
	  inodc(3,5) = 47
	  inodc(4,5) = 60
	  inodc(5,5) = 61
	  inodc(6,5) = 62
	  inodc(7,5) = 53
	  inodc(8,5) = 36
c
	  inodc(1,6) = 25
	  inodc(2,6) = 14
	  inodc(3,6) = 46
	  inodc(4,6) = 59
	  inodc(5,6) = 64
	  inodc(6,6) = 63
	  inodc(7,6) = 54
	  inodc(8,6) = 37
c
	  inodc(1,7) = 24
	  inodc(2,7) = 15
	  inodc(3,7) = 45
	  inodc(4,7) = 58
	  inodc(5,7) = 57
	  inodc(6,7) = 56
	  inodc(7,7) = 55
	  inodc(8,7) = 38
c
	  inodc(1,8) = 23
	  inodc(2,8) = 16
	  inodc(3,8) = 44
	  inodc(4,8) = 43
	  inodc(5,8) = 42
	  inodc(6,8) = 41
	  inodc(7,8) = 40
	  inodc(8,8) = 39
c	  
	end if
c
      do 100 l = 1, nenrx
         r = raone(l)
c
        if(nenx.eq.1) then
         shlone(1,l) = zero
         go to 100
        endif
c
        do 50 i = 1, nenx
         aa = one
         bb = one
         aax = zero
         do 40 j =1, nenx
          daj = one
          if (i .ne. j)then
          aa = aa * ( r - xaone(j))
          bb = bb * ( xaone(i) - xaone(j))
          endif
   40    continue
        shlone(i,l) = aa/bb
   50  continue
c
  100  continue      
c
      do ly=1,nenrx
      do lx=1,nenry
        l = inodc(lx,ly)
	 do iy=1,neny
	 do ix=1,nenx
         j = inod(ix,iy)
	   shl(j,l) = shlone(ix,lx)*shlone(iy,ly)
	 end do
	 end do
	end do
	end do
c
      return
      end
c**** new **********************************************************************
c**** new ********************************************************************** 
      subroutine shapesd(shl,nen,nenc)
c
c
      implicit real*8 (a-h,o-z)
c
c.... remove above card for single precision operation
c
      dimension raone(8),xaone(8)
      dimension shl(8,8)
      common /consts/ zero,pt1667,pt25,pt5,one,two,three,four,five,six
c
      if (nen.eq.1) xaone(1) = zero
      if (nenc.eq.1) raone(1) = zero
c
c
      if (nen.eq.2) then
         xaone(1) = -one
         xaone(2) =  one
      endif
c
      if (nenc.eq.2) then
         raone(1) = -one
         raone(2) =  one
      endif
c
c 
      if(nen.eq.3) then
         xaone(1)= -one
         xaone(2)= one
         xaone(3)= zero
       endif
c 
      if(nenc.eq.3) then
         raone(1)= -one
         raone(2)= one
         raone(3)= zero
       endif
c
c
      if (nen.eq.4) then
         xaone(1) = -one
         xaone(2) = one
         xaone(3) = -.333333333333333
         xaone(4) =  .333333333333333
       endif
c
      if (nenc.eq.4) then
         raone(1) = -one
         raone(2) = one
         raone(3) = -.333333333333333
         raone(4) =  .333333333333333
       endif
c
c
       if(nen.eq.5) then
         xaone(1)= -one 
         xaone(2)=  one 
         xaone(3)= -pt5
         xaone(4)= zero
         xaone(5)= pt5
        endif
c
       if(nenc.eq.5) then
         raone(1)= -one 
         raone(2)=  one 
         raone(3)= -pt5
         raone(4)= zero
         raone(5)= pt5
       endif
c
c
        if(nen.eq.6) then
         xaone(1) = -one
         xaone(2) =  one
         xaone(3) = -.600000000000000
         xaone(4) = -.200000000000000
         xaone(5) =  .200000000000000
         xaone(6) =  .600000000000000
        endif
c
        if(nenc.eq.6) then
         raone(1) = -one
         raone(2) =  one
         raone(3) = -.600000000000000
         raone(4) = -.200000000000000
         raone(5) =  .200000000000000
         raone(6) =  .600000000000000
        endif
c
c
        if(nen.eq.7) then
         xaone(1) = -one
         xaone(2) =  one
         xaone(3) = -.666666666666666
         xaone(4) = -.333333333333333
         xaone(5) = zero
         xaone(6) =  .333333333333333
         xaone(7) =  .666666666666666
        endif
c
        if(nenc.eq.7) then
         raone(1) = -one
         raone(2) =  one
         raone(3) = -.666666666666666
         raone(4) = -.333333333333333
         raone(5) = zero
         raone(6) =  .333333333333333
         raone(7) =  .666666666666666
        endif
c
c
        if(nen.eq.8) then
         xaone(1) = -one
         xaone(2) =  one
         xaone(3) = -0.71428571428571
         xaone(4) = -0.42857142857143
         xaone(5) = -0.14285714285714
         xaone(6) = 0.14285714285714
         xaone(7) = 0.42857142857143
         xaone(8) = 0.71428571428571
        endif
c
        if(nenc.eq.8) then
         raone(1) = -one
         raone(2) =  one
         raone(3) = -0.71428571428571
         raone(4) = -0.42857142857143
         raone(5) = -0.14285714285714
         raone(6) = 0.14285714285714
         raone(7) = 0.42857142857143
         raone(8) = 0.71428571428571
        endif
c
c
c
      do 100 l = 1, nenc
         r = raone(l)
c
        if(nen.eq.1) then
         shl(1,l) = zero
         go to 100
        endif
c
        do 50 i = 1, nen
         aa = one
         bb = one
         aax = zero
         do 40 j =1, nen
          daj = one
          if (i .ne. j)then
          aa = aa * ( r - xaone(j))
          bb = bb * ( xaone(i) - xaone(j))
          endif
   40    continue
        shl(i,l) = aa/bb
   50  continue
c
  100  continue      
c
      return
      end
c**** new **********************************************************************
c******* nova subrotina ***************************************************************
c                                              
      subroutine shltpbk(shl,nen,nside,nnods,nints)
c
c     Objetivo: calcular pesos, funcoes de interpolacao e derivadas locais
c               para elementos triangulares
c
c----------------------------------------------------------------------
c
c 
c.... program to calculate integration-rule weights, shape functions 
c        and local derivatives for a triangular element 
c 
c        c1, c2, c3 = local element coordinates ("l1", "l2", "l3".) 
c        shl(j,i,l) = local ("j") derivative of shape function 
c        shl(3,i,l) = local  shape function 
c              w(l) = integration-rule weight 
c                 i = local node number 
c                 l = integration point number 
c              nint = number of integration points, eq. 1 or 4 
c 
      implicit real*8 (a-h,o-z) 
c 
c.... remove above card for single precision operation 
c 
      dimension shl(3,nen,*),cl1(18),cl2(18),cl3(18),csi(6),ra(6)
      data   zero,pt1667,pt25,pt5 
     &      /0.0d0,0.1666666666666667d0,0.25d0,0.5d0/, 
     &       one,two,three,four,five,six 
     &      /1.0d0,2.0d0,3.0d0,4.0d0,5.0d0,6.0d0/ 
      data r1/0.33333333333333333333d00/,w1/1.d00/, 
     &     r2/0.5d00                   /,w2/0.3333333333333333333d00/, 
     &     r3a/0.3333333333333333333d00/,w3a/-0.5625d00/, 
     &     r3b1/0.6d00                 /,w3b/0.520833333333333333d00/, 
     &     r3b2/0.2d00                 /, 
     &     r7a/0.3333333333333333333d00/,w7a/0.225d00/, 
     &     r7b/0.0597158717d00         /,w7b/0.1323941527d00/, 
     &     r7c/0.4701420641d00         /, 
     &     r7d/0.7974269853d00         /,w7d/0.1259391805d00/, 
     &     r7e/0.1012865073d00         /, 
     &     r13a/0.3333333333333333333d00/,w13a/-0.1495700444676d0/, 
     &     r13b/0.479308067841923d0         /,w13b/0.175615257433204d0/, 
     &     r13c/0.260345966079038d0         /,w13c/0.053347235608839d0/ 
     &     r13d/0.869739794195568d0         /,w13d/0.077113760890257d0/, 
     &     r13e/0.065130102902216d0        /, 
     &     r13f/0.638444188568809d0         /, 
     &     r13g/0.312865496004875d0         / 
     &     r13h/0.048690315425316d0         / 
c 
      if (nints.eq.1) then
         ra(1) = zero
      endif
c
      if (nints.eq.2) then
         ra(1)=-.577350269189626
         ra(2)= .577350269189625
      endif
c
      if (nints.eq.3) then
         ra(1)=-.774596669241483
         ra(2)= zero
         ra(3)= .774596669241483
      endif
c
      if (nints.eq.4) then
         ra(1)=-.861136311594053
         ra(2)=-.339981043584856
         ra(3)= .339981043584856
         ra(4)= .861136311594053
      endif
c
       if(nints.eq.5) then
        ra(1)=-.906179845938664
        ra(2)=-.538469310105683
        ra(3)= zero
        ra(4)= .538469310105683
        ra(5)= .906179845938664
       endif
c
       if(nints.eq.6) then
         ra(1)=-.932469514203152
         ra(2)=-.661209386466265
         ra(3)=-.238619186083197
         ra(4)= .238619186083197
         ra(5)= .661209386466365
         ra(6)= .932469514203152
        endif
c
c
      lb = 0
      do 100 ns=1,3
      do 200 ls=1,nints
c
      lb = lb + 1
c
      if(ns.eq.1) then
        cl3(ls) = 0.d00
        cl1(ls) = (1.d00 - ra(ls))/two
        cl2(ls) = (1.d00 + ra(ls))/two
      end if 
c
      if(ns.eq.2) then
         cl1(ls) = 0.d00
         cl2(ls) = (1.d00 - ra(ls))/two
         cl3(ls) = (1.d00 + ra(ls))/two
      end if 
c
      if(ns.eq.3) then
         cl2(ls) = 0.d00
         cl3(ls) = (1.d00 - ra(ls))/two
         cl1(ls) = (1.d00 + ra(ls))/two
      end if 
c 
            c1 = cl1(ls) 
            c2 = cl2(ls) 
            c3 = cl3(ls) 
            shl(1,1,lb)= one 
            shl(2,1,lb)= zero 
            shl(3,1,lb)= c1 
            shl(1,2,lb)= zero 
            shl(2,2,lb)= one 
            shl(3,2,lb)= c2 
            shl(1,3,lb)=-one 
            shl(2,3,lb)=-one 
            shl(3,3,lb)= c3 
            if(nen.eq.6) then 
                  shl(1,4,lb)= four * c2 
                  shl(2,4,lb)= four * c1 
                  shl(3,4,lb)= four * c1 * c2 
                  shl(1,5,lb)=-four * c2 
                  shl(2,5,lb)= four * (c3 - c2) 
                  shl(3,5,lb)= four * c2 * c3 
                  shl(1,6,lb)= four * (c3 - c1) 
                  shl(2,6,lb)=-four * c1 
                  shl(3,6,lb)= four * c3 * c1 
c 
                  do 70 i=1,3 
                        shl(i,1,lb)=shl(i,1,lb) 
     &                  -pt5*(shl(i,4,lb)+shl(i,6,lb)) 
                        shl(i,2,lb)=shl(i,2,lb) 
     &                  -pt5*(shl(i,4,lb)+shl(i,5,lb)) 
                        shl(i,3,lb)=shl(i,3,lb) 
     &                  -pt5*(shl(i,5,lb)+shl(i,6,lb)) 
   70             continue 
            end if 
c
            if (nen.eq.10) then 
 
                  shl(1,1,lb)= pt5*((three*c1-one)*(three*c1-two) 
     &            +c1*three*(three*c1-two)+c1*three*(three*c1-one)) 
                  shl(2,1,lb)= zero 
                  shl(3,1,lb)= pt5*c1*(three*c1-one)*(three*c1-two) 
 
                  shl(1,2,lb)=zero 
                  shl(2,2,lb)=  pt5*((three*c2-one)*(three*c2-two) 
     &            +c2*three*(three*c2-two)+c2*three*(three*c2-one)) 
                  shl(3,2,lb)=  pt5*c2*(three*c2-one)*(three*c2-two) 
 
                  shl(1,3,lb)= -pt5*((three*c3-one)*(three*c3-two) 
     &            +c3*three*(three*c3-two)+c3*three*(three*c3-one)) 
                  shl(2,3,lb)=-pt5*((three*c3-one)*(three*c3-two) 
     &            +c3*three*(three*c3-two)+c3*three*(three*c3-one)) 
                  shl(3,3,lb)=  pt5*c3*(three*c3-one)*(three*c3-two) 
 
 
                  shl(1,4,lb) = (9.0d0/two)*(c2*(three*c1-one) 
     &            +c1*c2*three) 
                  shl(2,4,lb) = (9.0d0/two)*c1*(three*c1-one) 
                  shl(3,4,lb) = (9.0d0/two)*c1*c2*(three*c1-one) 
 
                  shl(1,5,lb) = (9.0d0/two)*c2*(three*c2-one) 
                  shl(2,5,lb) = (9.0d0/two)*(c1*(three*c2-one) 
     &            +c1*c2*three) 
                  shl(3,5,lb) = (9.0d0/two)*c1*c2*(three*c2-one) 
 
 
                 shl(1,6,lb) =- (9.0d0/two)*c2*(three*c2-one) 
                 shl(2,6,lb) =  (9.0d0/two)*(c3*(three*c2-one) 
     &           +c2*c3*three-c2*(three*c2-one)) 
                 shl(3,6,lb) = (9.0d0/two)*c3*c2*(three*c2-one) 
 
 
                  shl(1,7,lb) = (9.0d0/two)*(-c2*(three*c3-one) 
     &            -c2*c3*three) 
                  shl(2,7,lb) =  (9.0d0/two)*(c3*(three*c3-one) 
     &            -c2*c3*three -c2*(three*c3-one)) 
                  shl(3,7,lb) = (9.0d0/two)*c2*c3*(three*c3-one) 
 
                  shl(1,9,lb)= (9.0d0/two)*(-c1*(three*c1-one) 
     &            +c1*c3*three +c3*(three*c1-one)) 
                  shl(2,9,lb) =  (9.0d0/two)*(-c1*(three*c1-one)) 
                  shl(3,9,lb) = (9.0d0/two)*c3*c1*(three*c1-one) 
 
 
                  shl(1,8,lb) = (9.0d0/two)*(-c1*(three*c3-one) 
     &            -c1*c3*three +c3*(three*c3-one)) 
                  shl(2,8,lb) =  (9.0d0/two)*( 
     &            -c1*c3*three -c1*(three*c3-one)) 
                  shl(3,8,lb) = (9.0d0/two)*c1*c3*(three*c3-one) 
 
                 shl(1,10,lb) = 27.0d0*(c2*c3-c1*c2) 
                 shl(2,10,lb) = 27.0d0*(c1*c3-c1*c2) 
                 shl(3,10,lb) = 27.0d0*c1*c2*c3 
          endif 
  200     continue 
  100     continue
c
c 
           return 
           end
c**** new **********************************************************************
