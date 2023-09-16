      program SAWchains
c cccccccccccccccccccccccccccccccccccccccc
c     This program was modified from DNAslit.f that was in the DNAelectrophoresis
C     Try to turn this program to monitor Rg2, Rend, for a SAW chain only
c     Would allow to put bead-bead attractions (so we can model collapse if we want)
c     First attempted modifications: April 12, 2022
C***************************************************
c     version I, May 26, 2006
c     Authros: Yongmei Wang
c     version II: July 26, 2006: add monitoring of n_ads
c     also relax the systems first and then apply electrif field
c ccccccccccccccccccccccccccccccccccc
c     This is a dynamic MOnte Carlo simulation code, intended for the
c     study of DNA electrophorsis in a flat slab (a slit). This surface
c     can be easily removed, then you recover the DNA electrophoresis in
c     bulk free solution. The bulk runs can be chekced  for consistency.

c     The code will allow simulations of a solution of DNA with fixed length.
c     According to worm-like chain model, Kunh length (which is the step legnth
c     in the coarse-grained model), is twice of the persistent length. 
c     Double helix DNA has a persistent length around 50nm. Accordingly, the 
c     Kuhn length is 100nm. 
c     With simple cubic lattice, we will not be able to map exactly
c     the conditions used in Han and craighead study. Instead, we will 
c     assume that lattice spacing is about 50nm, and each segment on the 
c     lattice represnts 150bps of the DNA.
c*****************************************************************88
c  Input files needed
c     44rep.input --- simple text file that contained initial file in, 
c                      file out, nosnap, interval, noequil, and eletric field
c***********************************************************
c  Output files:
c     Monitor.dat --- monitor chains Rg2, R2, and fsuc, energy during the 
c                     simulation
c     chainscont.trj --- contained chains cooridnates during the run nosnap.
c     traj.info      --contained inforamtion for the trajectory. used
c                    by analytraj.f 
c ***********************************************************
c  How to compile:
c      when the lattice dimension is too large, on some of the computers
c      like nucleotide and peptide, you can use the following options
c      g77 -mcmodel=medium -O program.f -o program
c     the option -mcmodel=medium will allow you to compile a code with
c     memory grater than 2GB.
c*******************************************************************

c     Here are very crude mapping.
c     on simple cubic lattice, it has a C_inf=1.5. 
c     So we want the conditions: 1.5Nb^2 = nl_kuhn^2. and
c     Nb = n l_kuhb. This lead to b=l_kuhn/1.5 = 66.6 nm. l_kuhn is the 
c     kugh segment length of DNA which is 100nm. We will aproximate
c     b as b=50nm. Each segment on the lattice is about 150bps. 
c     So for our chain of N=100, this means DNA length of 15,000 bps. 
c     N=1000, that means DNA length of 150,000. Let's stick to this mapping.
c                                             z
c     lattice structure is: ----------------- |
c                                             |-->x
c                           ---------------- 
c     xlayer is labeled by LSX, ylayer is labeled by LSY, z layer is labeld
c     by lSZ. Their indices are given by xsite,ysite,zsite
c     ************************************************************

c here is the new notes; trying to make only a bulk lattice without surface.

      implicit none
      Integer LSXMAX, LSYMAX, LSZMAX, NCHMAX, LCHMAX, NSITEMAX
      parameter (LSXMAX=400, LSYMAX=400, LSZMAX=400)
      parameter (LCHMAX=1000, NCHMAX=100)
      integer chains(NCHMAX,LCHMAX,3)
      integer empty(LSXMAX,LSYMAX,LSZMAX)
      integer chainscont(NCHMAX,LCHMAX,3)
c      double precision VPOTENTIAL(LSXMAX,LSYMAX,LSZMAX),ex
      integer LSX, LSY, LSZ
      integer NCH, LCH, IX, IX1 
      integer nosnap, interval,noequil,iter,iter1,ichain,jbead,isnap

      common/e/nch,lch
      common/d/lsx,lsy,lsz
      COMMON/XXX/IX,IX1

      integer lsuc,ixseed,ixseed1,IFINE
      integer na1old,na2old,na1tes,na2tes,nasold,nastes
      double precision epas1,epas2, econf,epas,frep
      
      double precision rgx,rgy,rgz,rendx,rendy,rendz,rgsq,rendsq
      double precision rgxsum,rgysum,rgzsum,rendxsum,rendysum,rendzsum
      double precision rgsqsum,rendsqsum
      double precision avgE,avgE2,SArg,SArend
      double precision avgfsuc,fsuc
      double precision theta, bound
      integer nads
c     some local variables that do not bear importance
      integer i,j,k,l,m,n
      real tarray(2),CPU
      double precision rn
      logical success
      
c     I only retained three common blocks. In the future, I will eliminate common/xxx
c     the other two blcoks are never modified during after the initial reading
C     ***************************************************
      CHARACTER*20 FILEIN,FILEOUT
      open(unit=10,file='44rep.input')
      read(10,'(A)') FILEIN
      read(10,'(A)') FILEOUT
      read(10,*) nosnap
      read(10,*) interval
      read(10,*) noequil
      close(unit=10)
c      read(10,*) ex
c     delete those related to votage difference --April 2022
c     ex has the units of (qDV0/Lk_BT).  
c     the noequil is the first snapshots will be omitted before collecting the data


      OPEN (UNIT=11,FILE=FILEIN,status='old')
      read(11,*) epas,nch,lch
      read(11,*) lsx,lsy,lsz

c  I need to set epas, (A-bead and solvent interaction) to allow 
c  for possible collapse study. In this code, we will monitor A-Solvent
c  contacts;--

      if(lsx.gt.lsxmax) then
         write(*,*) 'x dimension exceeds the limit',lsxmax
         stop
      end if
      if(lsy.gt.lsymax) then
         write(*,*) 'y dimension exceeds the limit',lsymax
         stop
      end if
      if(lsz.gt.lszmax) then
         write(*,*) 'z dimension exceeds the limit',lszmax
         stop
      end if
      if(nch.gt.nchmax) then
         write(*,*) 'nch is greater than allowed',nchmax
         stop
      end if
      if(lch.gt.lchmax) then
         write(*,*) 'lch is greater than allowed',lchmax
         stop
      end if


      do I=1,nch
         do J=1,lch
            read(11,*) (chains(i,j,k),k=1,3)
         end do
      end do
      close(unit=11)

c     epas1 is interaction with the first surface at z=1 and epas2 is interaction
c     with the second surface at z=lsz
c     surface 1 site is labled by empty=5, surface 2 sites are labeled by empty=6
cc     ix=time() works on gcc, or gfortran, 
      ix=time()
      ix1=ix+666578
      ixseed=ix
      ixseed1=ix1
      frep=0.8
c     temporarily leave the seed as fixed here
c     these are bad random generators, need to check it later
c     do the same trick as entropic trap. First relax the system with
c     noequil of interation but with ex=0

      call set_lattice_bulk(lsx,lsy,lsz,empty,lsxmax,lsymax,lszmax)
c     modified to set a bulk lattice

      do 3072 i=1,nch
         do 3073 j=1,lch
            empty(chains(i,j,1),chains(i,j,2),chains(i,j,3))=1
 3073    continue
 3072 continue

c      call set_potential(0.0d0,lsx,lsy,lsz,VPOTENTIAL,lsxmax,lsymax,
c     $ lszmax)
c     set_potential will set the V for every grid points. 
c     for this simple application, this is uneccessary. but I code this way
c     so when we apply the study to Entropic trap, we don't need to rewrite the
c     code. We only need to rewrite the set_potential

c     At the begining of simulations, , get total A-S contacts 
      call CALNAB (nastes,empty,chains,lsxmax,lsymax,lszmax,
     $ nchmax,lchmax)      
      nasold=nastes

      call gregar(chains,chainscont,nchmax,lchmax)
c     Monitor.dat will contain iter, econf,nads,rg2,rend2 ...
c      CALL CHECKS(chains,chainscont,empty,nchmax,lchmax,lsxmax,
c     $ lsymax,lszmax,IFINE)
c      write(11,*) 'Ifine at the begining Ifine=',Ifine

c temporarily insert this here to make sure chains are connected. 

      open(unit=14,file='Monitor.dat')
      write(14,'(a)')'#iter,econf,rgx,rgy,rgz,rendx,rendy,rendz'


      rgxsum=0.0
      rgysum=0.0
      rgzsum=0.0
      rendxsum=0.0
      rendysum=0.0
      rendzsum=0.0
      rgsqsum=0.0
      rendsqsum=0.0
      avgE=0.0
      avgE2=0.0
      avgfsuc=0.0
      isnap=0


      do iter=1,noequil+nosnap
         do iter1=1,interval
            do ichain=1,nch

c     insert reptation move here--it is faster. 

               call rnd(rn)
               if(rn.le.frep) then
                  goto 401
               else
                  goto 402
               endif

 401           call repmove(ichain,chains,
     $              empty,nasold,epas,success,nchmax,
     $              lchmax,lsxmax,lsymax,lszmax)

c     finish the reptation move here ,and move to 405
               goto 501
c     temporarily insert goto to skip internal moves
 402           do jbead=1,lch
c     with equal probablity, we will try one-bead and two bead move
                  call rnd(rn)
                  if(rn.le.0.5) then
                     call onebeadmove(ichain,jbead,chains,chainscont,
     $ empty,nasold,epas,success,nchmax,
     $                    lchmax,lsxmax,lsymax,lszmax)
                  else if( rn.gt.0.5) then
                     call twobeadmove(ichain,jbead,chains,chainscont,
     $ empty,nasold,epas,success,nchmax,
     $ lchmax,lsxmax,lsymax,lszmax)
                  end if
               end do

c     if we did reptatation move, skip over the internal move
 501           continue
c     the other two end do ends the iter=1,nch and iter=1,isnap 
            end do
         end do
         econf=nasold*epas

c     in this code, chainscont is continually updated during the move---after inserting
c     reptation move, chainscont is no longer updated, so call it here to get chainscont

         call gregar(chains,chainscont,nchmax,lchmax)
         call gyration (chainscont,rgx,rgy,rgz,rgsq,rendx,
     $        rendy,rendz,rendsq,nchmax,lchmax)

         write(14,*) iter*interval,econf,rgx,rgy,rgz,rendx,rendy,rendz
c     all the quantity output here, rgx, etc are "squared terms". rg is actually rg^2.

         if(iter.gt.noequil) then
            isnap=isnap+1
            rgxsum=rgxsum+rgx
            rgysum=rgysum+rgy
            rgzsum=rgzsum+rgz
            rendxsum=rendxsum+rendx
            rendysum=rendysum+rendy
            rendzsum=rendzsum+rendz
            rgsqsum=rgsqsum+rgsq
            rendsqsum=rendsqsum+rendsq
            avgE=avgE+econf
            avgE2=avgE2+econf*econf
         end if
      end do
      
C end of big do loop of iter
 135   format(10(1x,I6))

       open(unit=12,file=fileout)
       write(12,*) epas,nch,lch
       write(12,*) lsx,lsy,lsz
       do I=1,nch
          do J=1,lch
             write(12,*) (chains(i,j,k),k=1,3)
          end do
       end do

c     check to make sure the update of nas contacts are correct
      call CALNAB (nastes,empty,chains,lsxmax,lsymax,lszmax,
     $ nchmax,lchmax)      
      write(12,*) 'nasold',nasold,'nastes',nastes
      write(12,*) 'nosnap',nosnap,'interval',interval,'noequil',noequil

      CALL CHECKS(chains,chainscont,empty,nchmax,lchmax,lsxmax,
     $ lsymax,lszmax,IFINE)
      WRITE(12,*) 'IFINE=1 mean good. CHeck your IFINE=',IFINE

      close (unit=12)

      open(unit=13, file='summary.dat')
      write(13,*) 'program is SAWchain version of  APril 2022'
      write(13,*) 'The systems contain nch',nch, 'length of chain',lch
      write(13,*) 'lattice size used',lsx,lsy,lsz
      write(13,*) 'noequil',noequil,'interval',interval,'nosnap',nosnap
      write(13,*) 'total iteration=nosnap+noequil',iter

      write(13,*) 'seed used',ixseed,ixseed1
      write(13,*) '                '
      
      rgzsum=rgzsum/float(isnap)
      rgysum=rgysum/float(isnap)
      rgxsum=rgxsum/float(isnap)
      rgsqsum=rgsqsum/float(isnap)

      SArg=rgsqsum-(rgzsum+rgxsum+rgysum)**2

      rendxsum=rendxsum/float(isnap)
      rendysum=rendysum/float(isnap)
      rendzsum=rendzsum/float(isnap)
      rendsqsum=rendsqsum/float(isnap)
      SArend=rendsqsum-(rendxsum+rendysum+rendzsum)**2

      write(13,*) 'summary for radius of guration'
      write(13,*) 'rgx2,rgy2,rgz2,rg2',rgxsum,rgysum,rgzsum,
     $ rgxsum+rgysum+rgzsum
      write(13,*) 'variance of rg2',SArg
      write(13,*) 'error bars for rg2',1.95*sqrt(SArg)
     $ /sqrt(float(nch*isnap))

      write(13,*) '               '
      write(13,*) 'summary for end to end distance'
      write(13,*) 'rendx2,rendy2,rendz2,rend2',rendxsum,rendysum,
     $ rendzsum,rendxsum+rendysum+rendzsum
      write(13,*) 'variance of rend2',SArend
      write(13,*) 'error bars for rend2',1.95*sqrt(SArend)
     $ /sqrt(float(nch*isnap))      

      write(13,*) '               '
      avgE=avgE/float(isnap)
      avgE2=avgE2/float(isnap)
      SArg=avgE2-avgE*avgE
      write(13,*) 'average Energy',avgE
      write(13,*) 'standard deivation or variance',SArg
      write(13,*) 'error bar',1.95*sqrt(SArg)/sqrt(float(isnap))

      CPU=etime(tarray)
      write(13,*) 'total CPU used',CPU/60.0,'minutes'
      close(unit=13)
      
      STOP
      END

c***********************
C     ***************************************************
      SUBROUTINE RND(RN)
      implicit none
C     THIS SUBROUTINE GENERATES A RANDOM NUMBER VIZ. RN
C     BETWEEN 0 AND 1
      integer d,m,IX,IX1
      double precision RN
      parameter (m=923456789, d=2147483647)
      COMMON/XXX/IX,IX1
      IX=IX*m
      rn=mod(abs(ix),d)
      rn=rn/float(d)
      RETURN
      END
C        ********************************************************
      SUBROUTINE RANDU(X,IRN)
C
C     THIS SUBROUTINE GENERATES A RANDOM NUMBER
C     BETWEEN 1 AND INT(X)
C
C     ***********************
      REAL RN,DUM
      integer d,m,IX,IX1,IRN
      parameter(m=992631493, d=2147483647)
      COMMON/XXX/IX,IX1

      IX1=IX1*m
      rn=mod(abs(ix1),d)
      rn=rn/float(d)
      IRN=IFIX(RN*X)+1
      IF (IRN.GT.IFIX(X)) THEN
      IRN=IRN-1
      END IF
      RETURN
      END
C      **********************************************
c     ************************************************************
c     no longer needed
      function kar(xsite,ysite,zsite)
      implicit none
      integer kar,xsite,ysite,zsite
      integer lsx,lsy,lsz
      common/d/lsx,lsy,lsz
c     through this function call, I avoid the declaration of kar
c                                             z
c     lattice structure is: ----------------- |
c                                             |-->x
c                           ---------------- 

      kar=(zsite-1)*lsx*lsy+(ysite-1)*lsx+xsite
      return
      end 

c     *********************************
      subroutine set_lattice_bulk(lsx,lsy,lsz,empty,
     $ lsxmax,lsymax,lszmax)
      implicit none
      integer lsx,lsy,lsz,lsxmax,lsymax,lszmax
      integer empty(lsxmax,lsymax,lszmax)
      integer i,j,k

c     set lattice for the slit. if runing a bulk simulation, then
c     just need to set a different type lattice
 
      do i=1,lsx
         do j=1,lsy
            do k=1,lsz
               empty(i,j,k)=0
c               if(k.eq.1) then
c                  empty(i,j,k)=5
c               end if
c               if(k.eq.lsz) then
c                  empty(i,j,k)=6
c               end if

            end do
         end do
      end do
      
      return
      end

c     **************************************************************
      SUBROUTINE CALNAB (na1,empty,chains,lsxmax,lsymax,lszmax,
     $ nchmax,lchmax)

c modified: this sub calcualte number A-solvent contants, A-A contants
c  na1 is A-solvent, 

      implicit none
      integer na1
      integer nchmax,lchmax,lsxmax,lsymax,lszmax
      integer nch,lch
      integer lsx,lsy,lsz
      common/e/nch,lch
      common/d/lsx,lsy,lsz
      INTEGER CHAINS(nchmax,lchmax,3),empty(lsxmax,lsymax,lszmax)
      integer i,j,k,isite
      integer nnn(6,3),nb(3)

      na1=0
      
      DO i=1,NCH
         DO j=1,lch
            do k=1,3
               nb(k)=chains(i,j,k)
            end do
            call n4mod(nb,nnn)
            do k=1,6
               if (empty(nnn(k,1),nnn(k,2),nnn(k,3)).eq.0) na1=na1+1
            end do
         end do
      end do
      
      RETURN
      END
c     **********************************************************************
      SUBROUTINE N4MOD (NB,NN)
      implicit none
      
      integer nb(3), nn(6,3),vector(6,3)
      integer lsx,lsy,lsz
      common/d/lsx,lsy,lsz
      integer vectors(6,3),ivector
      data vectors/1,0,0, -1,0,0, 0,1,0, 0,-1,0, 0,0,1, 0,0,-1/

C     ***************************************************
c     This subroutine gives six nearest neightbor sites. ndum array 
c     temporary array, which stores the i,j,k three coordinate 
c     varaible.
c     *******************************************************

         do ivector=1,6
            nn(ivector,1)=nb(1)+vectors(ivector,1)
            if(nn(ivector,1).lt.1) nn(ivector,1)=lsx
            if(nn(ivector,1).gt.lsx) nn(ivector,1)=1

            nn(ivector,2)=nb(2)+vectors(ivector,2)
            if(nn(ivector,2).lt.1) nn(ivector,2)=lsy
            if(nn(ivector,2).gt.lsy) nn(ivector,2)=1

            nn(ivector,3)=nb(3)+vectors(ivector,3)
            if(nn(ivector,3).lt.1) nn(ivector,3)=lsz
            if(nn(ivector,3).gt.lsz) nn(ivector,3)=1            
         end do
 
      RETURN
      END
c     ****************************************************************
 
         SUBROUTINE METRO (NANEW,NAOLD,epas, IFLAG)
         implicit none
         integer nanew,naold
         double precision epas,deltaE,cond
         integer iflag
         double precision rn1

         IFLAG=1
         deltaE= (NAnew-NAold)*epas
 
         IF (DELTAE.GT.0) THEN
             COND=EXP(-DELTAE)
             CALL RND (RN1)
             IF (RN1.GT.COND) IFLAG=0
         END IF
 
         RETURN
         END
c     ****************************************************************

      subroutine repmove(ichain,chains,empty,
     $ nasold,epas,success,nchmax,lchmax,
     $ lsxmax,lsymax,lszmax)

      implicit none
      integer ichain,jbead,nasold,nasnew
      integer lsxmax,lsymax,lszmax,nchmax,lchmax,nch,lch
      integer lsx,lsy,lsz
      common/e/nch,lch
      common/d/lsx,lsy,lsz
      integer empty(lsxmax,lsymax,lszmax),chains(nchmax,lchmax,3)
      integer chainscont(nchmax,lchmax,3),chain(lchmax,3)
      integer newchain(lchmax,3)
c      double precision vpotential(lsxmax,lsymax,lszmax),deltaV,ex
      double precision epas1,epas2,epas,rn
      logical success
      integer iflag,ibflag,irev
      integer enpnew(3),enpold(3)
      integer i,j,k,l,m,nc,l1,l2,l3,isite,isitehead,isitetail
      integer nb(3),nn(6,3),nbplus(3),nbminus(3),ntail(3)

c     perform reptation move here within the subroutine, perform on the ichain 


      success=.false.
      do j=1,lch
         do k=1,3
            chain(j,k)=chains(ichain,j,k)
         end do
      end do
c     first decide move head or ends

      call rnd(rn)
      if(rn.lt.0.5) then
         irev=0
      else if(rn.gt.0.5) then
         call rev(lch,chain,lchmax)
         irev=1
      end if

      do k=1,3
         nb(k)=chain(1,k)
         nbminus(k)=chain(2,k)
         ntail(k)=chain(lch,k)
      end do

      call endfind(nb,nbminus,enpnew)
c      Hossain modified; commented below line      
c      if(empty(enpnew(1),enpnew(2),enpnew(3)).eq.0) then
      nasnew=nasold

c     now you can try to move the chain to the new head at enpnew
         call n4mod(enpnew,nn)
         do 10 j=1,6
            l1=nn(j,1)
            l2=nn(j,2)
            l3=nn(j,3)
            isite=lsx*lsy*l3+lsy*l2+l1
            if(empty(l1,l2,l3).eq.0) nasnew=nasnew+1
            if(empty(l1,l2,l3).eq.1) nasnew=nasnew-1
            if(isite.eq.isitetail) nasnew=nasnew+1
 10      continue

c     tail is emptied out
         call n4mod(ntail,nn)
         do 20 k=1,6
            l1=nn(k,1)
            l2=nn(k,2)
            l3=nn(k,3)
            isite=lsx*lsy*l3+lsy*l2+l1
            if(empty(l1,l2,l3).eq.0) nasnew=nasnew-1
            if(empty(l1,l2,l3).eq.1) nasnew=nasnew+1
c            if(isite.eq.isitehead) then
c             nadd=nadd-1
c             nsub=nsub-1
 20      continue
         call METRO (NASNEW,NASOLD,EPAS,iflag)

         if(iflag.eq.1) then

c     move the chain
            do k=1,3
               newchain(1,k)=enpnew(k)
            end do
            do I=2,lch
               do k=1,3
                  newchain(I,k)=chain(I-1,k)
               end do
            end do
c     newchain stores the position of temporarily
            do I=1,lch
               do k=1,3
                  chain(I,k)=newchain(I,k)
               end do
            end do
c     now update empty
            empty(ntail(1),ntail(2),ntail(3))=0
            empty(enpnew(1),enpnew(2),enpnew(3))=1
            nasold=nasnew
            success=.true.
         endif
      end if
      if(irev.eq.1) then
         call rev(lch,chain,lchmax)
      end if
c     need to update chains and chainscont
      do j=1,lch
         do k=1,3
            chains(ichain,j,k)=chain(j,k)
         end do
      end do
      return
      end
 
c     *******************************************************************
      subroutine onebeadmove(ichain,jbead,chains,chainscont,empty,
     $ nasold,epas,success,nchmax,lchmax,
     $ lsxmax,lsymax,lszmax)
      implicit none
      integer ichain,jbead,na1old,na2old,na1new,na2new,nasold,nasnew
      integer lsxmax,lsymax,lszmax,nchmax,lchmax,nch,lch
      integer lsx,lsy,lsz
      common/e/nch,lch
      common/d/lsx,lsy,lsz
      integer empty(lsxmax,lsymax,lszmax),chains(nchmax,lchmax,3)
      integer chainscont(nchmax,lchmax,3)
c      double precision vpotential(lsxmax,lsymax,lszmax),deltaV,ex
      double precision epas1,epas2,epas
      logical success
      integer iflag,bflag
c     local variables
      integer enpnew(3),enpold(3)
      integer i,j,k,l,m,nc
      integer nb(3),nn(6,3),nbplus(3),nbminus(3)
      integer dx(3),lsl(3)

c     lsl(3) defined temorarily for lsx,lsy,lsz
      lsl(1)=lsx
      lsl(2)=lsy
      lsl(3)=lsz
      

      success=.false.

c ith chain, jth bead will be moved. as one bead move
c: nbplus, nbminus is rememerbing neighboring beads along the chain. 
c if jbead=1, jbminus does not exist. if jbead=lch, nbplus does not exist.      

      do k=1,3
         nb(k)=chains(ichain,jbead,k)
         if(jbead.eq.1) then
            nbplus(k)=chains(ichain,jbead+1,k)
            nbminus(k)=0
         else if (jbead.eq.lch) then
            nbminus(k)=chains(ichain,jbead-1,k)
            nbplus(k)=0
         else
            nbminus(k)=chains(ichain,jbead-1,k)
            nbplus(k)=chains(ichain,jbead+1,k)
         end if
      end do

      call n4mod(nb,nn)
      if(jbead.eq.1) then
ccheck ebfindnew, bfind for random walk
         call ebfindnew(nb,nbplus,enpnew,bflag)
      else if(jbead.eq.lch) then
         call ebfindnew(nb,nbminus,enpnew,bflag)
      else
         call bfind(nb,nbminus,nbplus,enpnew,bflag)
      end if

c     use bflag to distinguish the case when enpnew is not assigned
c     ibflag will tell us if we have a "place to move to"

      if(bflag.eq.1) then
         do k=1,3
            enpold(k)=nb(k)
         end do
         if(empty(enpnew(1),enpnew(2),enpnew(3)).eq.0) then
c            na1new=na1old
c            na2new=na2old
            nasnew=nasold
c     record the nas counts

            call n4mod(enpold,nn)
            do k=1,6
               if(empty(nn(k,1),nn(k,2),nn(k,3)).eq.0) nasnew=nasnew-1
               if(empty(nn(k,1),nn(k,2),nn(k,3)).eq.1) nasnew=nasnew+1
            end do

            call n4mod(enpnew,nn)
            do k=1,6
               if(empty(nn(k,1),nn(k,2),nn(k,3)).eq.0) nasnew=nasnew+1
               if(empty(nn(k,1),nn(k,2),nn(k,3)).eq.1) nasnew=nasnew-1
            end do

            call METRO (NASNEW,NASOLD,EPAS,iflag)
            if(iflag.eq.1) then
               success=.true.
c     now move
               do k=1,3
                  chains(ichain,jbead,k)=enpnew(k)
               end do
c     now update chainscont

               do k=1,3
                  dx(k)=enpnew(k)-enpold(k)
                  if(dx(k).eq.(1-lsl(k))) dx(k)=+1
                  if(dx(k).eq.(lsl(k)-1)) dx(k)=-1
                  chainscont(ichain,jbead,k)=chainscont(ichain,jbead,k)
     $                 +dx(k)
               end do
               empty(nb(1),nb(2),nb(3))=0
               empty(enpnew(1),enpnew(2),enpnew(3))=1
               nasold=nasnew
            end if
         end if
      end if
      return
      end

c     ********************************
      subroutine endfind(enpold,enpold2,enpnew)

c     this is called during reptation move. Pick a site near the end to move into
c     enpold is the first bead of the chain, enpew is the new site to move into
c     you need to randomly pick one out of five neighbor, 
c     empty is not passed here. 
      implicit none
      integer enpold(3),enpnew(3),iflag,enpold2(3)
      integer nd(3),md(3),nn(6,3),mm(5,3)
c     local variables
      common/d/lsx,lsy,lsz
      integer lsx,lsy,lsz
      integer isite,ipk
      integer i,j,k,isum,isiteold,isiteold2,isitenew

c     this iflag is always set to 1

c      iflag=1

      isiteold=lsx*lsy*enpold(3)+lsx*enpold(2)+enpold(1)
      isiteold2=lsx*lsy*enpold2(3)+lsx*enpold2(2)+enpold2(1)
c      do k=1,3
c         write(*,*) 'enpold',enpold(1),enpold(2),enpold(3)
c         write(*,*) 'enpold2',enpold2(1),enpold2(2),enpold2(3)
c      end do

      call n4mod(enpold,nn)
c      isite=0
c      do j=1,6
c         isitenew=lsx*lsy*nn(j,3)+lsx*nn(j,2)+nn(j,1)
cHossain modified; commented below line
c         if(isitenew.ne.isiteold2) then
c            isite=isite+1
c            do k=1,3
c               mm(isite,k)=nn(j,k)
c            end do
c         end if
c      end do
c      if(isite.eq.6) then
c         write(*,*) 'isite',isite,'isiteold2',isiteold2,
c     $ 'isiteold',isiteold
c         do k=1,3
c            write(*,*) 'enpold',enpold(1),enpold(2),enpold(3)
c            write(*,*) 'enpold2',enpold2(1),enpold2(2),enpold2(3)
c         end do
c         do j=1,6
c            isitenew=lsx*lsy*nn(j,3)+lsx*nn(j,2)+nn(j,1)
c            write(*,*) 'isitenew',isitenew
c         end do
c      end if
c     now mm contained 6 sites, the enpold2 is eliminated

      call randu(6.0, ipk)
      do k=1,3
         enpnew(k)=nn(ipk,k)
      end do
c     when exit from this routine, enpnew is assigned, and iflag

      return
      end
c     ******************************************************


         

c     ********************************
      subroutine ebfindnew(enpold,enpold2,enpnew,bflag)
c     this is called when the picked bead is an end bead.
c     this will select the four positions that form 
c     diagonal to the bond bewteen enpold-enpold2
      implicit none
      integer enpold(3),enpold2(3),enpnew(3),bflag
      integer lsx,lsy,lsz
      integer nd(3),md(3),nn(6,3),mm(4,3)
      common/d/lsx,lsy,lsz
c     local variables
      integer isite,ipk
      integer i,j,k,isum

c     this bflag is always set to 1
      bflag=1
      do 10 i=1,3
         nd(i)=enpold(i)-enpold2(i)
 10   continue
c     I don't worry about the periodic boudnary problem here. Since 
c     all I care is the cross product should be zero.

      call n4mod(enpold2,nn)
      isite=0
      do j=1,6
         do i=1,3
            md(i)=nn(j,i)-enpold2(i)
         end do
c     now determine the cross product
         isum=0
         do i=1,3
            isum=isum+md(i)*nd(i)
         end do
         if(isum.eq.0) then
            isite=isite+1
            do i=1,3
               mm(isite,i)=nn(j,i)
            end do
         end if
      end do
c     now mm contained 4 sites.

      call randu(4.0, ipk)
      do k=1,3
         enpnew(k)=mm(ipk,k)
      end do
      return
      end
c     ******************************************************



c     ********************************
      subroutine bfind(nb,nbminus,nbplus,enpnew,bflag)
c     this is called when the picked bead is an internal bead
c     and we try to flip the bead to its diagonal position.
c      ...
c         |
c         .....
c     The movement is only possible when nbminus-nb and nb-nbplus bonds
c     form 90 degree. One of the three coordinates for the three beads are the
c     same

      implicit none
      integer nb(3),nbminus(3),nbplus(3),enpnew(3)    
      integer lsx,lsy,lsz
      integer bflag
      common/d/lsx,lsy,lsz
c     local variables
      integer i,j,k,irem,iequal

c     find the common ccordinates bewteen nb(3), nbminus(3),nbplus(3)

      iequal=0
      do i=1,3
         if((nb(i).eq.nbminus(i)).and.(nb(i).eq.nbplus(i))) then
            iequal=iequal+1
            irem=i
         end if
      end do
      if(iequal.eq.0) then
         write(*,*) 'iequal is zero'
         return
      end if
      if(iequal.gt.1) then
         bflag=0
         enpnew(1)=0
         enpnew(2)=0
         enpnew(3)=0
      else if(iequal.eq.1) then
         bflag=1
         enpnew(irem)=nb(irem)
         do k=1,3
            if(nb(k).ne.nbminus(k)) enpnew(k)=nbminus(k)
            if(nb(k).ne.nbplus(k)) enpnew(k)=nbplus(k)
         end do
      end if

      return
      end
c**************************************************************************
      subroutine twobeadmove(ichain,jbead,chains,chainscont,empty,
     $ nasold,epas,success,nchmax,lchmax,
     $ lsxmax,lsymax,lszmax)
      implicit none
      integer ichain,jbead,na1old,na2old,na1new,na2new,nasold,nasnew
      integer lsxmax,lsymax,lszmax,nchmax,lchmax,nch,lch
      integer lsx,lsy,lsz
      common/e/nch,lch
      common/d/lsx,lsy,lsz
      integer empty(lsxmax,lsymax,lszmax),chains(nchmax,lchmax,3)
      integer chainscont(nchmax,lchmax,3)
c      double precision vpotential(lsxmax,lsymax,lszmax),deltaV,ex
      double precision epas1,epas2,epas
      logical success
      integer iflag,nmove,bflag
c     local variables
      integer enpnew1(3),enpnew2(3),enpold1(3),enpold2(3),enpnew(3)
      integer head(3),tail(3)
      integer i,j,k,l,m,nc
      integer nb(3),nn(6,3),nbplus(3),nbminus(3)
      integer lsl(3),dx(3)
      
c     in this subroutnie I used a feature "return". Return should exit out of the
c     subroutine, but I am not 100% sure yet. Need to test it.

      success=.false.
      if((jbead.eq.1).or.(jbead.eq.lch)) return

      do k=1,3
         nb(k)=chains(ichain,jbead,k)
         nbplus(k)=chains(ichain,jbead+1,k)
         nbminus(k)=chains(ichain,jbead-1,k)
      end do
      
      call bfind(nb,nbminus,nbplus,enpnew,bflag)
      if(bflag.eq.0) return
      if(bflag.eq.1) then
         nmove=0
         if(jbead.lt.(lch-1)) then
            iflag=1
            do k=1,3
               if(enpnew(k).ne.chains(ichain,jbead+2,k))iflag=0
            end do
            if(iflag.eq.1) then
               nmove=1
            end if
         end if
         if(jbead.gt.2) then
            iflag=1
            do k=1,3
               if(enpnew(k).ne.chains(ichain,jbead-2,k))iflag=0
            end do
            if(iflag.eq.1) then
               nmove=-1
            end if
         end if
      end if

c     here the assumption is that the crank-type geometry can either
c     involve +1 bead, or -1 bead

c     I always name the picked bead as enpold1, and the bead connected with
c     enpold1 as the head. enpold2 will change depending on the nmove sign
c     the one connected with enpold2 is called tail. 
      if(nmove.eq.0) return
      if(nmove.eq.1) then
         do k=1,3
            enpold1(k)=nb(k)
            enpold2(k)=nbplus(k)
            head(k)=nbminus(k)
            tail(k)=chains(ichain,jbead+2,k)
         end do
      else if (nmove.eq.-1) then
         do k=1,3
            enpold1(k)=nb(k)
            enpold2(k)=nbminus(k)
            tail(k)=chains(ichain,jbead-2,k)
            head(k)=nbplus(k)
         end do
      end if

      call crank(enpold1,enpold2,head,tail,enpnew1,enpnew2)
c      if(empty(enpnew1(1),enpnew1(2),enpnew1(3)).ne.0) return
c      if(empty(enpnew2(1),enpnew2(2),enpnew2(3)).ne.0) return

      if((empty(enpnew1(1),enpnew1(2),enpnew1(3)).eq.0).and.
     $ (empty(enpnew2(1),enpnew2(2),enpnew2(3)).eq.0)) then
      
         nasnew=nasold

         call n4mod(enpold1,nn)
         do k=1,6
            if(empty(nn(k,1),nn(k,2),nn(k,3)).eq.0) nasnew=nasnew-1
            if(empty(nn(k,1),nn(k,2),nn(k,3)).eq.1) nasnew=nasnew+1
         end do
         nasnew=nasnew-2
c     I think enpold1 neighbor would have two beads along the chain that should not be counted
         call n4mod(enpnew1,nn)
         do k=1,6
            if(empty(nn(k,1),nn(k,2),nn(k,3)).eq.0) nasnew=nasnew+1
            if(empty(nn(k,1),nn(k,2),nn(k,3)).eq.1) nasnew=nasnew-1
         end do
         nasnew=nasnew+1 

         call n4mod(enpold2,nn)
         do k=1,6
            if(empty(nn(k,1),nn(k,2),nn(k,3)).eq.0) nasnew=nasnew-1
            if(empty(nn(k,1),nn(k,2),nn(k,3)).eq.1) nasnew=nasnew+1
         end do
         nasnew=nasnew-2

         call n4mod(enpnew2,nn)
         do k=1,6
            if(empty(nn(k,1),nn(k,2),nn(k,3)).eq.0) nasnew=nasnew+1
            if(empty(nn(k,1),nn(k,2),nn(k,3)).eq.1) nasnew=nasnew-1
         end do
         nasnew=nasnew-1

         call METRO (NASNEW,NASOLD,EPAS, IFLAG)

         if(iflag.eq.1) then
            success=.true.
            lsl(1)=lsx
            lsl(2)=lsy
            lsl(3)=lsz
            if(nmove.eq.1) then
               do k=1,3
                  chains(ichain,jbead,k)=enpnew1(k)
                  chains(ichain,jbead+1,k)=enpnew2(k)
               end do
c     now update chainscont
               do k=1,3
                  dx(k)=enpnew1(k)-enpold1(k)
                  if(dx(k).eq.(1-lsl(k))) dx(k)=+1
                  if(dx(k).eq.(lsl(k)-1)) dx(k)=-1
                  chainscont(ichain,jbead,k)=chainscont(ichain,jbead,k)
     $                 +dx(k)
               end do
               do k=1,3
                  dx(k)=enpnew2(k)-enpold2(k)
                  if(dx(k).eq.(1-lsl(k))) dx(k)=+1
                  if(dx(k).eq.(lsl(k)-1)) dx(k)=-1
                  chainscont(ichain,jbead+1,k)=chainscont(ichain,
     $                 jbead+1,k)+dx(k)
               end do
            else if(nmove.eq.-1) then
               do k=1,3
                  chains(ichain,jbead,k)=enpnew1(k)
                  chains(ichain,jbead-1,k)=enpnew2(k)
               end do
c     now update the chainscont
               do k=1,3
                  dx(k)=enpnew1(k)-enpold1(k)
                  if(dx(k).eq.(1-lsl(k))) dx(k)=+1
                  if(dx(k).eq.(lsl(k)-1)) dx(k)=-1
                  chainscont(ichain,jbead,k)=chainscont(ichain,jbead,k)
     $                 +dx(k)
               end do
               do k=1,3
                  dx(k)=enpnew2(k)-enpold2(k)
                  if(dx(k).eq.(1-lsl(k))) dx(k)=+1
                  if(dx(k).eq.(lsl(k)-1)) dx(k)=-1
                  chainscont(ichain,jbead-1,k)=chainscont(ichain,
     $                 jbead-1,k)+dx(k)
               end do
            end if
            empty(enpold1(1),enpold1(2),enpold1(3))=0
            empty(enpold2(1),enpold2(2),enpold2(3))=0
            empty(enpnew1(1),enpnew1(2),enpnew1(3))=1
            empty(enpnew2(1),enpnew2(2),enpnew2(3))=1
            nasold=nasnew
         end if
      end if
      return
      end
c     I corrected the problem in two bead move when I don't use "return" when 
c     empty(enpnew).ne.0) this maybe something weird about gfortan 

c     *******************************************************************
      subroutine crank(enpold1,enpold2,head,tail,enpnew1,enpnew2)

      implicit none
      integer enpold1(3),enpold2(3),head(3),tail(3)
      integer enpnew1(3),enpnew2(3)
      integer i,j,k,l,krem,new
      integer iequal
      integer lsl(3)
      integer lsx,lsy,lsz
      common/d/lsx,lsy,lsz
      double precision rn1
c     this will determine the enpnew1,enpnew2, based on the four
c     beads position in crank-motion
c     xxxx      
c     |  |
c     cccccc 

      lsl(1)=lsx
      lsl(2)=lsy
      lsl(3)=lsz
      
      iequal=0
      do k=1,3
         if((head(k).eq.enpold2(k)).and.(head(k).eq.enpold1(k))) then
            iequal=iequal+1
            krem=k
         end if
      end do

      if(iequal.ne.1) then
         write(*,*) 'crank has errors'
         stop
      end if
      
      call rnd(rn1)
      if(rn1.gt.0.5) then
         new=enpold1(krem)+1
         if(new.gt.lsl(krem)) new=1
      else if (rn1.le.0.5) then
         new=enpold1(krem)-1
         if(new.eq.0) new=lsl(krem)
      end if
      
      enpnew1(krem)=new
      enpnew2(krem)=new

      do k=1,3
         if(k.ne.krem) then
            enpnew1(k)=head(k)
            enpnew2(k)=tail(k)
         end if
      end do
      return
      end
c     ****************************************************

      SUBROUTINE GREGAR (chains,chainscont,nchmax,lchmax)

C     this sub eliminates the periodic boundary conditions for each
c     chains. So we can calculate the gyration directly.
c     for each chain, we start off the position of the first bead, and
c     trace out the chainscont
C     ***********************************************************
      implicit none

      integer nchmax,lchmax
      integer chains(nchmax,lchmax,3),chainscont(nchmax,lchmax,3)
      integer nch,lch,lsx,lsy,lsz
      common/e/nch,lch
      common/d/lsx,lsy,lsz

      integer i,j,k,I1,I2,J1,J2,K1,K2

      do i=1,nch
         do k=1,3
            chainscont(i,1,k)=chains(i,1,k)
         end do
      end do

      do i=1,nch
         I1=chains(i,1,1)
         J1=chains(i,1,2)
         K1=chains(i,1,3)
         do j=2,lch
            I2=chains(i,j,1)
            J2=chains(i,j,2)
            K2=chains(i,j,3)
            
         IF ((I2-I1).EQ.0) chainscont(i,j,1)=chainscont(i,j-1,1)
         IF ((I2-I1).EQ.1) chainscont(i,j,1)=chainscont(i,j-1,1)+1
         IF ((I2-I1).EQ.-1) chainscont(i,j,1)=chainscont(i,j-1,1)-1
         IF ((I2-I1).EQ.(LSX-1)) chainscont(i,j,1)=chainscont(i,j-1,1)-1
         IF ((I2-I1).EQ.(1-LSX)) chainscont(i,j,1)=chainscont(i,j-1,1)+1

         IF ((J2-J1).EQ.0) chainscont(i,j,2)=chainscont(i,j-1,2)
         IF ((J2-J1).EQ.1) chainscont(i,j,2)=chainscont(i,j-1,2)+1
         IF ((J2-J1).EQ.-1) chainscont(i,j,2)=chainscont(i,j-1,2)-1
         IF ((J2-J1).EQ.(LSY-1)) chainscont(i,j,2)=chainscont(i,j-1,2)-1
         IF ((J2-J1).EQ.(1-LSY)) chainscont(i,j,2)=chainscont(i,j-1,2)+1


         IF ((K2-K1).EQ.0) chainscont(i,j,3)=chainscont(i,j-1,3)
         IF ((K2-K1).EQ.1) chainscont(i,j,3)=chainscont(i,j-1,3)+1
         IF ((K2-K1).EQ.-1) chainscont(i,j,3)=chainscont(i,j-1,3)-1
         IF ((K2-K1).EQ.(LSZ-1)) chainscont(i,j,3)=chainscont(i,j-1,3)-1
         IF ((K2-K1).EQ.(1-LSZ)) chainscont(i,j,3)=chainscont(i,j-1,3)+1
         I1=I2
         J1=J2
         K1=K2
         end do
      end do
      
      RETURN
      END

c     ********************************************************
      SUBROUTINE gyration (chainscont,rgx,rgy,rgz,rgsq,rendx,
     $ rendy,rendz,rendsq,nchmax,lchmax)

c     this subroutine will calcaluate the gyration and end-to-end
c     distance using chainscont cooridnates
c     all rg,rend is squared value

      implicit none
      integer nchmax,lchmax,nch,lch
      common/e/nch,lch
      integer chainscont(nchmax,lchmax,3)
      double precision chncm(nchmax,3)
      double precision rgx,rgy,rgz,rgsq,rendx,rendy,rendz,rendsq
      double precision rgxi,rgyi,rgzi,rendxi,rendyi,rendzi
      integer i,j, k

c     pass out the rg2 and r4=(rg2)^2 averages.
c     pass out the R2 and R4 averages

      do i=1,nch
         do k=1,3
            chncm(i,k)=0.0
         end do
      end do

      do i=1,nch
         do k=1,3
            do j=1,lch
               chncm(i,k)=chncm(i,k)+chainscont(i,j,k)
            end do
            chncm(i,k)=chncm(i,k)/float(lch)
         end do
      end do

      rgx=0.0
      rgy=0.0
      rgz=0.0
      rgsq=0.0
      rendx=0.0
      rendy=0.0
      rendz=0.0
      rendsq=0.0

      do i=1,nch
         rgzi=0.0
         rgxi=0.0
         rgyi=0.0
         do j=1,lch
            rgzi=rgzi+(chainscont(i,j,3)-chncm(i,3))**2
            rgyi=rgyi+(chainscont(i,j,2)-chncm(i,2))**2
            rgxi=rgxi+(chainscont(i,j,1)-chncm(i,1))**2
         end do
         rgzi=rgzi/float(lch)
         rgyi=rgyi/float(lch)
         rgxi=rgxi/float(lch)
         rgz=rgz+rgzi
         rgy=rgy+rgyi
         rgx=rgx+rgxi
         rgsq=rgsq+(rgzi+rgyi+rgxi)**2

         rendzi=(chainscont(i,lch,3)-chainscont(i,1,3))**2
         rendyi=(chainscont(i,lch,2)-chainscont(i,1,2))**2
         rendxi=(chainscont(i,lch,1)-chainscont(i,1,1))**2
         rendx=rendx+rendxi
         rendy=rendy+rendyi
         rendz=rendz+rendzi
         rendsq=rendsq+(rendzi+rendyi+rendxi)**2
      end do

      rgz=rgz/float(nch)
      rgy=rgy/float(nch)
      rgx=rgx/float(nch)
      rgsq=rgsq/float(nch)
      rendx=rendx/float(nch)
      rendy=rendy/float(nch)
      rendz=rendz/float(nch)
      rendsq=rendsq/float(nch)

      RETURN
      END


C     **********************************************
 
       SUBROUTINE CHECKS (chains,chainscont, empty,nchmax,lchmax,
     $ LSXMAX,LSYMAX,lSZMAX,IFINE)

C       THIS SUB CHECKS THE FOLLOWING: CHAIN CONNECTIVITY,NO 2 BEADS
C       OCCUPY ON THE SAME LATTICE POINT AND AN OCCUPIED LATTICE POINT
C       IS NOT EMPTY IN THE EMPTY ARRAY. IF EVERYTHING IS OKAY 
c       Ifine=1
       
       implicit none
       INTEGER CHAINS,EMPTY,IFINE,chainscont
       integer nchmax,lchmax,lsxmax,lsymax,lszmax
       integer nch,lch,lsx,lsy,lsz
       common/e/nch,lch
       common/d/lsx,lsy,lsz
       DIMENSION CHAINS(NCHMAX,LCHMAX,3),EMPTY(LSXMAX,LSYMAX,LSZMAX)
       dimension chainscont(nchmax,lchmax,3)
       integer nn(6,3),nb(3),nb2(3)
       integer i,j,k,isite1,isite2,isum
       integer kar
       
       IFINE=1
       DO I=1,NCH
          DO J=1,LCH
             do K=1,3
                nb(k)=chains(i,j,k)
             end do
             if(empty(nb(1),nb(2),nb(3)).ne.1) then
                ifine=0
                write(12,*) nb(1),nb(2),nb(3),'its EMPTY value',
     $               empty(nb(1),nb(2),nb(3)),'ith chain jbead',i,j
                return
             end if
          end do
       end do

       do i=1,nch
          do j=1,lch-1
             do K=1,3
                nb(k)=chains(i,j,k)
                nb2(k)=chains(i,j+1,k)
             end do
             call n4mod(nb,nn)
             isite1=kar(nb2(1),nb2(2),nb2(3))
             isum=0
             do k=1,6
                isite2=kar(nn(k,1),nn(k,2),nn(k,3))
                if(isite2.eq.isite1) isum=isum+1
             end do
             if(isum.ne.1) then
                ifine=0
                write(12,*) 'chain is not connected'
                write(12,*) 'ith j chain',i,j 
                return
             end if
          end do
       end do

c     check chainscont

       do i=1,nch
          do j=1,lch-1
             do K=1,3
                nb(k)=chainscont(i,j,k)
                nb2(k)=chainscont(i,j+1,k)
             end do
             isum=0
             do k=1,3
                isum=isum+(nb2(k)-nb(k))
             end do
             if(abs(isum).ne.1) then
                ifine=0
                write(12,*) 'chainscont is not connected'
                write(12,*) 'ith jth chain',i, j
                return
             end if
          end do
       end do
              
       
       RETURN
       END

c  ****************************************************8
      subroutine adsorption(chains,nchmax,lchmax,
     $ nads,bound,theta)

       implicit none
       integer nchmax,lchmax
       INTEGER CHAINS(nchmax,lchmax,3)
       integer nch,lch,lsx,lsy,lsz
       common/e/nch,lch
       common/d/lsx,lsy,lsz
       integer nads
       double precision bound, theta,bound_each
       integer i,j
       logical flag_ads
c      adsorption is considered whenever a chain has a bead in layer z=2
c     so this can notbeen directly used for entropic trap systems

       theta=0.0
       bound=0.0
       nads=0.0
       do i=1,nch
          flag_ads=.false.
          bound_each=0.0
          do j=1,lch
             if(chains(i,j,3).eq.2) then
                flag_ads=.true.
                theta=theta+1.0
                bound_each=bound_each+1.0
             end if
          end do
          if(flag_ads) then
             nads=nads+1
             bound_each=bound_each/float(lch)
             bound=bound+bound_each
          end if
       end do

       if(nads.ne.0) then
          bound=bound/float(nads)
       end if
       theta=theta/float(lsx*lsy)
       
       return
       end

       
c   ******************************************
C     THIS SUB JUST REVERSES THE CHAIN AND STORES IT BACK IN CHAIN
      SUBROUTINE REV(LCH,chain,lchmax)
      implicit none
      integer chain(lchmax,3),ITEMP(lchMAX,3)
      integer i,lch,idum,lchmax,k

      DO 100 I=1,LCH
         IDUM=LCH+1-I
         do k=1,3
            ITEMP(IDUM,k)=CHAIN(I,k)
         end do
 100  CONTINUE
 
      DO 101 I=1,LCH
         do k=1,3
            CHAIN(I,k)=ITEMP(I,k)
         end do
 101  CONTINUE
 
      RETURN
      END

c ***************************************
