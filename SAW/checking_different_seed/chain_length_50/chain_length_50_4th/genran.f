c     Version of Dec 7, 1999. Corrected two bugs (one in main,one in ebfindm)
c     version of May 25, 2006. Just modified slightly for generating
c     chains in the slit. 

      Integer chains(1000,2000),kar(500,500,200),empty(500,500,200)
      integer enpnew,ifail(2000),weight(2000)
      real etime_struct(2),elapsed
      integer IX,IX1
      COMMON/D1/LSR/D2/LSC/D3/LSL
      COMMON/B/KAR
      common/c/empty
      common/xxx/ix,ix1
      parameter (ifailmax=10)
      External Time
      open(unit=11,file='initial.dat')

c       call time(ix)
      ix=556066
      write(*,*)'ix',ix      

      write(6,*) 'what is the lsx, lsy, maximum is 250'
      read(5,*) lsr,lsc
      write(6,*) 'what is lsl, lsl is the layers in z maximum 200'
      read(5,*) lsl

      write(6,*) 'what is the number of  chains maximum 1000'
      read(5,*) nch
      write(6,*) 'what is the length of the chain A,maximum 2000'
      read(5,*) lch

      itrial=float((nch*lch))/float(lsr*lsc*lsl)*1000000

        do 1001 i=1,lsr
          do 1001 j=1,lsc
             do 1002 k=1,lsl
               kar(i,j,k)=(k-1)*lsc*lsr+(j-1)*lsr+i
               empty(i,j,k)=0
               if(k.eq.1) empty(i,j,k)=5
               if(k.eq.lsl) empty(i,j,k)=5
 1002          continue
 1001          continue


       ix1=ix+50000
       ixseed=ix
       ix1seed=ix1
       write(6,*) 'the seed is',ixseed,ix1seed

       do 10 i=1,nch

          do 11 k=1,lch
             weight(k)=0
             ifail(k)=0
 11       continue

          ifail2=0
 20       j=1
          call rndu(float(lsr),ipk1)
          call rndu(float(lsc),ipk2)
          call rndu(float(lsl),ipk3)
          if(empty(ipk1,ipk2,ipk3).eq.0) then
             chains(i,1)=kar(ipk1,ipk2,ipk3)
             empty(ipk1,ipk2,ipk3)=1
          else
             ifail2=ifail2+1
             if(ifail2.gt.itrial) then
                write(6,*) 'too many failed times'
                goto 100
             else   
                goto 20
             end if
          end if

 40       j=j+1
          if(j.gt.lch) goto 10

          call endfindm(chains(i,j-1),enpnew,iempty)
          if(enpnew.ne.0) then
             weight(j)=iempty
             chains(i,j)=enpnew
             call findm (enpnew,l1,l2,l3)
             empty(l1,l2,l3)=1
             goto 40
          else
           ifail(j)=ifail(j)+1
           if(j.eq.2) goto 43
           if(ifail(j).gt.ifailmax) goto 43
           jsub=1
 41        if(jsub.ge.(j-1)) goto 43
           if(weight(j-jsub).le.1) then
              jsub=jsub+1
              goto 41
           else if(weight(j-jsub).gt.1) then
              do 42 k=j-jsub,j-1
                 ifail(k)=ifail(k)+1
                 call findm(chains(i,j-1),l1,l2,l3)
                 empty(l1,l2,l3)=0
                 weight(k)=0
 42           continue
              j=j-jsub-1
              goto 40
           end if
        end if
 43     ifail2=ifail2+1
        if(ifail2.gt.itrial) then
           write(6,*) 'failed in inserting chain',i
           goto 100
        else
           do 44 k=1,j-1
              call findm(chains(i,k),l1,l2,l3)
              empty(l1,l2,l3)=0
              weight(k)=0
              ifail(k)=0
 44        continue
           goto 20
        end if      
 10   continue



 200          write(6,*) 'successful'
              write(6,*) 'what surface interaction you want'
              read(5,*) eaw
              write(11,*) eaw,nch,lch
              write(11,*) lsr,lsc, lsl
              do I=1,nch
                 do j=1,lch
                    call findm(chains(i,j),I1,I2,I3)
                    write(11,*) I1,I2,I3
                 end do
              end do
              goto 111

 100              write(6,*) 'failed too many times','ifail2',ifail2
                  write(6,*) 'i,j',i,j
 111        elapsed=etime(etime_struct)
            write(6,*) 'total CPU used',elapsed,'seconds'

          stop
                  end



      SUBROUTINE RNDU(X,IRN)
C
C     THIS SUBROUTINE GENERATES A RANDOM NUMBER
C     BETWEEN 1 AND INT(X)
C     ***********************
      REAL RN
      integer d
      parameter(m=992631493, d=2147483647)
      COMMON/XXX/IX,IX1
      IX1=IX1*m
      rn=mod(abs(ix1),d)
      rn=rn/float(d)
      IRN=IFIX(RN*X)+1
      if(irn.gt.ifix(x)) then
      IRN=IRN-1
      END IF
      RETURN
      END

C     ***************************************************
      SUBROUTINE RND(RN)

C     THIS SUBROUTINE GENERATES A RANDOM NUMBER VIZ. RN
C     BETWEEN 0 AND 1
c     *******************************
      REAL RN
      integer d
      parameter(m=923456789, d=2147483647)
      COMMON/XXX/IX,IX1
      IX=IX*m
      rn=mod(abs(ix),d)
      rn=rn/float(d)
      RETURN
      END

 
      SUBROUTINE FINDM (NB,IR,IC,IL)
C     ***********************************************
C     THIS MODIFIED SUBROUTINE FINDS THE COORDINATES OF THE
C     POINT NB BY DIVIDING IT WITH LSR. THE DIVISOR IS THE
C     COLUMN AND THE REMAINDER IS THE ROW, EXCEPT FOR THE
C     ROW LSR.
 
      COMMON/D1/LSR/D2/LSC/D3/LSL

C     FIRST FIND THE LAYER #
 
                   LPROD=LSR*LSC
 
         IF (MOD(NB,LPROD).EQ.0) THEN
            IL=NB/LPROD
         ELSE
            IL=(NB/LPROD)+1
         END IF
 
C        NOW PROJECT NB ON TO THE FIRST LAYER AND PROCEED AS IN 2-D
                   NBB=NB-(LPROD*(IL-1))
 
          IF (MOD(NBB,LSR).EQ.0) THEN
             IR=LSR
             IC=NBB/LSR
          RETURN
          ELSE
              IC=NBB/LSR
              IC=IC+1
              IR=MOD(NBB,LSR)
          END IF
 
      RETURN
      END
C     ***********************


 
       SUBROUTINE N4MOD (NB,NN)
       DIMENSION NN(6),ndum(3,2),kar(500,500,200) 
       COMMON/B/KAR
       COMMON/D1/LSR/D2/LSC/D3/LSL
C     ***************************************************
c     This subroutine gives six nearest neightbor sites. ndum array 
c     temporary array, which stores the i,j,k three coordinate 
c     varaible.

         call findm(nb,i,j,k)

          ndum(1,1)=i-1
          ndum(1,2)=i+1
          ndum(2,1)=j-1
          ndum(2,2)=j+1
          ndum(3,1)=k-1
          ndum(3,2)=k+1

          if(ndum(2,1).eq.0) ndum(2,1)=lsc
          if(ndum(2,2).eq.lsc+1) ndum(2,2)=1

          if(ndum(1,1).eq.0) ndum(1,1)=lsr
          if(ndum(1,2).eq.lsr+1) ndum(1,2)=1

              if(ndum(3,1).eq.0) ndum(3,1)=lsl
              if(ndum(3,2).eq.lsl+1) ndum(3,2)=1

        nn(1)=kar(ndum(1,1),j,k)
        nn(2)=kar(ndum(1,2),j,k)
        nn(3)=kar(i,ndum(2,2),k)
        nn(4)=kar(i,ndum(2,1),k)
        nn(5)=kar(i,j,ndum(3,1))
        nn(6)=kar(i,j,ndum(3,2))
 
      RETURN
      END


       SUBROUTINE ENDFINDM (ENPOLD,ENPNEW,iempty)
C           GIVEN ENPOLD THIS SUB RANDOMLY CHOOSES ONE OF THE VOIDS
C     AROUND IT I.E ENPNEW
 
C  **************************************
 
      INTEGER ENPOLD,ENPNEW,EMPTY(500,500,200)
      DIMENSION NN(6),NDK(6)
      COMMON/C/EMPTY

       call n4mod(enpold,nn)
       call findm(enpold,l1,l2,l3)

       iempty=0
       do 151 k=1,6
          call findm(nn(k),n1,n2,n3)
          if(empty(n1,n2,n3).eq.0) then
             iempty=iempty+1
             ndk(iempty)=nn(k)
           end if
 151       continue

          if(iempty.ne.0) then
          call rndu(float(iempty),ipk)
          enpnew=ndk(ipk)
          else
             enpnew=0
           end if  

             return
             end
 
 
