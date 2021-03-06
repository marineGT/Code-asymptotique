c****************************************************************
c       
c 	Gauss adaptative avec 24 - points pour 
c 	l'int�gration d'une fonction 'F' dans 
c 	l'intervalle [a,b]         
c     
c
c     Appel�e par le sous programme intepoly pour le calcul
c     de l'integrale de la force
c
c     Programmeur: M.Guemas; Fevrier 2014.
c****************************************************************

      function gaussint1dim2(verifsigma0zzB,a,b,Nt,c,zetap,
     &                       alpha,alpha0,Bn,Dn,prec)
c     ============================================================  
c     Declaration des variables 
c     -------------------------
      implicit double precision (a-h,o-z)
      integer Nt,i
      character name*(*)
      parameter (name = 'gaussint1dim2')
      dimension wter(12),x(12)
      double precision x,c1,c2,prec      

c$$$c     Declaration des interfaces
c$$$c     ---------------------------
c$$$      interface    
c$$$        function verifsigma0zzB(teta,Nt,c,zetap,alpha,
c$$$     &                          alpha0,Bn,Dn)
c$$$         integer, intent(in) :: Nt      
c$$$         double precision, intent(in) :: teta
c$$$         double precision, intent(in) :: c, zetap
c$$$         double precision, dimension(Nt) :: alpha
c$$$         double precision, intent(in) :: alpha0 
c$$$         double precision, dimension(Nt+1) :: Bn,Dn
c$$$         double precision ::  verifsigma0zzB
c$$$        end function verifsigma0zzB
c$$$      end interface
      double precision :: verifsigma0zzB
      double precision :: c, zetap      
      double precision, dimension(Nt) :: alpha
      double precision :: alpha0 
      double precision, dimension(Nt+1) :: Bn,Dn
      external  verifsigma0zzB

c     Declaration des parametres
c     -------------------------
      parameter (z1ter=1.d0,hf=5.d-1,cst=5.d-3)
      data x
     #        /0.96028 98564 97536 23168 35608 68569 47D0,
     #         0.79666 64774 13626 73959 15539 36475 83D0,
     #         0.52553 24099 16328 98581 77390 49189 25D0,
     #         0.18343 46424 95649 80493 94761 42360 18D0,
     #         0.98940 09349 91649 93259 61541 73450 33D0,
     #         0.94457 50230 73232 57607 79884 15534 61D0,
     #         0.86563 12023 87831 74388 04678 97712 39D0,
     #         0.75540 44083 55003 03389 51011 94847 44D0,
     #         0.61787 62444 02643 74844 66717 64048 79D0,
     #         0.45801 67776 57227 38634 24194 42983 58D0,
     #         0.28160 35507 79258 91323 04605 01460 50D0,
     #         0.95012 50983 76374 40185 31933 54249 58D-1/
c2345678912345678912345678912345678912345678912345678912345678912
      data wter
     #        /0.10122 85362 90376 25915 25313 54309 96D0,
     #         0.22238 10344 53374 47054 43559 94426 24D0,
     #         0.31370 66458 77887 28733 79622 01986 60D0,
     #         0.36268 37833 78361 98296 51504 49277 20D0,
     #         0.27152 45941 17540 94851 78057 24560 18D-1,
     #         0.62253 52393 86478 92862 84383 69943 78D-1,
     #         0.95158 51168 24927 84809 92510 76022 46D-1,
     #         0.12462 89712 55533 87205 24762 82192 02D0,
     #         0.14959 59888 16576 73208 15017 30547 48D0,
     #         0.16915 65193 95002 53818 93120 79030 36D0,
     #         0.18260 34150 44923 58886 67636 67969 22D0,
     #         0.18945 06104 55068 49628 53967 23208 28D0/
c2345678912345678912345678912345678912345678912345678912345678912
c     ==================================================================  
      eps=prec
c      eps=1.d-14
      h=0
c      WRITE(*,*) 'debut gaussint2, c= ', c
      if(b .eq. a) go to 99
      const=cst/abs(b-a)
      bb=a
 1    aa=bb
      bb=b
 2    c1=hf*(bb+aa)
      c2=hf*(bb-aa)      
      s8=0
      do 3 i = 1,4
      u=c2*x(i)
      s8=s8+wter(i)*(verifsigma0zzB(c1+u,Nt,c,zetap,alpha,
     & alpha0,Bn,Dn)+verifsigma0zzB(c1-u,Nt,c,zetap,alpha,
     & alpha0,Bn,Dn))
 3    continue
      s16=0
      do 4 i = 5,12
      u=c2*x(i)
      s16=s16+wter(i)*(verifsigma0zzB(c1+u,Nt,c,zetap,alpha,
     & alpha0,Bn,Dn)+verifsigma0zzB(c1-u,Nt,c,zetap,alpha,
     & alpha0,Bn,Dn))
 4    continue
      s16=c2*s16
      if(abs(s16-c2*s8) .le. eps*(1+abs(s16))) then
c       write (*,*) 'passage'
       h=h+S16       
       if(bb .ne. b) go to 1
      else
       bb=c1
       if(1+const*abs(c2) .ne. 1) go to 2
       h=0
       write(*,*) name,'la precision demandee est trop grande'
       go to 99
      end if
 99   gaussint1dim2=h
c      WRITE(*,*) 'FIN gaussint, c =', c
      return
      END
c****************************************************************
c2345678912345678912345678912345678912345678912345678912345678912
