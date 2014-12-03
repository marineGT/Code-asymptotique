      subroutine resolutionsystLU(Nt,Nx,Ng,precision,lsa,c,c2,dx,
     &           zetap,lambda0,gammahat,tabsB,tabsB0,alpha,
     &           alpha0,Bn,Dn,tabx,zeta1B,Rn,R0,tabsigma0zz)
c     ==================================================================
c     Aucune declaration implicite
c     ----------------------------
      implicit none
c     ==================================================================
c     Declaration des parametres
c     --------------------------

c     Declaration des variables en entree
c     -----------------------------------
      integer, intent(in) :: Nt, Nx, Ng
      double precision, intent(in) :: precision, lsa
      double precision, intent(in) :: c, c2, dx, zetap 
      double precision, intent(in) :: lambda0, gammahat
      double precision, intent(in) :: tabsB0, alpha0,R0      
      double precision, dimension(Nt) :: alpha,Rn
      double precision, dimension(Nt+1) :: Bn,Dn      
      double precision, dimension(Nt+1) :: tabsB

c     Declaration des variables en sortie
c     -----------------------------------
      double precision, dimension(0:Nx) :: tabx      
      double precision, dimension(0:Nx) :: zeta1B
      double precision, dimension(Nx+1) :: tabsigma0zz

c     Declaration des variables locales
c     ---------------------------------
      integer :: i, n, m, methode1, methode2
      integer :: k, INFO
      double precision :: eta, deta, invdeta
      double precision :: teta
      double precision :: invdeta2, pi, ch, sh
      double precision :: sum,sum2,aux,aux2
      double precision :: Sn, xn, res, xi, A      
      double precision :: prec, e, f, F0verif
      double precision, dimension(Nt+3) :: Pn
      double precision, dimension(Nt+1) :: Un, dUn
      character*1 :: TRANS, NORM
      double precision,dimension (Nx+1) :: IPIV
      double precision,dimension (3) ::coef,IPIV1
      double precision,dimension (Nx+1,Nx+1) :: amat
      double precision,dimension (Nx+1,1) :: yse
      double precision,dimension(7) :: x,wter
      double precision,dimension (Nx+1) :: Am
      external dgetrf,dgetrs
 
c     Declaration des interfaces
c     --------------------------
      interface
         function secondmembreB(Nt,zetap,c,c2,precision,lambda0,
     &                gammahat,tabsB,tabsB0,xi,Rn,R0)
            integer,intent(in) :: Nt 
            double precision, intent(in) :: c, c2, precision
            double precision, intent(in) :: xi, zetap
            double precision, intent(in) :: lambda0, tabsB0
            double precision, intent(in) :: gammahat, R0
            double precision, dimension(Nt+1),intent(in) :: tabsB
            double precision, dimension(Nt),intent(in) :: Rn
	    double precision :: secondmembreB            
         end function secondmembreB

        function sigma0zzB(Nt,c,zetap,xi,alpha,alpha0,Pn,
     &                      Un,dUn) 
            integer, intent(in)  ::  Nt
            double precision, intent(in) :: c, zetap 
            double precision, intent(in) :: xi, alpha0
            double precision, dimension(Nt), intent(in) :: alpha           
            double precision, dimension(Nt+3), intent(in) :: Pn
            double precision, dimension(Nt+1), intent(in) :: Un  
            double precision, dimension(Nt+1), intent(in) :: dUn
            double precision :: sigma0zzB 
         end function sigma0zzB

        function sigma0zzBCherZap(Nt,c,zetap,xi,alpha,alpha0,
     &                          Pn,Un,dUn) 
            integer, intent(in)  ::  Nt
            double precision, intent(in) :: c, zetap 
            double precision, intent(in) :: xi, alpha0
            double precision, dimension(Nt), intent(in) :: alpha           
            double precision, dimension(Nt+3), intent(in) :: Pn
            double precision, dimension(Nt+1), intent(in) :: Un  
            double precision, dimension(Nt+1), intent(in) :: dUn
            double precision :: sigma0zzBCherZap
         end function sigma0zzBCherZap

        subroutine Legendre(n,x,resultat)
           integer :: n
           double precision :: x,resultat           
       end subroutine Legendre

       function verifsigma0zzB(teta,Nt,c,zetap,alpha,alpha0,
     &                         Bn,Dn)      
         integer, intent(in) ::  Nt      
         double precision, intent(in) :: teta
         double precision, intent(in) ::  c, zetap
         double precision, dimension(Nt) :: alpha
         double precision, intent(in) :: alpha0 
         double precision, dimension(Nt+1) :: Bn,Dn
         double precision :: verifsigma0zzB
        end function verifsigma0zzB

      end interface 

      double precision :: gaussint1dim2
      external gaussint1dim2

c2345678912345678912345678912345678912345678912345678912345678912   
c     ==================================================================
c     Soit l'equation differentielle a resoudre sur la bulle

      ! sin(eta)*(ch-x)*F" - [cos(eta)*(ch-x)+2*(1-x*ch)]*F'
      !     + 2*ch*sin(eta)*F = sin(eta)*( S(x) + (Bo/Ca)*z' )

      ! Condition aux limites

c      F' = 0 quand eta = 0 ou eta = pi
c     ==================================================================
c     calcul de pi
c     ------------
      pi=2.d0*dasin(1.d0)

c     Calcul de deta et de ses inverses 
c     ---------------------------------
      deta=pi/dble(Nx)
      invdeta=1.D0/deta
      invdeta2=invdeta*invdeta     
      ch=dcosh(zetap)
      sh=dsinh(zetap)

c     On construit les tableaux Pn, un et dUn pour la fonction sigma0zzB

c     Construction du tableau Pn,Un et dUn 
c     ------------------------------------

      !initialisation
      do 100 k=1,Nt+3
       Pn(k)=0.d0
 100   continue

      do 200 k=1,Nt+1  
       Un(k)=0.d0
       dUn(k)=0.d0
 200  continue

c     Construction du tableau Pn, Un et dUn 
c     -------------------------------------
c     On sait que Pn0=1.D0  
c      open(unit=11,file='yse_Sn-tot.dat')

      do 10 n=1,Nt+1 
        xn=dble(n)      
        Un(n)=+Bn(n)*dsinh(-(xn-5.d-1)*zetap)
     &     +Dn(n)*dsinh(-(xn+(3.d0/2.d0))*zetap)
        dUn(n)=+(xn-5.d-1)*Bn(n)*dcosh((xn-5.d-1)*zetap)
     &         +(xn+(3.d0/2.d0))*Dn(n)
     &         *dcosh((xn+(3.d0/2.d0))*zetap)    
c        write(16,*) n, Un(n), dUn(n)
c        dUn(n)=0.d0
 10   continue  
      
c     Constitution du second membre yse(n,1) pour le systeme
c     lineaire des Hn.

c     ATTENTION selon le cas/conditions choisi(e)s, on doit penser à
c     modifier la matrice amat aussi
      do 1 n=1,Nx+1
         yse(n,1)=0.d0
1     continue    

c     la première ligne de yse = 0 (CL)
c     --------------------------------
      yse(1,1)=0.d0
      write(*,*) 'yse(1,1)',yse(1,1)

c     valeur de yse pour n=2,Nx
c     -------------------------
      eta=deta
      do 2 n=2,Nx
c$$$        xi=-1.D0+dx*dble(n-1)        
c$$$        do i=1,Nt+3 
c$$$          xn=dble(i)       
c$$$          call Legendre(i,xi,res)
c$$$          Pn(i)=res
c$$$        enddo
c        Sn=secondmembreB(Nt,zetap,c,c2,precision,lambda0,
c     &                gammahat,tabsB,tabsB0,xi,Rn,R0)
c$$$        Sn=sigma0zzB(Nt,c,zetap,xi,alpha,alpha0,Pn,
c$$$     &              Un,dUn)
c$$$        Sn=sigma0zzBCherZap(Nt,c,zetap,xi,alpha,alpha0,
c$$$     &                      Pn,Un,dUn) 
c$$$       aux=3.D0*lambda0*(c*sh/(ch-xi)) 
c$$$
c$$$        yse(n,1)=-Sn*c+aux*c
c$$$        write(11,*) xi, yse(n,1),Sn*c,aux*c,
c$$$     &              (ch-xi),sh/(ch-xi)    
        eta=dble(n-1)*deta
        xi=dcos(eta)   
c        xi=(1.d0-dcos(eta)*ch)/(ch-cos(eta))
        do i=1,Nt+3 
          xn=dble(i)       
          call Legendre(i,xi,res)
          Pn(i)=res
c          write(16,*) xi, i, Pn(i)
        enddo
        tabsigma0zz(n)=sigma0zzB(Nt,c,zetap,xi,alpha,
     &                  alpha0,Pn,Un,dUn)
        aux=3.D0*lambda0*(1.d0-dcos(eta)*ch)
     &           /(ch-dcos(eta))
        yse(n,1)=dsin(eta)*c*(aux+tabsigma0zz(n))

c        write(11,*) eta, xi, yse(n,1),tabsigma0zz(n),
c     &              aux 
 2    continue  
c      close(11)
      
c     valeur de yse pour n=Nx+1
cc     -------------------------  
      yse(Nx+1,1)=0.d0

c     fin  constitution du second membre  
c     ___________________________________________________________
c     --------constitution de la matrice amat du systeme---------

c     initialisation de la matrice
c     ------------------------------
      do 3 n=1,Nx+1
       do 4 m=1,Nx+1
        amat(n,m)=0.d0
 4     continue
 3    continue

c     constitution de la premiere ligne de amat en eta=0
c     --------------------------------------------------
c      F' = 0 en eta = 0

       amat(1,1)=-3.d0/(2.d0*deta)
       amat(1,2)=2.d0*invdeta
       amat(1,3)=-(1.d0/2.d0)*invdeta     
c$$$       write(*,*)'amat(1,1)=',amat(1,1)
c$$$       write(*,*)'amat(1,2)=',amat(1,2)
c$$$       write(*,*)'amat(1,3)=',amat(1,3)         

c     constitution des lignes 2 à la ligne Nx de amat
c     -----------------------------------------------
c     On resout ici l'eq differentielle    
c     fn+1.[(an/h²) + (bn/2h)] + fn.[(-2an/h²) + cn] 
c       + fn-1.[(an/h²) - (bn/2h)] = dn   

      eta=deta
      do 6 n=2,Nx
       eta=dble(n-1)*deta
       xi=dcos(eta)
       amat(n,n-1)=dsin(eta)*(ch-xi)*invdeta2
     &           -(xi*(ch-xi)+2.d0*(1.d0-ch*xi))
     &            *invdeta/2.d0            
       amat(n,n)=-2.d0*dsin(eta)*(ch-xi)*invdeta2
     &           +2.d0*ch*dsin(eta)
       amat(n,n+1)=dsin(eta)*(ch-xi)*invdeta2
     &           +(xi*(ch-xi)+2.d0*(1.d0-ch*xi))
     &            *invdeta/2.d0    
6     continue

c     constitution de la derniere ligne de amat en eta=pi
c     --------------------------------------------------
c           F' = 0 en eta = pi

       amat(Nx+1,Nx-1)=(1.d0/2.d0)*invdeta
       amat(Nx+1,Nx)=-2.d0*invdeta      
       amat(Nx+1,Nx+1)=3.d0/(2.d0*deta)
c$$$       write(*,*)
c$$$       write(*,*)'amat(Nx+1,Nx+1)=',amat(Nx+1,Nx+1)
c$$$       write(*,*)'amat(Nx+1,Nx)=',amat(Nx+1,Nx)
c$$$       write(*,*)'amat(Nx+1,Nx-1)=',amat(Nx+1,Nx-1)


c$$$       do n=2,Nx
c$$$          aux=amat(n,n-1)+amat(n,n)+amat(n,n+1)
c$$$          write(17,*) n, aux
c$$$       enddo


c$$$c     Test du second membre
c$$$         xi=-1.d0
c$$$         sum1= yse(1,1)/((ch-xi)**3)
c$$$         do 63 n=2,Nx
c$$$            xi=-1.D0+dx*dble(n-1)
c$$$            sum1=sum1+2.d0* yse(n,1)/((ch-xi)**3)
c$$$63      continue
c$$$        xi=1.d0
c$$$        sum1=sum1+ yse(Nx,1)/((ch-xi)**3)
c$$$        sum1=sum1*dx/(2.d0)
c$$$        write(*,*) 'test intégrale second membre', sum1
c$$$c     Fin du test du second membre

c     fin constitution de la matrice amat du systeme

c     ___________________________________________________________
c     ----factorisation LU de la matrice amat--------------------

c     definition des parametres pour utiliser les subsections 
c     -----------------------------------------------------
      INFO=12
      TRANS='N'
      NORM='1'

c     Appel du sous programme dgetrf
c     ------------------------------
      call dgetrf(Nx+1,Nx+1,amat,Nx+1,IPIV,INFO)
c      write(*,*) 'INFO',INFO

c     resolution du systeme, Appel du sous programme dgetrs 
c     ----------------------------------------------------
      call dgetrs(TRANS,Nx+1,1,amat,Nx+1,IPIV,yse,Nx+1,INFO)
c      write(*,*) 'INFO',INFO

c     Calcul du Determinant
c     ---------------------
      aux=0.d0
      do 777 n=1,Nx+1
         aux=aux+amat(n,n)         
 777     continue
c         write(*,*) 'det matrice=',aux

c     ecriture des resultats dans le tableau zeta1B(n)
c     ---------------------------------------------
      zeta1B(0)=yse(1,1)  
      do n=1,Nx
       zeta1B(n)=yse(n+1,1)
       eta=dble(n)*deta       
       xi=dcos(eta)
      enddo
      
c      write (*,*) 'zeta1B(0)=',zeta1B(0)
c      write (*,*) 'zeta1B(Nx)=',zeta1B(Nx)  

c     ==============================================
      !Test sur la matrice


c     Test de la conservation du volume
c     ---------------------------------
c     methode des trapezes      
         sum2=0.d0
         do n=2,Nx+1
          eta=dble(n-1)*deta
          xi=dcos(eta)
          aux2=(ch-xi)**3.d0          
          sum2=sum2+(1.d0-xi*ch)*dsin(eta)/aux2 
       enddo
        sum2=sum2*deta
c        write(*,*) 'test conservation volume', sum2 

c     Calcul conservation barycentre pour la solution homogene
c     ---------------------------------------------------------
      !Soit int_{-1}^{1} {(1-x*ch)*sin(eta)*deta/(ch-x)}

      sum2=0.d0
      do n=2,Nx+1 
        eta=dble(n-1)*deta
        xi=dcos(eta)
        aux2=(ch-xi)**4.d0
        sum2=sum2+(1.d0-xi*ch)*dsin(eta)/aux2
      enddo
      sum2=sum2*deta
c      write(*,*) 'integrale CM pour Sol homogene', sum2 
c      write(*,*)
c      write(*,*) 'res analytique = ',
c     &               -2.d0/(3.d0*(sh**4.d0)) 

c     Calcul du coeff A de la solution homogene
c     ------------------------------------------
      !Soit A=3*(sh**4)/2*int_{-1}^{1} {Fp*sin(eta)*deta/(ch-x)}

      sum2=0.d0
      do n=2,Nx+1 
        eta=dble(n-1)*deta
        xi=dcos(eta)
        aux2=(ch-xi)**4.d0
        sum2=sum2+zeta1B(n-1)*dsin(eta)/aux2
      enddo
      sum2=sum2*deta
c      write(*,*) 'integrale volume pour Fp', sum2 
      A=(3.d0/2.d0)*(sh**4.d0)*sum2
c      write(*,*) 'le coeff A =', A  

c     Test du calcul de coeff A
c     -------------------------
c     methode des trapezes       
         sum2=0.d0
         do n=2,Nx 
           eta=dble(n-1)*deta
           xi=dcos(eta)
           aux2=(ch-xi)**4.d0
           sum2=sum2+A*(1.d0-xi*ch)*dsin(eta)/aux2
         enddo
        sum2=sum2*deta
c$$$        write(*,*)
c$$$        write(*,*) 'test integral avec A', sum2 
c$$$        aux=-(2.d0*A)/(3.d0*(sh**4.d0))
c$$$        write(*,*) 'integral analytique', aux

c     Calcul de la solution totale F(eta)
c     -----------------------------------
      ! F(eta)= Fp + A*[1 - cos(eta)*ch]

        do  n=1,Nx+1 
            eta=dble(n-1)*deta
            xi=dcos(eta)            
            yse(n,1)=zeta1B(n-1)+A*(1.d0-xi*ch)                             
        enddo      

c     Test de la conservation du centre de masse
c     ------------------------------------------
c     methode des trapezes       
        sum2=0.d0
         do n=2,Nx 
           eta=dble(n-1)*deta
           xi=dcos(eta)
           aux2=(ch-xi)**4.d0
           sum2=sum2+yse(n,1)*dsin(eta)/aux2
         enddo
        sum2=sum2*deta
c        write(*,*)
c        write(*,*) 'test barycentre', sum2 

c     Calcul de F0 avec sigma0zz en spherique
c     ---------------------------------------
c$$$      sum2=0.d0
c$$$       do n=2,Nx 
c$$$         eta=dble(n-1)*deta
c$$$         aux2=tabsigma0zz(n)*dcos(eta)
c$$$c         aux2=dcos(eta)*dcos(eta) 
c$$$         sum2=sum2+aux2*dsin(eta)
c$$$       enddo
c$$$       sum2=sum2*deta
c$$$       write(*,*)
c$$$       write(*,*) 'Calcul de F0 avec sigma0zz', sum2 

c     On calcule l'integrale de sigma0zzB
c     -----------------------------------  
c      prec=1.d-5 ! precision de l'integrale
c      e=0.d0
c      f=pi
c      F0verif=gaussint1dim2(verifsigma0zzB,e,f,
c     &   Nt,c,zetap,alpha,alpha0,Bn,Dn,prec)     

c$$$      teta=pi/4.d0
c$$$      F0verif=verifsigma0zzB(teta,Nt,c,zetap,
c$$$     &            alpha,alpha0,Bn,Dn)
      
c      write(*,*) 'Calcul de F0verif', F0verif
c      write(*,*) '2*lambda0', 2.d0*lambda0
      
c     Solution totale zeta1B(n)
c     -------------------------
       do  n=1,Nx+1 
            eta=dble(n-1)*deta
            xi=dcos(eta)            
            zeta1B(n-1)=yse(n,1)                        
       enddo         

c     =======================================
c     Fin du sous-programme diagonalessysteme1
c     =======================================
      end subroutine resolutionsystLU
