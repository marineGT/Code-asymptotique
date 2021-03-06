      subroutine resolutionsystLU(Nt,Nx,Ng,precision,lsa,c,c2,dx,
     &            zetap,lambda0,gammahat,tabsB,tabsB0,alpha,alpha0,
     &            Bn,Dn,tabx,zeta1B,Rn,R0)
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
      double precision,dimension(Nt) :: alpha,Rn
      double precision, dimension(Nt+1) :: Bn,Dn      
      double precision, dimension(Nt+1) :: tabsB

c     Declaration des variables en sortie
c     -----------------------------------
      double precision, dimension(0:Nx) :: tabx      
      double precision, dimension(0:Nx) :: zeta1B

c     Declaration des variables locales
c     ---------------------------------
      integer :: i, n, m, methode1, methode2, INFO
      double precision :: invdx,invdx2,dx2,xi
      double precision :: ch,grxi,sum,aux,sum1
      double precision :: Sn, sh, xn, res
      integer :: k
      double precision, dimension(Nt+3) :: Pn
      double precision, dimension(Nt+1) :: Un, dUn    
      double precision,dimension(7) :: x,wter
      character*1 :: TRANS, NORM
      double precision,dimension (Nx+1) :: IPIV
      double precision,dimension (3) ::coef,IPIV1
      double precision,dimension (Nx+1,Nx+1) :: amat
      double precision,dimension (Nx+1,1) :: yse
      double precision,dimension (Nx+1) :: Am
      external dgetrf,dgetrs
 
c     Declaration des interfaces
c     --------------------------
      interface
         function secondmembreB(Nt,zetap,c,c2,precision,lambda0,
     &                gammahat,tabsB,tabsB0,xi,Rn,R0)
            integer,intent(in) :: Nt 
            double precision, intent(in) :: c, c2, precision, xi,zetap 
            double precision, intent(in) :: lambda0, tabsB0, gammahat
            double precision, intent(in) :: R0
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

      end interface     
 
c      if (Ng.eq.7) then      
      data x
     #        / 0.0 D0,
     #          0.40584 51513 77397 16690 D0, 
     #         -0.40584 51513 77397 16690 D0,
     #          0.74153 11855 99394 43986 D0,
     #         -0.74153 11855 99394 43986 D0,
     #          0.93246 95142 03152 02781 D0, 
     #         -0.93246 95142 03152 02781 D0/
c2345678912345678912345678912345678912345678912345678912345678912
      data wter
     #        / 0.41795 91836 73469 38775 D0,
     #          0.38183 00505 05118 94495 D0,
     #          0.38183 00505 05118 94495 D0,
     #          0.27970 53914 89276 66790 D0,
     #          0.27970 53914 89276 66790 D0, 
     #          0.12948 49661 68869 69327 D0,  
     #          0.12948 49661 68869 69327 D0/
c     if (Ng.eq.4) then     
c$$$      data x
c$$$     #        / 0.33998 10435 84856 26480 D0,
c$$$     #         -0.33998 10435 84856 26480 D0, 
c$$$     #          0.86113 63115 94052 57522 D0,    
c$$$     #         -0.86113 63115 94052 57522 D0/
c$$$c2345678912345678912345678912345678912345678912345678912345678912
c$$$      data wter
c$$$     #        / 0.65214 51548 62546 14262 D0,
c$$$     #          0.65214 51548 62546 14262 D0,
c$$$     #          0.34785 48451 37453 85737 D0,  
c$$$     #          0.34785 48451 37453 85737 D0/
      
c     
c     ==================================================================

c     Calcul des inverses de dx et dx2
c     --------------------------------
      invdx=1.D0/dx
      invdx2=invdx*invdx
      ch=dcosh(zetap)
      sh=-dsinh(zetap)

c     Condition limites sur l'ensemble de la fonction zeta1
c     -----------------------------------------------------
c     le premiere methode est celle des trapezes :
c        int(-1,1){f(x)dx}=h*sum(m=1,M){f(m-1)+f(m)}      
 
c     la deuxieme est celle des points de Gauss:
c         int(-1,1){f(x)dx}=sum(i=1,N){w(i)*f(ni)}
      
c     Pour nous f(x)=zetaB1/[ch(zetap)-xi]
      methode1=0
      methode2=0
      write(*,*) 'methode1=',methode1
      write(*,*) 'methode2=',methode2
c     ******************************************************************

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
      open(unit=11,file='yse_n-1_test-zetap.dat')

      do 10 n=1,Nt+1 
        xn=dble(n)      
        Un(n)=+Bn(n)*dsinh(-(xn-5.d-1)*zetap)
     &     +Dn(n)*dsinh(-(xn+(3.d0/2.d0))*zetap)
        dUn(n)=+(xn-5.d-1)*Bn(n)*dcosh((xn-5.d-1)*zetap)
     &         +(xn+(3.d0/2.d0))*Dn(n)
     &         *dcosh((xn+(3.d0/2.d0))*zetap)        
c        dUn(n)=0.d0
 10   continue  
      
c     Constitution du second membre yse(n,1) pour le systeme
c     lineaire des Hn.

c     ATTENTION selon le cas/conditions choisi(e)s, on doit penser �
c     modifier la matrice amat aussi
      do 1 n=1,Nx+1
         yse(n,1)=0.d0
1     continue    

c     la premi�re ligne de yse = 0 (CL)
c     --------------------------------
      yse(1,1)=0.d0
      write(*,*) 'yse(1,1)',yse(1,1)

c     valeur de yse pour n=2,Nx
c     -------------------------
      
      do 2 n=2,Nx
        xi=-1.D0+dx*dble(n-1)        
        do i=1,Nt+3 
          xn=dble(i)       
          call Legendre(i,xi,res)
          Pn(i)=res
        enddo
c        Sn=secondmembreB(Nt,zetap,c,c2,precision,lambda0,
c     &                gammahat,tabsB,tabsB0,xi,Rn,R0)
        Sn=sigma0zzB(Nt,c,zetap,xi,alpha,alpha0,Pn,
     &              Un,dUn)
c        Sn=sigma0zzBCherZap(Nt,c,zetap,xi,alpha,alpha0,
c     &                      Pn,Un,dUn) 
        aux=3.D0*lambda0*(c*sh/(ch-xi)) 

        yse(n,1)=-Sn*c+aux*c

        write(11,*) xi, yse(n,1),Sn*c,aux*c,
     &              (ch-xi),sh/(ch-xi)                
  
 2    continue  
      close(11)
      
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

c     constitution de la premiere ligne de amat qui est 
c     la condition limite sur le volume 
c     -----------------------------------------

c     methode des trapezes
      if (methode1.eq.1) then
         dx2=dx*dx
         amat(1,1)=dx2/((1.d0+1.d0/ch)**3)
         do 5 n=2,Nx 
            xi=-1.D0+dx*dble(n-1)
            tabx(n-1)=xi            
            amat(1,n)=2.d0*dx2/((1.d0-xi/ch)**3)
 5      continue
        amat(1,Nx+1)=dx2/((1.d0-1.d0/ch)**3)
      endif
c     Fin methode des trap�zes

c     Methode des points de gauss
c     Calcul des coefficiant Am(m)    
         do 50 n=1,Nx+1
            Am(n)=0.d0
 50      continue 
         do 51 m=1,Nx
            xi=-1.D0+dx*dble(m-1)
            tabx(m-1)=xi
            do 511 i=1,Ng              
              if(m.ne.1) then 
                grxi=(x(i)+1.d0)/2.d0
                Am(m-1)=Am(m-1)-wter(i)*grxi*(1.d0-grxi)/2.d0
                Am(m)=Am(m)+wter(i)*(1.d0-(grxi**2.d0))
                Am(m+1)=Am(m+1)+wter(i)*grxi*(1.d0+grxi)/2.d0
              else
                grxi=(-x(i)+1.d0)/2.d0
                Am(m)=Am(m)+wter(i)*grxi*(1.d0+grxi)/2.d0 
                Am(m+1)=Am(m+1)+wter(i)*(1.d0-(grxi**2.d0))
                Am(m+2)=Am(m+2)-wter(i)*grxi*(1.d0-grxi)/2.d0  
              endif                               
 511        continue 
 51      continue  
c     Test des coefficiants Am(n)
        do 522 m=1,Nx+1
             sum=sum+Am(m)*(dx/2.d0)
     &           *dexp(-1.d0+((dble(m)-1.d0)*dx))
 522     continue   
        write(*,*) 'Nx=',Nx, 'sum=',sum, 
     &            'err=',sum-dexp(1.d0)+1.d0/dexp(1.d0)

c     Calcul de la premiere ligne de la matrice amat avec Gauss
      if (methode2.eq.1) then
         do 52 n=1,Nx+1 
           xi=-1.D0+dx*dble(n-1)  
           amat(1,n)=Am(m)/((1.d0-xi/ch)**3)
c           write(13,*) xi, Am(m)
 52      continue
      endif 
c      Fin Methode des points de gauss

c     constitution de la premiere ligne de amat en Xi=-1
c     --------------------------------------------------
c       zeta1'(X) + zeta1*ch/(1+ch) = 0  en  X=-1 (tangente horizontale)
          
c$$$       amat(1,1)=ch/(1.d0+ch)-(3.d0/2.d0)*invdx 
c$$$       amat(1,2)=2.d0*invdx
c$$$       amat(1,3)=-(1.d0/2.d0)*invdx
c$$$       write(*,*)'amat(1,1)=',amat(1,1)
c$$$       write(*,*)'amat(1,2)=',amat(1,2)
c$$$       write(*,*)'amat(1,3)=',amat(1,3)
c     constitution des lignes 2 � la ligne Nx de amat
c     -----------------------------------------------
c     On resout ici l'eq differentielle    
c     fn+1.[(an/h�) + (bn/2h)] + fn.[(-2an/h�) + cn] 
c       + fn-1.[(an/h�) - (bn/2h)] = dn   

      do 6 n=2,Nx
       xi=-1.D0+dx*dble(n-1)
       tabx(n-1)=xi
       amat(n,n-1)=(1.d0-(xi**2))*((ch-xi)**2)*invdx2
     &             -(1.d0-ch*xi)*(ch-xi)*invdx            
       amat(n,n)=2.d0*ch*(ch-xi)-2.d0*(1.d0-(xi**2))
     &           *((ch-xi)**2)*invdx2       
       amat(n,n+1)=(1.d0-(xi**2))*((ch-xi)**2)*invdx2
     &             +(1.d0-ch*xi)*(ch-xi)*invdx     
6     continue

c     constitution de la ligne Nx+1 de amat
c     -------------------------------------
     
c      xNt=dble(Nt)
c      amat(Nx+1,Nx-1)=
c      amat(Nx+1,Nx)=
c      amat(Nx+1,Nx+1)=

c     constitution de la derniere ligne de amat en Xi=1
c     --------------------------------------------------
c     zeta1'(X) - zeta1*ch/(1-ch) = 0 en  X=1 (tangente horizontale) 

       amat(Nx+1,Nx-1)=(1.d0/2.d0)*invdx
       amat(Nx+1,Nx)=-2.d0*invdx ! valeur correcte
c       amat(Nx+1,Nx)=-4.d0*invdx ! test de la CL
       amat(Nx+1,Nx+1)=ch/(1.d0-ch)+(3.d0/2.d0)*invdx
       write(*,*)
       write(*,*)'amat(Nx+1,Nx+1)=',amat(Nx+1,Nx+1)
       write(*,*)'amat(Nx+1,Nx)=',amat(Nx+1,Nx)
       write(*,*)'amat(Nx+1,Nx-1)=',amat(Nx+1,Nx-1)
c     constitution de la ligne Nx+1 de amat
c     Condition au limite sur le centre de masse
c     ------------------------------------------
c     methode des trapezes
      if (methode1.eq.1) then
         dx2=dx*dx
         amat(Nx+1,1)=dx2/((1.d0+1.d0/ch)**4)
         do 61 n=2,Nx
            xi=-1.D0+dx*dble(n-1)
            amat(Nx+1,n)=2.d0*dx2/((1.d0-xi/ch)**4)
61      continue
        amat(Nx+1,Nx+1)=dx2/((1.d0-1.d0/ch)**4)
      endif
c     Fin methode des trap�zes

c     Test du second membre
         xi=-1.d0
         sum1= yse(1,1)/((ch-xi)**3)
         do 63 n=2,Nx
            xi=-1.D0+dx*dble(n-1)
            sum1=sum1+2.d0* yse(n,1)/((ch-xi)**3)
63      continue
        xi=1.d0
        sum1=sum1+ yse(Nx,1)/((ch-xi)**3)
        sum1=sum1*dx/(2.d0)
        write(*,*) 'test int�grale second membre', sum1
c     Fin du test du second membre

c     methode des points de gauss
      if (methode2.eq.1) then 
         do 62 n=1,Nx+1 
           xi=-1.D0+dx*dble(n-1) 
           amat(Nx+1,n)=Am(m)/((1.d0-xi/ch)**4)
 62      continue
      endif      
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

c     resolution du systeme, Appel du sous programme dgetrs 
c     ----------------------------------------------------
      call dgetrs(TRANS,Nx+1,1,amat,Nx+1,IPIV,yse,Nx+1,INFO)

c     Calcul du Determinant
c     ---------------------
      aux=0.d0
      do 777 n=1,Nx+1
         aux=aux*amat(n,n)         
 777     continue
         write(*,*) 'det matrice=',aux

c     ecriture des resultats dans le tableau zeta1B(n)
c     ---------------------------------------------
      
c     ATTENTION le tableau H_n commence par l'indice 1
c     il faut donc definir le tableau tabh0 qui est la
c     premiere ligne de la matrice
      do 66 n=1,Nx
       zeta1B(n)=yse(n+1,1)
 66    continue
      zeta1B(0)=yse(1,1)
      write (*,*) 'zeta1B(0)=',zeta1B(0)
      write (*,*) 'zeta1B(Nx)=',zeta1B(Nx)
      tabx(0)=-1.d0
      tabx(1)=tabx(0)+dx
      tabx(Nx)=1.d0
      close(13)
c     =======================================
c     Fin du sous-programme diagonalessysteme1
c     =======================================
      end subroutine resolutionsystLU
