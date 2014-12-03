      subroutine resolutionsystLU(Nt,Nx,Ng,precision,lsa,c,c2,dx,
     &            zetap,lambda0,gammahat,tabsB,tabsB0,
     &            R0,Rn,tabx,zeta1B)
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
      double precision, intent(in) :: tabsB0, R0
      double precision, dimension(Nt),intent(in) :: Rn
      double precision, dimension(Nt+1) :: tabsB

c     Declaration des variables en sortie
c     -----------------------------------
      double precision, dimension(0:Nx) :: tabx      
      double precision, dimension(0:Nx) :: zeta1B

c     Declaration des variables locales
c     ---------------------------------
      integer :: i, n, m, methode1, methode2, INFO
      double precision :: invdx,invdx2,xi,ch,x
      character*1 :: TRANS, NORM
      double precision,dimension (Nx+1) :: IPIV
      double precision,dimension (3) ::coef,IPIV1
      double precision,dimension (Nx+1,Nx+1) :: amat
      double precision,dimension (Nx+1,1) :: yse
      double precision,dimension (Nx-1) :: Am
      external dgetrf,dgetrs
 
c     Declaration des interfaces
c     --------------------------
      interface

         function secondmembreB(Nt,c,c2,precision,lambda0,
     &                gammahat,tabsB,tabsB0,R0,Rn,xi)
            integer,intent(in) :: Nt 
            double precision, intent(in) :: c, c2, precision, xi 
            double precision, intent(in) :: lambda0, gammahat
            double precision, intent(in) :: tabsB0, R0 
            double precision, dimension(Nt),intent(in) :: Rn
            double precision, dimension(Nt+1),intent(in) :: tabsB	   
	    double precision :: secondmembreB            
         end function secondmembreB
     
      end interface
      
      data x
      #        / 0.33998 10435 84856 26480 D0,
      #         -0.33998 10435 84856 26480 D0,
      #         -0.86113 63115 94052 57522 D0,
      #          0.86113 63115 94052 57522 D0/
c2345678912345678912345678912345678912345678912345678912345678912
      data wter
     #        / 0.65214 51548 62546 14262 D0,
     #         -0.65214 51548 62546 14262 D0,
     #         -0.34785 48451 37453 85737 D0,
     #          0.34785 48451 37453 85737 D0/
c     
c     ==================================================================

c     Calcul des inverses de dx et dx2
c     --------------------------------
      invdx=1.D0/dx
      invdx2=invdx*invdx
      ch=dcosh(zetap)

c     Condition limites sur l'ensemble de la fonction zeta1
c     -----------------------------------------------------
c     le premiere methode est celle des trapezes :
c        int(-1,1){f(x)dx}=h*sum(m=1,M){f(m-1)+f(m)}      
 
c     la deuxieme est celle des points de Gauss:
c         int(-1,1){f(x)dx}=sum(i=1,N){w(i)*f(ni)}
      
c     Pour nous f(x)=zetaB1/[ch(zetap)-xi]
      methode1=0
      methode2=1

c     ******************************************************************
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
c      open(unit=13,file='secondmembreB_res.dat')

      do 2 n=2,Nx
        xi=-1.D0+dx*dble(n) 
c        yse(n,1)=c*((1.d0-xi)**(3.d0/2.d0))*secondmembreB(Nt,
c     &         c,c2,precision,Bo,Ca,tabsB,tabsB0,R0,Rn,xi)  
        yse(n,1)=((1.d0-xi)**(3.d0/2.d0))*secondmembreB(Nt,
     &        c,c2,precision,lambda0,gammahat,
     &         tabsB,tabsB0,R0,Rn,xi)  
c      write(*,*) 'n,yse(n,1)',n,yse(n,1)
c        write(13,*) xi, yse(n,1)  
 2    continue
  
c      close(13)

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
         amat(1,1)=invdx2/((1.d0+1.d0/ch)**3)
         do 5 n=2,Nx 
            xi=-1.D0+dx*dble(n-1)
            tabx(n-1)=xi
            amat(1,n)=2.d0*invdx2/((1.d0-xi/ch)**3)
 5      continue
        amat(1,Nx+1)=invdx2/((1.d0-1.d0/ch)**3)
      endif

c     methode des points de gauss
      if (methode2.eq.1) then

         do 50 n=1,Nx-1
            Am(n)=0.d0
 50      continue 

         do 51 m=1,Nx-1
            xi=-1.D0+dx*dble(n-1)
            do 511 i=1,Ng              
              if(m.ne.1) then 
                grxi=(x(i)-xi)*invdx 
                Am(m-1)=Am(m-1)-wter(i)*grxi*(1.d0+grxi)/2.d0
                Am(m)=Am(m)+wter(i)*(1.d0-(grxi**2.d0))
                Am(m+1)=Am(m+1)+wter(i)*grxi*(1.d0+grxi)/2.d0
              else
                grxi=(xi-x(i))*invdx  
                Am(m)=Am(m)+wter(i)*(1.d0-(grxi**2.d0))
                Am(m+1)=Am(m+1)+wter(i)*grxi*(1.d0-grxi)/2.d0 
                Am(m+2)=Am(m+2)-wter(i)*grxi*(1.d0-grxi)/2.d0  
              endif                               
 511        continue 
 51      continue 

      endif
      

  
c     constitution des lignes 2 à la ligne Nx de amat
c     -----------------------------------------------
c     On resout ici l'eq differentielle    
c     fn+1.[(an/h²) + (bn/2h)] + fn.[(-2an/h²) + cn] 
c                              + fn-1.[(an/h²) - (bn/2h)] = dn      

      do 6 n=2,Nx    
       xi=-1.D0+dx*dble(n-1)
       amat(n,n-1)=(1.d0-xi**2)*((ch-xi)**2)*invdx2
     &             -(1.d0-ch*xi)*(ch-xi)*invdx            
       amat(n,n)=2.d0*ch*(ch-xi)-2.d0*(1.d0-xi**2)
     &           *((ch-xi)**2)*invdx2       
       amat(n,n+1)=(1.d0-xi**2)*((ch-xi)**2)*invdx2
     &             +(1.d0-ch*xi)*(ch-xi)*invdx     
6     continue


c     constitution de la ligne Nx+1 de amat
c     Condition au limite sur le centre de masse
c     ------------------------------------------

c     methode des trapezes
      if (methode1.eq.1) then
         amat(Nx+1,1)=invdx2/((1.d0+1.d0/ch)**4)
         do 61 n=2,Nx
            xi=-1.D0+dx*dble(n-1)
            amat(Nx+1,n)=2.d0*invdx2/((1.d0-xi/ch)**4)
61      continue
        amat(Nx+1,Nx+1)=invdx2/((1.d0-1.d0/ch)**4)
      endif


c     methode des points de gauss
      if (methode2.eq.1) then

         do 60 n=1,Nx-1
            Am(n)=0.d0
 60      continue 

         do 61 m=1,Nx-1
            xi=-1.D0+dx*dble(n-1)
            do 611 i=1,Ng
              if(m.ne.1) then 
                grxi=(x(i)-xi)*invdx 
                Am(m-1)=Am(m-1)-wter(i)*grxi*(1.d0+grxi)/2.d0
                Am(m)=Am(m)+wter(i)*(1.d0-(grxi**2.d0))
                Am(m+1)=Am(m+1)+wter(i)*grxi*(1.d0+grxi)/2.d0
              else
                grxi=(xi-x(i))*invdx 
                Am(m)=Am(m)+wter(i)*(1.d0-(grxi**2.d0))
                Am(m+1)=Am(m+1)+wter(i)*grxi*(1.d0-grxi)/2.d0 
                Am(m+2)=Am(m+2)-wter(i)*grxi*(1.d0-grxi)/2.d0  
              endif                               
 611        continue 
 61      continue 

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
c      write (*,*) 'INFO1=',INFO

c     resolution du systeme, Appel du sous programme dgetrs 
c     ----------------------------------------------------
      call dgetrs(TRANS,Nx+1,1,amat,Nx+1,IPIV,yse,Nx+1,INFO)
c      write (*,*) 'INFO2=',INFO

c$$$c     Calcul du Determinant
c$$$c     ---------------------
c$$$      aux=0.d0
c$$$      do 777 n=1,Nx+1
c$$$         aux=aux+amat(n,n)         
c$$$ 777     continue
c$$$         write(*,*) 'det matrice=',aux


c     ecriture des resultats dans le tableau zeta1B(n)
c     ---------------------------------------------
      
c     ATTENTION le tableau H_n commence par l'indice 1
c     il faut donc definir le tableau tabh0 qui est la
c     premiere ligne de la matrice
      do 66 n=1,Nx
       zeta1B(n)=yse(n+1,1)
 66    continue
      zeta1B(0)=yse(1,1)

      tabx(0)=-1.d0
      tabx(1)=tabx(0)+dx
      tabx(Nx)=1.d0
c     =======================================
c     Fin du sous-programme diagonalessysteme1
c     =======================================
      end subroutine resolutionsystLU
