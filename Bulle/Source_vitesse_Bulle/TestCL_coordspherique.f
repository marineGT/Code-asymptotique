      subroutine TestCL(Nt,Nx,c,zetap,lambda0,alpha,alpha0,Bn,
     &                  Dn,tabsigma0zz,phi,tabteta)
c     ==================================================================
c     Aucune declaration implicite
c     ----------------------------
      implicit none
c     ==================================================================
c     Declaration des parametres
c     --------------------------

c     Declaration des variables en entree
c     -----------------------------------
      integer, intent(in) :: Nt, Nx
      double precision, intent(in) :: c, zetap
      double precision, intent(in) :: lambda0, alpha0
      double precision,dimension(Nt) :: alpha
      double precision, dimension(Nt+1) :: Bn,Dn
      double precision, dimension(Nx+1) :: tabsigma0zz

c     Declaration des variables en sortie
c     -----------------------------------     
      double precision, dimension(0:Nx) :: phi
      double precision, dimension(0:Nx) :: tabteta

c     Declaration des variables locales
c     ---------------------------------
      integer :: i, n, m, k, INFO, Nx1
      double precision :: teta, dteta, invdteta
      double precision :: invdteta2, pi, Sn
      double precision :: grxi,sum,aux,sum1,sum2
      double precision :: aux1, aux2, aux3
      double precision :: prec, e, f
      double precision :: cstSn, cstCM
      double precision :: testf
      double precision,dimension(7) :: x,wter
      character*1 :: TRANS, NORM
      double precision,dimension (Nx+1) :: IPIV
      double precision,dimension (3) ::coef,IPIV1
      double precision,dimension (Nx+1,Nx+1) :: amat
      double precision,dimension (Nx+1,1) :: yse
      double precision,dimension (Nx+1) :: Am
      double precision, dimension(0:Nx) :: phipart
      double precision, dimension(Nx) :: b1, e1, d1
c      double precision, dimension(Nx+1) :: b2, e2, d2
      double precision, dimension(Nx) :: b2, e2, d2
      external dgetrf,dgetrs

c     Declaration des interfaces
c     --------------------------
      interface

       function calculsigma0zzB(teta,Nt,c,zetap,alpha,
     &                       alpha0,Bn,Dn)      
         integer, intent(in) ::  Nt      
         double precision, intent(in) :: teta
         double precision, intent(in) ::  c, zetap
         double precision, dimension(Nt) :: alpha
         double precision, intent(in) :: alpha0 
         double precision, dimension(Nt+1) :: Bn,Dn
         double precision :: calculsigma0zzB
        end function calculsigma0zzB

        function solparticuliere(b,teta,Nt,c,zetap,
     &              lambda0,alpha,alpha0,Bn,Dn,cstSn)
         integer, intent(in)  ::  Nt      
         double precision, intent(in)  :: b, teta
         double precision, intent(in) ::  c, zetap
         double precision, intent(in) :: lambda0
         double precision, dimension(Nt) :: alpha
         double precision, intent(in) :: alpha0 
         double precision, dimension(Nt+1) :: Bn,Dn
         double precision, intent(in) :: cstSn
         double precision :: solparticuliere
        end function solparticuliere

        function constanteSn(teta,Nt,c,zetap,
     &                 lambda0,alpha,alpha0,Bn,Dn)
         integer, intent(in)  ::  Nt      
         double precision, intent(in)  :: teta
         double precision, intent(in) ::  c, zetap
         double precision, intent(in) :: lambda0
         double precision, dimension(Nt) :: alpha
         double precision, intent(in) :: alpha0 
         double precision, dimension(Nt+1) :: Bn,Dn
         double precision :: constanteSn
        end function constanteSn

        function constanteCM(teta,Nx,phi,tabteta,
     &                       b1,e1,d1)
         integer, intent(in)  ::  Nx
         double precision, intent(in)  :: teta
         double precision, dimension(Nx) :: b1, e1, d1
         double precision, dimension(0:Nx) :: phi
         double precision, dimension(0:Nx) :: tabteta
         double precision :: constanteCM
        end function constanteCM

        function testvolume(teta,Nx,g,tabg,b,e,d)
          integer, intent(in)  ::  Nx
          double precision, intent(in)  :: teta
          double precision, dimension(Nx) :: b, e, d
          double precision, dimension(0:Nx) :: g
          double precision, dimension(0:Nx) :: tabg 
          double precision :: testvolume
        end function testvolume

        function testcentremasse(teta,Nx,f,tabf,b3,e3,d3)
          integer ::  Nx
          double precision, intent(in)  :: teta
          double precision, dimension(Nx) :: b3, e3, d3
          double precision, dimension(0:Nx) :: f
          double precision, dimension(0:Nx) :: tabf 
          double precision :: testcentremasse
        end function testcentremasse

        function solpartitest(teta)
           double precision, intent(in)  :: teta
           double precision :: solpartitest
        end function solpartitest

        function Sntot(teta,Nt,c,zetap,lambda0,alpha,
     &               alpha0,Bn,Dn)
          integer, intent(in)  ::  Nt      
          double precision, intent(in)  :: teta
          double precision, intent(in) ::  c, zetap
          double precision, intent(in) :: lambda0
          double precision, dimension(Nt) :: alpha
          double precision, intent(in) :: alpha0 
          double precision, dimension(Nt+1) :: Bn,Dn
          double precision :: Sntot        
        end function Sntot
      
        subroutine spline(x, y, b, c, d, n)
           integer  n
           double precision  x(n),y(n),b(n),c(n),d(n)
        end subroutine spline

       end interface 
 
       double precision gaussint1dim3
       external gaussint1dim3
       double precision gaussint1dim4
       external gaussint1dim4
       double precision gaussint1dim5
       external gaussint1dim5
       double precision gaussint1dim6
       external gaussint1dim6
       double precision gaussint1dim7
       external gaussint1dim7
       double precision gaussint1dimtest
       external gaussint1dimtest
       double precision gaussint1dimSn
       external gaussint1dimSn       

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
c     
c     ==================================================================
c     calcul de pi
c     ------------
      pi=2.d0*dasin(1.d0)

c     Calcul de dx et de ses inverses 
c     --------------------------------
      dteta=pi/dble(Nx)
c      write(*,*) dteta
      invdteta=1.D0/dteta
      invdteta2=invdteta*invdteta

c     Constitution du second membre yse(n,1) pour le systeme
c     lineaire des Hn.

c     ATTENTION selon le cas/conditions choisi(e)s, on doit penser à
c     modifier la matrice amat aussi
      do 1 n=1,Nx+1
         yse(n,1)=0.d0
         phipart(n-1)=0.d0
1     continue    

c     la première ligne de yse = 0 (CL)
c     --------------------------------
      yse(1,1)=0.d0
c      write(*,*) 'yse(1,1)',yse(1,1)

c     valeur de yse pour n=2,Nx
c     -------------------------
      teta=dteta
      do 2 n=2,Nx
        teta=dble(n-1)*dteta
        tabteta(n-1)=teta
c        yse(n,1)=-7.d0*dcos(3.d0*teta)*dsin(teta)
c     &           -3.d0*dcos(teta)*dsin(3.d0*teta) !test solution homogene
c         yse(n,1)=dcos(teta)*dsin(teta)
c          yse(n,1)=dsin(teta)
          Sn=calculsigma0zzB(teta,Nt,c,zetap,
     &            alpha,alpha0,Bn,Dn)

          ! si Bo=rho g a*a/gamma0
          aux=3.D0*lambda0*dcos(teta)

           yse(n,1)=-dsin(teta)*(Sn+aux)
c         yse(n,1)=-2.d0*dcos(2.d0*teta)*dsin(teta)
c     &           -2.d0*dcos(teta)*dsin(2.d0*teta)
c        yse(n,1)=0.d0
c        write(17,*) teta, yse(n,1),Sn,aux
 2    continue 

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

c     constitution de la premiere ligne de amat en teta=0
c     --------------------------------------------------
c       phi'(0) = 0  en teta=0 (tangente horizontale)

       !derivee premiere nulle en teta=0.d0  
c$$$       amat(1,1)=-3.d0/(2.d0*dteta)
c$$$       amat(1,2)=2.d0*invdteta
c$$$       amat(1,3)=-(1.d0/2.d0)*invdteta
       !derivee seconde nulle en teta=0.d0
c$$$       amat(1,1)=invdteta2
c$$$       amat(1,2)=-2.d0*invdteta2
c$$$       amat(1,3)=invdteta2
       amat(1,1)=1.d0 !solution particuliere

c$$$       write(*,*)'amat(1,1)=',amat(1,1)
c$$$       write(*,*)'amat(1,2)=',amat(1,2)
c$$$       write(*,*)'amat(1,3)=',amat(1,3)
c     constitution des lignes 2 à la ligne Nx de amat
c     -----------------------------------------------
c     On resout ici l'eq differentielle    
c     fn+1.[(an/h²) + (bn/2h)] + fn.[(-2an/h²) + cn] 
c       + fn-1.[(an/h²) - (bn/2h)] = dn   
       teta=dteta
      do 6 n=2,Nx    
       teta=dble(n-1)*dteta
       amat(n,n-1)=dsin(teta)*invdteta2
     &             -(1.d0/2.d0)*dcos(teta)*invdteta         
       amat(n,n)=-2.d0*dsin(teta)*invdteta2
     &           +2.d0*dsin(teta)
       amat(n,n+1)=dsin(teta)*invdteta2
     &             +(1.d0/2.d0)*dcos(teta)*invdteta
6     continue

c     constitution de la derniere ligne de amat en teta=pi
c     --------------------------------------------------
c     phi'(teta=pi) = 0 en  teta=pi (tangente horizontale) 

       amat(Nx+1,Nx-1)=(1.d0/2.d0)*invdteta
       amat(Nx+1,Nx)=-2.d0*invdteta      
       amat(Nx+1,Nx+1)=(3.d0/2.d0)*invdteta
c       write(*,*)
       write(*,*)'amat(Nx+1,Nx+1)=',amat(Nx+1,Nx+1)
       write(*,*)'amat(Nx+1,Nx)=',amat(Nx+1,Nx)
       write(*,*)'amat(Nx+1,Nx-1)=',amat(Nx+1,Nx-1)


c     Test du second membre
c     ---------------------
       sum1=0.d0
         teta=dteta        
         do 63 n=2,Nx
            teta=dble(n-1)*dteta
c            sum1=sum1+yse(n,1)*dcos(teta)
            sum1=sum1+yse(n,1)
63      continue
c        teta=pi
c        sum1=sum1+yse(Nx,1)
        sum1=-0.5d0*sum1*dteta
        write(*,*) 'constante ajoutee au second membre', 
     &               sum1

       prec=1.d-14 ! precision de l'integrale
       e=0.d0 
       f=pi
       cstSn=gaussint1dim4(constanteSn,e,f,Nt,c,zetap,
     &             lambda0,alpha,alpha0,Bn,Dn,prec)
       cstSn=-0.5d0*cstSn
        write(*,*) 'constante ajoutee au second membre', 
     &              cstSn

c     Modification du second membre !!!
c     ---------------------------------
      do n=2,Nx
        teta=dble(n-1)*dteta        
        yse(n,1)=yse(n,1)+cstSn*dsin(teta)
      enddo

c     Test du second membre
c     ---------------------
      sum1=0.d0
      teta=dteta        
      do n=2,Nx
       teta=dble(n-1)*dteta
       sum1=sum1+yse(n,1)*dcos(teta)
      enddo
      sum1=sum1*dteta
      write(*,*) 'test integrale second membre', sum1  

      prec=1.d-14 ! precision de l'integrale
      e=0.d0 
      f=pi
      sum1=gaussint1dimSn(Sntot,e,f,Nt,c,zetap,lambda0,
     &                    alpha,alpha0,Bn,Dn,prec)
       
      write(*,*) 'test integrale second membre', sum1  


c$$$c     Test solution phi=cos(3teta)
c$$$c     ---------------------------  
c$$$        teta=dteta        
c$$$c        do n=1,Nx+1
c$$$        do n=2,Nx
c$$$           aux=0.d0
c$$$           do m=2,Nx
c$$$              teta=dble(m-1)*dteta
c$$$              aux1=(1.d0+dcos(teta)
c$$$     &             *dlog(dtan(teta/2.d0)))
c$$$c              aux2=-(1.d0/6.d0)*dcos(teta)*(dlog(2.d0)
c$$$c     &             +4.d0*dlog(dsin(teta)))
c$$$c              aux2=-(1.d0/3.d0)*dlog(dsin(teta))*dcos(teta) 
c$$$c     &             +(1.d0/3.d0)*aux1
c$$$              aux2=-dcos(teta)*dlog(dtan(teta/2.d0))/2.d0
c$$$              aux=aux+amat(n,m)*aux2
c$$$c              aux=aux+amat(n,m)*dcos(teta) 
c$$$           enddo
c$$$c           write(*,*) 'n=',n,'num=',aux,'theo=',yse(n,1)            
c$$$c           write(*,*) 'n=',n, 'theo-num=', yse(n,1)-aux 
c$$$        enddo        
c$$$
c$$$        do  m=2,Nx         
c$$$            teta=dble(m-1)*dteta
c$$$            aux1=(1.d0+dcos(teta)
c$$$     &           *dlog(dtan(teta/2.d0)))
c$$$            aux2=-(1.d0/3.d0)*dlog(dsin(teta))*dcos(teta) 
c$$$     &             +(1.d0/3.d0)*aux1
c$$$            aux3=(1.d0/3.d0)*(dcos(teta)*(1.d0-dcos(teta))
c$$$     &            /dsin(teta)+dsin(teta)*(dlog(dsin(teta))
c$$$     &            -dlog(dtan(teta/2.d0))))
c$$$c            write(15,*) teta, aux2, aux3
c$$$        enddo

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
      write(*,*) 'INFO',INFO

c     resolution du systeme, Appel du sous programme dgetrs 
c     ----------------------------------------------------
      call dgetrs(TRANS,Nx+1,1,amat,Nx+1,IPIV,yse,Nx+1,INFO)
      write(*,*) 'INFO',INFO

c     Calcul du Determinant
c     ---------------------
      aux=0.d0
      do 777 n=1,Nx+1
         aux=aux+amat(n,n)         
 777     continue
c         write(*,*) 'det matrice=',aux

c     ecriture des resultats dans le tableau phi(n)
c     --------------------------------------------- 
      phi(0)=yse(1,1)
      write(15,*) 0.d0, phi(0), 1.d0
      teta=dteta
      do 66 n=1,Nx
       phi(n)=yse(n+1,1)
       teta=dble(n)*dteta
c       write(15,*) teta, phi(n), dcos(2.d0*teta),
c     &             dcos(teta)   
        write(15,*) teta, phi(n), dcos(teta)   
 66    continue
       write(15,*) pi, phi(Nx), -1.d0
      write (*,*) 'phi(0)=',phi(0)
      write (*,*) 'phi(Nx)=',phi(Nx)
      tabteta(0)=0.d0
      tabteta(Nx)=pi

c     Test de la conservation du volume
c     ---------------------------------
c     methode des trapezes      
        sum2=(phi(0)+phi(Nx))/2.d0
        do 5 n=2,Nx
            teta=dble(n-1)*dteta
            sum2=sum2+phi(n-1)*dsin(teta)
 5      continue
        sum2=sum2*dteta
        write(*,*) 'test conservation volume', sum2 
c     Fin methode des trapèzes

c     Test sur le centre de masse
c     ---------------------------
      sum2=(phi(0)+phi(Nx))/2.d0
       do n=2,Nx
         teta=dble(n-1)*dteta
         sum2=sum2+phi(n-1)*dsin(teta)
     &       *dcos(teta)   
       enddo
        sum2=sum2*dteta
        write(*,*) 'constante conservation CM', sum2 
      prec=1.d-14
      e=0.d0
      f=pi
      call spline(tabteta,phi,b1,e1,d1,Nx) 
      
      cstCM=gaussint1dim5(constanteCM,e,f,Nx,phi,
     &                   tabteta,b1,e1,d1,prec)

      write(*,*) 'constante conservation CM', cstCM

c     Solution totale fp + Acos(teta)
c     -------------------------------        
      do  n=1,Nx+1 
        teta=dble(n-1)*dteta
        yse(n,1)=yse(n,1)-3.d0*cstCM*dcos(teta)/2.d0
c        yse(n,1)=yse(n,1)-3.d0*sum2*dcos(teta)/2.d0
        write(17,*) teta, yse(n,1)    
      enddo

c     Calcul de la solution particuliere totale
c     ------------------------------------------
      !Phip(teta)= int_{0}^{teta} { sin(u)S(u)
      !            [(1+cos(teta)*log(tan(teta/2)))cos(u)
      !             -cos(teta)*(1+cos(u)*log(tan(u/2)))] du 
      !            }
        
c     Avec S(u) = -[sigma0zzB(u) + 3*lambda0*cos(u)]
       prec=1.d-5 ! precision de l'integrale
       e=0.d0 
       teta=dteta
      do n=2,Nx+1 
        teta=dble(n-1)*dteta
        f=teta
        phipart(n-1)=gaussint1dim3(solparticuliere,e,
     &           f,teta,Nt,c,zetap,lambda0,alpha,
     &           alpha0,Bn,Dn,cstSn,prec)
        write(18,*) teta, phipart(n-1)
      enddo
      f=pi
      teta=pi
      sum2=0.d0
      sum2=gaussint1dim3(solparticuliere,e,f,teta,Nt,c,
     &     zetap,lambda0,alpha,alpha0,Bn,Dn,cstSn,prec)
      write(18,*) pi, sum2

      sum2=0.d0
      do n=1,Nx+1
        teta=dble(n-1)*dteta
        sum2=dcos(3.d0*teta)-dcos(teta)
        write(20,*) teta, sum2 
      enddo

      !test volume avec f=cos(3*teta)
      
c       e=0.d0 
c       f=pi
c       prec=1.d-12
c       testf=gaussint1dimtest(solpartitest,e,f,prec)

c       write(*,*) 'test volume testf', testf
       

c      write(*,*) 'Solution particuliere en pi', sum2            
      
c$$$c     Test de la conservation du volume
c$$$c     ---------------------------------
c$$$c     methode des trapezes      
c$$$         sum2=(phipart(Nx))/2.d0
c$$$         do n=2,Nx
c$$$            teta=dble(n-1)*dteta
c$$$            sum2=sum2+phipart(n-1)*dsin(teta)
c$$$         enddo
c$$$        sum2=sum2*dteta
c$$$        write(*,*) 'test volume sol theo', sum2 
c$$$c     Fin methode des trapèzes
c$$$
c$$$      call spline(tabteta,phipart,b1,e1,d1,Nx)      
c$$$      prec=1.d-12
c$$$      sum2=gaussint1dim6(testvolume,e,f,Nx,phipart,
c$$$     &                   tabteta,b1,e1,d1,prec)
c$$$
c$$$      write(*,*) 'test volume sol theo', sum2 


c     Solution totale fp + Acos(teta)
c     ------------------------------------------
        do n=1,Nx+1 
          teta=dble(n-1)*dteta 
          phi(n-1)=yse(n,1)
c          write(15,*) teta, phi(n-1)
        enddo   

c     Test de la conservation du centre de masse
c     ------------------------------------------
c     methode des trapezes      
         sum2=(yse(1,1)+yse(Nx+1,1))/2.d0
         do n=2,Nx 
            teta=dble(n-1)*dteta
            sum2=sum2+yse(n,1)*dsin(teta)
     &           *dcos(teta)   
         enddo
        sum2=sum2*dteta
        write(*,*) 'test barycentre', sum2 
        write(*,*) 'phi(Nx)=', phi(Nx)
        
        do n=1,Nx
           b1(n)=0.d0
           e1(n)=0.d0
           d1(n)=0.d0
        enddo
      sum2=0.d0
      e=0.d0 
      f=pi
      prec=1.d-12       
      call spline(tabteta,phi,b1,e1,d1,Nx)      
      
      sum2=gaussint1dim7(testcentremasse,e,f,Nx,
     &         phi,tabteta,b1,e1,d1,prec)

      write(*,*) 'test barycentre', sum2 
      
c$$$      call spline(tabteta,phipart,b1,e1,d1,Nx) 
c$$$      
c$$$      sum2=gaussint1dim5(constanteCM,e,f,Nx,phipart,
c$$$     &                   tabteta,b1,e1,d1,prec)
c$$$
c$$$      write(*,*) 'CM avec sol theo', sum2 



c     =======================================
c     Fin du sous-programme diagonalessysteme1
c     =======================================
      end subroutine TestCL
