
c****************************************************************
c    Programme de constitution de la derivee de P_n(x) ou P 
c     est le polynome de Legendre d'ordre n et x sont argument 
c      avec |x|\leq 1 et n\geq.1
c
c    Notations:
c     n: ordrer du polynome (entier naturel positif ou nul)
c     x: argument
c     resultat: valeur de la derivee de P_n(x)
c     Programmeur: A.Sellier; mars 2004.
c****************************************************************
c2345678912345678912345678912345678912345678912345678912345678912
c     declaration des variables
      subroutine DLegendre(n,x,resultat)
      external LAegendre
      integer n,k,kmax
      double precision x,resultat,xn,aux
c     test
      if (n.le.0) then
c      write(*,*) 'probleme1 sur ordre du polynome de DLegrendre' 
       resultat=0.d0
       goto 123
      endif
      if (n.ge.3001) then 
      write(*,*) 'probleme2 sur ordre du polynome de DLegrendre'
      endif
      xn=(-1.d0)**n
      if (xn.eq.(1.d0)) then
       kmax=n/2-1
      else
       kmax=(n-1)/2
      endif
      resultat=0.d0
c      write (*,*) 'kmax=',kmax
      do 1 k=0,kmax
       call LAegendre(n-2*k-1,x,aux) 
c       write (*,*) 'k=',k,'aux=',aux
c       write (*,*) 'k=',k
       resultat=resultat+dble(2*n-4*k-1)*aux     
 1    continue
c      write (*,*) 'coucou2'
c      write (*,*) 'n=',n,'kmax=',kmax
 123  continue
      end
c****************************************************************
c2345678912345678912345678912345678912345678912345678912345678912
c     declaration des variables
      subroutine LAegendre(n,x,resultat)
      integer n,k,m
      double precision x,resultat,a,xk
      dimension a(3000)
c     test
      if (n.lt.0) then
      write(*,*) 'probleme1 sur ordre du polynome de Legrendre' 
      endif
      if (n.ge.3001) then 
      write(*,*) 'probleme2 sur ordre du polynome de Legrendre' 
      endif
      if (n.eq.0) resultat=1.d0
      if (n.eq.1) resultat=x
      if (n.eq.2) resultat=(3.d0*(x**2)-1.d0)/2.d0
      if (n.ge.3) then
c      initialisation des coefficients a
       do 1 m=1,3000
        a(m)=0.d0
 1     continue
        a(1)=x
        a(2)=(3.d0*(x**2)-1.d0)/2.d0
        do 2 k=2,n-1
         xk=dble(k)
         a(k+1)=-xk*a(k-1)/(xk+1.d0)
     +          +x*(2.d0*xk+1.d0)*a(k)/(xk+1.d0)
 2      continue
        resultat=a(n)       
      endif
      end
c****************************************************************
