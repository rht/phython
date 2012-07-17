from pylab import *

@vectorize
def cube(x): return x ** 3

def compose(f,g):
    return lambda x: f(g(x))

def f(x,y):
    return square(x) * cube(y)

def g(x,y):
    return (up (f x y) y))

def h(x, y):
    return f(f(x,y),y)

def f1(v):
    x = v[0]
    y = v[1]
    return square(x) * cube(y)

def g1(v):
    x = v[0]
    y = v[1]

    return (up (f x y) y))
    return (up (f1 v) y)))

def h1: return compose(f1, g1)

#a
(pe (((partial 0) f) 'x 'y)) ;(* 2 x (expt y 3))
(pe (((partial 1) f) 'x 'y)) ;(* 3 (expt x 2) (expt y 2))

;;a tuple version
 (pe (((partial 0) f1) (up 'x 'y))) ;(* 2 x (expt y 3))   SAME
 (pe (((partial 1) f1) (up 'x 'y))) ;(* 3 (expt x 2) (expt y 2))    SAME



;;b
;;(pe (((partial 0) (compose f g)) 'x 'y))
;;(pe (((partial 1) (compose f g)) 'x 'y))
(pe (((partial 0) f) (f 'x 'y) 'y))
;(* 2 (expt x 2) (expt y 6))
(pe (((partial 1) f) (f 'x 'y) 'y))
;(* 3 (expt x 4) (expt y 8))

;;b tuple version
(pe (((partial 0) f1) (g1 (up 'x 'y))))  ;can receive output of other function that outputs tuple SAME
;(* 2 (expt x 2) (expt y 6))
(pe (((partial 1) f1) (g1 (up 'x 'y))))  ;SAME
;(* 3 (expt x 4) (expt y 8))





;;c
(pe (((partial 0) g) 'x 'y)) 
;(up (* 2 x (expt y 3)) 0)
(pe (((partial 1) g) 'x 'y)) 
;(up (* 3 (expt x 2) (expt y 2)) 1)

;;c tuple version
(pe (((partial 0) g1) (up 'x 'y))) 
;(up (* 2 x (expt y 3)) 0)
(pe (((partial 1) g1) (up 'x 'y))) 
;(up (* 3 (expt x 2) (expt y 2)) 1)



;;d
(pe ((D f) 'a 'b))
;(down (* 2 a (expt b 3)) (* 3 (expt a 2) (expt b 2)))
(pe ((D g) 3 5))
;(down (up 750 0) (up 675 1))
(pe ((D h) (* 3 (square 'a)) (* 5 (cube 'b))))
;(down (* 210937500 (expt a 6) (expt b 27))
;      (* 284765625 (expt a 8) (expt b 24)))


;;d tuple version
(pe((D f1) (up 'a 'b)))
;(down (* 2 a (expt b 3)) (* 3 (expt a 2) (expt b 2)))
(pe((D g1) (up 3 5)))
;(down (up 750 0) (up 675 1))
(pe ((D h1) (up (* 3 (square 'a)) (* 5 (cube 'b)))))
;(down (* 210937500 (expt a 6) (expt b 27))
;      (* 284765625 (expt a 8) (expt b 24)))



