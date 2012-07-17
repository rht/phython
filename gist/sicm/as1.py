#ex 1.2-1.3
#a. Three Juggling pins      : 3[number] * (3[center of mass of each pin] + 3[orientation of each pin]) = 18
#b. Spherical pendulum       : [theta] + [phi] = 2
#c. Spherical double pendulum: 2[pendulum1's theta and phi] + 2[pendulum2's theta and phi] = 4
#d. Point mass on wire       : 1[distance travelled along the wire]
#e. attached symmetric top   : [yaw] + [pitch] + [roll] = 3
#f. attached asymmetric top  : [yaw] + [pitch] + [roll] = 3

#ex 1.4


 
#ex 1.5
# uh, python is awkward here
def L_harmonic(m, k):
    def Lh(local):
        q = coordinate(local)
        v = velocity(local)
        return 1/2 * m * square(v) - 1/2 * k * square(q)

(define win2 (frame 0. :pi/2 0. 1.2))

(define ((parametric-path-action Lagrangian t0 q0 t1 q1)
        intermediate-qs)
    (let ((path (make-path t0 q0 t1 q1 intermediate-qs)))
      ;; display path
      (graphics-clear win2)
      (plot-function win2 path t0 t1 (/ (- t1 t0) 100))
      ;; compute action
      (Lagrangian-action Lagrangian path t0 t1)))

(find-path (L-harmonic 1. 1.) 0. 1. :pi/2 0. 2)



