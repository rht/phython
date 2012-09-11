(define ((em-state-derivative foo) state)
  (let ((q (ref state 1))
        (v (ref state 2))
        (a (ref state 3)))
    (let ((r (sqrt (square q))))
      (up 1
          v
          a
          (up (*  20 (+ (* 1 (ref a 0)) (/ (ref q 0) (cube r))))
              (*  20 (+ (* 1 (ref a 1)) (/ (ref q 1) (cube r)))))))))

(define winem (frame -1000 1000 -1000 1000))
(define ((monitor win) state)
  (let ((q (coordinate state)))
    (plot-point win (ref q 0) (ref q 1))))

((evolve em-state-derivative 1)
 (up 0
     (up 50.  0.)
     (up 0. 90.)
     (up -1  0.))
 (monitor winem)
 0.001
 200
 1e-12)
