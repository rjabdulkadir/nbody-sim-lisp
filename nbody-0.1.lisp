;Simulates the gravitational interaction between bodies.

(ql:quickload :numcl)
(ql:quickload :nodgui)

(defun get-acc (pos mass gravitational-constant softening)
  (let*
      (
       ; x, y and z coordinates
       (x (numcl:aref pos t '(0 1)))
       (y (numcl:aref pos t '(1 2)))
       (z (numcl:aref pos t '(2 3)))
       ; differences between positions of each body
       (dx (numcl:- (numcl:transpose x) x)) 
       (dy (numcl:- (numcl:transpose y) y))
       (dz (numcl:- (numcl:transpose z) z))
       ; square of distances between object, softening value to avoid distance of 0
       (inv-r3-tmp (numcl:asarray (numcl:+       
           (numcl:expt dx 2) (numcl:expt dy 2) (numcl:expt dz 2) (numcl:expt softening 2))))
       ; matrix to store 1/r^3 for all separations
       (inv-r3 (numcl:expt inv-r3-tmp -1.5))
       ; intermediate values arrays (for cleaner code, maybe remove)
       (dx-tmp (numcl:asarray (numcl:* dx inv-r3)))
       (dy-tmp (numcl:asarray (numcl:* dy inv-r3)))
       (dz-tmp (numcl:asarray (numcl:* dz inv-r3)))
       ; multiply by mass and gravitational constant
       (ax (numcl:* gravitational-constant (numcl:matmul dx-tmp mass)))
       (ay (numcl:* gravitational-constant (numcl:matmul dy-tmp mass)))
       (az (numcl:* gravitational-constant (numcl:matmul dz-tmp mass))))
    ; pack into a single array and return
    (numcl:concatenate (list ax ay az) :axis 1)))

(defun get-energy (pos vel mass gravitational-constant)
  (let* (
         ; calculate kinetic energy
         (kinetic-energy (* 0.5 (numcl:sum (numcl:* mass (numcl:expt vel 2)))))
         ; x, y and z coordinates - position difference as above
         (x (numcl:aref pos t '(0 1)))
         (y (numcl:aref pos t '(1 2)))
         (z (numcl:aref pos t '(2 3)))
         (dx (numcl:- (numcl:transpose x) x))
         (dy (numcl:- (numcl:transpose y) y))
         (dz (numcl:- (numcl:transpose z) z))
         ; distances between each body according to pythagoras' law
         (inv-r-tmp (numcl:asarray 
                     (numcl:+ 
                      (numcl:expt dx 2) (numcl:expt dy 2) (numcl:expt dz 2))))
         (inv-r-tmp (numcl:expt inv-r-tmp 0.5))
         ; get 1/distances for those above zero
         (idx-above-zero (numcl:where (numcl:triu inv-r-tmp 1) (lambda (x) (> x 0))))
         (inv-r (numcl:/ 1.0 (numcl:asarray (numcl:take inv-r-tmp idx-above-zero))))
         ; calculate mass interactions for each interaction
         (mass-itxn (numcl:asarray (numcl:take (numcl:* mass (numcl:transpose mass)) idx-above-zero)))
         ; calculate potential energy according to P=GMm/d^2
         (potential-energy (* gravitational-constant (numcl:sum (numcl:* (numcl:- mass-itxn) inv-r)))))
    (list kinetic-energy potential-energy)))


(defun main (&key (save-data nil))
  "n-body simulatrion"
  (nodgui:with-nodgui ()
    (nodgui:wm-title nodgui:*tk* "n-body Simulation")
    (let* ((particle-count 10)
           (current-time 0)  
           (end-time 1.0)
           (time-step 0.01)
           (softening 0.1)
           (gravitational-constant 1.0)
           (canvas (make-instance 'nodgui:canvas :width 1000 :height 800 :background :white))
           (add-oval #'(lambda (x y w c)
                         (nodgui:configure (nodgui:make-oval canvas x y (+ x w) (+ y w))        ; 2
                                    :fill c
                                    :width 0
                                    :tag "currentparticle")))
           ; generate initial conditions - masses, positions and velocities 
           (particle-masses (numcl:ones (list particle-count 1) :type 'double-float))
           (mass (numcl:/ (numcl:* 20 particle-masses) particle-count))
           (pos (numcl:normal 0d0 1d0 (list particle-count 3)))
           (vel (numcl:normal 0d0 1d0 (list particle-count 3)))
           ; convert to center of mass frame
           (center-of-mass 
                 (numcl:/ (numcl:mean (numcl:* mass vel) :axes 0) (numcl:mean mass)))
           (vel (numcl:- vel center-of-mass))
           ; calculate initial acceleration
           (acc (get-acc pos mass gravitational-constant softening))
           ; calculate initial kinetic and potential energies
           (energies (get-energy pos vel mass gravitational-constant))
           ; number of steps
           (number-of-steps (ceiling (/ end-time time-step)))
           ; create arrays to save  positions and energies
           (saved-positions (numcl:zeros (list particle-count 3 (+ number-of-steps 1)) :type 'double-float))           
           (saved-kinetic-energies (numcl:zeros (+ number-of-steps 1)  :type 'double-float)) 
           (saved-potential-energies (numcl:zeros (+ number-of-steps 1) :type 'double-float)))
      ; save initial positions and energies
      (setf (numcl:aref saved-positions '(0 100) '(0 3) 0) pos)
      (setf (numcl:aref saved-kinetic-energies 0) (first energies))
      (setf (numcl:aref saved-potential-energies 0) (second energies))
      ; setup window
      (nodgui:grid canvas 0 0 :sticky "news")
      (nodgui:grid-columnconfigure nodgui:*tk* 0 :weight 1)
      (nodgui:grid-rowconfigure nodgui:*tk* 0 :weight 1)
      ; loop over simulation steps
      (dotimes (i number-of-steps)
        (nodgui:clear canvas)
        (dotimes (j particle-count)
	  (let (
             (x-pos (numcl:aref pos j 0))
             (y-pos (numcl:aref pos j 1))
             (xs (numcl:aref saved-positions t '(0 1) (list (max (- i 30) 0) (+ i 1))))
             (ys (numcl:aref saved-positions t '(1 2) (list (max (- i 30) 0) (+ i 1)))))
             ; loop over each body - display trail
             (loop for a across (numcl:aref xs j 0)
                for b across (numcl:aref ys j 0)
                do (funcall add-oval (+ 500 (* a 100)) (+ 400 (* b 100)) 5 :gray75))
             (funcall add-oval (+ 500 (* x-pos 100)) (+ 400 (* y-pos 100)) 7 :blue)))
        ; calculate next values for vel. pos, acc and energies
        (setf vel (numcl:+ vel (numcl:* acc (numcl:/ time-step 2.0))))
        (setf pos (numcl:+ pos (numcl:* vel time-step)))
        (setf acc (get-acc pos mass gravitational-constant softening))
        (setf current-time (+ current-time time-step))
        (setf energies (get-energy pos vel mass gravitational-constant))
        ; save positions, and energies
        (setf (numcl:aref saved-positions '(0 100) '(0 3) (+ i 1)) pos)
        (setf (numcl:aref saved-kinetic-energies (1+ i)) (first energies))
        (setf (numcl:aref saved-potential-energies (1+ i)) (second energies))
        (format t "~a, Energies: KE:~,vf PE:~,vf KE+PE:~,vf~%" i 3 (first energies) 3 (second energies) 3 (apply '+ energies)))
      (if save-data
	  (with-open-file (stream "/home/riyadh/nbody.dat" :direction :output)
	    (dotimes (i number-of-steps)
	      (format stream "~a~%" (numcl:aref saved-positions 0))))))))

(main :save-data t)
