(restrict-compiler-policy 'speed 3 3)
(restrict-compiler-policy 'debug 0 0)
(restrict-compiler-policy 'safety 0 0)
(setf *block-compile-default* t)

;; ;(setf *features* (delete :cl-mpm-pic *features*))
;; ;(asdf:compile-system :cl-mpm/examples/tpb :force T)

;; (ql:quickload :cl-mpm)
;; (ql:quickload :cl-mpm/setup)
;; (ql:quickload :cl-mpm/particle)
;; (ql:quickload :cl-mpm/mpi)
;; (ql:quickload "magicl")
;; (ql:quickload "cl-mpm/examples/lbar")

(in-package :cl-mpm/examples/lbar)

(defmethod cl-mpm::update-node-forces ((sim cl-mpm::mpm-sim))
  (with-accessors ((damping cl-mpm::sim-damping-factor)
                   (mass-scale cl-mpm::sim-mass-scale)
                   (mesh cl-mpm::sim-mesh)
                   (dt cl-mpm::sim-dt))
      sim
    (cl-mpm::iterate-over-nodes
     mesh
     (lambda (node)
       (cl-mpm::calculate-forces-cundall-conservative node damping dt mass-scale)))))

(defmacro rank-0-time (rank &rest body)
  `(if (= ,rank 0)
      (time
        (progn
          ,@body))
      (progn
        ,@body)))

(defun rectangle-sdf (position size)
  (lambda (pos)
    (let* ((pos (magicl:from-list (list
                                        (magicl:tref pos 0 0)
                                        (magicl:tref pos 1 0)) '(2 1) :type 'double-float))
           (position (magicl:from-list position '(2 1) :type 'double-float))
           (dist-vec (magicl:.- (magicl:map! #'abs (magicl:.- pos position))
                                (magicl:from-list size '(2 1) :type 'double-float))))

      (+ (sqrt (magicl::sum
                (magicl:map! (lambda (x) (* x x))
                             (magicl:map! (lambda (x) (max 0d0 x)) dist-vec))))
         (min (max (magicl:tref dist-vec 0 0)
                   (magicl:tref dist-vec 1 0)
                   ) 0d0)))))

(defun apply-pullout (sim load-mps push-rate)
  (with-accessors ((mps cl-mpm:sim-mps)
                   (mesh cl-mpm:sim-mesh))
      sim
    (loop for mp in load-mps
          do
             (cl-mpm::iterate-over-neighbours
              mesh
              mp
              (lambda (mesh mp node svp grad fsvp fgrad)
                (with-accessors ()
                    mp
                  (with-accessors ((pos cl-mpm/mesh::node-position)
                                   (vel cl-mpm/mesh::node-velocity)
                                   (active cl-mpm/mesh::node-active)
                                   (acc cl-mpm/mesh::node-acceleration)
                                   )
                      node
                    (when active
                      (setf (magicl:tref vel 1 0) push-rate)))))))))

(defparameter *current-load* 0d0)
(defun apply-force (sim load-mps push-rate)
  (with-accessors ((mps cl-mpm:sim-mps)
                   (mesh cl-mpm:sim-mesh))
      sim
    (loop for mp in load-mps
          do
             (setf (magicl:tref (cl-mpm/particle:mp-body-force mp) 1 0)
                   (* *current-load* (/ 1d0 (* (cl-mpm/particle:mp-volume mp) (length load-mps))))))
    ))

(defun energy-norm (sim)
  (/ (loop for mp across (cl-mpm:sim-mps *sim*)
          sum (magicl:norm (cl-mpm/particle:mp-velocity mp)))
     (length (cl-mpm:sim-mps *sim*))))

(defun get-disp (load-mps)
  (if (and load-mps
           (> (length load-mps) 0))
      (coerce  (- (/ 
                    (loop for mp in load-mps
                           sum (-
                                 (magicl:tref (cl-mpm/particle::mp-position mp) 1 0)
                                 (* 0.5d0 (magicl:tref (cl-mpm/particle::mp-domain-size mp) 1 0))
                                 )) 
                    (length load-mps))
                  *initial-surface*
                  ) 'double-float)
    0d0))

(defun get-force-mps (sim load-mps)
  (with-accessors ((mps cl-mpm:sim-mps)
                   (mesh cl-mpm:sim-mesh))
      sim
    (let ((force 0d0))
      (loop for mp in load-mps
            do
               (cl-mpm::iterate-over-neighbours
                mesh
                mp
                (lambda (mesh mp node svp &rest args)
                  (incf force
                        (* svp
                           (cl-mpm/fastmath::mag (cl-mpm/mesh::node-force node))
                           )
                        ))))
      (/ force (length load-mps))
      )))

(defun get-reaction-force (load-nodes)
  ;; (cl-mpm/fastmath::mag (cl-mpm/mesh::node-force (nth 0 load-nodes)))
  (loop for mp in load-nodes
        sum
        ;; (cl-mpm/fastmath::mag (cl-mpm/mesh::node-force mp))
        (- (magicl:tref (cl-mpm/mesh::node-force mp) 1 0))
        ))

(defparameter *target-displacement* 0d0)
(defun apply-disp-penalty (sim load-mps)
  (with-accessors ((mps cl-mpm:sim-mps)
                   (mesh cl-mpm:sim-mesh))
      sim
    (let* ((penalty 1d5)
           (displacement *target-displacement*)
           (pos
             (get-disp load-mps))
           (force (* penalty (- displacement pos))))
      (incf *current-load* force)
      (loop for mp in load-mps
            do
               (setf (magicl:tref (cl-mpm/particle:mp-body-force mp) 1 0)
                     (* force (/ 1d0 (* (cl-mpm/particle:mp-volume mp) (length load-mps))))
                     ))
      )))

;(defparameter *tip-velocity* -0.02d-3)
(defparameter *tip-velocity* -0.000d-3)

(defclass mpm-sim-quasi-static-damage (;cl-mpm::mpm-sim-quasi-static 
               ;cl-mpm/mpi::mpm-sim-usl-mpi-nodes-damage
               ;cl-mpm/mpi::mpm-sim-mpi-nodes-damage
               
               ) ())

(defun setup-test-column (size block-size offset &optional (e-scale 1) (mp-scale 1)
                                                   (kappa-scale 1d0)
                                                   (lc-scale 1d0)
                                                   )
  (let* ((sim (cl-mpm/setup::make-block
               (/ 1d0 e-scale)
               (mapcar (lambda (x) (* x e-scale)) size)
               ;; 'cl-mpm::mpm-sim-usf
               :sim-type 'cl-mpm/mpi::mpm-sim-usl-mpi-nodes-damage

               ;; 'cl-mpm/damage::mpm-sim-damage
               ;; 'mpm-sim-quasi-static-damage
               ;'cl-mpm/mpi::mpm-sim-usl-mpi-nodes-damage
               ;; 'cl-mpm/mpi::mpm-sim-mpi-nodes-damage
               ;; 'cl-mpm/mpi::mpm-sim-mpi-nodes
               ))
         (h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim)))
         (h-x (/ h 1d0))
         (h-y (/ h 1d0))
         (density 2.2d3)
         (elements (mapcar (lambda (s) (* e-scale (/ s 2))) size))
         )
    (declare (double-float h density))
    (progn
      (let* ((scaler (/ 1d0 (sqrt 1))))
        (setf (cl-mpm:sim-mps sim)
              (cl-mpm/setup::make-mps-from-list
               (append
                (cl-mpm/setup::make-block-mps-list
                 offset
                 block-size
                 (mapcar (lambda (e mp-s)
                           (round  (* e mp-s) h-x)
                           ) block-size
                             (list mp-scale 
                                   mp-scale 
                                   mp-scale 
                                   ;2
                                   )
                             )
                 density
                 'cl-mpm/particle::particle-limestone-delayed
                 :E 25.85d9
                 ;:E 21.00d9
                 :nu 0.18d0
                 ;; :elastic-approxmation :plane-stress
                 :fracture-energy 95d0
                 :initiation-stress (* 2.7d6 1d0)
                 :critical-damage 1.000d0
                 ;; :internal-length (* 25d-3 1d0 scaler)
                 :local-length (* 25d-3 scaler)
                 :local-length-damaged (* 25d-3 scaler)
                 ;:ductility 15.4d0;6.8d0
                 ;:ductility 12.6d0
                 ;:ductility 6.8d0
                 :ductility 26.85d0
                 
                 :compression-ratio 10d0
                 :delay-time 1.0d0
                 :enable-plasticity (= lc-scale 1d0)
                 :phi (* 35d0 (/ pi 180))
                 :psi 0d0
                 :c 3000d3
                                  
                 ;; :local-length-damaged 0.01d0
                 :gravity 0.0d0
                 :gravity-axis (cl-mpm/utils:vector-from-list '(0d0 0d0 0d0))
                 )
                ;; impactors
                )
               )))
      (format t "Plasticity enabled: ~A~%" (= lc-scale 1d0))
      (setf (cl-mpm:sim-allow-mp-split sim) nil)
      (setf (cl-mpm::sim-enable-damage sim) nil)
      (setf (cl-mpm::sim-nonlocal-damage sim) t)
      (setf (cl-mpm::sim-allow-mp-damage-removal sim) nil)
      (setf (cl-mpm::sim-mp-damage-removal-instant sim) nil)
      (setf (cl-mpm::sim-mass-filter sim) 1d-10)
      (let ((ms 1d0))
        (setf (cl-mpm::sim-mass-scale sim) ms)
        (setf (cl-mpm:sim-damping-factor sim)
              ;(* 1d-2 density ms)
              ;; 1d0
              ;0.05d0
              ;(* 8d0 density)
              ;; (* 4.0d0 density)
              ;(* 2d-3 (cl-mpm/setup::estimate-critical-damping sim))
              0.70d0

              ))

      (dotimes (i 0)
        (dolist (dir (list :x :y))
          (cl-mpm::split-mps-criteria
           sim
           (lambda (mp h)
             (when
                 (and
                  (> (magicl:tref (cl-mpm/particle:mp-position mp) 0 0)
                     (+ (first offset) (* (first block-size) 0.45))
                     )
                  (< (magicl:tref (cl-mpm/particle:mp-position mp) 0 0)
                     (+ (first offset) (* (first block-size) 0.55))
                     )
                  )
               dir
               )))))

      (let ((dt-scale 0.5d0))
        (setf
         (cl-mpm:sim-dt sim)
         (* dt-scale h
            (sqrt (cl-mpm::sim-mass-scale sim))
            (sqrt (/ density (cl-mpm/particle::mp-p-modulus (aref (cl-mpm:sim-mps sim) 0)))))))

      (format t "Estimated dt ~F~%" (cl-mpm:sim-dt sim))

      (let ((cut-size 0.25d0))
        (cl-mpm/setup::remove-sdf
         sim
         (rectangle-sdf
          (list
           (+ (first offset) (second block-size))
           (+ (second offset) 0d0))
          (list
           cut-size
           cut-size
           ))))

      (let* ((crack-width 0.00d0)
             (crack-pos ;; Max x pos
               (loop for mp across (cl-mpm:sim-mps sim)
                     when
                     (not (= (cl-mpm/particle::mp-index mp) 1))
                     maximize (magicl:tref (cl-mpm/particle:mp-position mp) 0 0)))
             (above-crack
               (loop for mp across (cl-mpm:sim-mps sim)
                     when
                     (and
                      (not (= (cl-mpm/particle::mp-index mp) 1))
                      (>= (magicl:tref (cl-mpm/particle:mp-position mp) 0 0)
                          (- crack-pos crack-width)))
                     collect mp)
               )
             (min-pos (loop for mp in above-crack
                            minimize (magicl:tref (cl-mpm/particle:mp-position mp) 1 0)))
             )
        (defparameter *terminus-mps*
          (loop for mp in above-crack
                when (= min-pos
                        (magicl:tref
                         (cl-mpm/particle:mp-position mp)
                         1 0))
                  collect mp)))

      (loop for mp in *terminus-mps*
            do (setf (cl-mpm/particle::mp-index mp) 1))

      (let ((left-node-pos
              (list
               (round (first offset) h-x)
               (round (second offset) h-x)
               0))
            (right-node-pos
              (list
               (round (+ (first offset) (first block-size)) h-x)
               (round (second offset) h-x)
               0)
              ))
        (defparameter *fixed-nodes* (mapcar (lambda (id) (cl-mpm/mesh:get-node (cl-mpm:sim-mesh sim)
                                                                                     id))
                                                  (list left-node-pos right-node-pos)
                                                  ))
        (format t "Fixed node ~A ~%" left-node-pos)
        (format t "Roller node ~A ~%" right-node-pos)
        (setf (cl-mpm:sim-bcs sim)
              (cl-mpm/bc::make-bcs-from-list
               (append
                (cl-mpm/bc::make-outside-bc-var-list
                 (cl-mpm:sim-mesh sim)
                 (lambda (i) (cl-mpm/bc::make-bc-fixed i '(0 nil nil)))
                 (lambda (i) (cl-mpm/bc::make-bc-fixed i '(0 nil nil)))
                 (lambda (i) (cl-mpm/bc::make-bc-fixed i '(nil 0 nil)))
                 (lambda (i) (cl-mpm/bc::make-bc-fixed i '(0 0 nil)))
                 (lambda (i) (cl-mpm/bc::make-bc-fixed i '(nil nil 0)))
                 (lambda (i) (cl-mpm/bc::make-bc-fixed i '(nil nil 0)))
                 )
                (list
                 ;(cl-mpm/bc::make-bc-fixed left-node-pos
                 ;                          '(0 0 nil))

                 ;(cl-mpm/bc::make-bc-fixed right-node-pos
                 ;                          '(nil 0 nil))
                 )
                ))))
      (defparameter *initial-surface*
        (loop for mp in *terminus-mps*
              when (= 1 (cl-mpm/particle::mp-index mp))
              minimizing (magicl:tref
                          (magicl:.- (cl-mpm/particle:mp-position mp)
                                   (magicl:scale (cl-mpm/particle::mp-domain-size mp) 0.5d0)
                                   )
                        1 0)))

      (format t "~A~%" h-x)
      (setf (cl-mpm::sim-bcs-force-list sim)
            (list
             (cl-mpm/bc:make-bcs-from-list
              (list
               ;; *floor-bc*
               (cl-mpm/bc::make-bc-closure
                nil
                (lambda ()
                  (with-accessors ((mesh cl-mpm:sim-mesh)
                                   (dt cl-mpm::sim-dt))
                      sim
                    (let ((datum (* 1d0 (+ *initial-surface* *target-displacement*)))
                          (normal (cl-mpm/utils:vector-from-list  '(0d0 1d0 0d0))))
                      (setf *terminus-mps*
                            (loop for mp across (cl-mpm:sim-mps *sim*)
                                  when (= (cl-mpm/particle::mp-index mp) 1)
                                    collect mp))
                      (cl-mpm/penalty::apply-displacement-control-mps 
                        mesh 
                        (coerce *terminus-mps* 'vector)
                        dt
                        normal
                        datum
                        ;(* density 1d5)
                        (* 25.85d9 1d0)
                        0d0)
                      ))))))))
      sim))
  )
(defparameter *sim* nil)
(defparameter *run-sim* t)
(defparameter *t* 0)
(defparameter *sim-step* 0)
(defparameter *refine* (/ 1d0 2d0))
(defun setup (&optional (kappa-scale 1d0)
                (lc-scale 1d0)
                (refine 1d0))
  (let* ((mesh-size (/ 0.025 refine))
         (mps-per-cell 6);; (floor (* 2 kappa-scale)) 
         (shelf-height 0.500d0)
         (shelf-length 0.500d0)
         (domain-length (+ shelf-length 0.2d0))
         (domain-height (+ shelf-length 0.1d0))
         (offset (list
                   0.1d0
                   0d0
                   ;0d0
                   )))
    (defparameter *sim*
      (setup-test-column (list domain-length
                               domain-height
                               ;(* 3 mesh-size)
                               )
                         (list shelf-length 
                               shelf-height
                               ;mesh-size
                               )
                         offset
                         (/ 1d0 mesh-size) mps-per-cell
                         kappa-scale
                         lc-scale
                         ))
    ;; (let ((cut-size 0.25d0))
    ;;   (cl-mpm/setup::remove-sdf
    ;;    *sim*
    ;;    (rectangle-sdf
    ;;     (list
    ;;      (+ (first offset) shelf-height)
    ;;      (+ (second offset) 0d0))
    ;;     (list
    ;;      cut-size
    ;;      cut-size
    ;;      ))))
    (defparameter *current-load* 0d0)
    ;; (loop for mp across (cl-mpm:sim-mps *sim*)
    ;;       do
    ;;          (setf (cl-mpm/particle:mp-damage mp) (random 0.1d0)))
    ;; (cl-mpm/setup::damage-sdf
    ;;  *sim*
    ;;  (lambda (p)
    ;;    (cl-mpm/setup::line-sdf (magicl:from-list (list (magicl:tref p 0 0)
    ;;                                                    (magicl:tref p 1 0))
    ;;                                              '(2 1))
    ;;                            (list (- shelf-length shelf-height) shelf-height)
    ;;                            (list shelf-length soil-boundary)
    ;;                            10d0
    ;;                            )) 0.8d0)
    ;(let ((sdf
    ;        (lambda (p)
    ;          (cl-mpm/setup::line-sdf (magicl:from-list (list (magicl:tref p 0 0)
    ;                                                          (magicl:tref p 1 0))
    ;                                                    '(2 1))
    ;                                  (list (- shelf-length shelf-height) shelf-height)
    ;                                  (list shelf-length 0d0)
    ;                                  20d0
    ;                                  ))
    ;        ))
    ;  (loop for mp across (cl-mpm:sim-mps *sim*)
    ;        do (with-accessors ((pos cl-mpm/particle:mp-position)
    ;                            (damage cl-mpm/particle:mp-damage)) mp
    ;             (when (>= 0 (funcall sdf pos))
    ;               (setf damage (min 1d0 (max 0d0 (coerce (* (funcall sdf pos) -0.1d0) 'double-float)))))
    ;             )))


    )
  (defparameter *target-displacement* 0d0)
  (format t "MPs: ~D~%" (length (cl-mpm:sim-mps *sim*)))
  (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./outframes/")) do (uiop:delete-file-if-exists f))
  (defparameter *run-sim* t)
  (defparameter *t* 0)
  (defparameter *sim-step* 0))


(defparameter *data-force* '())
(defparameter *data-displacement* '(0d0))
(defparameter *data-load* '(0d0))
(defparameter *data-time* '(0d0))
(defparameter *data-node-load* '(0d0))
(defparameter *target-displacement* 0d0)
(defparameter *data-averaged* t)

(defparameter *data-full-time* '(0d0))
(defparameter *data-full-load* '(0d0))

(defun run ()
  (cl-mpm/output:save-vtk-mesh (merge-pathnames "output/mesh.vtk")
                          *sim*)
  (cl-mpm/output::save-vtk-nodes (merge-pathnames (format nil "output/sim_nodes.vtk")) *sim*)
                       
  (defparameter *data-force* '())
  (defparameter *data-displacement* '(0d0))
  (defparameter *data-load* '(0d0))
  (defparameter *data-node-load* '(0d0))
  (defparameter *data-mp-load* '(0d0))
  (defparameter *data-time* '(0d0))
  (defparameter *target-displacement* 0d0)
  (defparameter *data-full-time* '(0d0))
  (defparameter *data-full-load* '(0d0))

  (with-open-file (stream (merge-pathnames "output/disp.csv") :direction :output :if-exists :supersede)
    (format stream "disp,load ~%")
    (format stream "~f,~f~%" 0d0 0d0)
    )

  (let* ((target-time 0.5d0)
         (dt (cl-mpm:sim-dt *sim*))
         (substeps (floor target-time dt))
         (dt-scale 0.1d0)
         (disp-step 0.008d-3)
         )

    (setf cl-mpm/penalty::*debug-force* 0d0)
    (setf cl-mpm/penalty::*debug-force-count* 0d0)
    (time (cl-mpm::update-sim *sim*))

    (format t "Calculating dt~%")
    (multiple-value-bind (dt-e substeps-e) (cl-mpm:calculate-adaptive-time *sim* target-time :dt-scale dt-scale)
                    (format t "CFL dt estimate: ~f~%" dt-e)
                    (format t "CFL step count estimate: ~D~%" substeps-e)
                    (setf substeps substeps-e))
    (format t "Substeps ~D~%" substeps)
    ;; (incf *target-displacement* -0.000d-3)
    ;(incf *target-displacement* disp-step)
    (setf *target-displacement* 0d0)
    (time (loop for steps from 0 to 100
                while *run-sim*
                do
                   (progn
                     (format t "Step ~d ~%" steps)
                     (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" *sim-step*)) *sim*)
                     (let ((average-force 0d0)
                           (average-disp 0d0)
                           (average-reaction 0d0))
                       (time
                        (dotimes (i substeps);)
                          (incf average-force (/
                                               cl-mpm/penalty::*debug-force*
                                               substeps
                                               ))
                          (incf average-reaction
                                (/ (get-reaction-force *fixed-nodes*) substeps)
                                )
                          (incf average-disp
                                (/ (get-disp *terminus-mps*) substeps)
                                )
                          (incf average-reaction (/ (get-reaction-force *fixed-nodes*) substeps))
                          (setf cl-mpm/penalty::*debug-force* 0d0)
                          (setf cl-mpm/penalty::*debug-force-count* 0d0)
                          (cl-mpm::update-sim *sim*)
                          (incf *target-displacement* (/ disp-step substeps))
                          (setf *t* (+ *t* (cl-mpm::sim-dt *sim*))))
                        )
                       ;; (incf *target-displacement* -0.01d-3)
                       (push
                        average-disp
                        *data-displacement*)
                       ;; (push
                       ;;  average-force
                       ;;  *data-mp-load*)
                       (push
                        average-reaction
                        *data-load*)

                       ;; (format t "Target load: ~f~%" (* *target-displacement* 20d9))
                       ;; (format t "Current load: ~f~%" (* (get-disp *terminus-mps*) 1d9))
                       ;; (format t "Node Load: ~f~%" (get-reaction-force *fixed-nodes*))
                       ;; (format t "Pen load: ~f~%" average-force)
                       (let ((mesh-size (cl-mpm/mesh::mesh-resolution (cl-mpm:sim-mesh *sim*))))
                         (with-open-file (stream (merge-pathnames "output/disp.csv") :direction :output :if-exists :append)
                           (format stream "~f,~f,~f~%"
                                   average-disp
                                   (/ average-force mesh-size)
                                   (/ average-reaction mesh-size)
                                   ))))


                     (format t "Target: ~f - Current: ~f Error: ~f - energy ~F~%"
                             *target-displacement*
                             (get-disp *terminus-mps*)
                             (* 100d0 (/ (- *target-displacement* (get-disp *terminus-mps*)) *target-displacement*))
                             (energy-norm *sim*))
                     ;; (format t "Debug load: ~f~%" cl-mpm/penalty::*debug-force*)
                     ;; (format t "Debug load count: ~f~%" cl-mpm/penalty::*debug-force-count*)

                     ;; (incf *target-displacement* -0.01d-3)

                     ;; (incf *current-load* (* -1d3 1d0))

                     ;; (loop for mp across (cl-mpm:sim-mps *sim*)
                     ;;       do
                     ;;          (when (>= (cl-mpm/particle:mp-damage mp) 1d0)
                     ;;            (let ((ms 1d2))
                     ;;              (setf (cl-mpm::sim-mass-scale *sim*) ms
                     ;;                    ;; target-time 1d0
                     ;;                    (cl-mpm:sim-damping-factor *sim*) (* 1d-8 ms)
                     ;;                    ))))

                     (multiple-value-bind (dt-e substeps-e) (cl-mpm:calculate-adaptive-time *sim* target-time :dt-scale dt-scale)
                       (format t "CFL dt estimate: ~f~%" dt-e)
                       (format t "CFL step count estimate: ~D~%" substeps-e)
                       (setf substeps substeps-e))
                     ;; (setf (cl-mpm:sim-damping-factor *sim*)
                     ;;       (* (cl-mpm:sim-damping-factor *sim*) (expt 1d-3 1/40)))

                     (incf *sim-step*)
                     (swank.live:update-swank)
                     ))))
  (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" *sim-step*)) *sim*))
(defun main-loop ()
  (let ((kappa (uiop:getenv "KAPPA"))
        (lc (uiop:getenv "lc"))
        (refine (uiop:getenv "REFINE")))
    (setf kappa (if kappa
                    (parse-float:parse-float kappa)
                    1d0))
    (setf lc (if lc
                 (parse-float:parse-float lc)
                 1d0))
    (setf refine (if refine
                     (parse-float:parse-float refine)
                     1d0))
    (setf *refine* refine)
    (format t "Kappa ~F - LC ~F - Refine ~F ~%" kappa lc refine)
    (setup kappa lc refine)
    (run-static (format nil "output-~f-~f-~f/" refine lc kappa))))
(defun mpi-loop ()
  (let ((rank (cl-mpi:mpi-comm-rank)))
    ;; (cl-mpi:mpi-init)
    (let ((rank (cl-mpi:mpi-comm-rank))
          (kappa (uiop:getenv "KAPPA"))
          (lc (uiop:getenv "lc"))
          (refine (uiop:getenv "REFINE")))
      ;; (cl-mpi:mpi-init)
      (setf kappa (if kappa 
                      (parse-float:parse-float kappa)
                      1d0))
      (setf lc (if lc 
                   (parse-float:parse-float lc)
                   1d0))
      (setf refine (if refine
                       (parse-float:parse-float refine)
                       1d0))
      (setf *refine* refine)
      (format t "Kappa ~F - LC ~F - Refine ~F ~%" kappa lc refine)
      ;(setup kappa lc refine)
      (setup kappa lc refine)
      (format t "Mesh size:~F~%" (cl-mpm/mesh::mesh-resolution (cl-mpm:sim-mesh *sim*)))
      (format t "Terminus mps:~F~%" (length *terminus-mps*))

      ;(let ((dsize (floor (sqrt (floor (cl-mpi:mpi-comm-size))))))
      ;  (setf (cl-mpm/mpi::mpm-sim-mpi-domain-count *sim*) (list dsize dsize 1)))
      (let ((dsize (floor (cl-mpi:mpi-comm-size))))
       (setf (cl-mpm/mpi::mpm-sim-mpi-domain-count *sim*) (list dsize 1 1)))

      ;; (let ((dsize (floor (sqrt (cl-mpi:mpi-comm-size)))))
      ;;   (setf (cl-mpm/mpi::mpm-sim-mpi-domain-count *sim*) (list dsize dsize 1)))

      (let ((dhalo-size (* 1 (cl-mpm/particle::mp-local-length (aref (cl-mpm:sim-mps *sim*) 0)))))
	      (setf (cl-mpm/mpi::mpm-sim-mpi-halo-damage-size *sim*) dhalo-size))
      (when (= rank 0)
        (format t "Damage halo size: ~f~%" (cl-mpm/mpi::mpm-sim-mpi-halo-damage-size *sim*))
        (format t "Sim MPs: ~a~%" (length (cl-mpm:sim-mps *sim*)))
        (format t "Decompose~%"))
      (cl-mpm/mpi::domain-decompose *sim*)
; 	  (cl-mpm/mpi::domain-decompose
; 	    *sim*
; 	    :domain-scaler
; 	    (lambda (domain)
; 	  	(destructuring-bind (x y z) domain
; 	  	  (let ((dnew (list (mapcar (lambda (p)
;                                      (if (and 
;                                            (> p 0d0)
;                                            (< p 1d0))
;                                          (* (+ p (/ 0.1 0.7)) 
;                                             )
;                                         p
;                                         ))
;                                       x)
; 	  						 y z)))
; 	  		(format t "Domain ~A ~A~%" (list x y z) dnew)
; 	  		dnew))))

      (format t "Rank: ~D - Sim MPs: ~a~%" rank (length (cl-mpm:sim-mps *sim*)))
 
	  ;(cl-mpm/mpi::domain-decompose *sim* :domain-scaler (lambda (domain)
	  ;													 (format t "~A~%" domain)
	  ;													 ;domain
	  ;													 (destructuring-bind (x y z) domain
	  ;													   (list
	  ;														 (mapcar
	  ;															  (lambda (p)
	  ;																(expt p 1.3)
	  ;																) x)
	  ;															 y z))))
      (when (= rank 0)
        (format t "Sim MPs: ~a~%" (length (cl-mpm:sim-mps *sim*))))
      (when (= rank 0)
        (format t "Run mpi~%"))

      (if (= rank 0)
          (time (cl-mpm::update-sim *sim*))
          (cl-mpm::update-sim *sim*))

      (run-mpi (format nil "output-~f-~f-~f/" refine lc kappa))
      ;; (run-mpi-static (format nil "output-~f-~f-~f/" refine lc kappa))
      (when (= rank 0)
        (format t "Done mpi~%")))
    )
  (cl-mpi:mpi-finalize)
  (sb-ext:quit)
  )
(defun mpi-average (value mp-count)
  (let ((sum 0d0))
    (static-vectors:with-static-vector (source 1 :element-type 'double-float :initial-element (coerce value 'double-float))
      (static-vectors:with-static-vector (dest 1 :element-type 'double-float :initial-element 0d0)
        (cl-mpi:mpi-allreduce source dest cl-mpi:+mpi-sum+)
        (setf sum (aref dest 0))
        ))

    (static-vectors:with-static-vector (source 1 :element-type 'double-float :initial-element (coerce mp-count 'double-float))
      (static-vectors:with-static-vector (dest 1 :element-type 'double-float :initial-element 0d0)
        (cl-mpi:mpi-allreduce source dest cl-mpi:+mpi-sum+)
        (setf sum (/ sum (max 1d0 (aref dest 0))))
        ))
    sum))
(defun mpi-sum (value)
  (let ((sum 0d0))
    (static-vectors:with-static-vector (source 1 :element-type 'double-float :initial-element value)
      (static-vectors:with-static-vector (dest 1 :element-type 'double-float :initial-element 0d0)
        (cl-mpi:mpi-allreduce source dest cl-mpi:+mpi-sum+)
        (setf sum (aref dest 0))))
    sum))
(defun run-mpi (&optional (output-folder "output/"))
  (ensure-directories-exist (merge-pathnames output-folder))
  (defparameter *data-force* '())
  (defparameter *data-displacement* '(0d0))
  (defparameter *data-load* '(0d0))
  (defparameter *data-node-load* '(0d0))
  (defparameter *data-mp-load* '(0d0))
  (defparameter *data-time* '(0d0))
  (defparameter *target-displacement* 0d0)
  (defparameter *data-full-time* '(0d0))
  (defparameter *data-full-load* '(0d0))

  (defparameter *terminus-mps*
    (loop for mp across (cl-mpm:sim-mps *sim*)
          when (= (cl-mpm/particle::mp-index mp) 1)
            collect mp))

  (let ((mass-scale 1d4)
        (h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh *sim*))))
    (setf (cl-mpm::sim-mass-scale *sim*) mass-scale)
    (setf (cl-mpm:sim-damping-factor *sim*)
          (* 1d2
             (expt h 2))))


  (setf (cl-mpm:sim-dt *sim*)
        (cl-mpm/setup:estimate-elastic-dt *sim* :dt-scale 0.8d0))

  (setf (cl-mpm::sim-enable-damage *sim*) t)
  ;; (setf (cl-mpm/damage::sim-damage-delocal-counter-max *sim*) 40)
  (setf (cl-mpm/damage::sim-damage-delocal-counter-max *sim*) (round (* 40 *refine*)))


  (let* ((target-time 0.1d0)
         (dt (cl-mpm:sim-dt *sim*))
         (dt-scale 0.8d0)
         (substeps (floor target-time dt))
         (rank (cl-mpi:mpi-comm-rank))
         (load-steps 100)
         (disp-total 0.8d-3)
         (disp-step (/ disp-total load-steps)))

    (when (= rank 0)
      (format t "Outputting mesh and load-disp graph")

      ;; (cl-mpm/output:save-vtk-mesh (merge-pathnames output-folder "mesh.vtk") *sim*)
      ;; (cl-mpm/output::save-vtk-nodes (merge-pathnames output-folder (format nil "sim_nodes.vtk")) *sim*)

      (with-open-file (stream (merge-pathnames output-folder "disp.csv") :direction :output :if-exists :supersede)
        (format stream "disp,load,reaction~%")
        (format stream "~f,~f,~f~%" 0d0 0d0 0d0)
        ))

    (setf cl-mpm/penalty::*debug-force* 0d0)
    (setf cl-mpm/penalty::*debug-force-count* 0d0)
    (when (= rank 0)
      (format t "Test run~%"))
    (cl-mpm::update-sim *sim*)
    ;(when (= rank 0)
    ;  (format t "Caluclating dt estimate~%"))
    ;(multiple-value-bind (dt-e substeps-e) (cl-mpm:calculate-adaptive-time *sim* target-time :dt-scale dt-scale)
    ;  (when (= rank 0)
    ;    (format t "CFL dt estimate: ~f~%" dt-e)
    ;    (format t "CFL step count estimate: ~D~%" substeps-e))
    ;  (setf substeps substeps-e))
    (when (= rank 0)
      (format t "Substeps ~D~%" substeps))
    (setf cl-mpm/penalty::*debug-force* 0d0)
    (time (loop for steps from 0 to 100
                while *run-sim*
                do
                (progn
                  (when (= rank 0)
                    (format t "Step ~d ~%" steps))
                  (cl-mpm/output:save-vtk (merge-pathnames output-folder (format nil "sim_~2,'0d_~5,'0d.vtk" rank *sim-step*)) *sim*)
                  ;; (cl-mpm/output::save-vtk-nodes (merge-pathnames output-folder (format nil "sim_nodes_~2,'0d_~5,'0d.vtk" rank *sim-step*)) *sim*)
                  (let ((average-force 0d0)
                        (average-disp 0d0)
                        (average-reaction 0d0))
                    (time
                      (dotimes (i substeps)
                        (defparameter *terminus-mps*
                          (loop for mp across (cl-mpm:sim-mps *sim*)
                                when (= (cl-mpm/particle::mp-index mp) 1)
                                collect mp))
                        (incf average-force (/
                                              cl-mpm/penalty::*debug-force*
                                              substeps
                                              ))
                        (incf average-reaction
                              (/ (get-reaction-force *fixed-nodes*) substeps)
                              )
                        (incf average-disp
                              (/ (get-disp *terminus-mps*) substeps)
                              )
                        (setf cl-mpm/penalty::*debug-force* 0d0)
                        (cl-mpm::update-sim *sim*)
                        (incf *target-displacement* (/ disp-step substeps))
                        (setf *t* (+ *t* (cl-mpm::sim-dt *sim*))))
                      )
                    (setf average-disp (cl-mpm/mpi::mpi-sum average-disp))
                    (setf average-force (cl-mpm/mpi::mpi-sum average-force))
                    (push
                      average-disp
                      *data-displacement*)
                    (push
                      average-force
                      *data-load*)

                    (when (= rank 0)
                      (let* ((mesh-size (cl-mpm/mesh::mesh-resolution (cl-mpm:sim-mesh *sim*)))
                             (mesh-size 1d0))
                        (with-open-file (stream (merge-pathnames output-folder "disp.csv") :direction :output :if-exists :append)
                          (format stream "~f,~f,~f~%"
                                  average-disp
                                  (/ average-force mesh-size)
                                  0d0
                                  ;(/ average-reaction mesh-size)
                                  ))))

                    (when (= rank 0)
                      (format t "Disp Target: ~E - Current: ~E~%"
                              *target-displacement*
                              average-disp
                              )
                      ))
                     (multiple-value-bind (dt-e substeps-e) (cl-mpm:calculate-adaptive-time *sim* target-time :dt-scale dt-scale)
                       (when (= rank 0)
                         (format t "CFL dt estimate: ~f~%" dt-e)
                         (format t "CFL step count estimate: ~D~%" substeps-e))
                      (setf substeps substeps-e))
                     (incf *sim-step*)
                     (swank.live:update-swank)
                     ))))
  ;(cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" *sim-step*)) *sim*)
  )



(defun run-static (&optional (output-folder "output/"))
  (ensure-directories-exist (merge-pathnames output-folder))
  (defparameter *data-force* '())
  (defparameter *data-displacement* '(0d0))
  (defparameter *data-load* '(0d0))
  (defparameter *data-node-load* '(0d0))
  (defparameter *data-mp-load* '(0d0))
  (defparameter *data-time* '(0d0))
  (defparameter *target-displacement* 0d0)
  (defparameter *data-full-time* '(0d0))
  (defparameter *data-full-load* '(0d0))

  (change-class *sim* 'cl-mpm/damage::mpm-sim-damage)
  (loop for mp across (cl-mpm:sim-mps *sim*)
        do (change-class mp 'cl-mpm/particle::particle-limestone))


  (let* ((target-time 0.2d0)
         (dt (cl-mpm:sim-dt *sim*))
         (substeps (floor target-time dt))
         (dt-scale 0.5d0)
         (load-steps 10)
         (disp-step (/ 0.8d-3 load-steps)))

    (format t "Outputting mesh and load-disp graph")
    ;; (cl-mpm/output:save-vtk-mesh (merge-pathnames output-folder "mesh.vtk") *sim*)
    (with-open-file (stream (merge-pathnames output-folder "disp.csv") :direction :output :if-exists :supersede)
      (format stream "disp,load,nload~%")
      (format stream "~f,~f,~f~%" 0d0 0d0 0d0)
      )

    (setf cl-mpm/penalty::*debug-force* 0d0)
    (setf cl-mpm/penalty::*debug-force-count* 0d0)

    (cl-mpm/output:save-vtk (merge-pathnames output-folder (format nil "sim_final.vtk")) *sim*)
    (format t "Substeps ~D~%" substeps)
    (time
      (loop for steps from 0 to load-steps
                while *run-sim*
                do
                (progn
                  (format t "Step ~d ~%" steps)
                  (let ((average-force 0d0)
                        (average-disp 0d0)
                        (average-reaction 0d0)
                        (conv-steps (floor (* 50 *refine*))))
                    (incf *target-displacement* disp-step)
                    (time
                     (progn
                       (setf (cl-mpm::sim-enable-damage *sim*) nil)
                       (dotimes (i conv-steps)
                         (setf cl-mpm/penalty::*debug-force* 0d0)
                         (cl-mpm:update-sim *sim*))
					   (format t "Estimated KE ~E - OOBF ~E~%" (cl-mpm/dynamic-relaxation::estimate-energy-norm *sim*)(cl-mpm/dynamic-relaxation::estimate-oobf *sim*))
                       ;(cl-mpm/dynamic-relaxation::converge-quasi-static
                       ; *sim*
                       ; :energy-crit 1d-3
                       ; :oobf-crit 1d-3
                       ; :dt-scale dt-scale
                       ; :substeps 50
                       ; :conv-steps conv-steps
                       ; :post-iter-step
                       ; (lambda (i fnorm oobf)
                       ;   (when (> i (* conv-steps 0.5))
                       ;     (setf 
                       ;       (cl-mpm:sim-damping-factor *sim*)
                       ;       0.9d0))))
                       ;(setf (cl-mpm::sim-enable-damage *sim*) t)
                       ;(cl-mpm/damage::calculate-damage *sim*)
                       ))
                    (incf average-disp (get-disp *terminus-mps*))
                    (incf average-force cl-mpm/penalty::*debug-force*)

                    (setf average-disp average-disp)
                    (setf average-force average-force)

                    (push
                     average-disp
                     *data-displacement*)
                    (push
                     average-force
                     *data-load*)

                    (let ((mesh-size (cl-mpm/mesh::mesh-resolution (cl-mpm:sim-mesh *sim*))))
                      (with-open-file (stream (merge-pathnames output-folder "disp.csv") :direction :output :if-exists :append)
                        (format stream "~f,~f,~f~%"
                                average-disp
                                (/ average-force mesh-size)
                                (/ average-reaction mesh-size))))

                    (format t "Target: ~E - Current: ~E ~%"
                            *target-displacement*
                            average-disp)
                    (format t "Load - ~E~%" average-force)
                    )
                  (incf *sim-step*))))
    (cl-mpm/output:save-vtk (merge-pathnames output-folder (format nil "sim_final.vtk")) *sim*)))

(defun run-mpi-static (&optional (output-folder "output/"))
  (ensure-directories-exist (merge-pathnames output-folder))
  (defparameter *data-force* '())
  (defparameter *data-displacement* '(0d0))
  (defparameter *data-load* '(0d0))
  (defparameter *data-node-load* '(0d0))
  (defparameter *data-mp-load* '(0d0))
  (defparameter *data-time* '(0d0))
  (defparameter *target-displacement* 0d0)
  (defparameter *data-full-time* '(0d0))
  (defparameter *data-full-load* '(0d0))

  (loop for mp across (cl-mpm:sim-mps *sim*)
        do (change-class mp 'cl-mpm/particle::particle-limestone))


  (let* ((target-time 0.2d0)
         (dt (cl-mpm:sim-dt *sim*))
         (substeps (floor target-time dt))
         (dt-scale 0.5d0)
         (load-steps 10)
         (disp-step (/ 0.8d-3 load-steps))
         (rank (cl-mpi:mpi-comm-rank)))

    (when (= rank 0)
      (format t "Outputting mesh and load-disp graph")
      (with-open-file (stream (merge-pathnames output-folder "disp.csv") :direction :output :if-exists :supersede)
        (format stream "disp,load,nload,energy,oobf~%")
        (format stream "~f,~f,~f,~f,~f~%" 0d0 0d0 0d0 0d0 0d0)))

    (setf (cl-mpm:sim-dt *sim*) (cl-mpm/setup::estimate-elastic-dt *sim* :dt-scale dt-scale))

    (setf cl-mpm/penalty::*debug-force* 0d0)
    (setf cl-mpm/penalty::*debug-force-count* 0d0)
    (cl-mpm/output:save-vtk (merge-pathnames output-folder (format nil "sim_final_~2,'0d.vtk" rank)) *sim*)
    (when (= rank 0)
      (format t "Substeps ~D~%" substeps))
    (rank-0-time
      rank
      (loop for steps from 0 to load-steps
                while *run-sim*
                do
                (progn
                  (when (= rank 0)
                    (format t "Step ~d ~%" steps))
                  (cl-mpm/output:save-vtk (merge-pathnames output-folder (format nil "sim_~2,'0d_~5,'0d.vtk" rank *sim-step*)) *sim*)
                  (cl-mpm/output::save-vtk-nodes (merge-pathnames output-folder (format nil "sim_nodes_~2,'0d_~5,'0d.vtk" rank *sim-step*)) *sim*)
                  ;(cl-mpm/output:save-vtk (merge-pathnames output-folder (format nil "sim_final_~2,'0d.vtk" rank)) *sim*) 
                  (let ((average-force 0d0)
                        (average-disp 0d0)
                        (average-reaction 0d0)
                        (conv-steps (floor (* 200 *refine*)))
                        (energy 0d0)
                        (oobf 0d0))
                    (incf *target-displacement* disp-step)
                    (rank-0-time
                      rank
                     (progn
                       (setf (cl-mpm::sim-enable-damage *sim*) nil)

                       (handler-bind ((cl-mpm/dynamic-relaxation::non-convergence-error
                                        (lambda (c)
                                          (format t "System failed to converge")
                                          (setf *run-sim* nil)
                                          (cl-mpm/output:save-vtk (merge-pathnames output-folder (format nil "sim_fail_~5,'0d.vtk" rank)) *sim*)
                                          (cl-mpm/output::save-vtk-nodes (merge-pathnames output-folder (format nil "sim_fail_nodes_~5,'0d.vtk" rank)) *sim*)
                                          (cl-mpi:mpi-barrier)
                                          (invoke-restart 'cl-mpm/dynamic-relaxation::continue)
                                          ;(error "Failed to converge")
                                          )))
                         (cl-mpm/dynamic-relaxation::converge-quasi-static
                           *sim*
                           :energy-crit 1d-3
                           :oobf-crit 1d-3
                           :dt-scale dt-scale
                           :substeps 50
                           :conv-steps conv-steps
                           :post-iter-step
                           (lambda (i fnorm oobf)
                             (setf *terminus-mps*
                                   (loop for mp across (cl-mpm:sim-mps *sim*)
                                         when (= (cl-mpm/particle::mp-index mp) 1)
                                         collect mp))
                             (let ((av-disp (cl-mpm/mpi::mpi-sum (get-disp *terminus-mps*))))
                               (when (= rank 0)
                                (format t "Target: ~E - Current: ~E ~%"
                                  *target-displacement*
                                  av-disp)))))
                         (setf (cl-mpm::sim-enable-damage *sim*) nil)
                         (cl-mpm/damage::calculate-damage *sim*)
                         )
                       ))
                    ;; (cl-mpi:mpi-waitall)
                    ;(format t "Rank ~D - Disp ~E~%" rank (get-disp *terminus-mps*))
                    ;(format t "Rank ~D - Disp ~E~%" rank cl-mpm/penalty::*debug-force*)
                    (incf average-disp (get-disp *terminus-mps*))
                    (incf average-force cl-mpm/penalty::*debug-force*)

                    (setf average-disp (cl-mpm/mpi::mpi-sum average-disp))
                    ;; (cl-mpi:mpi-waitall)
                    (setf average-force (cl-mpm/mpi::mpi-sum average-force))

                    (push
                     average-disp
                     *data-displacement*)
                    (push
                     average-force
                     *data-load*)

                    (when (= rank 0)
                      (let (;(mesh-size (cl-mpm/mesh::mesh-resolution (cl-mpm:sim-mesh *sim*)))
                            (mesh-size 1d0)
                            )
                        (with-open-file (stream (merge-pathnames output-folder "disp.csv") :direction :output :if-exists :append)
                          (format stream "~f,~f,~f,~f,~f~%"
                                  average-disp
                                  (/ average-force mesh-size)
                                  (/ average-reaction mesh-size)
                                  energy
                                  oobf
                                  ))))


                    (when (= rank 0)
                      (format t "Target: ~E - Current: ~E ~%"
                              *target-displacement*
                              average-disp)
                      (format t "Load - ~E~%" average-force)
                      ))
                  ;; (multiple-value-bind (dt-e substeps-e) (cl-mpm:calculate-adaptive-time *sim* target-time :dt-scale dt-scale)
                  ;;   (setf substeps substeps-e))
                  (incf *sim-step*)
                  (swank.live:update-swank)
                  )))
    (cl-mpm/output:save-vtk (merge-pathnames output-folder (format nil "sim_final_~2,'0d.vtk" rank)) *sim*))
  )

;(defun plot-load-disp ()
;  (let ((df (lisp-stat:read-csv
;	           (uiop:read-file-string #P"example_data/lbar/load-disp.csv"))))
;    (vgplot:xlabel "Displacement (mm)")
;    (vgplot:ylabel "Load (N)")
;    (let* ((x-model (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh *sim*)))
;           (x-experiment 0.1d0)
;           (x-scale (/ x-experiment x-model)))
;      (vgplot:plot
;       (lisp-stat:column df 'disp) (lisp-stat:column df 'load) "Data"
;       ;; (mapcar (lambda (x) (* x -1d3)) *data-displacement*) *data-node-load* "node"
;       (mapcar (lambda (x) (* x 1d3)) *data-displacement*) (mapcar (lambda (x) (* x x-scale)) *data-load*) "mpm-mps"
;       ;; (mapcar (lambda (x) (* x -1d3)) *data-displacement*) (mapcar (lambda (x) (* x -2d9)) *data-displacement*) "LE"
;       ))
;
;    ;; (vgplot:format-plot t "set xrange [~f:~f]"
;    ;;                     ;(* -1d3 (- 0.0000001d0 (reduce #'max *data-displacement*)))
;    ;;                     0d0
;    ;;                     (+ 1d-4 (* -1d3 (reduce #'min *data-displacement*)))
;    ;;                     )
;    ;; (vgplot:format-plot t "set yrange [~f:~f]"
;    ;;                     (reduce #'min (mapcar #'min *data-load* *data-mp-load*))
;    ;;                     (* 1.01 (reduce #'max (mapcar #'max *data-load* *data-mp-load*)))
;    ;;                     )
;    ;; (vgplot:axis (list nil nil nil nil))
;    )
;  )

(let ((threads (parse-integer (uiop:getenv "OMP_NUM_THREADS"))))
  (setf lparallel:*kernel* (lparallel:make-kernel threads :name "custom-kernel"))
  (format t "Thread count ~D~%" threads))
(defparameter *run-sim* nil)
;(main-loop)
(mpi-loop)

