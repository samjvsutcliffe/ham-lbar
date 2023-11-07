(restrict-compiler-policy 'speed 0 0)
(restrict-compiler-policy 'debug 3 3)
(restrict-compiler-policy 'safety 3 3)
;(setf *block-compile-default* t)
(setf *features* (delete :cl-mpm-pic *features*))
(ql:quickload "magicl")
(ql:quickload "cl-mpm")
(asdf:compile-system :cl-mpm :force T)
(ql:quickload "cl-mpm/examples/lbar")
(ql:quickload :cl-mpm/mpi)
;(asdf:compile-system :cl-mpm/examples/tpb :force T)
(in-package :cl-mpm/examples/lbar)

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
          sum (magicl:norm (cl-mpm/particle:mp-velocity mp))) (length (cl-mpm:sim-mps *sim*))))

(defun get-disp (load-mps)
  ;; (* *t* *tip-velocity*)
  (if (> (length load-mps) 0d0)
    (- (/ (loop for mp in load-mps
              sum (-
                   (magicl:tref (cl-mpm/particle::mp-position mp) 1 0)
                   (* 0.5d0 (magicl:tref (cl-mpm/particle::mp-domain-size mp) 1 0))
                   )) (length load-mps))
     *initial-surface*
     ) 
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

(defun setup-test-column (size block-size offset &optional (e-scale 1) (mp-scale 1))
  (let* ((sim (cl-mpm/setup::make-block
               (/ 1d0 e-scale)
               (mapcar (lambda (x) (* x e-scale)) size)
               #'cl-mpm/shape-function:make-shape-function-bspline
               ;; 'cl-mpm::mpm-sim-usf
               'cl-mpm/damage::mpm-sim-damage
               ))
         (h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim)))
         (h-x (/ h 1d0))
         (h-y (/ h 1d0))
         (density 2.2d3)
         (elements (mapcar (lambda (s) (* e-scale (/ s 2))) size))
         )
    (declare (double-float h density))
    (progn
      (let* ((impactor-size (list 10d-3
                                  (* 0.99 h-x))))
        (setf (cl-mpm:sim-mps sim)
              (cl-mpm/setup::make-mps-from-list
               (append
                (cl-mpm/setup::make-block-mps-list
                 offset
                 block-size
                 (mapcar (lambda (e)
                           (round  (* e mp-scale) h-x)
                           ) block-size)
                 density
                 'cl-mpm/particle::particle-limestone
                 :E 25.85d9
                 :nu 0.18d0
                 ;; :elastic-approxmation :plane-stress
                 :fracture-energy 95d0
                 :initiation-stress (* 2.7d6 1d0)
                 :critical-damage 1.000d0
                 :internal-length 25d-3
                 :local-length 25d-3
                 :local-length-damaged 25d-3
                 :compression-ratio 10d0
                 ;; :local-length-damaged 0.01d0
                 :gravity -0.0d0
                 :gravity-axis (cl-mpm/utils:vector-from-list '(0d0 1d0 0d0))
                 )
                ;; impactors
                )
               )))
      (setf (cl-mpm:sim-allow-mp-split sim) nil)
      (setf (cl-mpm::sim-enable-damage sim) t)
      (setf (cl-mpm::sim-nonlocal-damage sim) t)
      (setf (cl-mpm::sim-allow-mp-damage-removal sim) nil)
      (setf (cl-mpm::sim-mp-damage-removal-instant sim) nil)
      (setf (cl-mpm::sim-mass-filter sim) 0d0)
      (let ((ms 1d6))
        (setf (cl-mpm::sim-mass-scale sim) ms)
        (setf (cl-mpm:sim-damping-factor sim)
              (* 1d-3 density ms)
              ;; 1d0
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

      (let ((dt-scale 1d0))
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

      (let* ((crack-pos
               (loop for mp across (cl-mpm:sim-mps sim)
                     when
                     (not (= (cl-mpm/particle::mp-index mp) 1))
                     maximize (magicl:tref (cl-mpm/particle:mp-position mp) 0 0)))
             (above-crack
               (loop for mp across (cl-mpm:sim-mps sim)
                     when
                     (and
                      (not (= (cl-mpm/particle::mp-index mp) 1))
                      (= (magicl:tref (cl-mpm/particle:mp-position mp) 0 0) crack-pos))
                     collect mp)
               )
             (min-pos (loop for mp in above-crack
                            minimize (magicl:tref(cl-mpm/particle:mp-position mp) 1 0)))
             )
        (defparameter *terminus-mps*
          (loop for mp in above-crack
                when (= min-pos (magicl:tref
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
                 (cl-mpm/bc::make-bc-fixed left-node-pos
                                           '(0 0 nil))

                 (cl-mpm/bc::make-bc-fixed right-node-pos
                                           '(nil 0 nil)))
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
                '(0 0 0)
                (lambda ()
                  (with-accessors ((mesh cl-mpm:sim-mesh)
                                   (dt cl-mpm::sim-dt))
                      sim
                    (let ((datum (* 1d0 (+ *initial-surface* *target-displacement*)))
                          (normal (cl-mpm/utils:vector-from-list  '(0d0 1d0 0d0))))
                      (cl-mpm/penalty::apply-displacement-control-mps mesh (coerce *terminus-mps* 'vector)
                                                       dt
                                                       normal
                                                       datum
                                                       (* density 1d5)
                                                       0d0)
                      )
                    )))
               ))))

      ;; (let* ((terminus-size (second block-size))
      ;;        (ocean-y (* terminus-size 1.0d0)))
      ;;   (setf (cl-mpm::sim-bcs-force-list sim)
      ;;         (list
      ;;          (cl-mpm/bc:make-bcs-from-list
      ;;           (list
      ;;            (cl-mpm/buoyancy::make-bc-buoyancy-clip
      ;;             sim
      ;;             ocean-y
      ;;             1000d0
      ;;             (lambda (pos datum)
      ;;               t)
      ;;             ))))))

      sim))
  )
(defparameter *sim* nil)
(defparameter *run-sim* t)
(defparameter *t* 0)
(defparameter *sim-step* 0)
(defparameter *refine* (/ 1d0 2d0))
(let ((refine (uiop:getenv "REFINE")))
  (when refine
    (setf *refine* (parse-integer (uiop:getenv "REFINE")))
    ))
(defun setup (&key (undercut 0d0))
  ;; (let ((mps-per-dim 4))
  ;;   (defparameter *sim* (setup-test-column '(16 16) '(8 8)  '(0 0) *refine* mps-per-dim)))
  ;; (defparameter *sim* (setup-test-column '(1 1 1) '(1 1 1) 1 1))

  (let* ((mesh-size (/ 0.025 2.0d0))
         (mps-per-cell 2)
         (shelf-height 0.500d0)
         (shelf-length 0.500d0)
         ;; (shelf-length 0.225d0)
         (domain-length (+ shelf-length (* 8 mesh-size)))
         (offset (list
                  (* 2 mesh-size)
                  0d0
                  0d0))


         )
    (defparameter *sim*
      (setup-test-column (list domain-length
                               domain-length
                               (* 2 mesh-size))
                         (list shelf-length shelf-height
                               mesh-size)
                         offset
                         (/ 1d0 mesh-size) mps-per-cell))
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
    (format t "Total weight ~F~%"
            (loop for mp across (cl-mpm:sim-mps *sim*)
                  sum (* 9.8d0 (cl-mpm/particle:mp-mass mp))))

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
    (format stream "disp,load,nload~%")
    (format stream "~f,~f~%" 0d0 0d0 0d0)
    )

  (let* ((target-time 0.5d0)
         (dt (cl-mpm:sim-dt *sim*))
         (substeps (floor target-time dt))
         (dt-scale 1d0)
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
                       (with-open-file (stream (merge-pathnames "output/disp.csv") :direction :output :if-exists :append)
                         (format stream "~f,~f~%"
                                 average-disp
                                 (/ average-force mesh-size)
                                 (/ average-reaction mesh-size)
                                 )))


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
(defun mpi-loop ()
  (let ((rank (cl-mpi:mpi-comm-rank)))
    ;; (cl-mpi:mpi-init)
    (setup)

    (setf (cl-mpm/mpi::mpm-sim-mpi-domain-count *sim*) (list (cl-mpi:mpi-comm-size) 1 1))
    (when (= rank 0)
      (format t "Sim MPs: ~a~%" (length (cl-mpm:sim-mps *sim*)))
      (format t "Decompose~%"))
    (cl-mpm/mpi::domain-decompose *sim*)
    (when (= rank 0)
      (format t "Sim MPs: ~a~%" (length (cl-mpm:sim-mps *sim*))))
    (when (= rank 0)
      (format t "Run mpi~%"))
    ;; (time (cl-mpm/mpi::mpi-sync-momentum *sim*))
    ;(time (cl-mpm::update-sim *sim*))
    (if (= rank 0)
        (time (cl-mpm::update-sim *sim*))
        (cl-mpm::update-sim *sim*))
    (run-mpi)
    (when (= rank 0)
      (format t "Done mpi~%"))
    )
  (cl-mpi:mpi-finalize)
  (sb-ext:quit)
  )
(defun mpi-average (value mp-count)
  (let ((sum 0d0))
    (static-vectors:with-static-vector (source 1 :element-type 'double-float :initial-element value)
      (static-vectors:with-static-vector (dest 1 :element-type 'double-float :initial-element 0d0)
        (cl-mpi:mpi-allreduce source dest cl-mpi:+mpi-sum+)
        (setf sum (aref dest 0))
        ))
    (static-vectors:with-static-vector (source 1 :element-type 'double-float :initial-element (coerce mp-count 'double-float))
      (static-vectors:with-static-vector (dest 1 :element-type 'double-float :initial-element 0d0)
        (cl-mpi:mpi-allreduce source dest cl-mpi:+mpi-sum+)
        (setf sum (/ sum (min 1d0 (aref dest 0))))
        ))

    sum))
(defun run-mpi ()
  (defparameter *data-force* '())
  (defparameter *data-displacement* '(0d0))
  (defparameter *data-load* '(0d0))
  (defparameter *data-node-load* '(0d0))
  (defparameter *data-mp-load* '(0d0))
  (defparameter *data-time* '(0d0))
  (defparameter *target-displacement* 0d0)
  (defparameter *data-full-time* '(0d0))
  (defparameter *data-full-load* '(0d0))


  (let* ((target-time 2d0)
         (dt (cl-mpm:sim-dt *sim*))
         (substeps (floor target-time dt))
         (dt-scale 1d0)
         (disp-step -0.002d-3)
         (rank (cl-mpi:mpi-comm-rank))
         )

    (when (= rank 0)
      (format t "Outputting mesh and load-disp graph")
      (cl-mpm/output:save-vtk-mesh (merge-pathnames "output/mesh.vtk") *sim*)
      (with-open-file (stream (merge-pathnames "output/disp.csv") :direction :output :if-exists :supersede)
        (format stream "disp,load~%")
        (format stream "~f,~f~%" 0d0 0d0)
        ))

    (setf cl-mpm/penalty::*debug-force* 0d0)
    (setf cl-mpm/penalty::*debug-force-count* 0d0)
    (when (= rank 0)
      (format t "Test run~%"))
    (cl-mpm::update-sim *sim*)
    (when (= rank 0)
      (format t "Caluclating dt estimate~%"))
    (multiple-value-bind (dt-e substeps-e) (cl-mpm:calculate-adaptive-time *sim* target-time :dt-scale dt-scale)
      (when (= rank 0)
        (format t "CFL dt estimate: ~f~%" dt-e)
        (format t "CFL step count estimate: ~D~%" substeps-e))
      (setf substeps substeps-e))
    (when (= rank 0)
      (format t "Substeps ~D~%" substeps))
    (incf *target-displacement* disp-step)
    (time (loop for steps from 0 to 100
                while *run-sim*
                do
                   (progn
                     (when (= rank 0)
                       (format t "Step ~d ~%" steps))
                     (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~2,'0d_~5,'0d.vtk" rank *sim-step*)) *sim*)
                     (let ((average-force 0d0)
                           (average-disp 0d0)
                           (average-reaction 0d0))
                       (time
                        (dotimes (i substeps);)
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
                          (incf average-reaction (/ (get-reaction-force *fixed-nodes*) substeps))
                          (setf cl-mpm/penalty::*debug-force* 0d0)
                          (setf cl-mpm/penalty::*debug-force-count* 0d0)
                          (cl-mpm::update-sim *sim*)
                          (incf *target-displacement* (/ disp-step substeps))
                          (setf *t* (+ *t* (cl-mpm::sim-dt *sim*))))
                        )
                       ;; (incf *target-displacement* -0.01d-3)
                       (setf average-disp (mpi-average average-disp (length *terminus-mps*)))
                       (setf average-force (mpi-average average-force (length *terminus-mps*)))
                       (push
                        average-disp
                        *data-displacement*)
                       (push
                        average-force
                        *data-load*)

                       (when (= rank 0)
                         (with-open-file (stream (merge-pathnames "output/disp.csv") :direction :output :if-exists :append)
                           (format stream "~f,~f~%"
                                   average-disp
                                   average-force
                                   ))))


                     (when (= rank 0)
                       (format t "Target: ~f - Current: ~f Error: ~f - energy ~F~%"
                               *target-displacement*
                               (get-disp *terminus-mps*)
                               (* 100d0 (/ (- *target-displacement* (get-disp *terminus-mps*)) *target-displacement*))
                               (energy-norm *sim*)))
                     (multiple-value-bind (dt-e substeps-e) (cl-mpm:calculate-adaptive-time *sim* target-time :dt-scale dt-scale)
                       ;; (format t "CFL dt estimate: ~f~%" dt-e)
                       ;; (format t "CFL step count estimate: ~D~%" substeps-e)
                       (setf substeps substeps-e))
                     (incf *sim-step*)
                     (swank.live:update-swank)
                     ))))
  ;(cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" *sim-step*)) *sim*)
  )


(setf lparallel:*kernel* (lparallel:make-kernel 32 :name "custom-kernel"))
(defparameter *run-sim* nil)
(setup)
(format t "MP count:~D~%" (length (cl-mpm:sim-mps *sim*)))
(run)

;(mpi-loop)
