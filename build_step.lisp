;(restrict-compiler-policy 'speed 0 0)
;(restrict-compiler-policy 'debug 3 3)
;(restrict-compiler-policy 'safety 3 3)
(restrict-compiler-policy 'speed 3 3)
(restrict-compiler-policy 'debug 0 0)
(restrict-compiler-policy 'safety 0 0)
(setf *block-compile-default* t)

(declaim (sb-ext:muffle-conditions sb-ext:compiler-note))
(pushnew :cl-mpm-pic *features*)
(ql:quickload :cl-mpm-worker)
(in-package :cl-mpm-worker)
(ql:quickload :cl-mpm)
(ql:quickload :cl-mpm/setup)
(ql:quickload :cl-mpm/particle)
(ql:quickload :cl-mpm/mpi)
(ql:quickload "magicl")
(asdf:compile-system :cl-mpm/fastmath :force t)
(asdf:compile-system :cl-mpm/ext :force t)
(asdf:compile-system :cl-mpm/forces :force t)
(asdf:compile-system :cl-mpm/utils :force t)
(asdf:compile-system :cl-mpm :force t)
(asdf:compile-system :cl-mpm/setup :force t)
(asdf:compile-system :cl-mpm/particle :force t)
(asdf:compile-system :cl-mpm/mesh :force t)
(asdf:compile-system :cl-mpm/constitutive :force t)
(asdf:compile-system :cl-mpm/damage :force t)
(asdf:compile-system :cl-mpm/bc :force t)
(asdf:compile-system :cl-mpm/output :force t)
(asdf:compile-system :cl-mpm/mpi :force t)

(ql:quickload "cl-mpm/examples/lbar")
(ql:quickload :cl-mpm/mpi)

;(ql:quickload "cl-mpm/examples/slump")
;(require 'cl-mpm-worker)

;(defun cl-mpm-worker::primary-main ()
;  (format t "Running MPI with ~D jobs~%"  (cl-mpi::mpi-comm-size))
;  (cl-mpm/examples/slump::mpi-run (cl-mpi::mpi-comm-size)))


;(build)
;(asdf:make :cl-mpm-worker)
;lfarm-server 
;(ql:quickload :lfarm-server)
;(defun main (&rest args)
;  (print "hello")
;   (lfarm-server:start-server "127.0.0.1" 11111)
;  )
(sb-ext:save-lisp-and-die
   "mpi-worker"
    :executable t
    :toplevel #'main
    :save-runtime-options t)
(uiop:quit)
