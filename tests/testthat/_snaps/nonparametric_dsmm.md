# nonparametric_dsmm(); p and f are drifting

    Code
      obj_nonpar_model_1
    Output
      
      
      
      $dist$p_drift
      
      , , p_0
      
          AA  AC  CC
      AA 0.0 0.1 0.9
      AC 0.5 0.0 0.5
      CC 0.3 0.7 0.0
      
      , , p_(1/2)
      
          AA  AC  CC
      AA 0.0 0.6 0.4
      AC 0.7 0.0 0.3
      CC 0.6 0.4 0.0
      
      , , p_1
      
          AA  AC  CC
      AA 0.0 0.2 0.8
      AC 0.6 0.0 0.4
      CC 0.7 0.3 0.0
      
      
      
      
      $dist$f_drift
      
      , , l = 1, f_0
      
          AA  AC  CC
      AA 0.0 0.2 0.7
      AC 0.3 0.0 0.4
      CC 0.2 0.8 0.0
      
      , , l = 2, f_0
      
          AA   AC  CC
      AA 0.0 0.30 0.2
      AC 0.2 0.00 0.5
      CC 0.1 0.15 0.0
      
      , , l = 3, f_0
      
          AA   AC  CC
      AA 0.0 0.50 0.1
      AC 0.5 0.00 0.1
      CC 0.7 0.05 0.0
      
      , , l = 1, f_(1/2)
      
          AA        AC  CC
      AA 0.0 0.3333333 0.4
      AC 0.3 0.0000000 0.4
      CC 0.2 0.1000000 0.0
      
      , , l = 2, f_(1/2)
      
          AA        AC  CC
      AA 0.0 0.3333333 0.4
      AC 0.4 0.0000000 0.2
      CC 0.3 0.4000000 0.0
      
      , , l = 3, f_(1/2)
      
          AA        AC  CC
      AA 0.0 0.3333333 0.2
      AC 0.3 0.0000000 0.4
      CC 0.5 0.5000000 0.0
      
      , , l = 1, f_1
      
           AA  AC  CC
      AA 0.00 0.3 0.3
      AC 0.30 0.0 0.5
      CC 0.05 0.1 0.0
      
      , , l = 2, f_1
      
          AA  AC   CC
      AA 0.0 0.2 0.60
      AC 0.3 0.0 0.35
      CC 0.9 0.2 0.00
      
      , , l = 3, f_1
      
           AA  AC   CC
      AA 0.00 0.5 0.10
      AC 0.40 0.0 0.15
      CC 0.05 0.7 0.00
      
      
      $initial_dist
       AA  AC  CC 
      0.3 0.5 0.2 
      
      $states
      AA AC CC 
      
      $s
      3 
      
      $degree
      2 
      
      $k_max
      3 
      
      $model_size
      8000 
      
      $f_is_drifting
      TRUE 
      
      $p_is_drifting
      TRUE 
      
      $Model
      Model_1 
      
      
      Class: dsmm_nonparametric, dsmm

# nonparametric_dsmm(); p is drifting, f is not drifting.

    Code
      obj_nonpar_model_2
    Output
      
      
      
      $dist$p_drift
      
      , , p_0
      
          AA  AC  CC
      AA 0.0 0.1 0.9
      AC 0.5 0.0 0.5
      CC 0.3 0.7 0.0
      
      , , p_(1/2)
      
          AA  AC  CC
      AA 0.0 0.6 0.4
      AC 0.7 0.0 0.3
      CC 0.6 0.4 0.0
      
      , , p_1
      
          AA  AC  CC
      AA 0.0 0.2 0.8
      AC 0.6 0.0 0.4
      CC 0.7 0.3 0.0
      
      
      
      
      $dist$f_notdrift
      
      , , l = 1
      
          AA        AC  CC
      AA 0.0 0.3333333 0.4
      AC 0.3 0.0000000 0.4
      CC 0.2 0.1000000 0.0
      
      , , l = 2
      
          AA        AC  CC
      AA 0.0 0.3333333 0.4
      AC 0.4 0.0000000 0.2
      CC 0.3 0.4000000 0.0
      
      , , l = 3
      
          AA        AC  CC
      AA 0.0 0.3333333 0.2
      AC 0.3 0.0000000 0.4
      CC 0.5 0.5000000 0.0
      
      
      $initial_dist
       AA  AC  CC 
      0.7 0.1 0.2 
      
      $states
      AA AC CC 
      
      $s
      3 
      
      $degree
      2 
      
      $k_max
      3 
      
      $model_size
      10000 
      
      $f_is_drifting
      FALSE 
      
      $p_is_drifting
      TRUE 
      
      $Model
      Model_2 
      
      
      Class: dsmm_nonparametric, dsmm

# nonparametric_dsmm(); f is drifting, p is not drifting.

    Code
      obj_nonpar_model_3
    Output
      
      
      
      $dist$p_notdrift
      
          AA  AC  CC
      AA 0.0 0.2 0.8
      AC 0.6 0.0 0.4
      CC 0.7 0.3 0.0
      
      
      
      $dist$f_drift
      
      , , l = 1, f_0
      
          AA  AC  CC
      AA 0.0 0.2 0.7
      AC 0.3 0.0 0.4
      CC 0.2 0.8 0.0
      
      , , l = 2, f_0
      
          AA   AC  CC
      AA 0.0 0.30 0.2
      AC 0.2 0.00 0.5
      CC 0.1 0.15 0.0
      
      , , l = 3, f_0
      
          AA   AC  CC
      AA 0.0 0.50 0.1
      AC 0.5 0.00 0.1
      CC 0.7 0.05 0.0
      
      , , l = 1, f_(1/2)
      
          AA        AC  CC
      AA 0.0 0.3333333 0.4
      AC 0.3 0.0000000 0.4
      CC 0.2 0.1000000 0.0
      
      , , l = 2, f_(1/2)
      
          AA        AC  CC
      AA 0.0 0.3333333 0.4
      AC 0.4 0.0000000 0.2
      CC 0.3 0.4000000 0.0
      
      , , l = 3, f_(1/2)
      
          AA        AC  CC
      AA 0.0 0.3333333 0.2
      AC 0.3 0.0000000 0.4
      CC 0.5 0.5000000 0.0
      
      , , l = 1, f_1
      
           AA  AC  CC
      AA 0.00 0.3 0.3
      AC 0.30 0.0 0.5
      CC 0.05 0.1 0.0
      
      , , l = 2, f_1
      
          AA  AC   CC
      AA 0.0 0.2 0.60
      AC 0.3 0.0 0.35
      CC 0.9 0.2 0.00
      
      , , l = 3, f_1
      
           AA  AC   CC
      AA 0.00 0.5 0.10
      AC 0.40 0.0 0.15
      CC 0.05 0.7 0.00
      
      
      $initial_dist
       AA  AC  CC 
      0.3 0.4 0.3 
      
      $states
      AA AC CC 
      
      $s
      3 
      
      $degree
      2 
      
      $k_max
      3 
      
      $model_size
      10000 
      
      $f_is_drifting
      TRUE 
      
      $p_is_drifting
      FALSE 
      
      $Model
      Model_3 
      
      
      Class: dsmm_nonparametric, dsmm

