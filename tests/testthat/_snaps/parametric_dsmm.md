# parametric_dsmm(); p and f are drifting.

    Code
      obj_par_model_1
    Output
      
      
      
      $dist$p_drift
      
      , , p_0
      
                Dollar $  /1'2'3/   Z E T A  O_M_E_G_A
      Dollar $       0.0       0.1       0.4       0.5
       /1'2'3/       0.5       0.0       0.3       0.2
       Z E T A       0.3       0.4       0.0       0.3
      O_M_E_G_A      0.8       0.1       0.1       0.0
      
      , , p_1
      
                Dollar $  /1'2'3/   Z E T A  O_M_E_G_A
      Dollar $       0.0       0.3       0.6       0.1
       /1'2'3/       0.3       0.0       0.4       0.3
       Z E T A       0.5       0.3       0.0       0.2
      O_M_E_G_A      0.2       0.3       0.5       0.0
      
      
      
      
      $dist$f_drift_parametric
      
      , , f_0
      
                Dollar $    /1'2'3/   Z E T A   O_M_E_G_A 
      Dollar $  NA         "unif"    "dweibull" "nbinom"  
       /1'2'3/  "geom"     NA        "pois"     "dweibull"
       Z E T A  "dweibull" "pois"    NA         "geom"    
      O_M_E_G_A "pois"     NA        "geom"     NA        
      
      , , f_1
      
                Dollar $  /1'2'3/   Z E T A  O_M_E_G_A 
      Dollar $  NA       "pois"    "geom"    "nbinom"  
       /1'2'3/  "geom"   NA        "pois"    "dweibull"
       Z E T A  "unif"   "geom"    NA        "geom"    
      O_M_E_G_A "pois"   "pois"    "geom"    NA        
      
      
      
      
      $dist$f_drift_parameters
      
      , , 1, fpars_0
      
                Dollar $  /1'2'3/   Z E T A  O_M_E_G_A
      Dollar $        NA         5       0.4       4.0
       /1'2'3/       0.7        NA       5.0       0.6
       Z E T A       0.2         3        NA       0.6
      O_M_E_G_A      4.0        NA       0.4        NA
      
      , , 2, fpars_0
      
                Dollar $  /1'2'3/   Z E T A  O_M_E_G_A
      Dollar $        NA        NA       0.2       0.6
       /1'2'3/        NA        NA        NA       0.8
       Z E T A       0.6        NA        NA        NA
      O_M_E_G_A       NA        NA        NA        NA
      
      , , 1, fpars_1
      
                Dollar $  /1'2'3/   Z E T A  O_M_E_G_A
      Dollar $        NA       6.0       0.4       3.0
       /1'2'3/       0.7        NA       2.0       0.5
       Z E T A       3.0       0.6        NA       0.7
      O_M_E_G_A      6.0       0.2       0.7        NA
      
      , , 2, fpars_1
      
                Dollar $  /1'2'3/   Z E T A  O_M_E_G_A
      Dollar $        NA        NA        NA       0.6
       /1'2'3/        NA        NA        NA       0.8
       Z E T A        NA        NA        NA        NA
      O_M_E_G_A       NA        NA        NA        NA
      
      
      $initial_dist
       Dollar $  /1'2'3/   Z E T A  O_M_E_G_A 
            0.8       0.1       0.1       0.0 
      
      $states
      Dollar $  /1'2'3/   Z E T A  O_M_E_G_A 
      
      $s
      4 
      
      $degree
      1 
      
      $model_size
      10000 
      
      $f_is_drifting
      TRUE 
      
      $p_is_drifting
      TRUE 
      
      $Model
      Model_1 
      
      
      Class: dsmm_parametric, dsmm

# parametric_dsmm(); p is drifting, f is not drifting.

    Code
      obj_par_model_2
    Output
      
      
      
      $dist$p_drift
      
      , , p_0
      
                Dollar $  /1'2'3/   Z E T A  O_M_E_G_A
      Dollar $       0.0       0.1       0.4       0.5
       /1'2'3/       0.5       0.0       0.3       0.2
       Z E T A       0.3       0.4       0.0       0.3
      O_M_E_G_A      0.8       0.1       0.1       0.0
      
      , , p_1
      
                Dollar $  /1'2'3/   Z E T A  O_M_E_G_A
      Dollar $       0.0       0.3       0.6       0.1
       /1'2'3/       0.3       0.0       0.4       0.3
       Z E T A       0.5       0.3       0.0       0.2
      O_M_E_G_A      0.2       0.3       0.5       0.0
      
      
      
      
      $dist$f_notdrift_parametric
      
                Dollar $  /1'2'3/   Z E T A   O_M_E_G_A 
      Dollar $  NA       "pois"    NA         "nbinom"  
       /1'2'3/  "geom"   NA        "geom"     "dweibull"
       Z E T A  "unif"   "geom"    NA         "geom"    
      O_M_E_G_A "nbinom" "unif"    "dweibull" NA        
      
      
      
      $dist$f_notdrift_parameters
      
      , , 1
      
                Dollar $  /1'2'3/   Z E T A  O_M_E_G_A
      Dollar $        NA       0.2        NA       3.0
       /1'2'3/       0.2        NA       0.2       0.5
       Z E T A       3.0       0.4        NA       0.7
      O_M_E_G_A      2.0       3.0       0.7        NA
      
      , , 2
      
                Dollar $  /1'2'3/   Z E T A  O_M_E_G_A
      Dollar $        NA        NA        NA       0.6
       /1'2'3/        NA        NA        NA       0.8
       Z E T A        NA        NA        NA        NA
      O_M_E_G_A      0.2        NA       0.3        NA
      
      
      $initial_dist
       Dollar $  /1'2'3/   Z E T A  O_M_E_G_A 
            0.8       0.1       0.1       0.0 
      
      $states
      Dollar $  /1'2'3/   Z E T A  O_M_E_G_A 
      
      $s
      4 
      
      $degree
      1 
      
      $model_size
      10000 
      
      $f_is_drifting
      FALSE 
      
      $p_is_drifting
      TRUE 
      
      $Model
      Model_2 
      
      
      Class: dsmm_parametric, dsmm

# parametric_dsmm(); p is not drifting, f is drifting.

    Code
      obj_par_model_3
    Output
      
      
      
      $dist$p_notdrift
      
                Dollar $  /1'2'3/   Z E T A  O_M_E_G_A
      Dollar $       0.0      0.10      0.30       0.6
       /1'2'3/       0.4      0.00      0.10       0.5
       Z E T A       0.4      0.30      0.00       0.3
      O_M_E_G_A      0.9      0.01      0.09       0.0
      
      
      
      $dist$f_drift_parametric
      
      , , f_0
      
                Dollar $    /1'2'3/   Z E T A   O_M_E_G_A 
      Dollar $  NA         "unif"    "dweibull" "nbinom"  
       /1'2'3/  "geom"     NA        "pois"     "dweibull"
       Z E T A  "dweibull" "pois"    NA         "geom"    
      O_M_E_G_A "pois"     NA        "geom"     NA        
      
      , , f_1
      
                Dollar $  /1'2'3/   Z E T A  O_M_E_G_A 
      Dollar $  NA       "pois"    "geom"    "nbinom"  
       /1'2'3/  "geom"   NA        "pois"    "dweibull"
       Z E T A  "unif"   "geom"    NA        "geom"    
      O_M_E_G_A "pois"   "pois"    "geom"    NA        
      
      
      
      
      $dist$f_drift_parameters
      
      , , 1, fpars_0
      
                Dollar $  /1'2'3/   Z E T A  O_M_E_G_A
      Dollar $        NA         5       0.4       4.0
       /1'2'3/       0.7        NA       5.0       0.6
       Z E T A       0.2         3        NA       0.6
      O_M_E_G_A      4.0        NA       0.4        NA
      
      , , 2, fpars_0
      
                Dollar $  /1'2'3/   Z E T A  O_M_E_G_A
      Dollar $        NA        NA       0.2       0.6
       /1'2'3/        NA        NA        NA       0.8
       Z E T A       0.6        NA        NA        NA
      O_M_E_G_A       NA        NA        NA        NA
      
      , , 1, fpars_1
      
                Dollar $  /1'2'3/   Z E T A  O_M_E_G_A
      Dollar $        NA       6.0       0.4       3.0
       /1'2'3/       0.7        NA       2.0       0.5
       Z E T A       3.0       0.6        NA       0.7
      O_M_E_G_A      6.0       0.2       0.7        NA
      
      , , 2, fpars_1
      
                Dollar $  /1'2'3/   Z E T A  O_M_E_G_A
      Dollar $        NA        NA        NA       0.6
       /1'2'3/        NA        NA        NA       0.8
       Z E T A        NA        NA        NA        NA
      O_M_E_G_A       NA        NA        NA        NA
      
      
      $initial_dist
       Dollar $  /1'2'3/   Z E T A  O_M_E_G_A 
            0.3       0.2       0.2       0.3 
      
      $states
      Dollar $  /1'2'3/   Z E T A  O_M_E_G_A 
      
      $s
      4 
      
      $degree
      1 
      
      $model_size
      10000 
      
      $f_is_drifting
      TRUE 
      
      $p_is_drifting
      FALSE 
      
      $Model
      Model_3 
      
      
      Class: dsmm_parametric, dsmm

