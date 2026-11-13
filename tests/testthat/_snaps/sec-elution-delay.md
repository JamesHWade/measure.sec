# current bake rows validate selected columns and exact grids

    Code
      recipes::bake(prepared, new_data = mismatch)
    Condition
      Error in `purrr::map2()`:
      i In index: 1.
      Caused by error in `.f()`:
      ! Reference and target detector locations must align exactly for each current bake row; resampling is not performed.
