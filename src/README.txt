5-4-09
multiple fixes and additions
1 numerical_grad_lap_jas had a bug and did not calculate the gradient or laplacian correctly (fixed)
2 implemented generic function passing for gradient etc.
3 decoupled loop in VPI_main for evaluationg U_0 and U_1.  when the U_1 loop was terminated due to a collision,
the U_0 array became 'corrupted' with zeros with zeros. (fixed)
4 implemented cascading swap move based on algorithm proposaed (I think) bu krauth.  
