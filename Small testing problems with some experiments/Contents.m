% bidiag_gk.m   ... Golub-Kahan iterative bidiagonalization
% ==============================================================================
% Contents.m    ... this file
% ==============================================================================
% Experiments:
%
% ex_alpha_beta_estimate_rho.m  plots alphas, betas, and estimate rho for two
%                               problems (shaw, ilaplace) with and without 
%                               reorthogonalization
%                               (Figs. 11, 18)
%                               (Fig. 12 after some modificanion of code) 
% ex_distrib_function.m         plots the distribution function \omega and 
%                               the smallest Ritz values with the corresponding 
%                               weights of the matrix AA^T (i.e. eigenvalues of
%                               the matrix L_k L_k^T) growing with k
%                               (Fig 14(c) and (d))
% ex_distrib_function2.m        plots the distribution function \omega and 
%                               its approximation \omega^{(k)} 
%                               (not used in the paper)
% ex_dpc_ilaplace.m             illustrates satisfying the discrete Picard 
%                               condition for the ilaplace problem
%                               (not used in the paper)
% ex_dpc_shaw.m                 illustrates satisfying the discrete Picard 
%                               condition for the shaw problem
%                               (Fig. 1)
% ex_exact_noise_parts.m        computes the exact and the noise part of the 
%                               vectors s_j and plots and compares s_j^{exact} 
%                               and s_k^{noise} for the given index k (not
%                               used in the paper, similar to Figs. 8, 9)
% ex_heat_problem.m             plots the discrete Picard condition and the 
%                               smalest Ritz values (the noise revealing 
%                               indicator) for the problem heat with different
%                               values of parameter kappa (the problem changes
%                               from well-posed to the ill-posed)
%                               (not used in the paper) 
% ex_orthogonalization.m        animates the orthogonalization proces of left
%                               vectors from the bidiagonalization
%                               (Fig. 6)
% ex_ritz_values.m              plots ritz values of (L_k L_k^T) and 
%                               eigenvalues of AA^T 
%                               (Fig. 3)
% ex_signal_noise_space.m       plots vectors in signal and noise spaces 
%                               (Fig. 5)
% ex_sing_vecs.m                plots left singular vectors of the given problem
%                               (Fig. 2)
% ex_smalest_ritz_3pbs.m        plots the weights corresponding to the smalest 
%                               Ritz values for three problems (shaw, ilaplace, 
%                               and phillips) and compares results computed with
%                               and without reorthogonalization
%                               (not used in the paper) 
% ==============================================================================
% Small testing problems taken from [P. C. Hansen -- Regularization Tools]:
% 
% heat.m        
% ilaplace.m    
% phillips.m    
% shaw.m        
% shaw_vpa.m    ... problem shaw computed in various precision arithmetic 