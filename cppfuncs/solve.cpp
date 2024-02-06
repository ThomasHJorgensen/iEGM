#define MAIN
#include "myheader.h"

// include these again here to ensure that they are automatically compiled by consav
#ifndef MAIN
#include "precompute.cpp"
#endif

/////////////
// 5. MAIN //
/////////////

EXPORT void solve(sol_struct *sol, par_struct *par){
    
    // pre-compute intra-temporal optimal allocation
    precompute::precompute(sol,par);

    // loop backwards
    for (int t = par->T-1; t >= 0; t--){
        single::solve_single_to_single(t,sol,par); 
        single::solve_couple_to_single(t,sol,par); 
        couple::solve_couple(t,sol,par);
        couple::solve_single_to_couple(t,sol,par);
        single::expected_value_start_single(t,sol,par);
    }
}


EXPORT void simulate(sim_struct *sim, sol_struct *sol, par_struct *par){
    
    sim::model(sim,sol,par);

}


EXPORT void compute_margEV(sol_struct* sol, par_struct* par){
    for (int t = 0; t < par->T; t++){
        single::calc_marginal_value_single(t, woman, sol, par);
        single::calc_marginal_value_single(t, man, sol, par);

        for (int iP=0; iP<par->num_power; iP++){
            for (int iL=0; iL<par->num_love; iL++){
                int idx = index::couple(t,iP,iL,0,par);
                double* Vw = &sol->Vw_start_as_couple[idx];
                double* Vm = &sol->Vm_start_as_couple[idx];
                double* margV = &sol->margV_start_as_couple[idx];
                couple::calc_marginal_value_couple(t, iP, iL, Vw, Vm, margV, sol, par);
            }
        }
    }
}