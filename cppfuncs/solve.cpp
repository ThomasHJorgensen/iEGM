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

EXPORT void check_participation(int t, int iL, int iP, int iA, sol_struct* sol, par_struct* par){

    // 1. Setup
    /// a. lists
    int num = 5;
    double** list_start_as_couple = new double*[num]; 
    double** list_couple_to_couple = new double*[num];
    double* list_couple_to_single = new double[num];             

    // b. temporary arrays
    double* Sw = new double[par->num_power];
    double* Sm = new double[par->num_power];

    // c. index struct to pass to bargaining algorithm
    index::index_couple_struct* idx_couple = new index::index_couple_struct;

    // i. Get indices
    int idx_single = index::single(t,iA,par);
    idx_couple->t = t;
    idx_couple->iL = iL;
    idx_couple->iA = iA;
    idx_couple->par = par;

    // ii Calculate marital surplus
    for (int iP=0; iP<par->num_power; iP++){
        int idx_tmp = index::couple(t,iP,iL,iA,par);
        Sw[iP] = couple::calc_marital_surplus(sol->Vw_couple_to_couple[idx_tmp],sol->Vw_couple_to_single[idx_single],par);
        Sm[iP] = couple::calc_marital_surplus(sol->Vm_couple_to_couple[idx_tmp],sol->Vm_couple_to_single[idx_single],par);
    }

    // iii. setup relevant lists 
    int i = 0;
    list_start_as_couple[i] = sol->Vw_start_as_couple; i++;
    list_start_as_couple[i] = sol->Vm_start_as_couple; i++;
    list_start_as_couple[i] = sol->Cw_priv_start_as_couple; i++;
    list_start_as_couple[i] = sol->Cm_priv_start_as_couple; i++;
    list_start_as_couple[i] = sol->C_pub_start_as_couple; i++; //consider having two of these, one for each spouse
    i = 0;
    list_couple_to_couple[i] = sol->Vw_couple_to_couple; i++;
    list_couple_to_couple[i] = sol->Vm_couple_to_couple; i++;
    list_couple_to_couple[i] = sol->Cw_priv_couple_to_couple; i++;
    list_couple_to_couple[i] = sol->Cm_priv_couple_to_couple; i++;
    list_couple_to_couple[i] = sol->C_pub_couple_to_couple; i++; //consider having two of these, one for each spouse
    i = 0;
    list_couple_to_single[i] = sol->Vw_couple_to_single[idx_single]; i++;
    list_couple_to_single[i] = sol->Vm_couple_to_single[idx_single]; i++;
    list_couple_to_single[i] = sol->Cw_priv_couple_to_single[idx_single]; i++;
    list_couple_to_single[i] = sol->Cm_priv_couple_to_single[idx_single]; i++;
    list_couple_to_single[i] = sol->Cw_pub_couple_to_single[idx_single]; i++; //consider having two of these, one for each spouse

    logs::write("barg_log.txt",0, "Checking participation constraints for couple at time %d with love shock %d and power shock %d and asset level %d \n",t,iP,iL,iA);

    // 2. Check participation constraints
    bargaining::check_participation_constraints_verbose(sol->power_idx, sol->power, Sw, Sm, idx_couple, list_start_as_couple, list_couple_to_couple, list_couple_to_single, num, par, true);



}