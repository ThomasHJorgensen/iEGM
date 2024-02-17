typedef struct sol_struct
{
 double* Vw_single_to_single;
 double* Vm_single_to_single;
 double* Cw_priv_single_to_single;
 double* Cm_priv_single_to_single;
 double* Cw_pub_single_to_single;
 double* Cm_pub_single_to_single;
 double* Cw_tot_single_to_single;
 double* Cm_tot_single_to_single;
 double* EmargUw_single_to_single_pd;
 double* C_totw_single_to_single_pd;
 double* Mw_single_to_single_pd;
 double* EmargUm_single_to_single_pd;
 double* C_totm_single_to_single_pd;
 double* Mm_single_to_single_pd;
 double* Vw_couple_to_single;
 double* Vm_couple_to_single;
 double* Cw_priv_couple_to_single;
 double* Cm_priv_couple_to_single;
 double* Cw_pub_couple_to_single;
 double* Cm_pub_couple_to_single;
 double* Cw_tot_couple_to_single;
 double* Cm_tot_couple_to_single;
 double* EVw_start_as_single;
 double* EVm_start_as_single;
 double* EmargVw_start_as_single;
 double* EmargVm_start_as_single;
 double* Vw_couple_to_couple;
 double* Vm_couple_to_couple;
 double* V_couple_to_couple;
 double* Cw_priv_couple_to_couple;
 double* Cm_priv_couple_to_couple;
 double* C_pub_couple_to_couple;
 double* C_tot_couple_to_couple;
 double* Sw;
 double* Sm;
 int* power_idx;
 double* power;
 double* EmargU_pd;
 double* C_tot_pd;
 double* M_pd;
 double* V_couple_to_couple_pd;
 double* Vw_single_to_couple;
 double* Vm_single_to_couple;
 double* V_single_to_couple;
 double* Cw_priv_single_to_couple;
 double* Cm_priv_single_to_couple;
 double* C_pub_single_to_couple;
 double* Cw_tot_single_to_couple;
 double* Cm_tot_single_to_couple;
 double* initial_power;
 int* initial_power_idx;
 double* Vw_start_as_couple;
 double* Vm_start_as_couple;
 double* margV_start_as_couple;
 double* EVw_start_as_couple;
 double* EVm_start_as_couple;
 double* EmargV_start_as_couple;
 double* Cw_priv_start_as_couple;
 double* Cm_priv_start_as_couple;
 double* C_pub_start_as_couple;
 double* C_tot_start_as_couple;
 double* savings_vec;
 double* Vw_plus_vec;
 double* Vm_plus_vec;
 double* pre_Ctot_Cw_priv;
 double* pre_Ctot_Cm_priv;
 double* pre_Ctot_C_pub;
 double* solution_time;
} sol_struct;

double* get_double_p_sol_struct(sol_struct* x, char* name){

 if( strcmp(name,"Vw_single_to_single") == 0 ){ return x->Vw_single_to_single; }
 else if( strcmp(name,"Vm_single_to_single") == 0 ){ return x->Vm_single_to_single; }
 else if( strcmp(name,"Cw_priv_single_to_single") == 0 ){ return x->Cw_priv_single_to_single; }
 else if( strcmp(name,"Cm_priv_single_to_single") == 0 ){ return x->Cm_priv_single_to_single; }
 else if( strcmp(name,"Cw_pub_single_to_single") == 0 ){ return x->Cw_pub_single_to_single; }
 else if( strcmp(name,"Cm_pub_single_to_single") == 0 ){ return x->Cm_pub_single_to_single; }
 else if( strcmp(name,"Cw_tot_single_to_single") == 0 ){ return x->Cw_tot_single_to_single; }
 else if( strcmp(name,"Cm_tot_single_to_single") == 0 ){ return x->Cm_tot_single_to_single; }
 else if( strcmp(name,"EmargUw_single_to_single_pd") == 0 ){ return x->EmargUw_single_to_single_pd; }
 else if( strcmp(name,"C_totw_single_to_single_pd") == 0 ){ return x->C_totw_single_to_single_pd; }
 else if( strcmp(name,"Mw_single_to_single_pd") == 0 ){ return x->Mw_single_to_single_pd; }
 else if( strcmp(name,"EmargUm_single_to_single_pd") == 0 ){ return x->EmargUm_single_to_single_pd; }
 else if( strcmp(name,"C_totm_single_to_single_pd") == 0 ){ return x->C_totm_single_to_single_pd; }
 else if( strcmp(name,"Mm_single_to_single_pd") == 0 ){ return x->Mm_single_to_single_pd; }
 else if( strcmp(name,"Vw_couple_to_single") == 0 ){ return x->Vw_couple_to_single; }
 else if( strcmp(name,"Vm_couple_to_single") == 0 ){ return x->Vm_couple_to_single; }
 else if( strcmp(name,"Cw_priv_couple_to_single") == 0 ){ return x->Cw_priv_couple_to_single; }
 else if( strcmp(name,"Cm_priv_couple_to_single") == 0 ){ return x->Cm_priv_couple_to_single; }
 else if( strcmp(name,"Cw_pub_couple_to_single") == 0 ){ return x->Cw_pub_couple_to_single; }
 else if( strcmp(name,"Cm_pub_couple_to_single") == 0 ){ return x->Cm_pub_couple_to_single; }
 else if( strcmp(name,"Cw_tot_couple_to_single") == 0 ){ return x->Cw_tot_couple_to_single; }
 else if( strcmp(name,"Cm_tot_couple_to_single") == 0 ){ return x->Cm_tot_couple_to_single; }
 else if( strcmp(name,"EVw_start_as_single") == 0 ){ return x->EVw_start_as_single; }
 else if( strcmp(name,"EVm_start_as_single") == 0 ){ return x->EVm_start_as_single; }
 else if( strcmp(name,"EmargVw_start_as_single") == 0 ){ return x->EmargVw_start_as_single; }
 else if( strcmp(name,"EmargVm_start_as_single") == 0 ){ return x->EmargVm_start_as_single; }
 else if( strcmp(name,"Vw_couple_to_couple") == 0 ){ return x->Vw_couple_to_couple; }
 else if( strcmp(name,"Vm_couple_to_couple") == 0 ){ return x->Vm_couple_to_couple; }
 else if( strcmp(name,"V_couple_to_couple") == 0 ){ return x->V_couple_to_couple; }
 else if( strcmp(name,"Cw_priv_couple_to_couple") == 0 ){ return x->Cw_priv_couple_to_couple; }
 else if( strcmp(name,"Cm_priv_couple_to_couple") == 0 ){ return x->Cm_priv_couple_to_couple; }
 else if( strcmp(name,"C_pub_couple_to_couple") == 0 ){ return x->C_pub_couple_to_couple; }
 else if( strcmp(name,"C_tot_couple_to_couple") == 0 ){ return x->C_tot_couple_to_couple; }
 else if( strcmp(name,"Sw") == 0 ){ return x->Sw; }
 else if( strcmp(name,"Sm") == 0 ){ return x->Sm; }
 else if( strcmp(name,"power") == 0 ){ return x->power; }
 else if( strcmp(name,"EmargU_pd") == 0 ){ return x->EmargU_pd; }
 else if( strcmp(name,"C_tot_pd") == 0 ){ return x->C_tot_pd; }
 else if( strcmp(name,"M_pd") == 0 ){ return x->M_pd; }
 else if( strcmp(name,"V_couple_to_couple_pd") == 0 ){ return x->V_couple_to_couple_pd; }
 else if( strcmp(name,"Vw_single_to_couple") == 0 ){ return x->Vw_single_to_couple; }
 else if( strcmp(name,"Vm_single_to_couple") == 0 ){ return x->Vm_single_to_couple; }
 else if( strcmp(name,"V_single_to_couple") == 0 ){ return x->V_single_to_couple; }
 else if( strcmp(name,"Cw_priv_single_to_couple") == 0 ){ return x->Cw_priv_single_to_couple; }
 else if( strcmp(name,"Cm_priv_single_to_couple") == 0 ){ return x->Cm_priv_single_to_couple; }
 else if( strcmp(name,"C_pub_single_to_couple") == 0 ){ return x->C_pub_single_to_couple; }
 else if( strcmp(name,"Cw_tot_single_to_couple") == 0 ){ return x->Cw_tot_single_to_couple; }
 else if( strcmp(name,"Cm_tot_single_to_couple") == 0 ){ return x->Cm_tot_single_to_couple; }
 else if( strcmp(name,"initial_power") == 0 ){ return x->initial_power; }
 else if( strcmp(name,"Vw_start_as_couple") == 0 ){ return x->Vw_start_as_couple; }
 else if( strcmp(name,"Vm_start_as_couple") == 0 ){ return x->Vm_start_as_couple; }
 else if( strcmp(name,"margV_start_as_couple") == 0 ){ return x->margV_start_as_couple; }
 else if( strcmp(name,"EVw_start_as_couple") == 0 ){ return x->EVw_start_as_couple; }
 else if( strcmp(name,"EVm_start_as_couple") == 0 ){ return x->EVm_start_as_couple; }
 else if( strcmp(name,"EmargV_start_as_couple") == 0 ){ return x->EmargV_start_as_couple; }
 else if( strcmp(name,"Cw_priv_start_as_couple") == 0 ){ return x->Cw_priv_start_as_couple; }
 else if( strcmp(name,"Cm_priv_start_as_couple") == 0 ){ return x->Cm_priv_start_as_couple; }
 else if( strcmp(name,"C_pub_start_as_couple") == 0 ){ return x->C_pub_start_as_couple; }
 else if( strcmp(name,"C_tot_start_as_couple") == 0 ){ return x->C_tot_start_as_couple; }
 else if( strcmp(name,"savings_vec") == 0 ){ return x->savings_vec; }
 else if( strcmp(name,"Vw_plus_vec") == 0 ){ return x->Vw_plus_vec; }
 else if( strcmp(name,"Vm_plus_vec") == 0 ){ return x->Vm_plus_vec; }
 else if( strcmp(name,"pre_Ctot_Cw_priv") == 0 ){ return x->pre_Ctot_Cw_priv; }
 else if( strcmp(name,"pre_Ctot_Cm_priv") == 0 ){ return x->pre_Ctot_Cm_priv; }
 else if( strcmp(name,"pre_Ctot_C_pub") == 0 ){ return x->pre_Ctot_C_pub; }
 else if( strcmp(name,"solution_time") == 0 ){ return x->solution_time; }
 else {return NULL;}

}


int* get_int_p_sol_struct(sol_struct* x, char* name){

 if( strcmp(name,"power_idx") == 0 ){ return x->power_idx; }
 else if( strcmp(name,"initial_power_idx") == 0 ){ return x->initial_power_idx; }
 else {return NULL;}

}


