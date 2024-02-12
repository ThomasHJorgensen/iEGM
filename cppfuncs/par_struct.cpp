typedef struct par_struct
{
 double R;
 double beta;
 double div_A_share;
 double div_cost;
 double inc_w;
 double inc_m;
 double rho_w;
 double rho_m;
 double alpha1_w;
 double alpha1_m;
 double alpha2_w;
 double alpha2_m;
 double phi_w;
 double phi_m;
 int T;
 int num_A;
 double max_A;
 int num_power;
 int num_love;
 double max_love;
 double sigma_love;
 int num_shock_love;
 double p_meet;
 double* prob_partner_A_w;
 double* prob_partner_A_m;
 bool interp_inverse;
 char* interp_method;
 bool interp_intra_period_in_num_inverse;
 int num_Ctot;
 double max_Ctot;
 bool do_egm;
 bool analytic_marg_u_single;
 bool analytic_inv_marg_u_single;
 int num_A_pd;
 double max_A_pd;
 int num_marg_u;
 int seed;
 int simT;
 int simN;
 bool use_external_solution;
 int num_A_true;
 double max_A_true;
 int num_power_true;
 int num_love_true;
 double max_love_true;
 bool do_cpp;
 int threads;
 bool centered_gradient;
 bool interp_power;
 double* grid_A;
 double* grid_Aw;
 double* grid_Am;
 double* grid_power;
 double* grid_power_flip;
 double* grid_love;
 double* grid_shock_love;
 double* grid_weight_love;
 double* grid_Ctot;
 double* grid_marg_u;
 double* grid_marg_u_for_inv;
 double* grid_C_for_marg_u;
 double* grid_inv_marg_u;
 double* grid_marg_u_single_w;
 double* grid_marg_u_single_w_for_inv;
 double* grid_marg_u_single_m;
 double* grid_marg_u_single_m_for_inv;
 double* grid_A_pd;
 double* grid_Aw_pd;
 double* grid_Am_pd;
 double* prob_repartner;
 double* prob_partner_love;
 double* cdf_partner_Aw;
 double* cdf_partner_Am;
 double* grid_A_true;
 double* grid_Aw_true;
 double* grid_Am_true;
 double* grid_power_true;
 double* grid_love_true;
} par_struct;

double get_double_par_struct(par_struct* x, char* name){

 if( strcmp(name,"R") == 0 ){ return x->R; }
 else if( strcmp(name,"beta") == 0 ){ return x->beta; }
 else if( strcmp(name,"div_A_share") == 0 ){ return x->div_A_share; }
 else if( strcmp(name,"div_cost") == 0 ){ return x->div_cost; }
 else if( strcmp(name,"inc_w") == 0 ){ return x->inc_w; }
 else if( strcmp(name,"inc_m") == 0 ){ return x->inc_m; }
 else if( strcmp(name,"rho_w") == 0 ){ return x->rho_w; }
 else if( strcmp(name,"rho_m") == 0 ){ return x->rho_m; }
 else if( strcmp(name,"alpha1_w") == 0 ){ return x->alpha1_w; }
 else if( strcmp(name,"alpha1_m") == 0 ){ return x->alpha1_m; }
 else if( strcmp(name,"alpha2_w") == 0 ){ return x->alpha2_w; }
 else if( strcmp(name,"alpha2_m") == 0 ){ return x->alpha2_m; }
 else if( strcmp(name,"phi_w") == 0 ){ return x->phi_w; }
 else if( strcmp(name,"phi_m") == 0 ){ return x->phi_m; }
 else if( strcmp(name,"max_A") == 0 ){ return x->max_A; }
 else if( strcmp(name,"max_love") == 0 ){ return x->max_love; }
 else if( strcmp(name,"sigma_love") == 0 ){ return x->sigma_love; }
 else if( strcmp(name,"p_meet") == 0 ){ return x->p_meet; }
 else if( strcmp(name,"max_Ctot") == 0 ){ return x->max_Ctot; }
 else if( strcmp(name,"max_A_pd") == 0 ){ return x->max_A_pd; }
 else if( strcmp(name,"max_A_true") == 0 ){ return x->max_A_true; }
 else if( strcmp(name,"max_love_true") == 0 ){ return x->max_love_true; }
 else {return NAN;}

}


int get_int_par_struct(par_struct* x, char* name){

 if( strcmp(name,"T") == 0 ){ return x->T; }
 else if( strcmp(name,"num_A") == 0 ){ return x->num_A; }
 else if( strcmp(name,"num_power") == 0 ){ return x->num_power; }
 else if( strcmp(name,"num_love") == 0 ){ return x->num_love; }
 else if( strcmp(name,"num_shock_love") == 0 ){ return x->num_shock_love; }
 else if( strcmp(name,"num_Ctot") == 0 ){ return x->num_Ctot; }
 else if( strcmp(name,"num_A_pd") == 0 ){ return x->num_A_pd; }
 else if( strcmp(name,"num_marg_u") == 0 ){ return x->num_marg_u; }
 else if( strcmp(name,"seed") == 0 ){ return x->seed; }
 else if( strcmp(name,"simT") == 0 ){ return x->simT; }
 else if( strcmp(name,"simN") == 0 ){ return x->simN; }
 else if( strcmp(name,"num_A_true") == 0 ){ return x->num_A_true; }
 else if( strcmp(name,"num_power_true") == 0 ){ return x->num_power_true; }
 else if( strcmp(name,"num_love_true") == 0 ){ return x->num_love_true; }
 else if( strcmp(name,"threads") == 0 ){ return x->threads; }
 else {return -9999;}

}


double* get_double_p_par_struct(par_struct* x, char* name){

 if( strcmp(name,"prob_partner_A_w") == 0 ){ return x->prob_partner_A_w; }
 else if( strcmp(name,"prob_partner_A_m") == 0 ){ return x->prob_partner_A_m; }
 else if( strcmp(name,"grid_A") == 0 ){ return x->grid_A; }
 else if( strcmp(name,"grid_Aw") == 0 ){ return x->grid_Aw; }
 else if( strcmp(name,"grid_Am") == 0 ){ return x->grid_Am; }
 else if( strcmp(name,"grid_power") == 0 ){ return x->grid_power; }
 else if( strcmp(name,"grid_power_flip") == 0 ){ return x->grid_power_flip; }
 else if( strcmp(name,"grid_love") == 0 ){ return x->grid_love; }
 else if( strcmp(name,"grid_shock_love") == 0 ){ return x->grid_shock_love; }
 else if( strcmp(name,"grid_weight_love") == 0 ){ return x->grid_weight_love; }
 else if( strcmp(name,"grid_Ctot") == 0 ){ return x->grid_Ctot; }
 else if( strcmp(name,"grid_marg_u") == 0 ){ return x->grid_marg_u; }
 else if( strcmp(name,"grid_marg_u_for_inv") == 0 ){ return x->grid_marg_u_for_inv; }
 else if( strcmp(name,"grid_C_for_marg_u") == 0 ){ return x->grid_C_for_marg_u; }
 else if( strcmp(name,"grid_inv_marg_u") == 0 ){ return x->grid_inv_marg_u; }
 else if( strcmp(name,"grid_marg_u_single_w") == 0 ){ return x->grid_marg_u_single_w; }
 else if( strcmp(name,"grid_marg_u_single_w_for_inv") == 0 ){ return x->grid_marg_u_single_w_for_inv; }
 else if( strcmp(name,"grid_marg_u_single_m") == 0 ){ return x->grid_marg_u_single_m; }
 else if( strcmp(name,"grid_marg_u_single_m_for_inv") == 0 ){ return x->grid_marg_u_single_m_for_inv; }
 else if( strcmp(name,"grid_A_pd") == 0 ){ return x->grid_A_pd; }
 else if( strcmp(name,"grid_Aw_pd") == 0 ){ return x->grid_Aw_pd; }
 else if( strcmp(name,"grid_Am_pd") == 0 ){ return x->grid_Am_pd; }
 else if( strcmp(name,"prob_repartner") == 0 ){ return x->prob_repartner; }
 else if( strcmp(name,"prob_partner_love") == 0 ){ return x->prob_partner_love; }
 else if( strcmp(name,"cdf_partner_Aw") == 0 ){ return x->cdf_partner_Aw; }
 else if( strcmp(name,"cdf_partner_Am") == 0 ){ return x->cdf_partner_Am; }
 else if( strcmp(name,"grid_A_true") == 0 ){ return x->grid_A_true; }
 else if( strcmp(name,"grid_Aw_true") == 0 ){ return x->grid_Aw_true; }
 else if( strcmp(name,"grid_Am_true") == 0 ){ return x->grid_Am_true; }
 else if( strcmp(name,"grid_power_true") == 0 ){ return x->grid_power_true; }
 else if( strcmp(name,"grid_love_true") == 0 ){ return x->grid_love_true; }
 else {return NULL;}

}


bool get_bool_par_struct(par_struct* x, char* name){

 if( strcmp(name,"interp_inverse") == 0 ){ return x->interp_inverse; }
 else if( strcmp(name,"interp_intra_period_in_num_inverse") == 0 ){ return x->interp_intra_period_in_num_inverse; }
 else if( strcmp(name,"do_egm") == 0 ){ return x->do_egm; }
 else if( strcmp(name,"analytic_marg_u_single") == 0 ){ return x->analytic_marg_u_single; }
 else if( strcmp(name,"analytic_inv_marg_u_single") == 0 ){ return x->analytic_inv_marg_u_single; }
 else if( strcmp(name,"use_external_solution") == 0 ){ return x->use_external_solution; }
 else if( strcmp(name,"do_cpp") == 0 ){ return x->do_cpp; }
 else if( strcmp(name,"centered_gradient") == 0 ){ return x->centered_gradient; }
 else if( strcmp(name,"interp_power") == 0 ){ return x->interp_power; }
 else {return false;}

}


char* get_char_p_par_struct(par_struct* x, char* name){

 if( strcmp(name,"interp_method") == 0 ){ return x->interp_method; }
 else {return NULL;}

}


