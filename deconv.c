// Stuff that's weird:
//
// some unitialized variables 
// matrix mlriplication is wrong?
// gcc seems to dislike control_a?
//

#include "headers.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_fit.h>
#include <time.h>
#include <math.h>
#include <unistd.h>
//#define PI (3.1415)
//#define k_b (1.3806488*pow(10, -23))
#define ALPHA_GUARD (1.0)	/* In principle the value here should not matter. Zero is unfortunate for numerical reasons. */
#define P_GUARD (pow(10,-2))

// Alpha guard is to prevent division by zero
// Starting a new function file for all those damned different convolution things.
// Maximum entropy routines are supposed to go here.
void get_signal_with_var(double * signal, int len, double z_start, double z_final,
			 int traces, int samples, double * lambda_ar, double K,
			 double ** x_ar, int c, double * bin_centers, double beta, double * var){

  int i,k,l;
  double s;
  double * bin_edges = (double *) malloc (sizeof(double) * len+1);
  //double * bin_centers = (double *) malloc (sizeof(double) * len);
  double * signal_counter = (double *) malloc (sizeof(double) * len);
  double ** ext_work = (double **) malloc ( sizeof(double * )* traces );
  double * exp_avarage = (double *) malloc (sizeof(double) * samples);
  double * exp2_avarage = (double *) malloc (sizeof(double) * samples);
  double * var_temp = (double *) malloc (sizeof(double) * len);
  double d_bin = (z_final - z_start)/len;
  // bins
  for(i=0;i<len+1;i++){
    bin_edges[i] = z_start + i*d_bin;
  }
  for(i=0;i<len;i++){
    bin_centers[i] = (bin_edges[i+1] + bin_edges[i])*0.5;
  }
  for(l=0;l<len;l++){		/* Is not put to zero in previous file. Might be a mistake */
    signal[l] = 0;
    var_temp[l] = 0;
  }
  

  // ext work
  for(k=0;k<traces;k++){
    ext_work[k] = (double *) malloc (sizeof(double) * samples);
    ext_work[k][0] = 0;
    for(i=1;i<samples;i++){
      ext_work[k][i] = ext_work[k][i-1] - 0.5*K*((x_ar[k][i] - lambda_ar[i])
						 +(x_ar[k][i-1] - lambda_ar[i-1])
						 )*(lambda_ar[i]-lambda_ar[i-1]);
    }
    //*work_mean += ext_work[k][samples-1]/traces;
  }
  // weighted work, avaraging over samples (time)
  for(i=0;i<samples;i++){
    exp_avarage[i] = 0;
    exp2_avarage[i] = 0;
    for(k=0;k<traces;k++){
      exp_avarage[i] +=  exp(-beta*ext_work[k][i]);
      exp2_avarage[i] += exp(-2*beta*ext_work[k][i]);
    }
    exp_avarage[i] = exp_avarage[i]/traces;
    exp2_avarage[i] = exp2_avarage[i]/traces;
  }
  // Avareging over recon bins.
  for(i=0;i<samples;i++){
    for(l=0;l<len;l++){
      if(bin_edges[l] < lambda_ar[i] && bin_edges[l+1] > lambda_ar[i]){
	// in this case the support is in l'th bin
	// printf("Hist hit. in  %E < %E <%E \n", bin_edges[l],prot[i], bin_edges[l+1samples]);
	signal[l] += exp_avarage[i];
	var_temp[l] += exp2_avarage[i];
	signal_counter[l] += 1;
      }
    }
  }
  for(l=0;l<len;l++){
    if(signal_counter[l] != 0){
      signal[l] = signal[l]/signal_counter[l]; /* Normalize if there ARE hits. */
      var_temp[l] = var_temp[l]/signal_counter[l];
    }
    //var[l] = (var_temp[l] - pow(signal[l],2))/traces;
    var[l] = (var_temp[l] - pow(signal[l],2))/sqrt(traces);
    //var[l] = (var_temp[l] - pow(signal[l],2));
    //var[l] = (var_temp[l] - pow(signal[l],2))/pow(traces, 0.8);
    //printf("vartest %E \n", pow(traces, 0.5));
  }
  



  if(c){
    s = vec_sum(signal, len);
    for(i=0;i<len;i++){signal[i] = signal[i]/s;}
  }
  
  // might include option for normalization here.
}

void get_signal(double * signal, int len, double z_start, double z_final,
		int traces, int samples, double * lambda_ar, double K,
		double ** x_ar, int c, double * bin_centers, double beta){

  int i,k,l;
  double s;
  double * bin_edges = (double *) malloc (sizeof(double) * len+1);
  //double * bin_centers = (double *) malloc (sizeof(double) * len);
  double * signal_counter = (double *) malloc (sizeof(double) * len);
  double ** ext_work = (double **) malloc ( sizeof(double * )* traces );
  double * exp_avarage = (double *) malloc (sizeof(double) * samples);
  double d_bin = (z_final - z_start)/len;
  // bins
  for(i=0;i<len+1;i++){
    bin_edges[i] = z_start + i*d_bin;
  }
  for(i=0;i<len;i++){
    bin_centers[i] = (bin_edges[i+1] + bin_edges[i])*0.5;
  }


  // ext work
  for(k=0;k<traces;k++){
    ext_work[k] = (double *) malloc (sizeof(double) * samples);
    ext_work[k][0] = 0;
    for(i=1;i<samples;i++){
      ext_work[k][i] = ext_work[k][i-1] - 0.5*K*((x_ar[k][i] - lambda_ar[i])
						 +(x_ar[k][i-1] - lambda_ar[i-1])
						 )*(lambda_ar[i]-lambda_ar[i-1]);
    }
    //*work_mean += ext_work[k][samples-1]/traces;
  }
  // weighted work
  for(i=0;i<samples;i++){
    exp_avarage[i] = 0;
    for(k=0;k<traces;k++){
      exp_avarage[i] +=  exp(-beta*ext_work[k][i]);
    }
    exp_avarage[i] = exp_avarage[i]/traces;
  }
  for(i=0;i<samples;i++){
    for(l=0;l<len;l++){
      if(bin_edges[l] < lambda_ar[i] && bin_edges[l+1] > lambda_ar[i]){
	// in this case the support is in l'th bin
	// printf("Hist hit. in  %E < %E <%E \n", bin_edges[l],prot[i], bin_edges[l+1]);
	signal[l] += exp_avarage[i];
	signal_counter[l] += 1;
      }
    }
  }
  for(l=0;l<len;l++){
    if(signal_counter[l] != 0){
      signal[l] = signal[l]/signal_counter[l]; /* Normalize if there ARE hits. */
    }
  }

  if(c){
    s = vec_sum(signal, len);
    for(i=0;i<len;i++){signal[i] = signal[i]/s;}
  }
  
  // might include option for normalization here.
 }

void fill_transfer_matrix(double ** H, int len1, int len2, double beta, double K,
			  double * bin_centers, double * bin_centers_reduced){
  // we will actually have to normalize this bugger!
  // most inefficient code ever, hopefully normalized though
  int i,j;
  double * norm_vector = (double *) malloc(sizeof(double)* len2);
  for(j=0;j<len2;j++){norm_vector[j] = 0;} /* setting norm-array to zero */
  for(i=0;i<len1;i++){
    for(j=0;j<len2;j++){
      //H[i][j] = sqrt(2*M_PI/K*beta)*exp(-0.5*beta*K*pow(bin_centers[i] - bin_centers_reduced[j],2)); /* proportional to this */
      H[i][j] = exp(-0.5*beta*K*pow(bin_centers[i] - bin_centers_reduced[j],2)); /* proportional to this */
      norm_vector[j] += exp(-0.5*beta*K*pow(bin_centers[i] - bin_centers_reduced[j],2));
    }
  }
  for(i=0;i<len1;i++){
    for(j=0;j<len2;j++){
      H[i][j] = H[i][j]/norm_vector[j];
    }
  }  
}

void fill_transfer_matrix_svd(double ** H, int len1, int len2, double beta, double K,
			      double * bin_centers, double * bin_centers_reduced){
  
  int i,j;
  double dz = bin_centers[1] - bin_centers[0]; /* all equidistant so this is ok */
  for(i=0;i<len1;i++){
    for(j=0;j<len2;j++){
      //H[i][j] = sqrt(2*M_PI/K*beta)*exp(-0.5*beta*K*pow(bin_centers[i] - bin_centers_reduced[j],2)); /* proportional to this */
      H[i][j] = dz*exp(-0.5*beta*K*pow(bin_centers[i] - bin_centers_reduced[j],2)); /* proportional to this */
    }
  }
}

void right_matrix_product(double ** A, double *x, double * y,  int len1, int len2){
  // perfroms y = Hx,
  int i,j;
  for(i=0;i<len1;i++){
    y[i] = 0;
    for(j=0;j<len2;j++){
      y[i] += A[i][j]*x[j];
    }
  }
}

void vector_addition(double * x, double *z, int len){
  //z = x+y
  int i;
  for(i=0;i<len;i++){z[i] = x[i] + z[i];}
}
// double * y;

void matrix_product(double ** C, double ** A, double ** B, int len1, int len2){
  // Not sure if it is ok to put in the same pointer for A and B here.
  int i,j,k;
  for(i=0;i<len2;i++){
    for(j=0;j<len2;j++){
      C[i][j] = 0;
      for(k=0;k<len1;){
	    C[i][j] += A[k][i]*B[k][j];
	// loop over data index. (small index.)
      }
    }
  }
}
void matrix_product2(double ** C, double ** A, int len1, int len2){
  // Not sure if it is ok to put in the same pointer for A and B here.
  int i,j,k;
  for(i=0;i<len2;i++){
    for(j=0;j<len2;j++){
      C[i][j] = 0;
      for(k=0;k<len1;k++){
	C[i][j] += A[k][i]*A[k][j];
	// loop over data index. (small index.)
      }
    }
  }
}
void update_A(double ** A, double ** H, int len, int N_recon, double prior_sigma, double * P, double alpha){
  int k,j;
  double temp = pow(prior_sigma,-2);
  matrix_product2(A, H, len, N_recon);
  // multiply by inverse variance
  for(j=0;j<N_recon;j++){
    for(k=0;k<N_recon;k++){
      A[j][k] = A[j][k]*temp;
    }
  }
  // add to diagonal
  for(k=0;k<N_recon;k++){
    A[k][k] += alpha/P[k];
    //printf("A[k][k] = %E \n", A[k][k]);
  }
  
}

int step_A(double ** A, double ** H, int N_recon, double * P, double alpha, 
		double * u, double * f, double * s, double * precond_grad, 
		double * u_old, int len, double prior_sigma){
  
  int i,j,k;
  double top_temp = 0;
  double bottom_temp = 0;
  /* matrix_product(A, H, H, len, N_recon); /\* updates A *\/ */
  //matrix_product2(A, H, len, N_recon); /* updates A */
  update_A(A,H,len,N_recon,prior_sigma,P,alpha);
  for(j=0;j<N_recon;j++){
    //A[j][j] += alpha/P[j];
    precond_grad[j] = pow(A[j][j],-1)*(f[j] - alpha*s[j]); /* this is search direction in first iteration */
    if(A[j][j] == 0){printf("Warning! diagonal A[j][j] = 0, meaning precond grad is infinite \n");}
    u_old[j] = u[j];					   /* save this */
    //printf("from A: u[%d] = %E \n", j , u[j]);
  }
  // calculating factors in eq 25
  for(i=0;i<N_recon;i++){
    for(k=0;k<N_recon;k++){
      top_temp += u_old[i]*A[i][k]*precond_grad[k];
      bottom_temp += u_old[i]*A[i][k]*u_old[k];
      //printf("A[%d][%d] = %E \n", i,k, A[i][k]);
    }
  }
  //printf("BOTTOM_TEMP %E \n", bottom_temp);
  if(bottom_temp == 0){printf("Warning! direction vector u undefined. \n");}
  
  // Updating new search direction
  for(i=0;i<N_recon;i++){
    //printf("precond_grad[%d] = %E || precond_grad_prev[%d] = %E \n",i, precond_grad[i], i, u_old[i] );
    u[i] = precond_grad[i] - (top_temp/bottom_temp)*u_old[i];
    //printf("from A: u[%d] = %E \n", i , u[i]);
  }
  return 1;
}
// double * v;

int step_B(double * S_u, double * F_u, double *F_uu, double * v, double sigma_prior, int J, int I,
	   double *u, double *s, double * f, double ** H){
  int j;
  *S_u=0;*F_u=0;*F_uu=0;
  double u_size = 0;		/* just for debug purposes */
  double v_size = 0;
  //for(i=0;i<I;i++){printf("v[%d] = %E \n", i, v[i]);}
  //for(i=0;i<J;i++){printf("u[%d] = %E \n", i, u[i]);} /* this is zeros. */
  right_matrix_product(H, u, v,  I, J);
  for(j=0;j<J;j++){
    *S_u +=u[j]*s[j];
    *F_u +=u[j]*f[j];
    u_size += u[j]*u[j];
    v_size += v[j]*v[j];	/* v is dim I so this underestimates it. */
    //printf("%E \n",pow(v[j],2));//*pow(v[j],2));
    *F_uu += pow(sigma_prior,-2)*pow(v[j],2);
  }
  //printf("F_u = %E \n", *F_u);
  /* printf("step B checks ------- \n"); */
  /* printf("F_uu = %E \n", *F_uu); */
  /* printf("|u|^2 = %E \n", u_size); */
  /* printf("|v|^2 = %E \n", v_size); */
  /* printf("S_u = %E \n", *S_u); */
  /* printf("F_u = %E \n", *F_u); */
  
  if(*F_uu==0){return 1;}else{return 0;} /* 1 we break */
}

int step_C(double * xi_0, double F_u, double F_uu, double S_u,
	   double epsilon, double * P, double *u, int J, double epsilon_prime,
	   double alpha){

  int control = 0;		/* if we update we increment, so that we return True if something was changed */
  int i;
  //printf("C again. \n");
  *xi_0 = -(F_u - alpha*S_u)/F_uu;
  //printf("sign check xi_0 = %E \n", *xi_0);
  for(i=0;i<J;i++){
    if(P[i] + (*xi_0*u[i]) < epsilon*P[i]){
      control=1;
      //printf("u[%d] = %E ---> ", i, u[i]);
      u[i] = ((epsilon_prime*P[i] - P[i])/(*xi_0)); /* you should only have to do thi once? Why does it chew on this... */
      //printf("%E \n", u[i]);
    }
  }
  return control;		/* if return is 1, we changed u. And must go back to step B */
}

void step_D_init(double * w, double ** H, double *xi, int J, int I, double * v){
  int i,j;
  *xi = 0;			/* initialize. */
  for(j=0;j<J;j++){
    w[j] = 0;
    for(i=0;i<I;i++){
      w[j] += H[i][j]*v[i];
    }
  }
}

int step_D(double * w, double ** H, double *u, double ** A, double * xi,int J, int I,
	   double * G_uu, double * F_u, double *S_u, double * alpha, double * F, double * s,
	   double * M, double * f, double sigma_prior, double xi_threshold, double * P, double F_uu){
  int j,k;
  double d_xi;
  double temp1=0;double temp2=0;
  update_A(A,H, I, J, sigma_prior, P, *alpha);
  *G_uu = 0;
  for(j=0;j<J;j++){
    for(k=0;k<J;k++){
      *G_uu += A[j][k]*u[j]*u[k];
    }
  }
  // substep a;
  d_xi = (-(*F_u - (*alpha)*(*S_u))/(*G_uu));
  // substeb b,c:
  *xi += d_xi;			/* Warning! There must be some typo in the code here. this does not make sense */
  //*xi += fabs(d_xi);		/* WARNING WARNING WARNING AD HOC CHANGE NO BASIS FOR IT!!!!!! */
  //printf("d_xi = %E \n", d_xi);
  
  
  
  //printf("SEGTEST 1 \n");
  for(j=0;j<J;j++){
    P[j] += d_xi*u[j];
    if(P[j]<0){printf("Warning! Negative pixel value.\n");}
    //printf("From Update step: P[%d] = %E \n", j, P[j]);
    s[j] = -log(P[j]/(M[j])); /* to make flat distribution take maximum entropy value. */
    if(P[j]/M[j] <= 0){printf("Warning! Invaid entropy \n");}
    f[j] += d_xi*pow(sigma_prior, -2)*w[j];
    temp1 += P[j]*pow(f[j],2);
    temp2 += P[j]*pow(s[j],2);
  }
  
  *F += 2*d_xi*(*F_u) + pow(d_xi,2)*F_uu; /* this is the culprit! */
  *F_u += d_xi*F_uu;
  // uncomment this if things get weird
  //printf("S_u_old = %E -->", *S_u);
  *S_u = 0;
  for(j=0;j<J;j++){
    *S_u += s[j]*u[j];
  }
  //printf("%E \n", *S_u);
  // You have to recalculate s here. But you dont touch S_u ... Why? This must reasonably be updated too?
  //printf("SEGTEST 2 \n");
  
  // updating alpha
  if(temp2==0){*alpha=ALPHA_GUARD;}else{*alpha = sqrt((I*temp1)/((*F)*temp2));}
  //printf("ALPHA CHECK: %E \n", *alpha);
  //printf("F CHECK: %E \n", *F);
  
  // Run loop if d_xi is too big! that is if we return 1, return of zero means break of loop.
  // This should probably check that xi is in the interval xi_1 < xi < xi_2
  // where x_1 and x_2 is given by eq 34. Thats not what it says though. .. 
  if(fabs(d_xi) >= fabs(xi_threshold*(*xi))){
    //printf("xi = %E  | |d_xi/xi| = %E \n",*xi, fabs(d_xi)/fabs(*xi));
    return 1;
  }else{
    return 0;
  }
}
int step_E(double xi, double * u, int J, double P_threshold){
  // This is just to test for convergence. If dp is diminishing appropriately then its probably working.
  int j;
  double temp = 0;
  double temp2;
  for(j=0;j<J;j++){
    temp += u[j]*u[j];
  }
  // must figure out some stopping condition here...
  // this ought to take come counter. If
  temp2 = temp;
  temp = sqrt(temp2)*fabs(xi);
  printf("dp = %E \n", temp);
  if(temp > P_threshold){
    return 1;
  }
  else{
    return 0;
  }
}
//int iter_count

void max_ent(double * signal, int len, double z_start, double z_final,
	     int traces, int samples, double * lambda_ar, double K,
	     double ** x_ar, int c, double * bin_centers, int cutoff,
	     double * G, double beta, double * bin_centers_reduced, double prior_sigma,
	     int MAX_ITER, double U, double sigma, double * U_lc_surf){

  get_signal(signal,len,z_start,z_final,
	     traces,samples,lambda_ar,K,
	     x_ar,c, bin_centers,beta);
  
  // 0) STORAGE
  int i,j,k,l,N_recon;
  double temp1,temp2;
  N_recon = len-cutoff;

  int control_a,control_b,control_c,control_d,control_e;
  printf("R test %d \n", MAX_ITER);
  MAX_ITER=600;
  //MAX_ITER=1;
  double F,F_u,F_uu;
  double S_u;
  double G_uu;
  double alpha;
  double xi_0;
  double xi;
  double xi_threshold= 0.5*pow(10,-3);		/* There are no specifications for this. */
  double P_threshold=pow(10,-10);
  
  // WARNING! SET THESE BEFORE RUNNING, values taken from article recommendations.
  double epsilon = 0.1;		
  double epsilon_prime = 0.5;
  
  double * P = (double *) malloc (sizeof(double) * (N_recon));
  double * M = (double *) malloc (sizeof(double) * (N_recon));
  double * s = (double *) malloc (sizeof(double) * (N_recon));
//  double * r = (double *) malloc (sizeof(double) * (N_recon));
  double * u = (double *) malloc (sizeof(double) * (N_recon));
  double * u_old = (double *) malloc (sizeof(double) * (N_recon));
  double * v = (double *) malloc (sizeof(double) * (len));
  double * w = (double *) malloc (sizeof(double) * (N_recon));
  
  double * e = (double *) malloc (sizeof(double) * (len));
  double * f = (double *) malloc (sizeof(double) * (N_recon));
  double ** H = (double **) malloc (sizeof(double *) * (len));
  double ** A = (double **) malloc (sizeof(double *)*N_recon);
  double * precond_grad = (double *) malloc (sizeof(double) * (N_recon));

  for(l=0;l<len-cutoff;l++){
    bin_centers_reduced[l] = bin_centers[l+cutoff/2];
  }
  // uncomment for gaussian innputs
  /* for(l=0;l<len-cutoff;l++){ */
  /*   signal[l] = exp(-pow(bin_centers_reduced[l] -  bin_centers_reduced[l/2],2)); */
  /* } */
  
  
  for(i=0;i<len;i++){
    H[i] = (double *) malloc(sizeof(double)*N_recon);
  }
  fill_transfer_matrix(H, len, N_recon,beta, K,
		       bin_centers, bin_centers_reduced);
  
  // I) INITIALIZE PIXEL ARRAY AND COMPUTE ENTROPY GRADIENT. All of these steps are instantiations.
  for(i=0;i<(N_recon);i++){
    //P[i] = exp(-100*pow(bin_centers_reduced[i],2));
    P[i] = 1.0/(N_recon);      /* normalized to 1 */
    M[i] = 1.0/(N_recon);      /* flat prior. Needed later. */
    s[i] = -log(P[i]/(M[i]));    /* just for clarity. (Takes no time.) */
    A[i] = (double *) malloc(sizeof(double)*N_recon); /* instantiating matrix */
  }

  // II) THE LONGEST STEP
  F=0;
  for(i=0;i<len;i++){
    e[i] =  signal[i];
    for(j=0;j<N_recon;j++){
      e[i] += -H[i][j]*P[j];
      //printf("H[%d][%d] = %E \n", i,j,H[i][j]);
    }
    F += e[i]*e[i]*pow(prior_sigma,-2);
  }
  //
  for(j=0;j<N_recon;j++){
    f[j] = 0;
    for(i=0;i<len;i++){
      f[j] += -pow(prior_sigma,-2)*H[i][j]*e[i];
      //printf("f[%d] = %E \n", j, f[j]);
    }
  }
  
  //
  temp1=0;temp2=0;
  for(j=0;j<N_recon;j++){
    temp1+=P[j]*pow(f[j],2);
    temp2+=P[j]*pow(s[j],2);
    //printf("f[%d] = %E , s[%d] = %E \n", j,f[j], j, s[j]);
  }
  if(temp2==0){alpha=ALPHA_GUARD;}else{alpha = sqrt((len*temp1)/(F*temp2));} /* this is OK. If all s_i vanish then its multiplicator s does not matter */
  printf("ALPHA CHECK: %E \n", alpha);
  printf("F CHECK: %E \n", F);
  //alpha = sqrt((N_recon*temp1)/(F*temp2));
  
  // THIS IS WHERE YOU WERE WHEN YOU LEFT MONDAY, YOU NEEDED TO INSTANTIATE A,v and w and G_uu
  // instantiate G_uu,v and w  
  //matrix_product2(A, H, len, N_recon); /* updates A */
  update_A(A,H,len,N_recon,prior_sigma,P,alpha);
  //printf("alpha = %E \n", alpha);      /* heres your problem */
  for(j=0;j<N_recon;j++){
    //A[j][j] += alpha/P[j];
    //printf("P[%d] = %E \n", j, P[j]);
    //printf("A[%d][%d] = %E \n", j,j, A[j][j]);
    precond_grad[j] = pow(A[j][j],-1)*(f[j] - alpha*s[j]); /* this is search direction in first iteration */
    u[j] = precond_grad[j];
    //printf("u[%d] = %E \n", j, u[j]);
   }
  for(i=0;i<len;i++){
    v[i] = 0;
    for(j=0;j<N_recon;j++){
      v[i] += H[i][j]*u[j];
    }
  }
  G_uu=0;
  for(j=0;j<N_recon;j++){
    w[j] = 0;
    for(i=0;i<len;i++){
      w[j] += H[i][j]*v[i];
    }
    for(k=0;k<N_recon;k++){
      G_uu += A[j][k]*v[k]*v[j];
    }
  }

  // print some stuff and see
  /* for(i=0;i<len;i++){printf("v[%d] = %E \n", i, v[i]);} */
  write_1darray(f, N_recon,"/tmp/f_astrid.dat");
  write_2darray(H, len, N_recon, "/tmp/H_astrid.dat");
  
  
  // THIS IS WHERE THE BIG LOOP STARTS.
  control_b=0;
  control_c=1;
  control_d=1;
  control_a=1;	// why did these change place?
  control_e=1;
  int iter_count = 0;
  int c_iter = 0;
  int ITER_MIN=3;
  while(control_e){
    control_d=1;control_b=0;control_c=1;control_e=1;
    printf("A loop. \n");
    if(iter_count !=0){
      //printf("UPDATING CONJ-DIRECTIONS \n");
      control_a=step_A(A,H,N_recon,P,alpha,u,
		       f,s,precond_grad,u_old,len,prior_sigma);
    } /* dont do this first step. There are no new search directions! */
    //printf("A testing. \n");
    while(control_c){
      //printf("C loop. \n");
      // run while step_C returns 1, for some change in u.
      control_b = step_B(&S_u,&F_u,&F_uu,v,prior_sigma,N_recon, len,
			 u, s, f, H);
      if(control_b){
	    break;			/* case means F_uu=0 */
      }
      // THIS STEP BELOW CHANGES SEARCH DIRECTIONS TO PREVENT NEGATIVE PIXEL VALUES.
      c_iter++;
      //if(c_iter<20){printf("%d ,xi_0 = %E \n",c_iter, xi_0);}
      control_c = step_C(&xi_0,F_u,F_uu,S_u,
      			 epsilon,P,u, N_recon,epsilon_prime,alpha);
      /* control_c=0; */
    }

    //for(j=0;j<N_recon;j++){printf("After conv: u[%d] = %E \n", j, u[j]);}
    /* if(control_b){ */
    /*   break;			/\* if we break because F_uu, break it all. *\/ */
    /* } */

    step_D_init(w,H,&xi,N_recon,len,v);
    printf("Starting D loop. Time consuming part. \n");
    while(control_d){
      //printf("D_step. \n");
      control_d = step_D(w,H,u,A,&xi,N_recon,len,
      			 &G_uu,&F_u,&S_u,&alpha,&F,s,
      			 M,f,prior_sigma,xi_threshold,P, F_uu);
      //printf("D_step done. \n");
      
    }
    printf("ALPHA CHECK: %E \n", alpha);
    printf("F: %E \n", F);

    if(iter_count > ITER_MIN){
      control_e = step_E(xi,u,N_recon,P_threshold);
    }
    
    iter_count++;
    if(iter_count>MAX_ITER){control_e=0;}
  }

double g_const = U_surface(bin_centers_reduced[0],U, sigma) - (-pow(beta,-1)*log(P[0]));
  printf("g_const = %E \n", g_const);
  
  /* // when all is done. Update vector */
  for(l=0;l<N_recon;l++){
    U_lc_surf[l] = U_surface(bin_centers_reduced[l],U, sigma);
    //U_lc_surf[l] = 0;
    if(P[l] > 0){
      G[l] = -pow(beta,-1)*log(P[l]);
    }
    G[l] += g_const;
  }

  /* // We are supposedly done at this stage. Time to convert the probabilities. */
  /* for(j=0;j<N_recon;j++){ */
  /*   G[j] = -pow(beta,-1)*log(P[j]); */
  /*   printf("P[%d] = %E || G[%d] = %E \n", j, P[j], j, G[j]); */
  /* } */

  /* write_1darray(G, N_recon,"/tmp/G_astrid.dat"); */
  /* write_1darray(P, N_recon,"/tmp/P_astrid.dat"); */  
}

void normed_lc(double ** x_ar, double * lambda_ar,
	       int traces, int samples, double K, double beta, int recon_samples,
	       double z_start, double z_final, 
	       double epsilon, double sigma, double * U_lc_surf,
	   	   double * bin_centers, double *lc_bins, double * G_lc, 
		   int cutoff, int ITER_MAX){
  int i,j,l,I,J,m;
  double P_j_temp;
  I = recon_samples;
  J = recon_samples-cutoff;
  double * signal = (double *) malloc(sizeof(double) * I);
  double * P = (double *) malloc(sizeof(double) * J);
  double * c = (double *) malloc(sizeof(double) * I);
//  double * U = (double *) malloc(sizeof(double) * I);
//  double * U_tilde = (double *) malloc(sizeof(double) * I);

  
  get_signal(signal,recon_samples,z_start,z_final,
	     traces,samples,lambda_ar,K,
	     x_ar,1, bin_centers,beta); /* meaning normalized. */
  for(l=0;l<J;l++){
    lc_bins[l] = bin_centers[l+cutoff/2];
    P[l] = 1.0/(J);		/* initial guess. */
  }
  double ** H = (double **) malloc (sizeof(double *) * (I));
  for(i=0;i<I;i++){H[i] = (double *) malloc(sizeof(double)*J);}

  fill_transfer_matrix(H, I, J,beta, K,
		       bin_centers, lc_bins);
  
  // Lets begin. If you want to see the evolution you should let P be a two dimensinal array of dimensions (ITER_MAX,J)
  for(m=0;m<ITER_MAX;m++){  
    // first we need to fill the c_i
    for(i=0;i<I;i++){
      c[i] = 0;
      for(j=0;j<J;j++){
	c[i] += H[i][j]*P[j]; 
      }
      //U[i] = -2*pow(T_param,-2)*(signal[i]*log(c[i]/signal[i]) - c[i] + signal[i]);
      //U_tilde[i] = min_f(U[i],1);
    }
    // c_i filled
    // updating P[j]
    for(j=0;j<J;j++){
      P_j_temp = P[j];
      P[j] = 0; // Scale problems if you leave this in. But it reasonably SHOULD be there.. Doesnt seem like it works...
      for(i=0;i<I;i++){
	P[j] += P_j_temp*(signal[i]/c[i])*H[i][j];
      }
    }  
  }

  /* for(j=0;j<J;j++){ */
  /*   G_lc[j] = -pow(beta,-1)*log(P[j]); */
  /* } */
    //printf("P[%d] = %E || G[%d] = %E \n", j, P[j], j, G[j]);
  
  // OFFSETTING
  double g_const = U_surface(lc_bins[0],epsilon, sigma) - (-pow(beta,-1)*log(P[0]));
  printf("g_const = %E \n", g_const);
  
  /* // when all is done. Update vector */
  for(l=0;l<J;l++){
    U_lc_surf[l] = U_surface(lc_bins[l],epsilon, sigma);
    //U_lc_surf[l] = 0;
    if(P[l] > 0){
      G_lc[l] = -pow(beta,-1)*log(P[l]);
    }
    G_lc[l] += g_const;
  }
}
//double* t_ar, G, recon_splope, recon_msd, work_mean, work_std, exp2_average, average work, rms, samples_work_variance, prot;
//int steps
//int * invalid count;

void get_weierstrass_signal(double * signal, double * A_prim, 
		double * A_primprim, int len, double z_start, double z_final,
		int traces, int samples, double * lambda_ar, double K,
		double ** x_ar, double * bin_centers, double beta){

  int i,k,l;
  //double s;
  double * bin_edges = (double *) malloc (sizeof(double) * len+1);
  //double * bin_centers = (double *) malloc (sizeof(double) * len);
  double * signal_counter = (double *) malloc (sizeof(double) * len);
  double ** ext_work = (double **) malloc ( sizeof(double * )* traces );
  double * exp_avarage = (double *) malloc (sizeof(double) * samples);
  double d_bin = (z_final - z_start)/len;
  /* double A_prim = (double *) malloc (sizeof(double)*(len)); /\* see hummer2 for descriptions of these. Not really as named! *\/ */
  /* double A_primprim = (double *) malloc (sizeof(double)*(len)); */
  double * F2 = (double *) malloc (sizeof(double)*(len)); /* ensemble avarage of force**2. See hummer2! */
  
  // bins
  for(i=0;i<len+1;i++){
    bin_edges[i] = z_start + i*d_bin;
  }
  for(i=0;i<len;i++){
    bin_centers[i] = (bin_edges[i+1] + bin_edges[i])*0.5;

    F2[l] = 0;
    A_prim[l] = 0;		/* init to zero */
    A_primprim[l] = 0;
    signal[l] = 0;
  }

  // ext work
  for(k=0;k<traces;k++){
    ext_work[k] = (double *) malloc (sizeof(double) * samples);
    ext_work[k][0] = 0;
    for(i=1;i<samples;i++){
      ext_work[k][i] = ext_work[k][i-1] - 0.5*K*((x_ar[k][i] - lambda_ar[i])
						 +(x_ar[k][i-1] - lambda_ar[i-1])
						 )*(lambda_ar[i]-lambda_ar[i-1]);
    }
    //*work_mean += ext_work[k][samples-1]/traces;
  }
  // weighted work
  for(i=0;i<samples;i++){
    exp_avarage[i] = 0;
    for(k=0;k<traces;k++){
      exp_avarage[i] +=  exp(-beta*ext_work[k][i]);
    }
    exp_avarage[i] = exp_avarage[i]/traces;
  }
  for(i=0;i<samples;i++){
    for(l=0;l<len;l++){
      if(bin_edges[l] < lambda_ar[i] && bin_edges[l+1] > lambda_ar[i]){
	// in this case the support is in l'th bin
	// printf("Hist hit. in  %E < %E <%E \n", bin_edges[l],prot[i], bin_edges[l+1]);
	signal[l] += exp_avarage[i];
	signal_counter[l] += 1;
	for(k=0;k<traces;k++){	/* Very dubious! */
	  A_prim[l] += -K*exp(-beta*ext_work[k][i])*(x_ar[k][i] - bin_centers[l])/traces;
	  F2[l] += pow(K,2)*exp(-beta*ext_work[k][i])*pow(x_ar[k][i] - bin_centers[l],2)/traces;
	    }	
      }
    }
  }
  for(l=0;l<len;l++){
    if(signal_counter[l] != 0){
      signal[l] = signal[l]/signal_counter[l]; /* Normalize if there ARE hits. */
      A_prim[l] = A_prim[l]/signal_counter[l];
      F2[l] = F2[l]/signal_counter[l];
    }
  }
  for(l=0;l<len;l++){
    A_primprim[l] = beta*(F2[l] - pow(A_prim[l],2))/K;
  }

  // Skip normalization if its not to be used for image reconstruction.
  /* if(c){ */
  /*   s = vec_sum(signal, len); */
  /*   for(i=0;i<len;i++){signal[i] = signal[i]/s;} */
  /* } */
  
  // might include option for normalization here.
}
// int c;

void weierstrass_deconv(double * signal, int len, double z_start, double z_final,
			int traces, int samples, double * lambda_ar, double K,
			double ** x_ar, double * bin_centers, double beta,
			double * lc_bins, double * G_lc, double * U_lc_surf, double U, double sigma){
  
  int k,l;
  
  double * A_prim = (double *) malloc (sizeof(double)*len);
  double * A_primprim = (double *) malloc (sizeof(double)*len);
  //shifted_bin_centers = (double *) malloc (sizeof(double)*len);
  
  get_weierstrass_signal(signal,A_prim,A_primprim,
			 len,z_start,z_final,
			 traces,samples,lambda_ar,K,
			 x_ar,bin_centers,beta);
  for(l=0;l<len;l++){
    G_lc[l] = -pow(beta,-1)*log(signal[l]) - pow(A_prim[l],2)*pow(2*K,-1) + pow(2*beta,-1)*log(A_primprim[l]);
    lc_bins[l] = bin_centers[l] - A_prim[k]/K;
    U_lc_surf[l] = U_surface(lc_bins[l],U, sigma); 
  }
  double g_const = U_surface(lc_bins[0],U,sigma) - G_lc[0];
  for(l=0;l<len;l++){
    G_lc[l] += g_const;
  }
}
// int c

void wiener_deconvolution(double * signal, int len, double z_start, double z_final,
			  int traces, int samples, double * lambda_ar, double K,
			  double ** x_ar, int c, double * bin_centers, double beta,
			  double * lc_bins, double * G_lc, double * U_lc_surf, double U, double sigma,
			  double regularization_param, int cutoff){
  // VERY IMPORTANT. THE ARGUMENT cutoff MUST BE SCRIPT ARGUMENT DIVIDED BY 2 (THIS IS ALWYS EVEN SO ITS FINE)
  // for now start with this. Later we can try to find ways to estimate this regularization in some smarter way.
  int i,count;
  double dz,norm;			/* bin size */
  char temp_char[256];
  char line[256];
  get_signal(signal,len,z_start,z_final,
	     traces,samples,lambda_ar,K,
	     x_ar,c, bin_centers,beta);
  double * kernel = (double *) malloc (sizeof(double)*len);
  double * deconv = (double *) malloc (sizeof(double)*len);

  dz = (z_final - z_start)/len;
  printf("dz = %E \n", dz);		       /* just checking that it has correct units */
  norm = 0;
  for(i=0;i<len;i++){kernel[i]=0;deconv[i]=0;}                      /* for simplicity  */
  for(i=0;i<cutoff;i++){
    signal[i] = 0;                  		        /* pad right */
    signal[len-1-i] = 0;	                	/* pad left  */
    kernel[i] = dz*exp(-0.5*beta*K*pow(i*dz,2));               /* not normal */
    kernel[len-1-i] = dz*exp(-0.5*beta*K*pow(i*dz,2));
    norm += kernel[i] + kernel[len-1-i];
  }
  // Testing to normalize it. This doesnt seem to have any effect at all.
  for(i=0;i<len;i++){kernel[i] = kernel[i]/norm;}
  // must put these guys in a temporary file
  FILE *outfile;
  outfile = fopen("pytemp/temp_arrays.dat", "w+");
  for(i=0;i<len;i++){
    fprintf(outfile,"%E \t %E \n", signal[i], kernel[i]);
    printf("%E \t %E \n", signal[i], kernel[i]);
  }
  fclose(outfile);
  
  // Calling python script.
  sprintf(temp_char, "python wiener_filter.py %E", regularization_param);
  system(temp_char);

  // reading file.
  count=0;
  outfile = fopen("pytemp/deconvolved_temp.dat", "r");
  while(fgets(line, sizeof(line), outfile)){
    deconv[count] = atof(line);
    count++;
  }
  
  // Offseting and creating U_surf etc.
  for(i=0;i<len-2*cutoff;i++){
    lc_bins[i] = bin_centers[cutoff + i];
    G_lc[i] = - pow(beta,-1)*log(deconv[i+cutoff]);
    U_lc_surf[i] = U_surface(lc_bins[i],U, sigma);
  }
  double g_offset = U_lc_surf[0] - G_lc[0];
  for(i=0;i<len-2*cutoff;i++){G_lc[i] += g_offset;}
}
//int c;

void svd_deconvolution(double * signal, int len, double z_start, double z_final,
		       int traces, int samples, double * lambda_ar, double K,
		       double ** x_ar, int c, double * bin_centers, double beta,
		       double * lc_bins, double param, int cutoff){

  int i,l;
  char temp_char1[256];
  //char temp_char2[256];
  double ** H = (double **) malloc (sizeof(double *)* len);
  
  for(i=0;i<len;i++){H[i] = (double *) malloc(sizeof(double)*(len-cutoff));}

  get_signal(signal,len,z_start,z_final,
	     traces,samples,lambda_ar,K,
	     x_ar,c, bin_centers,beta);
  for(l=0;l<len-cutoff;l++){
    lc_bins[l] = bin_centers[l+cutoff/2];
  }
  fill_transfer_matrix_svd(H, len, len-cutoff,  beta, K,
			   bin_centers,lc_bins); /* normalized. */

  // At this point we need to pipe this guys to the python file
  write_2darray(H, len, len-cutoff, "pytemp/H.dat");
  write_1darray(signal, len, "pytemp/temp_array.dat");

  // Calling svd_deconv.py
  sprintf(temp_char1, "python svd_deconv.py %E", param); /* not sure what param is yet. */
  system(temp_char1);					 /* calling python */

  // Reading svd_deconv.py out
  //FILE * outfile;
  /* count=0; */
  /* outfile = fopen("pytemp/svd_deconv.dat", "r"); */
  /* while(fgets(temp_char2, sizeof(temp_char2), outfile)){ */
  /*   deconv[count] = atof(temp_char2); */
  /*   count++; */
  /* } */  
}
// double* G_lc, U_lc_sur, U, sigma;

void svd_deconvolution2(double * signal, int len, double z_start, double z_final,
		int traces, int samples, double * lambda_ar, double K,
			double ** x_ar, int c, double * bin_centers, double beta,
			double * lc_bins, double param, int cutoff){

  // A second version where the errors are somewhat taken into account.
  // there will be some correlations, so maybe one will have to
  // multiply by cholesky factor in order to reduce problem to one with diagonal covariance.
 
  int i,l;
  char temp_char1[256];
  //char temp_char2[256];
  double ** H = (double **) malloc (sizeof(double *)* len);
  double * var = (double *) malloc (sizeof(double) *len); /* assuming diagonal to start with. */
  
  for(i=0;i<len;i++){H[i] = (double *) malloc(sizeof(double)*(len-cutoff));}
  
  get_signal_with_var(signal,len,z_start,z_final,
    traces,samples,lambda_ar,K,x_ar,c, bin_centers,beta,var);
  
  for(l=0;l<len-cutoff;l++){
    lc_bins[l] = bin_centers[l+cutoff/2];
  }
  fill_transfer_matrix_svd(H, len, len-cutoff,  beta, K,
			   bin_centers,lc_bins); /* non-normalized. */

  // At this point we need to pipe this guys to the python file
  write_2darray(H, len, len-cutoff, "pytemp/H.dat");
  write_1darray(signal, len, "pytemp/temp_array.dat");
  write_1darray(var, len, "pytemp/var.dat");
  
  // Calling svd_deconv.py
  sprintf(temp_char1, "python svd_deconv2.py %E", param); /* not sure what param is yet. */
  system(temp_char1);					 /* calling python */

  // Reading svd_deconv.py out
  //FILE * outfile;
  /* count=0; */
  /* outfile = fopen("pytemp/svd_deconv.dat", "r"); */
  /* while(fgets(temp_char2, sizeof(temp_char2), outfile)){ */
  /*   deconv[count] = atof(temp_char2); */
  /*   count++; */
  /* } */  
}
// double* G_lc, U_lc_surf, 
// double U, sigma;

void twomey_deconvolution(double * signal, int len, double z_start, double z_final,
			  int traces, int samples, double * lambda_ar, double K,
			  double ** x_ar, int c, double * bin_centers, double beta,
			  double * lc_bins, double * U_lc_surf, double U, double sigma,
			  double param, int cutoff){
  
  int i,l;
  char temp_char1[256];
//  char temp_char2[256];
  double ** H = (double **) malloc (sizeof(double *)* len);
  double * U_surf = (double *) malloc(sizeof(double) * len);
  
  
  get_signal(signal,len,z_start,z_final,
	     traces,samples,lambda_ar,K,
	     x_ar,c, bin_centers,beta);

  for(i=0;i<len;i++){
    H[i] = (double *) malloc(sizeof(double)*(len-cutoff));
    U_surf[i] =U_surface(bin_centers[i], U, sigma);
  }
  for(l=0;l<len-cutoff;l++){
    lc_bins[l] = bin_centers[l+cutoff/2];
  }
  fill_transfer_matrix_svd(H, len, len-cutoff,  beta, K,
  			   bin_centers,lc_bins); /* normalized. */

  /* fill_transfer_matrix(H, len, len-cutoff,  beta, K, */
  /* 		       bin_centers,lc_bins); /\* normalized. *\/ */
  
  /* double F_0 = U_surface(lc_bins[0],U, sigma); /\* confused by this :S *\/ */
  for(l=0;l<len-cutoff;l++){U_lc_surf[l] = U_surface(lc_bins[l], U, sigma);}

  // At this point we need to pipe this guys to the python file
  write_2darray(H, len, len-cutoff, "pytemp/H.dat");
  write_1darray(signal, len, "pytemp/temp_array.dat");
  write_1darray(lc_bins, len-cutoff, "pytemp/temp_bins.dat");
  write_1darray(bin_centers, len, "pytemp/temp_bins_full.dat");
  write_1darray(U_lc_surf, len-cutoff, "pytemp/temp_U_surface.dat");
  write_1darray(U_surf, len, "pytemp/temp_U_surface_full.dat");
  
  // Calling svd_deconv.py
  sprintf(temp_char1, "python twomey.py %E", param); /* not sure what param is yet. */
  system(temp_char1);					 /* calling python */

  // Reading svd_deconv.py out
  //FILE * outfile;
  /* count=0; */
  /* outfile = fopen("pytemp/svd_deconv.dat", "r"); */
  /* while(fgets(temp_char2, sizeof(temp_char2), outfile)){ */
  /*   deconv[count] = atof(temp_char2); */
  /*   count++; */
  /* } */  
}
// double* G_lc 

void twomey_deconvolution2(double * signal, int len, double z_start, double z_final,
			   int traces, int samples, double * lambda_ar, double K,
			   double ** x_ar, int c, double * bin_centers, double beta,
			   double * lc_bins, double * G_lc, double * U_lc_surf, double U, double sigma,
			   double param, int cutoff, double * wham_recon_arr, double g_const){
  int i,l;
  char temp_char1[256];
  //char temp_char2[256];
  double * reduced_wham_recon = (double *) malloc (sizeof(double)*(len-cutoff));
  double ** H = (double **) malloc (sizeof(double *)* len);
  double * U_surf = (double *) malloc(sizeof(double) * (len-cutoff));

  get_signal(signal,len,z_start,z_final,
	     traces,samples,lambda_ar,K,
	     x_ar,c, bin_centers,beta);
  for(l=0;l<len-cutoff;l++){
    lc_bins[l] = bin_centers[l+cutoff/2];
    reduced_wham_recon[l] = exp(-beta*wham_recon_arr[l + cutoff/2]);
    //reduced_wham_recon[l] = exp(-beta*wham_recon_arr[l + cutoff/2])*exp(beta*g_const); // Using offset.
    U_surf[l] = U_surface(lc_bins[l], U, sigma);
  }
  // constant in integral expression is not relevant normally BUT in this case it is because the seed
  // constrains the solution
  for(i=0;i<len;i++){
    H[i] = (double *) malloc(sizeof(double)*(len-cutoff));
    signal[i] = signal[i]*exp(beta*g_const);;
  }
  fill_transfer_matrix_svd(H, len, len-cutoff,  beta, K,
  			   bin_centers,lc_bins);

  // writing some files to python
  write_2darray(H, len, len-cutoff, "pytemp/H.dat");
  write_1darray(signal, len, "pytemp/temp_array.dat");
  write_1darray(reduced_wham_recon, len-cutoff, "pytemp/seed.dat");
  write_1darray(lc_bins, len-cutoff, "pytemp/temp_bins.dat");
  write_1darray(bin_centers, len, "pytemp/temp_bins_full.dat");
  write_1darray(U_surf, len-cutoff, "pytemp/reduced_surface.dat");
  
  sprintf(temp_char1, "python twomey2.py %E", param);
  system(temp_char1);
  
  /* FILE * outfile; */
  /* count=0; */
  /* outfile = fopen("pytemp/f.dat", "r"); */
  /* while(fgets(temp_char2, sizeof(temp_char2), outfile)){ */
  /*   G_lc[count] = atof(temp_char2); */
  /*   count++; */
  /* } */
  /* for(i=0;i<len-cutoff;i++){ */
  /*   U_lc_surf[i] = U_surface(lc_bins[i],U, sigma); */
  /* } */
  
}
















































































/* void normed_lc_damped(double ** x_ar, double * t_ar, double * lambda_ar, */
/* 		      int traces, int samples, double K, double beta, int recon_samples, */
/* 		      double z_start, double z_final, double * G, */
/* 		      double epsilon, double sigma, double * U_lc_surf, */
/* 		      double * recon_slope, double * recon_msd, double * work_mean, double * work_std, */
/* 		      double * exp2_avarage, double * avarage_work, double *rms, int * invalid_count, */
/* 		      double * bin_centers, double * samples_work_variances, */
/* 		      double *lc_bins, double * G_lc, int cutoff, int ITER_MAX, double * prot, int steps){ */
/*   int i,j,k,l,I,J,m; */
/*   double P_j_temp; */
/*   I = recon_samples; */
/*   J = recon_samples-cutoff; */
/*   double * signal = (double *) malloc(sizeof(double) * I); */
/*   double * P = (double *) malloc(sizeof(double) * J); */
/*   double * c = (double *) malloc(sizeof(double) * I); */
/*   double * U = (double *) malloc(sizeof(double) * I); */
/*   double * U_tilde = (double *) malloc(sizeof(double) * I); */

  
/*   get_signal(signal,recon_samples,z_start,z_final, */
/* 	     traces,samples,lambda_ar,K, */
/* 	     x_ar,1, bin_centers,beta); /\* meaning normalized. *\/ */
/*   for(l=0;l<J;l++){ */
/*     lc_bins[l] = bin_centers[l+cutoff/2]; */
/*     P[l] = 1.0/(J);		/\* initial guess. *\/ */
/*   } */
/*   double ** H = (double **) malloc (sizeof(double *) * (I)); */
/*   for(i=0;i<I;i++){H[i] = (double *) malloc(sizeof(double)*J);} */

/*   fill_transfer_matrix(H, I, J,beta, K, */
/* 		       bin_centers, lc_bins); */
  
/*   // Lets begin. If you want to see the evolution you should let P be a two dimensinal array of dimensions (ITER_MAX,J) */
/*   for(m=0;m<ITER_MAX;m++){   */
/*     // first we need to fill the c_i */
/*     for(i=0;i<I;i++){ */
/*       c[i] = 0; */
/*       for(j=0;j<J;j++){ */
/* 	c[i] += H[i][j]*P[j];  */
/*       } */
/*       U[i] = -2*pow(T_param,-2)*(signal[i]*log(c[i]/signal[i]) - c[i] + signal[i]); */
/*       U_tilde[i] = min_f(U[i],1); */
/*     } */
/*     // c_i filled */
/*     // updating P[j] */
/*     for(j=0;j<J;j++){ */
/*       P_j_temp = P[j]; */
/*       for(i=0;i<I;i++){ */
/* 	P[j] += P_j_temp*(signal[i]/c[i])*H[i][j]; */
/*       } */
/*     } */
    
/*   } */

/*   /\* for(j=0;j<J;j++){ *\/ */
/*   /\*   G_lc[j] = -pow(beta,-1)*log(P[j]); *\/ */
/*   /\* } *\/ */
/*     //printf("P[%d] = %E || G[%d] = %E \n", j, P[j], j, G[j]); */
  
/*   // OFFSETTING */
/*   double g_const = U_surface(lc_bins[0],epsilon, sigma) - (-pow(beta,-1)*log(P[0])); */
/*   printf("g_const = %E \n", g_const); */
  
/*   /\* // when all is done. Update vector *\/ */
/*   for(l=0;l<J;l++){ */
/*     U_lc_surf[l] = U_surface(lc_bins[l],epsilon, sigma); */
/*     //U_lc_surf[l] = 0; */
/*     if(P[l] > 0){ */
/*       G_lc[l] = -pow(beta,-1)*log(P[l]); */
/*     } */
/*     G_lc[l] += g_const; */
/*   } */
/* } */
