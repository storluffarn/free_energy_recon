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
//#define PI 4*(atan(1))
//#define k_b (1.38064852*pow(10, -23))

// different protocols, chose one.

double h_step(double t){
  if(t>0){return 1;}else{return 0;}
  //return (double) (t>0);
};

// ASSYMETRIC TENT PROTOCOL
/* void calculate_lambda(double * lambda, double t, double T_final, */
/* 		      double z_final, double A, double-A)*(z_final - ((z_final)/(T_final - A))*(t-A)   ); */
/* }; */


// SINGLE WELL OSCILLATION PROTOCOL
/* void calculate_lambda(double * lambda, double t, double T_final, */
/* 		      double z_final, double A, double B){ */
/*   *lambda = (B*z_final)*sin(2*M_PI*t*A/T_final); */
/* }; */


// LINEAR RAMP PROTOCOL
void calculate_lambda(double * lambda, double t, double T_final,
		      double z_final){
  *lambda = (z_final/T_final)*t;
};
// double A, B;

// SANG PROTOCOL (PLATEAU PROTOCOL)
/* void calculate_lambda(double * lambda, double t, double T_final, */
/* 		      double z_final, double A, double B){ */
/*   if(t<B*T_final){ */
/*     *lambda = (A/B)*(z_final/T_final)*t; /\* initial slope *\/ */
/*   }else if(t>T_final*(1-B)){ */
/*     *lambda = A*z_final + (A/B)*(z_final/T_final)*(t-(1-B)*T_final); */
/*   }else{ */
/*     *lambda = A*z_final; */
/*   } */
/* }; */

// HARMONIC PROTOCOL
/* void calculate_lambda(double * lambda, double t, double T_final, */
/* 		      double z_final, double A, double B){ */
/*   *lambda = 0.5*z_final*(1 - cos(A*M_PI*t/T_final + B*M_PI)); */
/* }; */

// LOGARITHMIC PROTOCOL
/* void calculate_lambda(double * lambda, double t, double T_final, */
/* 		      double z_final, double A, double B){ */
/*   *lambda = z_final*log(1 + ((exp(1) - 1)/T_final)*t ); */
/* }; */

// TENT PROTOCOL
/* void calculate_lambda(double * lambda, double t, double T_final, */
/* 		      double z_final, double A, double B){ */
/*   double T_per = T_final/A; */
/*   *lambda = (z_final/T_per)*(T_per  - fabs(fmod(t, 2*T_per) - T_per)); */
/* }; */

// LINE-SINE PROTOCOL
/* void calculate_lambda(double * lambda, double t, double T_final, */
/* 		      double z_final, double A, double B){ */
/*   *lambda = (z_final/T_final)*t + B*z_final*sin(2*M_PI*t*A/T_final); */
/* }; */

// STEP PROTOCOL
/* void calculate_lambda(double * lambda, double t, double T_final, */
/* 		      double z_final, double A, double B){ */

/*   /\* printf("Hello? \n"); *\/ */
/*   if(t<(0.5*B*T_final)){ */
/*     *lambda = ((A*z_final)/(B*T_final))*t; */
/*     //printf("Triggered first \n"); */
/*   } */
/*   else if(t >= 0.5*B*T_final && t < (1-0.5*B)*T_final){ */
/*     *lambda = 0.5*A*z_final + (((1-A)*z_final)/((1-B)*T_final))*(t - 0.5*B*T_final); */
/*     //printf("Triggered mid \n"); */
/*   }else{ */
/*     *lambda = (1-0.5*A)*z_final + ((A*z_final)/(B*T_final))*(t - (1-0.5*B)*T_final); */
/*     //printf("Triggered last \n"); */
/*   } */
/* }; */


double force(double x, double lambda, double beta, double D,
	     double epsilon, double sigma, double K){
  
  double f = -beta*D*(epsilon*(2*M_PI/sigma)*sin(2*M_PI*x/sigma)
		  +K*( x - lambda));

  //printf("force is: %E \n", f);

  return f;
}
// double well force...
/* double force(double x, double lambda, double beta, double D, */
/* 	     double epsilon, double sigma, double K){ */

  /* return -beta*D*(epsilon*64*pow(sigma,-4)*x*(x- 0.5*sigma)*(x-sigma) */
/* 		  + K*(x - lambda) */
/* 		  ); */
/* } */

/* double force(double x, double lambda, double beta, double D, */
/* 	     double epsilon, double sigma, double K){ */

/*   return -beta*D*(epsilon*64*pow(sigma,-4)*(pow(x,3) - 1.5*sigma*pow(x,2) + 0.5*pow(sigma,2)*x) */
/* 		  + K*(x - lambda) */
/* 		  ); */
/* } */
double ermak_force(double x, double lambda,
		   double epsilon, double sigma, double K){
  
  return -K*(x-lambda) -epsilon*sin(2*M_PI*x/sigma)*(2*M_PI/sigma);
  
}
// double beta, D

void ricci(double * x, double * p, double lambda1, double lambda2,
	   double K, double W, double epsilon, double sigma,
	   double t_m, double t2_2m, 
	   double e1, double e2, double e4, double e34,
	   double W_amp3, double x_rand_coeff, double h_2){
  
  /*  LOTS of speed issues to consider here later. */
  double x_temp;
  double p_rand = W*W_amp3;
  double x_rand = x_rand_coeff*p_rand; /* this is the smart thing about the algorithm, you only need 1 random */

  x_temp = *x + t_m*e2*(*p) + t2_2m*ermak_force(*x, lambda1, epsilon, sigma, K) + x_rand;
  /* now x(t+dt) is stored in x_temp */
  *p = e1*(*p) + h_2*e4*ermak_force(x_temp, lambda2, epsilon, sigma, K) + h_2*ermak_force(*x, lambda1, epsilon, sigma, K)*e34 + p_rand;
  *x = x_temp;
  /* Hopefully this should actually be it. */

}
// double W_amp, h, M, eta, beta, D

void ermak(double * x, double * p, double lambda1, double lambda2, double x_rand, double p_rand,
	   double eta_coeff, double eta_coeff2, double t2_2m, double exp_eta,
	   double sigma, double epsilon, double K, double m_i){
  double x_temp;
  x_temp = *x + (*p)*eta_coeff*m_i + t2_2m*ermak_force(*x, lambda1, epsilon, sigma, K) + x_rand;
  *p = exp_eta*(*p) + eta_coeff2*(ermak_force(*x, lambda1, epsilon, sigma, K)+ermak_force(x_temp, lambda2, epsilon, sigma, K) + p_rand);
  *x = x_temp; 
}
// double beta, D

double U_surface(double x, double epsilon, double sigma){
  return epsilon*(1 - cos(2*M_PI*x/sigma));
  // whatch out for the factor of 2, changes def of epsilon.
}

// Safety Runge-Kutta
double rk_force(double x, double lambda, double epsilon, double K,
		double sigma, double beta, double rand, double D){
  
  return beta*D*(rand +
		 -epsilon*0.5*(2*M_PI/sigma)*sin(2*M_PI*x/sigma)
		 -K*( x - lambda )
		 );
    };
void RK4(double * x, double t, double * rk,
	 double epsilon, double K,
	 double sigma, double beta, double rand,
	 double D, double h, double * lambda,
	 double T_final, double z_final){
  // this means to calls to calculate lambda because we have to calculate
  // lambda(t + 0.5h) and lambda(t + h) as midsteps.
  // we could actually just keep the last lambda value and not calculate it
  // again because the last call sets this here. A word of warning however:
  // this is not the case for other solvers, like heun for instance.
  rk[0] = rk_force(*x, *lambda, epsilon,K, sigma,
		   beta, rand,D );  
  calculate_lambda(lambda,t + 0.5*h ,T_final,z_final);
  rk[1] = rk_force(*x + 0.5*h*rk[0], *lambda, epsilon,K, sigma,
		   beta, rand,D );
  rk[2] = rk_force(*x + 0.5*h*rk[1], *lambda, epsilon,K, sigma,
		   beta, rand,D );
  calculate_lambda(lambda,t + h ,T_final,z_final);
  rk[3] = rk_force(*x + h*rk[2], *lambda, epsilon,K, sigma,
		   beta, rand,D );
  *x = *x + h*(1./6)*(rk[0] + 2*rk[1] + 2*rk[2] + rk[3]); 
};
//double h_i

double heun(double x, double lambda, double K, double W,
	    double W_amp, double h, double epsilon,
	    double sigma, double beta, double D){
  // This simplifies alot sense diffusion coefficient is constant.
  return x + force(x, lambda, beta, D, epsilon, sigma, K)*h + W_amp*W;
}

void write_2darray(double ** x_ar, int len1, int len2, char *x_file){
  int i,k;
  FILE * outfile;
  outfile = fopen(x_file, "w+");
  for(i=0;i<len1;i++){
    for(k=0; k<len2; k++){
      fprintf(outfile, "%E \t", x_ar[k][i]);
    }
    fprintf(outfile,"\n");
  }
  fclose(outfile);
}

void write_1darray(double * ar, int len, char *filename){
  int i;
  FILE *outfile;
  outfile = fopen(filename, "w+");
  for(i=0;i<len;i++){
    fprintf(outfile, "%E \n", ar[i]);
  }
  fclose(outfile);
}

double RMS(double * ar1, double * ar2, int len){
  double rms = 0;
  int i;
  for(i=0; i < len; i++){
    rms += pow(ar2[i] - ar1[i],2);
  }
  return sqrt(rms/len);
}
void write_paramfile(char * filename, double z_final, double sigma, double D, double T,
		     double T_final, unsigned long int steps, unsigned int traces,
		     int recon_samples, double z_start, char * information_string, 
		     double epsilon, double h, double recon_slope, double recon_msd,
		     double A, double B, double work_mean, double work_std, double z_min, double z_max){
  FILE *outfile;
  outfile = fopen(filename, "w+");
  fprintf(outfile, "General info: \n  %s \n \n", information_string);
  fprintf(outfile, "T_final = %E \n", T_final);
  fprintf(outfile, "z_final [sigma] = %E \n", z_final/sigma);
  fprintf(outfile, "z_start [sigma]= %E \n", z_start/sigma); /* should usually be zero */
  fprintf(outfile, "v_eff =z_final/T_final = %E \n", z_final/T_final);
  fprintf(outfile, "T = %E \n", T);
  fprintf(outfile, "sigma = %E \n", sigma);
  fprintf(outfile, "epsilon [k_BT] = %E \n", epsilon/(k_b*T));
  fprintf(outfile, "D = %E \n", D);
  fprintf(outfile, "steps = %lu \n", steps);
  fprintf(outfile, "traces = %d \n", traces);
  fprintf(outfile, "h = %E \n", h);
  fprintf(outfile, "recon_samples = %d \n", recon_samples);
  fprintf(outfile, "recon_slope[k_bT/sigma] = %E \n", recon_slope);
  fprintf(outfile, "recon chi2 [(k_b*T)^2] = %E \n", recon_msd);
  fprintf(outfile, "protocol param A = %E \n", A);
  fprintf(outfile, "protocol param B = %E \n", B);
  fprintf(outfile, "work_mean[k_b*T] = %E \n", work_mean);
  fprintf(outfile, "work_std[k_b*T]  = %E \n", work_std);
  fprintf(outfile, "z_min[sigma]  =    %E \n", z_min/sigma);
  fprintf(outfile, "z_max[sigma]  = %E \n", z_max/sigma);
  fprintf(outfile, "scan_length[sigma]  = %E \n", (z_max-z_min)/sigma);

  fclose(outfile);
}

int get_args(int argc, char *argv[],
	     long unsigned int *steps,  char * outfile,
	     double * epsilon_d, unsigned int * traces,
	     double * z_final_d, int * recon_samples, double * T_final,
	     char * info_file, char * outdir, double * A, double * B, int * p, double * protocol_constant,
	     char * paste_file, int * R){

  // Suppliment this function every time you need more commandline args!
  // Bad practice i know ... but multiarg functions are a hassle!
  // basically just add a switch and add [argname]: to the argstring.
  //int index;
  int c;
  opterr = 0;
  while( (c = getopt(argc, argv, "M:N:O:E:Z:r:t:i:d:A:B:p:C:P:R:")) != -1 ){
    switch (c){
    case 'M':
      *traces = (int) atof( optarg  ); // optarg should be string im pretty sure
      break;
    case 'r':
      *recon_samples = (int) atof( optarg  ); // optarg should be string im pretty sure
      break;
    case 'i':
      sprintf(info_file, "%s", optarg);
      break;
    case 'R':
      *R = (int) atof( optarg );
      break;
    case 'd':
      sprintf(outdir, "%s", optarg);
      break;
    case 'O':
      sprintf(outfile, "%s", optarg);
      break;
    case 'P':
      sprintf(paste_file, "%s", optarg);
      break;
    case 'Z':
      *z_final_d = (double) atof (optarg);
      break;
    case 't':
      *T_final = (double) atof (optarg);
      break;
    case 'N':
      *steps = *steps * ( (int) atof (optarg));
      break;
    case 'p':
      *p =  ((int) atof (optarg));
      break;
    case 'E':
      *epsilon_d = (double) atof (optarg);
      break;
    case 'A':
      *A = (double) atof (optarg);
      break;
    case 'B':
      *B = (double) atof (optarg);
      break;
    case 'C':
      *protocol_constant = (double) atof (optarg);
      break;
    case '?':
	if(optopt == 'c'){
	  fprintf(stderr, "Opt -%c requires arg \n", optopt);
	}else if(isprint (optopt)){
	  fprintf(stderr, "Unknown option -%c ... \n", optopt);
	}else{
	  fprintf(stderr, "Unknown opt char '\\x%x' \n", optopt);
	}
      return 1;
    default:
      abort();
    }
  }
  return 1;
};

void reconstruction_accuracy(double * G, double * U_surf, double * bin_centers , double beta, int recon_samples, double * recon_slope, double * recon_msd,
			     double sigma){
  int l;
  double * residuals = (double *) malloc (sizeof (double)*recon_samples);
  for(l=0;l<recon_samples;l++){
    residuals[l] = (G[l] - U_surf[l]); /* defined so that slope normally is positive. */
  }
  
  double c0,c1, cov00, cov01, cov11 , chisq;
  gsl_fit_linear(bin_centers, 1, residuals, 1, recon_samples, &c0, &c1, &cov00,&cov01,&cov11, &chisq);
  *recon_slope = c1*beta*sigma;
  *recon_msd = chisq*(pow(beta,2));
  
}

void get_work_dist(double ** x_ar, double * lambda_ar,
		   int traces, int samples, double K, double * work_vector){
  int i,k;
  double ** ext_work = (double **) malloc (sizeof(double *) * traces);

  for(k=0;k<traces;k++){
    ext_work[k] = (double *) malloc (sizeof(double) * samples);
    ext_work[k][0] = 0;
    for(i=1;i<samples;i++){
      ext_work[k][i] = ext_work[k][i-1] - 0.5*K*((x_ar[k][i] - lambda_ar[i])
						 +(x_ar[k][i-1] - lambda_ar[i-1])
						 )*(lambda_ar[i]-lambda_ar[i-1]);
    }
    work_vector[k] = ext_work[k][samples-1];
  }
}
// double  t_ar, beta, signma; int j, l

void get_third_cumulant(double ** ext_work, int samples, int traces,
			double * weight_work, double beta){
  int i,k;
  /* double * mean_work = (double *) malloc (sizeof(double ) * samples); */
  /* double * var_work = (double *) malloc (sizeof(double ) * samples); */

  double mom1,mom2,mom3;
  
  for(i=0;i<samples;i++){
    mom1 = 0; mom2=0; mom3=0;
    
    // collecting mean
    for(k=0;k<traces;k++){
      mom1 += ext_work[k][i]; /* make sure order of indices is correct. */
      mom2 += pow(ext_work[k][i],2);
      mom3 += pow(ext_work[k][i],3);
    }

    weight_work[i] = mom1 - 0.5*beta*(mom2 - pow(mom1,2)) + 0.25*pow(beta,2)*(2*pow(mom1,3) -3*mom1*mom2 + mom3);
  }
}
// double control; int j

void berkovic_works(double ** ext_work, double ** x_ar, double *lambda_ar, int samples, int traces,
		    double slip_distance, double K){
  int i,k;
  int slip_control=0;
  for(k=0;k<traces;k++){
    ext_work[k] = (double *) malloc (sizeof(double) * samples);
    ext_work[k][0] = 0;
    for(i=1;i<samples;i++){
      if(slip_control==0){
      ext_work[k][i] = ext_work[k][i-1] - 0.5*K*((x_ar[k][i] - lambda_ar[i])
						 +(x_ar[k][i-1] - lambda_ar[i-1])
						 )*(lambda_ar[i]-lambda_ar[i-1]);
      }else{
      ext_work[k][i] = 0.5*K*((x_ar[k][i] - lambda_ar[i])
						 +(x_ar[k][i-1] - lambda_ar[i-1])
						 )*(lambda_ar[i]-lambda_ar[i-1]);
      slip_control=0;
      }
      if(sqrt(pow(x_ar[k][i-1]-x_ar[k][i],2)) > slip_distance){
	//slip_control==1;
	printf("Slip at (%d, %d)! Reseting works to zero. \n", i ,k);
      }
    }
  }
}

void get_weight_work(double ** ext_work, int samples, int traces,
		     double * weight_work, double beta){
  // control 0 means normal exponentiated avarages, True valued variable gives second order cumulant
  // makes sure the sample variance is the unbiased one.
  int i,k;
  double * mean_work = (double *) malloc (sizeof(double ) * samples);
  double * var_work = (double *) malloc (sizeof(double ) * samples);
  
  for(i=0;i<samples;i++){
    mean_work[i] = 0; // init
    var_work[i] = 0; // init
    weight_work[i] = 0;
    // collecting mean
    for(k=0;k<traces;k++){
      mean_work[i] += ext_work[k][i]; /* make sure order of indices is correct. */
    }
    // using mean to collect unbiased sample var
    for(k=0;k<traces;k++){
      var_work[i] += pow(ext_work[k][i] - mean_work[i],2);
    }

    weight_work[i] = mean_work[i] - 0.5*beta*var_work[i]/(traces-1);
  } 
}
//double control; int j;

void make_protocol(int reverse, double * prot, double T_final, double z_final, double * lambda, double h,
		   int steps){
  int i;
  //k=0;
  /* double z_max = protocol_constant; */
  /* double z_min = protocol_constant; */
  if(reverse==0){
    for(i=0;i<steps;i++){
      calculate_lambda(lambda, i*h, T_final, z_final);
      prot[i] = *lambda;
      //printf("prot[%d] = %E \n", i, prot[i]);
    }    
  }else{
    for(i=0;i<steps;i++){
      calculate_lambda(lambda, h*(steps-i), T_final, z_final);
      prot[i] = *lambda;
    }    
  }
}
// int samples

void wham_bias2(int len, int bias_len ,double z_start, double z_final,
		int traces, int samples, double ** lambda_ar, double K,
		double ** x_ar, double * bin_centers, double beta,
		double z_f, double * G, int * invalid_count,
		double epsilon, double sigma, double * U_surf){
  // Now! we have to create a vector of the biases (that is, the differences between the support and the tip.)
  // you might want to also put this into potential for easier histogramming.
  printf("RECONSTRUCTING USING WHAM_BIAS2, MINH PRESCRIPTION \n");
  
  int i,k,l,m;
  //int j;
  char temp_char2[256];
  //double s;
  double * bin_edges = (double *) malloc (sizeof(double) * len+1);
  //double * bin_centers = (double *) malloc (sizeof(double) * len);
  double ** ext_work = (double **) malloc ( sizeof(double * )* traces );
  //double * exp_avarage = (double *) malloc (sizeof(double) * samples);

  double d_bin = (z_final - z_start)/len;
  //double d_bin_bias = (z_f - z_0)/bias_len;
  
  double * bias_works = (double *) malloc (sizeof(double) * bias_len);
  double ** hist = (double **) malloc (sizeof(double) * bias_len);
  double ** hist_counter = (double **) malloc (sizeof(double) * bias_len);
  //double ** u_bias = (double *) malloc (sizeof(double) * bias_len);

  
  double * bias_centers = (double *) malloc (sizeof(double) * bias_len);
  double * bias_edges = (double *) malloc (sizeof(double) * (bias_len+1));
  double * bias_counter = (double *) malloc (sizeof(double) * bias_len);
  double ** bias_array = (double **) malloc (sizeof(double *) * traces);
  
  // the bin edges we know, but the bias bins we do not yet know. We can only really know the length of them
  // we can find out the maximum bias. But its reasonable to assume its not going to be too big.
  // lets say the maximum bias is v(scan_lengt)
  
  for(i=0;i<len+1;i++){
    bin_edges[i] = z_start + i*d_bin;
  }
  for(i=0;i<len;i++){
    bin_centers[i] = (bin_edges[i+1] + bin_edges[i])*0.5;  
  }
  // collecting works, collect bias too
  double max_bias = 0;
  double max_displacement = 0;
  for(k=0;k<traces;k++){
    ext_work[k] = (double *) malloc (sizeof(double) * samples);
    ext_work[k][0] = 0;
    bias_array[k] = (double *) malloc (sizeof(double) * samples);
    bias_array[k][0] = 0.5*K*pow(x_ar[k][0] - lambda_ar[0][k],2);
    
    for(i=1;i<samples;i++){
      ext_work[k][i] = ext_work[k][i-1] - 0.5*K*((x_ar[k][i] - lambda_ar[i][k])
						 +(x_ar[k][i-1] - lambda_ar[i-1][k])
						 )*(lambda_ar[i][k]-lambda_ar[i-1][k]);
      bias_array[k][i] = 0.5*K*pow(x_ar[k][i] - lambda_ar[i][k],2);
      if(bias_array[k][i] > max_bias){max_bias = bias_array[k][i]; }
      if(fabs(x_ar[k][i]-lambda_ar[i][k]) > max_displacement){max_displacement = fabs(x_ar[k][i]-lambda_ar[i][k]);}
    }
  }
  
  // We need to set up the biasing bins.
  printf("max bias was %E which compares to 0.5K(z_f)** as %E \n", max_bias, max_bias/(0.5*K*pow(z_f,2)));
  printf("max displacement: %E \n", max_displacement/sigma);
  //double d_bias = max_bias/bias_len;
  double d_displacement = max_displacement/bias_len;
  for(i=0;i<bias_len+1;i++){
    //bias_edges[i] = i*d_bias;
    bias_edges[i] = 0.5*K*pow(i*d_displacement,2); /* nonlinear binning. */
    
  }
  for(i=0;i<bias_len;i++){
    bias_centers[i] = (bias_edges[i+1] + bias_edges[i])*0.5;
    hist[i] = (double *) malloc (sizeof(double) * len);
    hist_counter[i] = (double *) malloc (sizeof(double) * len);
  }
  // INITIALIZE CONTAINERS
  for(l=0;l<bias_len;l++){
    bias_works[l] = 0; bias_counter[l] = 0;
    for(m=0;m<len;m++){
      hist[l][m] = 0; hist_counter[l][m] = 0;
    }
  }
  
  // BIG LOOP,REALLY NASTY, COULD PROBABLY BE DONE MUCH FASTER IF THIS IS TOO SLOW
  for(i=0;i<samples;i++){
    for(k=0;k<traces;k++){
      for(l=0;l<bias_len;l++){
	
	if(bias_edges[l] < bias_array[k][i] && bias_edges[l+1] > bias_array[k][i]){
	  // Ok we have landed at some particular bias. That means we increment the work and the counter.
	  bias_works[l] += exp(-beta*ext_work[k][i]);
	  bias_counter[l] += 1;

	  for(m=0;m<len;m++){
	    if(bin_edges[m] < x_ar[k][i] && bin_edges[m+1] > x_ar[k][i]){
	      hist[l][m] += exp(-beta*ext_work[k][i]);
	      hist_counter[l][m] += 1;
	    }
	    
	  }
	}
      }
    }
  }
  //printf("DEBUG \n");
  for(l=0;l<bias_len;l++){
    //printf("bias_works[%d] = %E \n", l, bias_works[l]);
    if(bias_counter[l] != 0){
      bias_works[l] = bias_works[l]/bias_counter[l];
      //printf("bias_works[%d] = %E \n", l, bias_works[l]);
    }
    for(m=0;m<len;m++){
      if(hist_counter[l][m] != 0){
	hist[l][m] = hist[l][m]/ hist_counter[l][m];
      }
    }
  }
  // now you have to loop over the biases and see if they fall in some particular bin.
  double top = 0; double bottom = 0;
  double g_const;
  int count = 0;
  int count2 = 0;
  
  for(m=0;m<len;m++){
    top = 0; bottom=0;
    for(l=0;l<bias_len;l++){
      if(bias_works[l] != 0){
	top += hist[l][m]/bias_works[l];
	bottom += exp(-beta*bias_centers[l])/bias_works[l];
      }
    }
    if(top==0){
      printf("Encountered invalid value in reconstruction, setting to zero \n");
      G[m] = 0;			/* Putting to zero all those values which are undefined basically. */
      count++;
    }
    else{
      //printf("exp(G[%d]) ~ %E \n", m,  top/bottom);
      G[m] = -pow(beta,-1)*log(top/bottom);
    }
  }
  double * G_valid = (double *) malloc (sizeof(double) * (len-count));
  double * bin_valid = (double *) malloc (sizeof(double) * (len-count));
  //double * U_surf_valid = (double *) malloc (sizeof(double) * (len-count));
  *invalid_count = count;

  
  // putting valid things into new arrays.
  for(l=0;l<len;l++){
    //printf("G[l] = %E \n", G[l]);
    if(G[l] != 0){
      //printf("Valid point,  G[l] = %E \n", G[l]);
      G_valid[count2] = G[l];
      bin_valid[count2] = bin_centers[l];
      count2++;
    }
  }

  // forgot to pipe to python
  printf("THERE WAS %d INVALID COUNTS, INTERPOLATING. \n", count);
  write_1darray(G_valid, len-count, "pytemp/g_valid.dat");
  write_1darray(bin_valid, len-count, "pytemp/z_valid.dat");
  write_1darray(bin_centers, len, "pytemp/z_full.dat");

  //system("python interpolate_wham.py"); // Cubic-Splines.
  //system("python interpolate_wham2.py");  // Linear-interp.
  system("python interpolate_wham_lin.py");
  
  FILE * outfile;
  count=0;
  outfile = fopen("pytemp/g_interp.dat", "r");
  while(fgets(temp_char2, sizeof(temp_char2), outfile)){
    G[count] = atof(temp_char2);
    count++;
  }

  g_const = U_surface(bin_centers[0], epsilon,sigma) -G[0];
  
  // Setting proper gauge
  count2 = 0;
  for(l=0;l<len;l++){
    G[l] = G[l] + g_const;
    U_surf[l] = U_surface(bin_centers[l], epsilon,sigma);    
  } 
}
//double z_0;
// int c, k

void wham_bias(int len, int bias_len ,double z_start, double z_final,
	       int traces, int samples, double ** lambda_ar, double K,
	       double ** x_ar, double * bin_centers, double beta,
	       double z_0, double z_f, double * G, int * invalid_count,
	       double epsilon, double sigma, double * U_surf){
  
  // There is a big difference between the sets of z_start/finish etc.
  // the second pair sets the bias bin, and the first sets the recon bins.
  int i,k,l,m;
  //int j
  char temp_char2[256];
  //double s;
  double * bin_edges = (double *) malloc (sizeof(double) * len+1);
  //double * bin_centers = (double *) malloc (sizeof(double) * len);
  double ** ext_work = (double **) malloc ( sizeof(double * )* traces );
  //double * exp_avarage = (double *) malloc (sizeof(double) * samples);

  double d_bin = (z_final - z_start)/len;
  double d_bin_bias = (z_f - z_0)/bias_len;
  
  double * bias_works = (double *) malloc (sizeof(double) * bias_len);
  double ** hist =      (double **) malloc (sizeof(double) * bias_len);
  double ** hist_counter = (double **) malloc (sizeof(double) * bias_len);
  double ** u_bias = (double **) malloc (sizeof(double) * bias_len);
  
  
  double * bias_centers = (double *) malloc (sizeof(double) * bias_len);
  double * bias_edges = (double *) malloc (sizeof(double) * (bias_len+1));
  double * bias_counter = (double *) malloc (sizeof(double) * bias_len);

  
  
  // the recon bins should generally be larger than the bias bins probably.

  
  // bins
  
  for(i=0;i<bias_len+1;i++){
    bias_edges[i] = z_0 + i*d_bin_bias;
  }
  for(i=0;i<bias_len;i++){
    bias_centers[i] = (bias_edges[i+1] + bias_edges[i])*0.5;
    hist[i] = (double *) malloc(sizeof(double) * len);		/* borrowing for instantiation */
    hist_counter[i] = (double *) malloc(sizeof(double) * len);
    u_bias[i] = (double *) malloc(sizeof(double) * len);
    for(l=0;l<len;l++){
      hist[i][l] = 0;		/* instantiation */
      hist_counter[i][l] = 0;
    }
  }

  
  for(i=0;i<len+1;i++){
    bin_edges[i] = z_start + i*d_bin;
  }
  for(i=0;i<len;i++){
    bin_centers[i] = (bin_edges[i+1] + bin_edges[i])*0.5;  
  }
  
  
  
  // lambda is now trace dependent. We must gather all appropriate lambdas, not just one
  for(k=0;k<traces;k++){
    ext_work[k] = (double *) malloc (sizeof(double) * samples);
    ext_work[k][0] = 0;
  
    for(i=1;i<samples;i++){
      ext_work[k][i] = ext_work[k][i-1] - 0.5*K*((x_ar[k][i] - lambda_ar[i][k])
						 +(x_ar[k][i-1] - lambda_ar[i-1][k])
						 )*(lambda_ar[i][k]-lambda_ar[i-1][k]);
    }
  }
  
      
  for(l=0;l<bias_len;l++){
    bias_works[l] = 0;bias_counter[l] =0;
    
    for(i=0;i<samples;i++){
      for(k=0;k<traces;k++){
	// here we have to check if lambda is with its bounds. Then we can increment the work. You should histogram here too.
	if(bias_edges[l] < lambda_ar[i][k] && bias_edges[l+1] > lambda_ar[i][k]){
	  
	  bias_works[l] += exp(-beta*ext_work[k][i]);
	  bias_counter[l] += 1;
	  // now we must do ONE MORE loop, over histogramming index. lets call this m. This does NOT have to be
	  // the same bins as what we histogrammed previously. But its simple obviously. but probably not smart...
	  for(m=0;m<len;m++){
	    if(bin_edges[m] < x_ar[k][i] && bin_edges[m+1] > x_ar[k][i]){
	      hist[l][m] += exp(-beta*ext_work[k][i]);
	      hist_counter[l][m] += 1;
	    }
     	  }
	}
      }
    }
    if(bias_counter[l] != 0){bias_works[l] = bias_works[l]/bias_counter[l];}
  }
  printf("THIS IS FROM WHAM BIAS  %d \n", k);
  write_1darray(bias_works, bias_len, "bias_works_testing.dat");
  for(l=0;l<bias_len;l++){
    for(m=0;m<len;m++){
      u_bias[l][m] = 0.5*K*pow(bin_centers[m] - bias_centers[l],2);
      if(hist_counter[l][m] != 0){
	hist[l][m] = hist[l][m]/hist_counter[l][m];
      }
    }
  }
  double top = 0; double bottom = 0;
  double g_const;
  int count = 0;
  int count2 = 0;
  
  for(m=0;m<len;m++){
    top = 0; bottom=0;
    for(l=0;l<bias_len;l++){
      top += hist[l][m]/bias_works[l];
      bottom += exp(-beta*u_bias[l][m])/bias_works[l];
    }
    if(top==0){
      printf("Encountered invalid value in reconstruction, setting to zero \n");
      G[l] = 0;			/* Putting to zero all those values which are undefined basically. */
      count++;
    }
    else{
      //printf("tep = %E \n", top);
      //printf("bottom = %E \n", bottom);
      G[m] = -pow(beta,-1)*log(top/bottom);
      //printf("G[%d] = %E \n", l, G[m]);
    }
  }
  
  double * G_valid = (double *) malloc (sizeof(double) * (len-count));
  double * bin_valid = (double *) malloc (sizeof(double) * (len-count));
  //double * U_surf_valid = (double *) malloc (sizeof(double) * (len-count));
  *invalid_count = count;

  
  // putting valid things into new arrays.
  for(l=0;l<len;l++){
    printf("G[l] = %E \n", G[l]);
    if(G[l] != 0){
      printf("Valid point,  G[l] = %E \n", G[l]);
      G_valid[count2] = G[l];
      bin_valid[count2] = bin_centers[l];
      count2++;
    }
  }

  // forgot to pipe to python
  printf("THERE WAS %d INVALID COUNTS, INTERPOLATING. \n", count);
  write_1darray(G_valid, len-count, "pytemp/g_valid.dat");
  write_1darray(bin_valid, len-count, "pytemp/z_valid.dat");
  write_1darray(bin_centers, len, "pytemp/z_full.dat");
  system("python interpolate_wham.py");
  
  FILE * outfile;
  count=0;
  outfile = fopen("pytemp/g_interp.dat", "r");
  while(fgets(temp_char2, sizeof(temp_char2), outfile)){
    G[count] = atof(temp_char2);
    count++;
  }

  g_const = U_surface(bin_centers[0], epsilon,sigma) -G[0];
  
  // Setting proper gauge
  count2 = 0;
  for(l=0;l<len;l++){
    G[l] = G[l] + g_const;
    U_surf[l] = U_surface(bin_centers[l], epsilon,sigma);    
  }

  
}
// int c

void wham_reconstruct_interp(double ** x_ar, double * lambda_ar,
			     int traces, int samples, double K, double beta, int recon_samples,
			     double z_start, double z_final, double * G,
			     double epsilon, double sigma, double * U_surf,
			     double * recon_slope, double * recon_msd, 
				 double * work_mean, double * work_std,
			     double * exp2_avarage, double * avarage_work, 
				 double *rms, int * invalid_count,
			     double * bin_centers, double * samples_work_variances){
  int i,k,l;
  char temp_char2[256];
  *work_mean = 0; *work_std=0;
  double ** ext_work = (double **) malloc ( sizeof(double * )* traces );
  double * exp_avarage = (double *) malloc (sizeof(double) * samples);
  //double * bin_centers = (double *) malloc (sizeof (double) * recon_samples);
  double * bin_edges =  (double *) malloc (sizeof (double) * (recon_samples + 1) );
  double ** u_bias = (double **) malloc (sizeof (double) * recon_samples);
  for(i=0;i<recon_samples;i++){u_bias[i] = (double *) malloc (sizeof (double) * samples);}
  double d_bin = ((z_final - z_start)/recon_samples);
  double ** histogram = (double **) malloc (sizeof (double *) * samples);
  for(i=0;i<samples;i++){histogram[i] = (double *) malloc (sizeof (double) * recon_samples);}
  double w2_temp;
  // I) ext work loop
  
  for(k=0;k<traces;k++){
    ext_work[k] = (double *) malloc (sizeof(double) * samples);
    ext_work[k][0] = 0;
    for(i=1;i<samples;i++){
      ext_work[k][i] = ext_work[k][i-1] - 0.5*K*((x_ar[k][i] - lambda_ar[i])
  						 +(x_ar[k][i-1] - lambda_ar[i-1])
  						 )*(lambda_ar[i]-lambda_ar[i-1]);
    }
    *work_mean += ext_work[k][samples-1]/traces;
  }

  /* berkovic_works(ext_work,x_ar,lambda_ar,samples, traces, */
  /* 		 0.4*sigma, K); */
 
  // II) exp avarage loop
  for(i=0;i<samples;i++){
    exp_avarage[i] = 0;
    exp2_avarage[i] = 0;
    avarage_work[i] = 0;
    w2_temp = 0;
    for(k=0;k<traces;k++){
      exp_avarage[i] +=  exp(-beta*ext_work[k][i]);
      exp2_avarage[i] += exp(-beta*ext_work[k][i]); /* CHANGED! NAME IS MESSED UP NOW */
      avarage_work[i] += ext_work[k][i];
      w2_temp += pow(ext_work[k][i],2);
    }
    
    w2_temp = w2_temp/traces;
    exp_avarage[i] = exp_avarage[i]/traces;
    exp2_avarage[i] = exp2_avarage[i]/traces;
    avarage_work[i] = (avarage_work[i]/traces);
    samples_work_variances[i] = sqrt(  w2_temp - pow(avarage_work[i],2)    )*beta; /* in units of kbt */
    avarage_work[i] = (avarage_work[i])*beta; /* in units of kbt */
  }

  // No good!
  /* get_weight_work(ext_work,samples,traces, */
  /* 		  exp_avarage,1, beta); */

  // No good either!
  /* get_third_cumulant(ext_work,samples,traces, */
  /* 		     exp_avarage,1, beta); */
  
  
  // III) defining bin centers.
  for(i=0;i<recon_samples+1;i++){
    bin_edges[i] = z_start + i*d_bin;
  }
  for(i=0;i<recon_samples;i++){
    bin_centers[i] = (bin_edges[i+1] + bin_edges[i])*0.5;
  }
  // IV) Calculating deflection bias.
  for(i=0;i<samples;i++){
    for(k=0;k<recon_samples;k++){
      u_bias[k][i] = 0.5*K*pow(    bin_centers[k] - lambda_ar[i] , 2 );
    }
  }
  
  // V) histogramming (This is the tough part)
  for(l=0;l<recon_samples;l++){
    for(i=0;i<samples;i++){
      histogram[i][l] = 0;
      for(k=0;k<traces;k++){
	if(x_ar[k][i] > bin_edges[l] && x_ar[k][i] < bin_edges[l+1]){
	  histogram[i][l] += exp(-beta*ext_work[k][i])/traces;
	}
      }
    }
  }
  
  // VI) Compiling everything.
  double top = 0; double bottom = 0;
  double g_const;
  int count = 0;
  int count2 = 0;

  // CURRENT STRATEGY FOR REMOVING INFS.
  // A) check i the value that goes into the log is zero, if it is put that recon point to 0 and increment invalid counter (count) with 1
  // B) Loop over all recon points and put the ones that are valid into a new list.
  // C) loop over list and find a good value to set the zero level of the energy.
  // D) Once this is done, swap the arrays and make sure to only use the recon_samples-count of the array.

  for(l=0;l<recon_samples;l++){
    top = 0; bottom = 0;
    for(i=0;i<samples;i++){
      top += histogram[i][l]/exp_avarage[i];
      bottom += exp(-beta*u_bias[l][i])/exp_avarage[i];
    }
    if(top==0){
      printf("Encountered invalid value in reconstruction, setting to zero \n");
      G[l] = 0;			/* Putting to zero all those values which are undefined basically. */
      count++;
    }
    else{
    G[l] = -pow(beta,-1)*log(top/bottom);
    }
  }
  // at this stage, identical zero means a miss basically. We do not want these values.
  double * G_valid = (double *) malloc (sizeof(double) * (recon_samples-count));
  double * bin_valid = (double *) malloc (sizeof(double) * (recon_samples-count));
  *invalid_count = count;

  // putting valid things into new arrays.
  for(l=0;l<recon_samples;l++){
    if(G[l] != 0){
      G_valid[count2] = G[l];
      bin_valid[count2] = bin_centers[l];
      count2++;
    }
  }

  printf("THERE WAS %d INVALID COUNTS, INTERPOLATING. \n", count);
  // it is at this stage we should probably pipe G_valid, bin_valid etc to some python script
  // and then just use a simper interpolation to find the missing pieces.
  write_1darray(G_valid, recon_samples-count, "pytemp/g_valid.dat");
  write_1darray(bin_valid, recon_samples-count, "pytemp/z_valid.dat");
  write_1darray(bin_centers, recon_samples, "pytemp/z_full.dat");
  //system("python interpolate_wham.py"); // Cubic-splines. 
  //system("python interpolate_wham2.py");  // Linear interpolation.
  system("python interpolate_wham_lin.py");
  
  FILE * outfile;
  count=0;
  outfile = fopen("pytemp/g_interp.dat", "r");
  while(fgets(temp_char2, sizeof(temp_char2), outfile)){
    G[count] = atof(temp_char2);
    count++;
  }
  
  g_const = U_surface(bin_centers[0], epsilon,sigma) -G[0];
  
  // Setting proper gauge
  count2 = 0;
  for(l=0;l<recon_samples;l++){
    G[l] = G[l] + g_const;
    U_surf[l] = U_surface(bin_centers[l], epsilon,sigma);    
  }
  
  /* // swapping, always remembering that the last count of these arrays will be junk now. */
  /* for(l=0; l<recon_samples-count;l++){ */
  /*   G[l] = G_valid[l]; */
  /*   bin_centers[l] = bin_valid[l]; */
  /*   U_surf[l] = U_surf_valid[l]; */
  /* } */

  count=0;

  for(k=0;k<traces;k++){
    *work_std += pow(ext_work[k][samples-1] - *work_mean,2)/traces;
  }
  *work_std = sqrt(*work_std)*beta;
  *work_mean = *work_mean*beta;
  *rms = RMS(G, U_surf, recon_samples-count)*beta;
  reconstruction_accuracy(G,U_surf,bin_centers ,beta,recon_samples-count,recon_slope,recon_msd, sigma);
  printf("work_mean[k_b*T] = %E , work_std[(k_bT)] = %E \n", *work_mean, *work_std);
  printf("recon_slope [k_bT/sigma] = %E , recon_msd[(k_bT)^2] = %E \n",*recon_slope, *recon_msd);
}
// double * t_ar; int j;
//double * U_surf_valid = (double *) malloc (sizeof(double) * (recon_samples-count));

void wham_reconstruct(double ** x_ar, double * lambda_ar,
		      int traces, int samples, double K, double beta, int recon_samples,
		      double z_start, double z_final, double * G,
		      double epsilon, double sigma, double * U_surf,
		      double * recon_slope, double * recon_msd, 
			  double * work_mean, double * work_std,
		      double * exp2_avarage, double * avarage_work, 
			  double *rms, int * invalid_count,
		      double * bin_centers, double * samples_work_variances,
		      double * G_interp, double * U_interp, double * bin_interp,
		      double * bin_full, double * U_full){
  int i,k,l;
  //int j;
  char temp_char2[256];
  *work_mean = 0; *work_std=0;
  double ** ext_work = (double **) malloc ( sizeof(double * )* traces );
  double * exp_avarage = (double *) malloc (sizeof(double) * samples);
  //double * bin_centers = (double *) malloc (sizeof (double) * recon_samples);
  double * bin_edges =  (double *) malloc (sizeof (double) * (recon_samples + 1) );
  double ** u_bias = (double **) malloc (sizeof (double) * recon_samples);
  for(i=0;i<recon_samples;i++){u_bias[i] = (double *) malloc (sizeof (double) * samples);}
  double d_bin = ((z_final - z_start)/recon_samples);
  double ** histogram = (double **) malloc (sizeof (double *) * samples);
  for(i=0;i<samples;i++){histogram[i] = (double *) malloc (sizeof (double) * recon_samples);}
  double w2_temp;
  // I) ext work loop
  
  for(k=0;k<traces;k++){
    ext_work[k] = (double *) malloc (sizeof(double) * samples);
    ext_work[k][0] = 0;
    for(i=1;i<samples;i++){
      ext_work[k][i] = ext_work[k][i-1] - 0.5*K*((x_ar[k][i] - lambda_ar[i])
						 +(x_ar[k][i-1] - lambda_ar[i-1])
						 )*(lambda_ar[i]-lambda_ar[i-1]);
    }
    *work_mean += ext_work[k][samples-1]/traces;
  }
 
  // II) exp avarage loop
  for(i=0;i<samples;i++){
    exp_avarage[i] = 0;
    exp2_avarage[i] = 0;
    avarage_work[i] = 0;
    w2_temp = 0;
    for(k=0;k<traces;k++){
      exp_avarage[i] +=  exp(-beta*ext_work[k][i]);
      exp2_avarage[i] += exp(-beta*ext_work[k][i]); /* CHANGED! NAME IS MESSED UP NOW */
      avarage_work[i] += ext_work[k][i];
      w2_temp += pow(ext_work[k][i],2);
    }
    
    w2_temp = w2_temp/traces;
    exp_avarage[i] = exp_avarage[i]/traces;
    exp2_avarage[i] = exp2_avarage[i]/traces;
    avarage_work[i] = (avarage_work[i]/traces);
    samples_work_variances[i] = sqrt(  w2_temp - pow(avarage_work[i],2)    )*beta; /* in units of kbt */
    avarage_work[i] = (avarage_work[i])*beta; /* in units of kbt */
  }
   
  // III) defining bin centers.
  for(i=0;i<recon_samples+1;i++){
    bin_edges[i] = z_start + i*d_bin;
  }
  for(i=0;i<recon_samples;i++){
    bin_centers[i] = (bin_edges[i+1] + bin_edges[i])*0.5;
  }
  // IV) Calculating deflection bias.
  for(i=0;i<samples;i++){
    for(k=0;k<recon_samples;k++){
      u_bias[k][i] = 0.5*K*pow(    bin_centers[k] - lambda_ar[i] , 2 );
    }
  }
  
  // V) histogramming (This is the tough part)
  for(l=0;l<recon_samples;l++){
    for(i=0;i<samples;i++){
      histogram[i][l] = 0;
      for(k=0;k<traces;k++){
	if(x_ar[k][i] > bin_edges[l] && x_ar[k][i] < bin_edges[l+1]){
	  histogram[i][l] += exp(-beta*ext_work[k][i])/traces;
	}
      }
    }
  }
  
  // VI) Compiling everything.
  double top = 0; double bottom = 0;
  double g_const;
  int count = 0;
  int count2 = 0;

  // CURRENT STRATEGY FOR REMOVING INFS.
  // A) check i the value that goes into the log is zero, if it is put that recon point to 0 and increment invalid counter (count) with 1
  // B) Loop over all recon points and put the ones that are valid into a new list.
  // C) loop over list and find a good value to set the zero level of the energy.
  // D) Once this is done, swap the arrays and make sure to only use the recon_samples-count of the array.

  for(l=0;l<recon_samples;l++){
    top = 0; bottom = 0;
    for(i=0;i<samples;i++){
      top += histogram[i][l]/exp_avarage[i];
      bottom += exp(-beta*u_bias[l][i])/exp_avarage[i];
    }
    if(top==0){
      printf("Encountered invalid value in reconstruction, setting to zero \n");
      G[l] = 0;			/* Putting to zero all those values which are undefined basically. */
      count++;
    }
    else{
    G[l] = -pow(beta,-1)*log(top/bottom);
    }
  }
  // at this stage, identical zero means a miss basically. We do not want these values.
  double * G_valid = (double *) malloc (sizeof(double) * (recon_samples-count));
  double * bin_valid = (double *) malloc (sizeof(double) * (recon_samples-count));
  double * U_surf_valid = (double *) malloc (sizeof(double) * (recon_samples-count));

  double * G_interp_temp = (double *) malloc (sizeof(double) * (count));
  double * bin_interp_temp = (double *) malloc (sizeof(double) * (count));
//  double * U_interp_temp = (double *) malloc (sizeof(double) * (count));


  *invalid_count = count;
  int interp_count=0;
  // putting valid things into new arrays.
  for(l=0;l<recon_samples;l++){
    if(G[l] != 0){
      G_valid[count2] = G[l];
      bin_valid[count2] = bin_centers[l];
      count2++;
    }else{
      bin_interp_temp[interp_count] = bin_centers[l];
      interp_count++;
    }
  }

  /* printf("THERE WAS %d INVALID COUNTS, INTERPOLATING. \n", count); */
  /* // it is at this stage we should probably pipe G_valid, bin_valid etc to some python script */
  /* // and then just use a simper interpolation to find the missing pieces. */
  
  g_const = U_surface(bin_centers[0], epsilon,sigma) -G[0];
  
  // Setting proper gauge
  count2 = 0;
  for(l=0;l<recon_samples-count;l++){
    G_valid[l] = G_valid[l] + g_const;
    U_surf_valid[l] = U_surface(bin_valid[l], epsilon,sigma);    
  }
  if(*invalid_count !=0){
    // this is fine, we can pipe this to some interpolator. BUT ONLY KEEP THE HITS THAT WERE NOT RECONSTRUCTED
    write_1darray(G_valid, recon_samples-count, "pytemp/g_valid.dat");
    write_1darray(bin_valid, recon_samples-count, "pytemp/z_valid.dat");
    write_1darray(bin_interp_temp, count, "pytemp/z_interp.dat");
  
    system("python interpolate_wham2.py");
    
    FILE * outfile;
    count=0;
    outfile = fopen("pytemp/g_interp.dat", "r");
    while(fgets(temp_char2, sizeof(temp_char2), outfile)){
      G_interp_temp[count] = atof(temp_char2);
      count++;
    }
    for(l=0;l<*invalid_count;l++){
      G_interp[l] = G_interp_temp[l];
      bin_interp[l] = bin_interp_temp[l];
      U_interp[l] = U_surface(bin_interp_temp[l], epsilon,sigma);
    }
  }

  for(l=0;l<recon_samples;l++){
    bin_full[l] = bin_centers[l];
    U_full[l] = U_surface(bin_centers[l], epsilon,sigma);
  }
 
  // swapping, always remembering that the last count of these arrays will be junk now.
  for(l=0; l<recon_samples-count;l++){
    G[l] = G_valid[l];
    bin_centers[l] = bin_valid[l];
    U_surf[l] = U_surf_valid[l];
  }

  count=0;

  for(k=0;k<traces;k++){
    *work_std += pow(ext_work[k][samples-1] - *work_mean,2)/traces;
  }
  *work_std = sqrt(*work_std)*beta;
  *work_mean = *work_mean*beta;
  *rms = RMS(G, U_surf, recon_samples-count)*beta;
  reconstruction_accuracy(G,U_surf,bin_centers ,beta,recon_samples-count,recon_slope,recon_msd, sigma);
  printf("work_mean[k_b*T] = %E , work_std[(k_bT)] = %E \n", *work_mean, *work_std);
  printf("recon_slope [k_bT/sigma] = %E , recon_msd[(k_bT)^2] = %E \n",*recon_slope, *recon_msd);
}
// double t_ar;

void check_if_exists(char * paste_file ){
  char  temp_char[200];
  if(access(paste_file, F_OK) == -1){
    sprintf(temp_char, "touch %s", paste_file);
    system(temp_char);
  }
}
void paste_column(char * ar_file, char * paste_file){
  char  temp_char[200];
  check_if_exists(paste_file);
  

  sprintf(temp_char, "paste %s %s > temp.dat", paste_file, ar_file);
  printf("%s \n", temp_char);
  system(temp_char);
  sprintf(temp_char, "cat temp.dat > %s", paste_file);
  printf("%s \n", temp_char);
  system(temp_char);
}

/* void interactive_plot_recon(){ */
/*   FILE * gnuplot_pipe = popen("gnuplot -persistant", "w"); */
/*   fprintf(gnuplot_pipe, "plot 'G_out.dat' using 1,\'U_surface.dat' using 1 with line \n"); */
/*   fflush(gnuplot_pipe); */
/* } */

void interactive_plot_recon(){
  FILE * gnuplot_pipe = popen("gnuplot -persistant", "w");
  fprintf(gnuplot_pipe, "plot '<paste bin_out.dat G_out.dat' using 1:2,\'<paste bin_out.dat U_surface.dat' using 1:2 with line \n");
  fflush(gnuplot_pipe);
}

void interactive_plot_exp2(){
  FILE * gnuplot_pipe = popen("gnuplot -persistant", "w");
  fprintf(gnuplot_pipe, "plot 'exp2_avarage.dat' using 1 with line \n");
  fflush(gnuplot_pipe);
}
void interactive_plot_work(){
  FILE * gnuplot_pipe = popen("gnuplot -persistant", "w");
  fprintf(gnuplot_pipe, "plot 'avarage_work.dat' using 1 with line \n");
  fflush(gnuplot_pipe);
}
void set_sang_protocol(double *prot, double A, double B,
		       double z_final, double t_final, double h,
		       int N, double sigma){
  // nasty long one here!
  double temp_x,v_s,tau,t_s,last_bit;
  int count,n,temp_steps,k,m;
  
  // needed quantities.
  v_s = z_final/(B*t_final);
  n = (int) floor((z_final-A*sigma)/sigma);
  // last_bit = (z_final-A*sigma)%sigma;
  last_bit = fmod(z_final-A*sigma, sigma);
  tau = t_final*(1-B)/((float) (n+1)); 
  t_s = B*t_final/( 1 + n*(1.0/A) + ((last_bit)/(A*sigma))   );

  printf("number of sigma switches = %d \n", n);
  printf("last little bit %E [sigma]", last_bit/sigma);



  
  count=0;temp_x=0;
  // First step
  temp_steps=(int) floor((t_s/h));
  for(k=0;k<temp_steps;k++){
    prot[k + count] = temp_x + v_s*h*k;
  }
  temp_x += A*sigma;
  count  += temp_steps;
  
  // mid steps n plateaus and n sigma switches
  for(m=0;m<n; m++){
    temp_steps = (int) floor((tau/h));
    for(k=0;k<temp_steps;k++){
      prot[count + k] = temp_x;
    }  
    count += temp_steps;
    temp_steps = (int) floor(t_s*(1.0/A)/h);
    for(k=0;k<temp_steps;k++){
      prot[count + k] = temp_x + v_s*h*k;
    }
    temp_x += sigma;
    count += temp_steps;  
  }
  // Last plateau 
  temp_steps = (int) floor((tau/h));
  for(k=0;k<temp_steps;k++){
    prot[count + k] = temp_x;
  }
  count += temp_steps;
  // Last switch
  temp_steps = N-count;
  for(k=0;k<temp_steps;k++){
    prot[count + k] = temp_x +k*v_s*h;
  }
};
double min_f(double x, double y){
  // returns the maximum of x and y
  if(x<=y){return x;}
  else{return y;}
}

void lc_deconvolution(double ** x_ar, double * lambda_ar,
		      int traces, int samples, double K, double beta, int recon_samples,
		      double z_start, double z_final,
		      double epsilon, double sigma, double * U_lc_surf,
		      double * work_mean,
		      double * bin_centers,
		      double *lc_bins, double * G_lc, int cutoff, int ITER_MAX){
  
  double * bin_edges =  (double *) malloc (sizeof (double) * (recon_samples + 1) );
  double d_bin = ((z_final - z_start)/recon_samples);
  double ** ext_work = (double **) malloc ( sizeof(double * )* traces );
  double * exp_avarage = (double *) malloc (sizeof(double) * samples);
  double * W = (double *) malloc (sizeof(double) * samples);
  double * W2 = (double *) malloc (sizeof(double) * samples);
  
  int i,k,l,m;
  //int cutoff=4;			/* should be an even number. */
  double s;
  double F_temp;
  double * F_mes = (double *) malloc(sizeof(double)*recon_samples);
  int * F_mes_counter = (int *) malloc(sizeof(int) *recon_samples); /* might change to double. */
  double * F_deconv = (double *) malloc(sizeof(double)*(recon_samples-cutoff));
  double * deconv_bins = (double *) malloc(sizeof(double)*(recon_samples-cutoff)); /* might have to be an argument */
  double * lc_norm = (double *) malloc(sizeof(double)*recon_samples); /* vector for c_i */
  double gauss_norm = sqrt(2*M_PI/(beta*K));
  
  for(i=0;i<recon_samples+1;i++){
    bin_edges[i] = z_start + i*d_bin;
  }
  for(i=0;i<recon_samples;i++){
    bin_centers[i] = (bin_edges[i+1] + bin_edges[i])*0.5;
  }
  for(l=0;l<recon_samples;l++){F_mes[l]=0;F_mes_counter[l]=0;}
  for(l=0;l<recon_samples-cutoff;l++){
    F_deconv[l]=1;
    deconv_bins[l] = bin_centers[l+cutoff/2];
    lc_bins[l] = bin_centers[l+cutoff/2];
  } /* seed */
  for(k=0;k<traces;k++){
    ext_work[k] = (double *) malloc (sizeof(double) * samples);
    ext_work[k][0] = 0;
    for(i=1;i<samples;i++){
      ext_work[k][i] = ext_work[k][i-1] - 0.5*K*((x_ar[k][i] - lambda_ar[i])
						 +(x_ar[k][i-1] - lambda_ar[i-1])
						 )*(lambda_ar[i]-lambda_ar[i-1]);
    }
    *work_mean += ext_work[k][samples-1]/traces;
  }
  for(i=0;i<samples;i++){
    exp_avarage[i] = 0;
    W[i] = 0;
    W2[i] = 0;
    for(k=0;k<traces;k++){
      exp_avarage[i] +=  exp(-beta*ext_work[k][i]);
      W[i] += ext_work[k][i];
      W2[i] += pow(ext_work[k][i],2);
    }
    /* W[i] =  W[i]/traces; */
    /* W2[i] = W2[i]/traces; */
    /* exp_avarage[i] = exp(-beta*W[i] + 0.5*pow(beta,2)*(W2[i] - pow(W[i],2))); */
    exp_avarage[i] = exp_avarage[i]/traces;
  }
 
  // Now when we have bins and exponential works, we need to avarage these guys over the bins. In this case though, we need to avarage wi
  // th respect to the support position and not the tip position.
  // now... this loop incomming might be very demanding indeed.
  // loop over edges in some clever way.
  
  for(i=0;i<samples;i++){
    for(l=0;l<recon_samples;l++){
      if(bin_edges[l] < lambda_ar[i] && bin_edges[l+1] > lambda_ar[i]){
	// in this case the support is in l'th bin
	// printf("Hist hit. in  %E < %E <%E \n", bin_edges[l],prot[i], bin_edges[l+1]);
	F_mes[l] += exp_avarage[i];
	F_mes_counter[l] += 1; 
      }
    }
  }
  for(l=0;l<recon_samples;l++){
    if(F_mes_counter[l] != 0){
      F_mes[l] = F_mes[l]/F_mes_counter[l]; /* Normalize if there ARE hits. */
    }
    printf("F_mes[%d] = %E \n", l,F_mes[l]);
  }
  s = vec_sum(F_mes, recon_samples);
  for(i=0;i<recon_samples;i++){F_mes[i] = F_mes[i]/s;}

  for(m=0;m<ITER_MAX;m++){
    
    for(k=0;k<recon_samples;k++){
      lc_norm[k] = 0;		/* temp */
      // lc_sum[k] = 0;		/* temp */
      for(l=0;l<recon_samples-cutoff;l++){
	//lc_norm[k] += gauss_norm*d_bin*exp(-0.5*K*beta*(pow(bin_centers[k]-deconv_bins[l],2)))*F_deconv[l];
	lc_norm[k] += gauss_norm*exp(-0.5*K*beta*(pow(bin_centers[k]-deconv_bins[l],2)))*F_deconv[l];

      }
    }
    for(l=0;l<recon_samples-cutoff;l++){
      F_temp = F_deconv[l]; 	/* save old */
      F_deconv[l] = 0;		/* set new to zero before adding. */
      for(k=0;k<recon_samples;k++){
	//F_deconv[l]+=gauss_norm*d_bin*F_mes[k]*exp(-0.5*beta*K*pow(bin_centers[k]- deconv_bins[l],2))/lc_norm[k];
	F_deconv[l]+=gauss_norm*F_mes[k]*exp(-0.5*beta*K*pow(bin_centers[k]- deconv_bins[l],2))/lc_norm[k];
      }
      F_deconv[l] = F_temp*F_deconv[l];
      //printf("iter %d :F_deconv[%d] = %E \n",m, l, F_deconv[l]);
    }
  }
  double g_const = U_surface(deconv_bins[0],epsilon, sigma) -(-pow(beta,-1)*log(F_deconv[0]));
  printf("g_const = %E \n", g_const);
  // when all is done. Update vector
  for(l=0;l<recon_samples-cutoff;l++){
    U_lc_surf[l] = U_surface(deconv_bins[l],epsilon, sigma);
    if(F_deconv[l] > 0){
      G_lc[l] = -pow(beta,-1)*log(F_deconv[l]);
    }
    G_lc[l] += g_const;
  } 
}
// double* t_at, recon_slope, recon_msd, work_std, exp2_average. average_work, rsm, G, sambples_work_variances, prot; 
// int* invalid_count; // int steps;

void damped_lc_deconvolution(double ** x_ar, double * lambda_ar,
			     int traces, int samples, double K, double beta, 
				 int recon_samples, double z_start, double z_final, 
				 double epsilon, double sigma, double * U_lc_surf,
			     double * work_mean, double * bin_centers, 
				 double *lc_bins, double * G_lc, int cutoff, 
				 int ITER_MAX, double T_param){
  
  double * bin_edges =  (double *) malloc (sizeof (double) * (recon_samples + 1) );
  double d_bin = ((z_final - z_start)/recon_samples);
  double ** ext_work = (double **) malloc ( sizeof(double * )* traces );
  double * exp_avarage = (double *) malloc (sizeof(double) * samples);
  int N = 10;			/* for N=1 should be normal LR */
  double q=1.0;			/* ffor faster conv. q=1 is normal LR scheme */
  int i,k,l,m;
  //int cutoff=4;			/* should be an even number. */
  double F_temp;
  double * F_mes = (double *) malloc(sizeof(double)*recon_samples);
  int * F_mes_counter = (int *) malloc(sizeof(int) *recon_samples); /* might change to double. */
  double * F_deconv = (double *) malloc(sizeof(double)*(recon_samples-cutoff));
  double * top_temp = (double *) malloc(sizeof(double)*(recon_samples-cutoff));
  double * transfer_sum = (double *) malloc(sizeof(double)*(recon_samples-cutoff));
  double * deconv_bins = (double *) malloc(sizeof(double)*(recon_samples-cutoff)); /* might have to be an argument */
  double * lc_norm = (double *) malloc(sizeof(double)*recon_samples); /* vector for c_i */
  double * U = (double *) malloc(sizeof(double)*recon_samples);
  double * U_tilde = (double *) malloc(sizeof(double)*recon_samples);
  double gauss_norm = sqrt(2*M_PI/(beta*K));
  double transfer_temp;
  
  for(i=0;i<recon_samples+1;i++){
    bin_edges[i] = z_start + i*d_bin;
  }
  for(i=0;i<recon_samples;i++){
    bin_centers[i] = (bin_edges[i+1] + bin_edges[i])*0.5;
  }
  for(l=0;l<recon_samples;l++){F_mes[l]=0;F_mes_counter[l]=0;}
  for(l=0;l<recon_samples-cutoff;l++){
    F_deconv[l]=1;
    deconv_bins[l] = bin_centers[l+cutoff/2];
    lc_bins[l] = bin_centers[l+cutoff/2];
  } /* seed */
  for(k=0;k<traces;k++){
    ext_work[k] = (double *) malloc (sizeof(double) * samples);
    ext_work[k][0] = 0;
    for(i=1;i<samples;i++){
      ext_work[k][i] = ext_work[k][i-1] - 0.5*K*((x_ar[k][i] - lambda_ar[i])
						 +(x_ar[k][i-1] - lambda_ar[i-1])
						 )*(lambda_ar[i]-lambda_ar[i-1]);
    }
    *work_mean += ext_work[k][samples-1]/traces;
  }
  for(i=0;i<samples;i++){
    exp_avarage[i] = 0;
    for(k=0;k<traces;k++){
      exp_avarage[i] +=  exp(-beta*ext_work[k][i]);
    }  
    exp_avarage[i] = exp_avarage[i]/traces;
  }
 
  // Now when we have bins and exponential works, we need to avarage these guys over the bins. In this case though, we need to avarage wi
  // th respect to the support position and not the tip position.
  // now... this loop incomming might be very demanding indeed.
  // loop over edges in some clever way.
  
  for(i=0;i<samples;i++){
    for(l=0;l<recon_samples;l++){
      if(bin_edges[l] < lambda_ar[i] && bin_edges[l+1] > lambda_ar[i]){
	// in this case the support is in l'th bin
	// printf("Hist hit. in  %E < %E <%E \n", bin_edges[l],prot[i], bin_edges[l+1]);
	F_mes[l] += exp_avarage[i];
	F_mes_counter[l] += 1; 
      }
    }
  }
  for(l=0;l<recon_samples;l++){
    if(F_mes_counter[l] != 0){
      F_mes[l] = F_mes[l]/F_mes_counter[l]; /* Normalize if there ARE hits. */
    }
    printf("F_mes[%d] = %E \n", l,F_mes[l]);
  }

  // lets see now...
  for(m=0;m<ITER_MAX;m++){
    
    for(k=0;k<recon_samples;k++){
      lc_norm[k] = 0;		/* temp */
      // lc_sum[k] = 0;		/* temp */
      for(l=0;l<recon_samples-cutoff;l++){
	//lc_norm[k] += gauss_norm*d_bin*exp(-0.5*K*beta*(pow(bin_centers[k]-deconv_bins[l],2)))*F_deconv[l];
	lc_norm[k] += gauss_norm*exp(-0.5*K*beta*(pow(bin_centers[k]-deconv_bins[l],2)))*F_deconv[l];
      }
      U[k] = -2*pow(T_param,-2)*(F_mes[k]*log(lc_norm[k]/F_mes[k]) - lc_norm[k] + F_mes[k]);
      U_tilde[k] = min_f(U[k],1); 
    }
    for(l=0;l<recon_samples-cutoff;l++){
      F_temp = F_deconv[l]; 	/* save old */
      F_deconv[l] = 0;		/* set new to zero before adding. */
      top_temp[l] = 0;
      transfer_sum[l] = 0;
      for(k=0;k<recon_samples;k++){
	transfer_temp = gauss_norm*exp(-0.5*beta*K*pow(bin_centers[k]- deconv_bins[l],2));
	top_temp[l] += transfer_temp*(1 + pow(U_tilde[k],N)*(N - (N-1)*U_tilde[k])*( (F_mes[k] - lc_norm[k])/lc_norm[k]  ));
	transfer_sum[l] += transfer_temp;
	//F_deconv[l]+=gauss_norm*d_bin*F_mes[k]*exp(-0.5*beta*K*pow(bin_centers[k]- deconv_bins[l],2))/lc_norm[k];
	//F_deconv[l]+=gauss_norm*F_mes[k]*exp(-0.5*beta*K*pow(bin_centers[k]- deconv_bins[l],2))/lc_norm[k];
      }
      //F_deconv[l] = F_temp*F_deconv[l];
      F_deconv[l] = F_temp*pow(top_temp[l]/transfer_sum[l],q);
      
      //printf("iter %d :F_deconv[%d] = %E \n",m, l, F_deconv[l]);
    }
  }
  double g_const = U_surface(deconv_bins[0],epsilon, sigma) -(-pow(beta,-1)*log(F_deconv[0]));
  // when all is done. Update vector
  for(l=0;l<recon_samples-cutoff;l++){
    U_lc_surf[l] = U_surface(deconv_bins[l],epsilon, sigma);
    if(F_deconv[l] > 0){
      G_lc[l] = -pow(beta,-1)*log(F_deconv[l]);
    }
    G_lc[l] += g_const;
  } 
}
// double* t_at, recon_slope, recon_msd, work_std, exp2_average. average_work, rsm, G, sambples_work_variances, prot, K; 
// int* invalid_count; // int steps, samples;

double L2_norm(double * x, int len){
  // computes standard norm.
  int i;
  double ret = 0;
  for(i=0;i<len;i++){
    ret += x[i]*x[i];
  }
  return sqrt(ret);
}

void vec_diff(double * diff, double * ar, int len, double dx){
  // returns difference vector of length len-1
  int i;
  for(i=0; i<len-1; i++){
    diff[i] = (ar[i+1] - ar[i])/dx; 
  }
}

void vec_copy(double * from_ar, double * to_ar, int len){
  int i;
  for(i=0;i<len;i++){
    to_ar[i] = from_ar[i];
  }
}

double vec_sum(double * ar, int len){
  int i;
  double s;
  for(i=0;i<len;i++){s += ar[i];}
  return s;
}

void acc_damped_lc_deconvolution(double ** x_ar, double * lambda_ar,
				 int traces, int samples, double K, double beta, 
				 int recon_samples, double z_start, double z_final, 
				 double epsilon, double sigma, double * U_lc_surf,
				 double * work_mean, double * bin_centers, 
				 double *lc_bins, double * G_lc, int cutoff, 
				 int ITER_MAX, double T_param, double * signal){

  double * bin_edges =  (double *) malloc (sizeof (double) * (recon_samples + 1) );
  double d_bin = ((z_final - z_start)/recon_samples);
  double ** ext_work = (double **) malloc ( sizeof(double * )* traces );
  double * exp_avarage = (double *) malloc (sizeof(double) * samples);
  double * W = (double *) malloc (sizeof(double) * samples);
  double * W2 = (double *) malloc (sizeof(double) * samples);
  


  int N = 4;			/* for N=1 should be normal LR */
  double q=1.0;			/* ffor faster conv. q=1 is normal LR scheme */
  double q_old = 0;
  double dq_lim = 0;
  int i,k,l,m;
  double s;
  //int cutoff=4;			/* should be an even number. */
  double F_temp;
  double * F_mes = (double *) malloc(sizeof(double)*recon_samples);
  double ** grad = (double **) malloc (sizeof(double *) *2); /* only need two of these */
  for(i=0;i<2;i++){grad[i] = (double *) malloc(sizeof(double *)*(recon_samples - cutoff - 1)); }
  double grad_quote;
  
  int * F_mes_counter = (int *) malloc(sizeof(int) *recon_samples); /* might change to double. */
  double * F_deconv = (double *) malloc(sizeof(double)*(recon_samples-cutoff));
  double * top_temp = (double *) malloc(sizeof(double)*(recon_samples-cutoff));
  double * transfer_sum = (double *) malloc(sizeof(double)*(recon_samples-cutoff));
  double * deconv_bins = (double *) malloc(sizeof(double)*(recon_samples-cutoff)); /* might have to be an argument */
  double * lc_norm = (double *) malloc(sizeof(double)*recon_samples); /* vector for c_i */
  double * U = (double *) malloc(sizeof(double)*recon_samples);
  double * U_tilde = (double *) malloc(sizeof(double)*recon_samples);
  double gauss_norm = sqrt(2*M_PI/(beta*K));
  double transfer_temp;
  
  for(i=0;i<recon_samples+1;i++){
    bin_edges[i] = z_start + i*d_bin;
  }
  for(i=0;i<recon_samples;i++){
    bin_centers[i] = (bin_edges[i+1] + bin_edges[i])*0.5;
  }
  for(l=0;l<recon_samples;l++){F_mes[l]=0;F_mes_counter[l]=0;}
  for(l=0;l<recon_samples-cutoff;l++){
    F_deconv[l]=1;
    deconv_bins[l] = bin_centers[l+cutoff/2];
    lc_bins[l] = bin_centers[l+cutoff/2];
  } /* seed */
  for(k=0;k<traces;k++){
    ext_work[k] = (double *) malloc (sizeof(double) * samples);
    ext_work[k][0] = 0;
    for(i=1;i<samples;i++){
      ext_work[k][i] = ext_work[k][i-1] - 0.5*K*((x_ar[k][i] - lambda_ar[i])
						 +(x_ar[k][i-1] - lambda_ar[i-1])
						 )*(lambda_ar[i]-lambda_ar[i-1]);
    }
    *work_mean += ext_work[k][samples-1]/traces;
  }
  for(i=0;i<samples;i++){
    exp_avarage[i] = 0;
    W[i] = 0;
    W2[i] = 0;
    for(k=0;k<traces;k++){
      exp_avarage[i] +=  exp(-beta*ext_work[k][i]);
      W[i] += ext_work[k][i];
      W2[i] += pow(ext_work[k][i],2);
    }
    W[i] =  W[i]/traces;
    W2[i] = W2[i]/traces;
    exp_avarage[i] = exp(-beta*W[i] + 0.5*pow(beta,2)*(W2[i] - pow(W[i],2)));
    //exp_avarage[i] = exp_avarage[i]/traces;
    
  }
 
  // Now when we have bins and exponential works, we need to avarage these guys over the bins. In this case though, we need to avarage wi
  // th respect to the support position and not the tip position.
  // now... this loop incomming might be very demanding indeed.
  // loop over edges in some clever way.
  
  for(i=0;i<samples;i++){
    for(l=0;l<recon_samples;l++){
      if(bin_edges[l] < lambda_ar[i] && bin_edges[l+1] > lambda_ar[i]){
	// in this case the support is in l'th bin
	// printf("Hist hit. in  %E < %E <%E \n", bin_edges[l],prot[i], bin_edges[l+1]);
	F_mes[l] += exp_avarage[i];
	F_mes_counter[l] += 1; 
      }
    }
  }
  for(l=0;l<recon_samples;l++){
    if(F_mes_counter[l] != 0){
      F_mes[l] = F_mes[l]/F_mes_counter[l]; /* Normalize if there ARE hits. */
    }
    printf("F_mes[%d] = %E \n", l,F_mes[l]);
  }

  s = vec_sum(F_mes, recon_samples);
  vec_copy(F_mes, signal, recon_samples);
  for(i=0;i<recon_samples;i++){F_mes[i] = F_mes[i]/s;}

  // Up to this point is actually just getting the signal.

  // For the accellerated version we must start by doing two iterations and saving the data.
  
  for(m=0;m<2;m++){
    
    for(k=0;k<recon_samples;k++){
      lc_norm[k] = 0;		/* temp */
      // lc_sum[k] = 0;		/* temp */
      for(l=0;l<recon_samples-cutoff;l++){
	//lc_norm[k] += gauss_norm*d_bin*exp(-0.5*K*beta*(pow(bin_centers[k]-deconv_bins[l],2)))*F_deconv[l];
	lc_norm[k] += gauss_norm*exp(-0.5*K*beta*(pow(bin_centers[k]-deconv_bins[l],2)))*F_deconv[l];
      }
      U[k] = -2*pow(T_param,-2)*(F_mes[k]*log(lc_norm[k]/F_mes[k]) - lc_norm[k] + F_mes[k]);
      U_tilde[k] = min_f(U[k],1); 
    }
    for(l=0;l<recon_samples-cutoff;l++){
      F_temp = F_deconv[l]; 	/* save old */
      F_deconv[l] = 0;		/* set new to zero before adding. */
      top_temp[l] = 0;
      transfer_sum[l] = 0;
      for(k=0;k<recon_samples;k++){
	transfer_temp = gauss_norm*exp(-0.5*beta*K*pow(bin_centers[k]- deconv_bins[l],2));
	top_temp[l] += transfer_temp*(1 + pow(U_tilde[k],N)*(N - (N-1)*U_tilde[k])*( (F_mes[k] - lc_norm[k])/lc_norm[k]  ));
	transfer_sum[l] += transfer_temp;
	//F_deconv[l]+=gauss_norm*d_bin*F_mes[k]*exp(-0.5*beta*K*pow(bin_centers[k]- deconv_bins[l],2))/lc_norm[k];
	//F_deconv[l]+=gauss_norm*F_mes[k]*exp(-0.5*beta*K*pow(bin_centers[k]- deconv_bins[l],2))/lc_norm[k];
      }
      //F_deconv[l] = F_temp*F_deconv[l];
      F_deconv[l] = F_temp*pow(top_temp[l]/transfer_sum[l],q);
      
      //printf("iter %d :F_deconv[%d] = %E \n",m, l, F_deconv[l]);
    }
    vec_diff(grad[m], F_deconv, recon_samples-cutoff, d_bin );
    //vec_diff(grad[m], F_deconv, recon_samples-cutoff, 1 );
    
  }

  grad_quote = L2_norm(grad[1],recon_samples-cutoff) / L2_norm(grad[0],recon_samples-cutoff); 
  q = exp(grad_quote) - grad_quote;
  printf("First q = %E \n", q);

  // Iteration starting
  for(m=0;m<ITER_MAX;m++){
    
    for(k=0;k<recon_samples;k++){
      lc_norm[k] = 0;		/* temp */
      // lc_sum[k] = 0;		/* temp */
      for(l=0;l<recon_samples-cutoff;l++){
	//lc_norm[k] += gauss_norm*d_bin*exp(-0.5*K*beta*(pow(bin_centers[k]-deconv_bins[l],2)))*F_deconv[l];
	lc_norm[k] += gauss_norm*exp(-0.5*K*beta*(pow(bin_centers[k]-deconv_bins[l],2)))*F_deconv[l];
      }
      U[k] = -2*pow(T_param,-2)*(F_mes[k]*log(lc_norm[k]/F_mes[k]) - lc_norm[k] + F_mes[k]);
      U_tilde[k] = min_f(U[k],1); 
    }
    for(l=0;l<recon_samples-cutoff;l++){
      F_temp = F_deconv[l]; 	/* save old */
      F_deconv[l] = 0;		/* set new to zero before adding. */
      top_temp[l] = 0;
      transfer_sum[l] = 0;
      for(k=0;k<recon_samples;k++){
	transfer_temp = gauss_norm*exp(-0.5*beta*K*pow(bin_centers[k]- deconv_bins[l],2));
	top_temp[l] += transfer_temp*(1 + pow(U_tilde[k],N)*(N - (N-1)*U_tilde[k])*( (F_mes[k] - lc_norm[k])/lc_norm[k]  ));
	transfer_sum[l] += transfer_temp;
	//F_deconv[l]+=gauss_norm*d_bin*F_mes[k]*exp(-0.5*beta*K*pow(bin_centers[k]- deconv_bins[l],2))/lc_norm[k];
	//F_deconv[l]+=gauss_norm*F_mes[k]*exp(-0.5*beta*K*pow(bin_centers[k]- deconv_bins[l],2))/lc_norm[k];
      }
      //F_deconv[l] = F_temp*F_deconv[l];
      F_deconv[l] = F_temp*pow(top_temp[l]/transfer_sum[l],q);
      
      //printf("iter %d :F_deconv[%d] = %E \n",m, l, F_deconv[l]);
    }
    // We need to update q here.
    vec_copy(grad[1],grad[0], recon_samples-cutoff);
    vec_diff(grad[1], F_deconv, recon_samples-cutoff, d_bin );
    //vec_diff(grad[1], F_deconv, recon_samples-cutoff, 1 );
    q_old = q;
    q = exp(L2_norm(grad[1],recon_samples-cutoff) / L2_norm(grad[0],recon_samples-cutoff)) - grad_quote;
    printf("q[%d] = %E \n",m, q);
    if(fabs(q-q_old) < dq_lim){printf("%E < dq_lim(%E) breaking at %d %d %d \n",fabs(q-q_old),dq_lim, m,m,m);break;}
  }
  
  // Iterations and all that done at this point, just converting to free energy and setting zero level.
  double g_const = U_surface(deconv_bins[0],epsilon, sigma) -(-pow(beta,-1)*log(F_deconv[0]));
  // when all is done. Update vector
  for(l=0;l<recon_samples-cutoff;l++){
    U_lc_surf[l] = U_surface(deconv_bins[l],epsilon, sigma);
    if(F_deconv[l] > 0){
      G_lc[l] = -pow(beta,-1)*log(F_deconv[l]);
    }
    G_lc[l] += g_const;
  } 
}
// double* t_at, recon_slope, recon_msd, work_std, exp2_average, average_work, rsm, G, sambples_work_variances, prot; 
// int* invalid_count; // int steps;


