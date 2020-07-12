#include "headers.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_fit.h>
#include <time.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
//#define PI 4*(atan(1))

// TODO
// A) One danger remains. When you have few traces its not at all unlikely that you will get some infs.
//    and when you define the constant in G(z) = U(z) + constant. It will basically fizzle if G[0] is inf.
//    You could insert safeguards for this if necessary.

// B) Implement some function that quantifies reconstruction error. This is a little subtle, since in general it ought to depend
//    on the number of reconstruction samples.


int main(int argc, char * argv[]){

  gsl_rng * rng;
  rng = gsl_rng_alloc(gsl_rng_mt19937);
  int random_seed = (int)time(NULL);
  gsl_rng_set(rng, random_seed);
  struct stat st = {}; 	/* for checking if input dir exists */
 
  
  unsigned int i,j,k,m;
  // unsigned int l,n;
  //char temp_char[100];
  //char temp_char2[100];
  //int invalid_count;
  unsigned long steps = 1000000;
  unsigned int traces = 2;
  unsigned int equil_steps = 100000;
  int sample_freq;
  /* HARDCODED SAMPLE FREQUENCY! WATCH Out! */
  //int samples = 5000*2*1;
  int samples = 5000;
  
  int recon_samples = 40;
  double W,h,z_max,z_min;	/* So that you may start the protocol wherever you damn well please basically */
  z_max=0;z_min=0;
  double protocol_constant = 0;

  double T_final=0.5;
 // double z_start;
  double lambda;
  double x;
  //double t;
  int p = 0;
  DECLARE_CONSTANTS();
  //printf(" testing outdir %s", outdir);
  get_args(argc, argv,&steps,outfile,&epsilon_d, &traces,
	   &z_final_d, &recon_samples, &T_final, info_file, outdir,
	   &A, &B, &p, &protocol_constant, paste_file, &R);
  //samples=steps;
  //D_d = pow(10,-4);
  //K=5.7;			/* maximum of tribolever */
  UPDATE_CONSTANTS();
  printf("sample size= %d, sample_freq = %d \n", samples, sample_freq);
  printf("Stepsize is currently %E \n", h);

  printf("Size of (double *) pointer is %zd \n", sizeof(double));
  // Fucking retard... Completely different spring stiffnesses
  
  //double recon_slope, recon_msd, work_mean, work_std, rms;                                   /* msd for mean square deviation.  */
  //double rk[4]; 
  //double gauss_sigma = sqrt(2*pow(k_b*T, 2)/(D*h));               		     	/* just in case of RK4 solving.    */
  double ** x_ar = (double **) malloc (sizeof(double *) * traces);
  double * lambda_ar = (double *) malloc (sizeof (double)* samples);
  
  double * t_ar = (double *) malloc (sizeof (double)* samples);
  //double * G = (double *) malloc (sizeof (double) * recon_samples);
  //double * U_surf = (double *) malloc (sizeof (double) * recon_samples);
  //double * bin_centers = (double *) malloc (sizeof (double) * recon_samples);
  double * prot = (double *) malloc (sizeof (double)*steps);             /* Should add safeguard here so that you dont overallocate.  */
  
  //double * exp2_avarage = (double *) malloc (sizeof (double) * samples);
  //double * avarage_work = (double *) malloc (sizeof (double) * samples);
  //double * sample_work_variances = (double *) malloc (sizeof (double) * samples); /* changed size from traces */
  //double * work_vector = (double *) malloc (sizeof (double) * traces);
  // DISCRETIZING PROTOCOL, MAINLY FOR SPEED, This routine by itself is very slow though. Lets implement a reverse protocol thingy...

  
  MAKE_PROTOCOL();
  printf("DEBUG CHECK \n");
  /* set_sang_protocol(prot,A, B,z_final,T_final, h, */
  /* 		    steps, sigma); */
  /* z_min=0; */
  /* z_max = z_final; */
  
  printf("discretizing protocol complete. z_max[sigma] = %E , z_min[sigma] = %E \n", z_max/sigma, z_min/sigma);
  printf("protocol constant is set to %E \n", protocol_constant);
  for(i=0;i<traces;i++){
    x_ar[i] = (double *) malloc (sizeof(double *) * samples); 
  }
  
  //W_amp = 0;
  //W_amp = sqrt(h)*(2*D);
  for(m=0; m<traces; m++){
    //t = 0; 
    //calculate_lambda(&lambda, t, T_final, z_final, A , B);
    //lambda = lambda + protocol_constant;
    //x = lambda;
    x = prot[0] + protocol_constant;
    k=0;
    for(j=0; j<equil_steps; j++){
      W = gsl_ran_gaussian_ziggurat(rng,1);
      /* x = heun(x,lambda,K,W,W_amp,h, epsilon, */
      /* 	        sigma,beta, D); */
      x = heun(x, prot[0] + protocol_constant,K,W,W_amp,h, epsilon,
	        sigma,beta, D);
      /* if(j<100){ */
      /* 	printf("x = %E \n", x);	/\* nuclear option for testing *\/ */
      /* } */
      
    }
    /* printf("x_equil[%d] = %E \n", m, x); */
    for(i=0;i<steps;i++){
      // HEUN SOLVING
      W = gsl_ran_gaussian_ziggurat(rng,1);
      /* x = heun(x,lambda,K,W,W_amp,h, epsilon, */
      /* 	        sigma,beta, D); */
      x = heun(x,prot[i]+protocol_constant,K,W,W_amp,h, epsilon,
      	        sigma,beta, D);
      
      // RK4 SOLVING.
      /* W = gsl_ran_gaussian_ziggurat(rng,1)*gauss_sigma; */
      /* RK4(&x,t,rk,epsilon,K, sigma,beta, W, */
      /* 	  D,h, h_i, &lambda,T_final,z_final); */
      
      //t = t + h;
      //calculate_lambda(&lambda, t, T_final, z_final, A, B);
      if(i%sample_freq==0){
	x_ar[m][k] = x;
	if(m==0){lambda_ar[k] = prot[i]+protocol_constant;}
	//printf("x[%d][%d] = %E \n", m,k,x);
	k++;
      }
      
    }
    /* if(m==traces/2){ */
    /*   B=1; */
    /*   MAKE_PROTOCOL(); // its not quite this simple. You will need new biasing things in the wham! How you do this can be read in minh paper. */
      
    /* } */
    printf("Trace %d done \n", m);
  }

  printf("SEGTEST \n");

  if (stat(outdir, &st) == -1) {
    mkdir(outdir, 0700);	/* Works, checks if directory exists, if it doesn't then it simply creates */
  }
  chdir(outdir);		/* This changes directory so that we write in the proper outdir, after program exits it goes back */
  printf("directory argument is %s \n", outdir);
  
  /* Feel free to remove the ones you don't need, this will take alot of space otherwise. */

  printf("in run.c the traces are %d \n", traces);
  write_2darray(x_ar,samples,traces, "x_out.dat");
  write_1darray(t_ar,samples, "t_out.dat");
  // write_1darray(work_vector, traces, "work_dist.dat");
  write_1darray(lambda_ar, samples, "lambda_out.dat");
  /* write_paramfile("params.dat", z_final,sigma, D, T, T_final, steps,traces, */
  /* 		  recon_samples-invalid_count,z_start, info_file,  epsilon, h, recon_slope, recon_msd, */
  /* 		  A, B, work_mean, work_std, z_min, z_max); */
  
  // UNCOMMENT THIS BLOCK IF YOU WANT RECON
  /* write_1darray(G_lc, recon_samples-cutoff, "G_lc.dat"); */
  /* write_1darray(G, recon_samples-invalid_count, "G_out.dat"); */
  /* write_1darray(G_interp, invalid_count, "G_interp.dat"); */
  
  /* write_1darray(bin_centers, recon_samples-invalid_count, "bin_out.dat"); */
  /* write_1darray(bin_interp, invalid_count, "bin_interp.dat"); */
  /* write_1darray(bin_full, recon_samples, "bin_full.dat"); */

  /* write_1darray( U_surf, recon_samples-invalid_count, "U_surface.dat"); */
  /* write_1darray( U_interp, invalid_count, "U_interp.dat"); */
  /* write_1darray( U_full, recon_samples, "U_full.dat"); */

  
  /* /\* for(i=0;i<samples;i++){exp2_avarage[i] = log(exp2_avarage[i]);} *\/ */
  /* /\* write_1darray(exp2_avarage, samples, "exp2_avarage.dat"); *\/ */
  /* /\* write_1darray(avarage_work, samples, "avarage_work.dat"); *\/ */
  /* /\* write_1darray(sample_work_variances, samples, "sample_work_variance.dat"); *\/ */

  /* // UNCOMMENT IF YOU WANT TO PASTE FiLES */
  /* sprintf(temp_char, "G_out.dat"); */
  /* sprintf(temp_char2, "%s_G.dat", paste_file); */
  /* paste_column(temp_char, temp_char2); */

  /* sprintf(temp_char, "G_interp.dat"); */
  /* sprintf(temp_char2, "%s_G_interp.dat", paste_file); */
  /* paste_column(temp_char, temp_char2); */

  /* sprintf(temp_char, "bin_out.dat"); */
  /* sprintf(temp_char2, "%s_bins.dat", paste_file); */
  /* paste_column(temp_char, temp_char2); */
  
  /* sprintf(temp_char, "bin_interp.dat"); */
  /* sprintf(temp_char2, "%s_bins_interp.dat", paste_file); */
  /* paste_column(temp_char, temp_char2); */
  
  
  /* /\* sprintf(temp_char, "bin_out.dat"); *\/ */
  /* /\* sprintf(temp_char2, "%s_bins.dat", paste_file); *\/ */
  /* /\* paste_column(temp_char, temp_char2); *\/ */
  
  /* /\* sprintf(temp_char, "U_surface.dat"); *\/ */
  /* /\* sprintf(temp_char2, "%s_U.dat", paste_file); *\/ */
  /* /\* paste_column(temp_char, temp_char2); *\/ */

  /* // You need to collect to more pieces of data, the work and the final support position. */
  /* /\* sprintf(temp_char2, "%s_work.dat", paste_file); *\/ */
  /* /\* FILE * append_file = fopen(temp_char2, "ab+"); *\/ */
  /* /\* fprintf(append_file,"%E \t %E \n", (prot[steps-1] + protocol_constant)/sigma, work_mean ); *\/ */
  /* /\* fclose(append_file); *\/ */
  

  /* // DO NOT REMOVE THIS PART, UNCOMMENT WHEN YOU WANT TO APPEND TO A DATA FILE. */
  /* // check_if_exists(outfile); */

  /* /\* double miss_bias = ((double) recon_samples - invalid_count )/((double) recon_samples); *\/ */
  /* /\* append_file = fopen(outfile, "ab+"); *\/ */
  /* /\* printf("Caveman debug \n"); *\/ */
  /* /\* fprintf(append_file, "%E \t %E \t %E \t %E \t %E \t %E \t %E \t %E \n", *\/ */
  /* /\* 	  A,B,protocol_constant/sigma, recon_msd, recon_slope, rms, rms/miss_bias, ((z_max-z_min)/sigma) ); *\/ */
  /* /\* fclose(append_file); *\/ */
  

  /* if(p){ */
  /*   if(p==1){interactive_plot_recon();} */
  /*   if(p==2){interactive_plot_exp2();} */
  /*   if(p==3){interactive_plot_work();} */
  /* } */


  
  
  return 0;
}



/* sprintf(temp_char, "exp2_avarage.dat"); */
  /* sprintf(temp_char2, "%s_exp2.dat", paste_file); */
  /* paste_column(temp_char, temp_char2); */
  
  /* sprintf(temp_char, "lambda_out.dat"); */
  /* sprintf(temp_char2, "%s_lambda.dat", paste_file); */
  /* paste_column(temp_char, temp_char2); */
  
  /* sprintf(temp_char, "avarage_work.dat"); */
  /* sprintf(temp_char2, "%s_work.dat", paste_file); */
  /* paste_column(temp_char, temp_char2); */
  
