#ifndef _orig_headers_h
#define _orig_headers_h
#define k_b (1.3806488*pow(10, -23))
void get_signal_with_var(double * signal, int len, double z_start, double z_final,
			 int traces, int samples, double * lambda_ar, double K,
			 double ** x_ar, int c, double * bin_centers, double beta, double * var);
void wham_reconstruct_interp_rare(double ** x_ar, double * t_ar, double * lambda_ar,
				  int traces, int samples, double K, double beta, int recon_samples,
				  double z_start, double z_final, double * G,
				  double epsilon, double sigma, double * U_surf,
				  double * recon_slope, double * recon_msd, double * work_mean, double * work_std,
				  double * exp2_avarage, double * avarage_work, double *rms, int * invalid_count,
				  double * bin_centers, double * samples_work_variances,
				  double * rare_traj, double * work_min);
void berkovic_works(double ** ext_work, double ** x_ar, double *lambda_ar, int samples, int traces,
		    double slip_distance, double K);
void get_weight_work(double ** ext_work, int samples, int traces,
		     double * weight_work, int control, double beta);
void get_third_cumulant(double ** ext_work, int samples, int traces,
			double * weight_work, int control, double beta);
void make_protocol(int reverse, int samples, double * prot, double T_final, double z_final, double A, double B, double * lambda, double h,
		   int steps);
void wham_bias(int len, int bias_len ,double z_start, double z_final,
	       int traces, int samples, double ** lambda_ar, double K,
	       double ** x_ar, int c, double * bin_centers, double beta,
	       double z_0, double z_f, double * G, int * invalid_count,
	       double epsilon, double sigma, double * U_surf);
void wham_bias2(int len, int bias_len ,double z_start, double z_final,
		int traces, int samples, double ** lambda_ar, double K,
		double ** x_ar, int c, double * bin_centers, double beta,
		double z_0, double z_f, double * G, int * invalid_count,
		double epsilon, double sigma, double * U_surf);

void wiener_deconvolution(double * signal, int len, double z_start, double z_final,
			  int traces, int samples, double * lambda_ar, double K,
			  double ** x_ar, int c, double * bin_centers, double beta,
			  double * lc_bins, double * G_lc, double * U_lc_surf, double U, double sigma,
			  double regularization_param, int cutoff);
void svd_deconvolution(double * signal, int len, double z_start, double z_final,
		       int traces, int samples, double * lambda_ar, double K,
		       double ** x_ar, int c, double * bin_centers, double beta,
		       double * lc_bins, double * G_lc, double * U_lc_surf, double U, double sigma,
		       double param, int cutoff);
void twomey_deconvolution(double * signal, int len, double z_start, double z_final,
			  int traces, int samples, double * lambda_ar, double K,
			  double ** x_ar, int c, double * bin_centers, double beta,
			  double * lc_bins, double * G_lc, double * U_lc_surf, double U, double sigma,
			  double param, int cutoff);
void twomey_deconvolution2(double * signal, int len, double z_start, double z_final,
			   int traces, int samples, double * lambda_ar, double K,
			   double ** x_ar, int c, double * bin_centers, double beta,
			   double * lc_bins, double * G_lc, double * U_lc_surf, double U, double sigma,
			   double param, int cutoff, double * wham_recon_arr, double g_const);

double force(double x, double lambda, double beta, double D,
	     double epsilon, double sigma, double K);
double heun(double x, double lambda, double K, double W,
	    double W_amp, double h, double epsilon,
	    double sigma, double beta, double D);
double rk_force(double x, double lambda, double epsilon, double K,
		double sigma, double beta, double rand, double D );
void ermak(double * x, double * p, double lambda1, double lambda2, double x_rand, double p_rand,
	   double eta_coeff, double eta_coeff2, double t2_2m, double beta, double exp_eta,
	   double sigma, double D, double epsilon, double K, double m_i);
void ricci(double * x, double * p, double lambda1, double lambda2, double K, double W,
	   double W_amp, double h, double epsilon, double M,
	   double sigma, double beta, double D, double eta,
	   double t_m, double t2_2m, double e1,
	   double e2, double e4, double e34, double W_amp3, double x_rand_coeff, double h_2);
void fill_transfer_matrix_svd(double ** H, int len1, int len2, double beta, double K,
			      double * bin_centers, double * bin_centers_reduced);
void RK4(double * x, double t, double * rk,
	 double epsilon, double K,
	 double sigma, double beta, double rand,
	 double D, double h, double h_i, double * lambda,
	 double T_final, double z_final, double A, double B);

double h_step(double t);
void calculate_lambda(double * lambda, double t, double T_final,
		      double z_final, double A, double B);
void write_1darray(double * ar, int len, char *filename);
void write_2darray(double ** x_ar, int len1, int len2, char *x_file);
double U_surface(double x, double epsilon, double sigma);

void reconstruction_accuracy(double * G, double * U_surf, double * bin_centers , double beta, int recon_samples, double * recon_slope, double * recon_msd,
			     double sigma);

void wham_reconstruct(double ** x_ar, double * t_ar, double * lambda_ar,
		      int traces, int samples, double K, double beta, int recon_samples,
		      double z_start, double z_final, double * G,
		      double epsilon, double sigma, double * U_surf, double * recon_slope, double * recon_msd,
		      double * work_mean, double * work_std, double * exp2_avarage,
		      double * avarage_work, double * rms, int * invalid_count, double * bin_centers,
		      double * sample_work_variances,double * G_interp, double * U_interp, double * bin_interp,
		      double * bin_full, double * U_full);
// Interpolates misses instead, can be trouble if first and last are not valid. Otherwise should work
// (Ofc it won't be good but at least it should be fair.)
void wham_reconstruct_interp(double ** x_ar, double * t_ar, double * lambda_ar,
			     int traces, int samples, double K, double beta, int recon_samples,
			     double z_start, double z_final, double * G,
			     double epsilon, double sigma, double * U_surf,
			     double * recon_slope, double * recon_msd, double * work_mean, double * work_std,
			     double * exp2_avarage, double * avarage_work, double *rms, int * invalid_count,
			     double * bin_centers, double * samples_work_variances);

void lc_deconvolution(double ** x_ar, double * t_ar, double * lambda_ar,
		      int traces, int samples, double K, double beta, int recon_samples,
		      double z_start, double z_final, double * G,
		      double epsilon, double sigma, double * U_lc_surf,
		      double * recon_slope, double * recon_msd, double * work_mean, double * work_std,
		      double * exp2_avarage, double * avarage_work, double *rms, int * invalid_count,
		      double * bin_centers, double * samples_work_variances,
		      double *lc_bins, double * G_lc, int cutoff, int ITER_MAX,
		      double * prot, int steps);
void damped_lc_deconvolution(double ** x_ar, double * t_ar, double * lambda_ar,
			     int traces, int samples, double K, double beta, int recon_samples,
			     double z_start, double z_final, double * G,
			     double epsilon, double sigma, double * U_lc_surf,
			     double * recon_slope, double * recon_msd, double * work_mean, double * work_std,
			     double * exp2_avarage, double * avarage_work, double *rms, int * invalid_count,
			     double * bin_centers, double * samples_work_variances,
			     double *lc_bins, double * G_lc, int cutoff, int ITER_MAX, double * prot, int steps,
			     double T_param);

void vec_diff(double * diff, double * ar, int len, double dx);
void vec_copy(double * from_ar, double * to_ar, int len);
double L2_norm(double * x, int len);
double vec_sum(double * ar, int len);
void acc_damped_lc_deconvolution(double ** x_ar, double * t_ar, double * lambda_ar,
				 int traces, int samples, double K, double beta, int recon_samples,
				 double z_start, double z_final, double * G,
				 double epsilon, double sigma, double * U_lc_surf,
				 double * recon_slope, double * recon_msd, double * work_mean, double * work_std,
				 double * exp2_avarage, double * avarage_work, double *rms, int * invalid_count,
				 double * bin_centers, double * samples_work_variances,
				 double *lc_bins, double * G_lc, int cutoff, int ITER_MAX, double * prot, int steps,
				 double T_param,double * signal);
void write_paramfile(char * filename, double z_final, double sigma, double D, double T,
		     double T_final, unsigned long int steps, unsigned int traces,
		     int recon_samples, double z_start, char * information_string, double epsilon, double h,
		     double recon_slope, double recon_msd, double A, double B, double work_mean, double work_std,
		     double z_min, double z_max);
void paste_column(char * ar_file,  char * paste_file);
void check_if_exists(char * paste_file );
double RMS(double * ar1, double * ar2, int len);
#define DECLARE_CONSTANTS()			\
  double epsilon; double epsilon_d;		\
  double v; double v_d;				\
  double sigma; double sigma_d;			\
  double eta; double eta_d;			\
  double T; double T_d;				\
  double K; double K_d;				\
  double D; double D_d;				\
  double A=0; double B=0;		        \
  int R = 0;					\
  double W_amp;					\
  double beta;					\
  double z_final; double z_final_d;		\
  char outfile[80] = "default_out.dat";		\
  char paste_file[80] = "paste_default";	\
  double h_i;					\
  char info_file[100] = "No specifics";	\
  char outdir[100] = "default_out";		\
  T = 300;T_d=1;				\
  epsilon= k_b*T;epsilon_d=1;			\
  v=pow(10,-9);v_d=1;				\
  sigma=3*pow(10,-10);sigma_d=1;		\
  K=2;K_d=1;					\
  D = 530 * pow(10, -18);D_d=1;			\
  z_final_d=1;
  
#define UPDATE_CONSTANTS()				\
  T=T*T_d;						\
  epsilon=epsilon*epsilon_d;				\
  v=v*v_d;						\
  sigma=sigma*sigma_d;					\
  K=K*K_d;						\
  D=D*D_d;						\
  beta=pow(T*k_b,-1);					\
  h = T_final/steps;					\
  h_i = pow(h,-1);					\
  W_amp=sqrt(h)*sqrt(2*D);				\
  z_final = sigma;			       \
  z_final = z_final*z_final_d;		\
  sample_freq = steps/samples;	       \
  protocol_constant = protocol_constant*sigma;\



int get_args(int argc, char *argv[],
	     long unsigned int *steps,  char * outfile,
	     double * epsilon_d, unsigned int * traces,
	     double *z_final_d, int * recon_samples,
	     double * T_final, char * info_file, char * outdir, 
	     double * A, double * B, int * p, double * protocol_constant,
	     char * paste_file, int * R);
void get_work_dist(double ** x_ar, double * t_ar, double * lambda_ar,
		   int traces, int samples, double K, double beta,
		   double sigma, double * work_vector);
void interactive_plot_recon();
void interactive_plot_exp2();
void interactive_plot_work();



#define MAKE_PROTOCOL()						\
  k=0;								\
  R=0;								\
  z_min=protocol_constant;					\
  z_max=protocol_constant;					\
  for(i=0;i<steps;i++){						\
    if(R==0){								\
      calculate_lambda(&lambda, i*h, T_final, z_final, A, B);		\
      prot[i] = lambda;							\
      if(i%sample_freq==0){lambda_ar[k] = lambda + protocol_constant; t_ar[k]=i*h;k++;} \
    }									\
    else{								\
      calculate_lambda(&lambda, h*(steps-i), T_final, z_final, A, B);	\
      prot[i] = lambda;							\
      if(i%sample_freq==0){lambda_ar[k] = lambda + protocol_constant; t_ar[k]=i*h;k++;} \
    }									\
    									\
    if(lambda + protocol_constant > z_max){z_max=lambda + protocol_constant;} \
    if(lambda + protocol_constant < z_min){z_min=lambda + protocol_constant;} \
  }									\
  
void set_sang_protocol(double *prot, double A, double B,
		       double z_final, double t_final, double h,
		       int N, double sigma);


// MAX ENTROPY FUNCTIONS;
void get_signal(double * signal, int len, double z_start, double z_final,
		int traces, int samples, double * lambda_ar, double K,
		double ** x_ar, int c, double * bin_centers, double beta);
void fill_transfer_matrix(double ** H, int len1, int len2, double beta, double K,
			  double * bin_centers, double * bin_centers_reduced);
void right_matrix_product(double ** A, double *x, double * y,  int len1, int len2);
void vector_addition(double * x, double *y, double *z, int len);
void matrix_product(double ** C, double ** A, double ** B, int len1, int len2);
void matrix_product2(double ** C, double ** A, int len1, int len2);
void update_A(double ** A, double ** H, int len, int N_recon, double prior_sigma, double * p, double alpha);
int step_A(double ** A, double ** H, int N_recon, double * P, double alpha, double * u,
	   double * v, double * f, double * s, double * precond_grad, double * u_old, int len,
	   double prior_sigma);

int step_B(double * S_u, double * F_u, double *F_uu, double * v, double sigma_prior, int J, int I,
	   double *u, double *s, double * f, double ** H);
int step_C(double * xi_0, double F_u, double F_uu, double S_u,
	   double epsilon, double * P, double *u, int J, double epsilon_prime, double alpha);
void step_D_init(double * w, double ** H, double *xi, int J, int I, double * v);
int step_D(double * w, double ** H, double *u, double ** A, double * xi,int J, int I,
	   double * G_uu, double * F_u, double * S_u, double * alpha, double * F, double * s,
	   double * M, double * f, double sigma_prior, double xi_threshold, double * P, double F_uu);
int step_E(double xi, double * u, int J, double P_threshold, int iter_count);
void max_ent(double * signal, int len, double z_start, double z_final,
	     int traces, int samples, double * lambda_ar, double K,
	     double ** x_ar, int c, double * bin_centers, int cutoff,
	     double * G, double beta, double * bin_centers_reduced, double prior_sigma,
	     int MAX_ITER,double U, double sigma, double * U_lc_surf);

void weierstrass_deconv(double * signal, int len, double z_start, double z_final,
			int traces, int samples, double * lambda_ar, double K,
			double ** x_ar, int c, double * bin_centers, double beta,
			double * lc_bins, double * G_lc, double * U_lc_surf, double U,double sigma);

void get_weierstrass_signal(double * signal, double * A_prim, double * A_primprim, int len, double z_start, double z_final,
			    int traces, int samples, double * lambda_ar, double K,
			    double ** x_ar, int c, double * bin_centers, double beta);


#endif
