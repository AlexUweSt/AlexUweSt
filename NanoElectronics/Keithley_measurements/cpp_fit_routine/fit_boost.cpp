
/*
Vectorization headers and libs for further acceleration of execution
*/

#if 0
#pragma GCC option("arch=native","tune=native","no-zero-upper") //Enable AVX
#pragma GCC target("avx2")  //Enable AVX
#include <x86intrin.h> //AVX/SSE Extensions
#endif

#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlinear.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_fit.h>
//#include <gsl/gsl_ieee_utils.h>

#include <boost/filesystem.hpp>

#include <fstream>
#include <iostream>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/thread.hpp>

#include <vector>
#include <array>
#include <string>
#include <utility>
#include<cmath>
#include <iomanip>
#include <numeric>
//#include <experimental/filesystem>
//#include <thread>
#include <chrono>  
#include <valarray>

#include <unistd.h>
#include <atomic> 
#include <functional>
#include <mutex>
#include <condition_variable> 
#include <time.h>
#include <chrono>
#include <algorithm>







using namespace std;
using namespace boost::iostreams;

#define N      800    /* number of data points to fit */

#define TMAX   (40.0) /* time variable in [0,TMAX] */
const unsigned short int no_of_robust_iterations = 5;
int max_V=800;/*max of applied voltage in mV*/
//WorkerPool variables
std::mutex Queue_Mutex;
std::condition_variable condition;
std::vector<std::function<void()>> Queue;
std::vector<boost::thread> Pool;



//atomic constants for parameter storage

double GOF_array_d  [10000];
double E0_d         [10000];
double Gamma1_d     [10000];
double Gamma2_d     [10000];
double GOF_array_u  [10000];
double E0_u         [10000];
double Gamma1_u     [10000];
double Gamma2_u     [10000];

//initialize previous arrays
//std::fill(E0_d, E0_d+10000, 0)
std::vector<std::array<double,800>> weivec_d(10000);
std::vector<std::array<double,800>> tvec_d(10000);
std::vector<std::array<double,800>> yvec_d(10000);
std::vector<std::array<double,800>> weivec_u(10000);
std::vector<std::array<double,800>> tvec_u(10000);
std::vector<std::array<double,800>> yvec_u(10000);
std::vector<std::array<double,800>> smvec_d(10000);
std::vector<std::array<double,800>> smvec_u(10000);
std::vector<bool>wasthere(10000);

std::atomic<unsigned int> progress (0);
 
// Class for dealing with parallelization execution of threads on workers. It manages creation, assigning and destruction.

class WorkerPool
{
  public:
      //int Num_Threads = thread::hardware_concurrency();
      int Num_Threads = boost::thread::hardware_concurrency();
      static int Active_Workers;
      bool StopSignal=0;

      
      WorkerPool()
      {      
          cout<<"Number of available system threads "<<Num_Threads<<endl;
          
      }

      bool StartPool(unsigned int no_of_workers)
      {   
          
          if (no_of_workers>(Num_Threads-Active_Workers)){return 0;}
          init_workers=no_of_workers;
          Active_Workers=Active_Workers+no_of_workers; 
          cout<<"Number of active threads "<<Active_Workers<<endl;

          for(int ii = 0; ii < no_of_workers; ii++)
          {
              //cout<<ii<<endl;
              Pool.push_back(std::move(boost::thread([this](){Infinite_loop_function();})));   
          }

          std::cout<<"Threads are working ";
          
      return(0);
      
          

      }

      bool StartPool(void)
      {   
          
          //if (no_of_workers>(Num_Threads-Active_Workers)){return 0;}
          init_workers=Num_Threads-2;
          Active_Workers=Active_Workers+init_workers; 
          cout<<"Number of active threads "<<Active_Workers<<endl;

          for(int ii = 0; ii < init_workers; ii++)
          {
              //cout<<ii<<endl;
              Pool.push_back(std::move(boost::thread([this](){Infinite_loop_function();})));   
          }

          std::cout<<"Threads are working ";
          
      return(0);
      
          

      }

      bool StopPool()
      {
          StopSignal=1;
          condition.notify_all();
          //// do not use init_workers
            for(int ii = 0; ii < init_workers; ii++)
          {
              //cout<<"inlop"<<ii<<endl;
              Pool.at(ii).join();
          }
          cout<<Active_Workers<<endl;
          WorkerPool::Active_Workers=0;
      }
      
      void Infinite_loop_function()
      {
      function<void()> Job;
      //cout<<"incool\n";
      while(true && !(StopSignal))
      {
          //std::this_thread::sleep_for(std::chrono::seconds(1));
          {
              unique_lock<mutex> locker(Queue_Mutex);
              //condition.wait(locker);
              //locker.lock();
              //cout<<"Thread ID "<<pthread_self()<<"\n";
              condition.wait(locker,[&]{return (!Queue.empty()|| StopSignal);});
              if (StopSignal){break;}
              //cout<<"break"<<endl;
              Job = Queue.back();
              Queue.pop_back();
            // locker.unlock();
            //cout<<"Step1";
          }
          //cout<<"Step2";
          
          Job(); // function<void()> type
      }
      
      }
      
      void Add_Job(function<void()> New_Job)
      {
          {
              unique_lock<mutex> lock(Queue_Mutex);
              Queue.push_back(New_Job);
              //cout<<"Number of Jobs "<<Queue.size()<<"\n";
          }
          condition.notify_one();
      }

      
      ~WorkerPool()
      {
      }
  private:
  int init_workers=0;

};

//Class for reading the data from the file and splitting into 4 arrays(I,V(down, up direction))

class MeasFile
{
  public:
      std::vector<std::string> m_header;
      std::array<std::vector<double>, 3> m_data;
      
      std::array<double,800> data_downsweep_V;
      std::array<double,800> data_downsweep_I;
      std::array<double,800> data_upsweep_V;
      std::array<double,800> data_upsweep_I;
      int sweeppoints=0;/////
      MeasFile(std::istream &filecont)
      {
        
        m_header.reserve(32);
          while(filecont)
          {
            auto a = filecont.peek();
            if(a != '#')
              break;
            std::string transfer;
              getline(filecont, transfer);
              m_header.push_back(std::move(transfer));
          }

          for(auto& a : m_data)
            a.reserve(2048);
          //std::string V1, I1, t1;
          std::string temp_string;
          double V,I,t;
          std::string::size_type sz; 
          //while(filecont >> V >> I >> t)
          const std::string WHITESPACE = " \n\r\t\f\v";
          size_t start ;
          //for (int i=0;i<2401;i++)
          while(getline(filecont, temp_string))
          {
            //getline(filecont, temp_string);
            //cout<<temp_string<<endl;
            std::size_t found = temp_string.find("nan");
            if (found!=std::string::npos)
            {
            std::cout << "first 'nan' found at: " << found << '\n';
            m_data[0].push_back(0);
            m_data[1].push_back(0);
            m_data[2].push_back(0);
            }
            else 
            {
            //temp_string >> V1 >> I1 >> t1;
            V=std::stod(temp_string,&sz);
            //cout<<"broj je"<<i <<" "<<temp_string<<endl;
            temp_string=temp_string.substr(sz);
            //start = temp_string.find_first_not_of(WHITESPACE);
            //temp_string=(start == std::string::npos) ? "" : temp_string.substr(start);
            //cout<<temp_string<<endl;
            I=std::stod(temp_string,&sz);
            temp_string=temp_string.substr(sz);
            t=std::stod(temp_string,&sz);
            {m_data[0].push_back(V);} 
            {m_data[1].push_back(I);} 
            {m_data[2].push_back(t);} 
            //cout<<V<<"   "<<I<<"   "<<t<<"   "<<endl;
            }
          }

          
        sweeppoints = (m_data.at(0)).size();
        //cout<<"sweep"<<sweeppoints;F
        const int sweeppart1a = ((sweeppoints-1)/6 + 1);
        const int sweeppart1b = ((sweeppoints-1)/2 + 1);
        const int sweeppart2a = ((sweeppoints-1)/2 + 1);
        const int sweeppart2b = (sweeppoints-(sweeppoints-1)/6);
        
        std::copy(&m_data[0].at(0)+sweeppart1a,&m_data[0].at(0)+sweeppart1b,data_downsweep_V.begin()); //more elegant way?
        std::copy(&m_data[1].at(0)+sweeppart1a,&m_data[1].at(0)+sweeppart1b,data_downsweep_I.begin());
        std::copy(&m_data[0].at(0)+sweeppart2a,&m_data[0].at(0)+sweeppart2b,data_upsweep_V.begin());
        std::copy(&m_data[1].at(0)+sweeppart2a,&m_data[1].at(0)+sweeppart2b,data_upsweep_I.begin());
        /*
        ofstream report;
        auto temp_file12="/home/filipk/Desktop/check.txt";
        report.open (temp_file12);
        //int curve_no=0;
        //report<<"Curve No."<< <<endl; instead use header...
        /*
        std::vector<double> I_fit_n(2000);
        I_fit_n.resize(number_of_points);
        std::copy(y+int((n-number_of_points)/2),y+n-int((n-number_of_points)/2),I_fit_n.begin());
        */
        /*
        for(int i=0;i<800;i++)
        {
        report<<fixed << setw( 11 )<<internal<< setprecision( 15 ) << data_downsweep_V.at(i)<<"       "<< setprecision( 15 ) << data_downsweep_I.at(i)<<endl;
        
        }
        report.close();
        */


      }
      ~MeasFile()
      {
      }
      void PrintHeader()
      {
      for(int i=0;i<32;i++)
      {
        cout<<m_header[i]<<'\n';
      }   

      }
  private:
};






















struct data {
  size_t n;
  double * t;
  double * y;
};

struct int_p{double E0; double G1; double G2; double V;};
struct fit_p{double E0; double G1; double G2;};

double Transmission4(double E,double E0_pn, double Gamma1_pn, double Gamma2_pn, double V)
{
  //Takes care so the parameters dont go into negative side....
double E0=std::abs(E0_pn);
double Gamma1=abs(Gamma1_pn);
double Gamma2=abs(Gamma2_pn);
//
double SumGammaSq=(Gamma1 + Gamma2)*(Gamma1 + Gamma2);
double E0split=(E - (E0 + ((Gamma1 - Gamma2) / (Gamma1 + Gamma2)) * copysign(1.0, V)*sqrt(abs(V))/(2)))*(E - (E0 + ((Gamma1 - Gamma2) / (Gamma1 + Gamma2)) * copysign(1.0, V)*sqrt(abs(V))/2));
double Transmission=(4* Gamma1*Gamma2)/(E0split + SumGammaSq);
return Transmission;
}

double Transmission5(double E,double E0_pn, double Gamma1_pn, double Gamma2_pn, double V)
{
  //Takes care so the parameters dont go into negative side....
double E0=std::abs(E0_pn);
double Gamma1=abs(Gamma1_pn);
double Gamma2=abs(Gamma2_pn);
//
double SumGammaSq=(Gamma1 + Gamma2)*(Gamma1 + Gamma2);
double E0split=(E - (E0 + ((Gamma1 - Gamma2) / (Gamma1 + Gamma2)) * V/2))*(E - (E0 + ((Gamma1 - Gamma2) / (Gamma1 + Gamma2)) * V/2));
double Transmission=(4* Gamma1*Gamma2)/(E0split + SumGammaSq);
return Transmission;
}


double Fermi(double E,double mu)
{
double b=40;
double Fermi = 1 / (1 + exp((E - mu) * b));
return Fermi;
}


double Transmission3(double E,double E0_pn, double Gamma1_pn, double Gamma2_pn, double V)
{
  //Takes care so the parameters dont go into negative side....
double E0=std::abs(E0_pn);
double Gamma1=abs(Gamma1_pn);
double Gamma2=abs(Gamma2_pn);
//
double SumGammaSq=pow(Gamma1,2);
double E0split=pow((E - (E0 + Gamma2 * abs(V*V*V))),2);
double Transmission=(pow(Gamma1,2))/(E0split + SumGammaSq);
return Transmission;
}


double Transmission(double E,double E0_pn, double Gamma1_pn, double Gamma2_pn, double V)
{
  //Takes care so the parameters dont go into negative side....2
double E0=std::abs(E0_pn);
double Gamma1=abs(Gamma1_pn);
double Gamma2=abs(Gamma2_pn);
//
double SumGammaSq=pow(Gamma1,2);
double E0split=pow((E - (E0 + Gamma2 * V)),2);
double Transmission=(pow(Gamma1,2))/(E0split + SumGammaSq);
return Transmission;
}


double TransmissionSym(double E,double E0_pn, double Gamma1_pn, double Gamma2_pn, double V)
{
  //Takes care so the parameters dont go into negative side....
double E0=std::abs(E0_pn);
double Gamma1=abs(Gamma1_pn);
double Gamma2=abs(Gamma2_pn);
//
double SumGammaSq=(Gamma1 + Gamma1)*(Gamma1 + Gamma1);
double E0split=(E - E0)*(E - E0);
double Transmission=(4* Gamma1*Gamma1)/(E0split + SumGammaSq);
return Transmission;
}


double Integrand(double x, void * params)
{
  struct int_p * parameters = (struct int_p*)params;
  double E0 = (parameters->E0);
  double G1 = (parameters->G1);
  double G2 = (parameters->G2);
  double V = (parameters->V);
  double Integrand = Transmission5(x,E0,G1,G2,V)*(Fermi(x,V/2)-Fermi(x,-V/2)); 
  return Integrand;
}




double Integrate(double const E0,double const G1, double const G2, double const V)
{
  gsl_integration_cquad_workspace * w
    = gsl_integration_cquad_workspace_alloc (650);//was 650

  double result, error;
  size_t num;
  struct int_p alpha = {E0,G1,G2,V};

  gsl_function F;
  F.function = &Integrand;
  F.params = &alpha;

  gsl_integration_cquad (&F, -2, 2, 0, 1e-9,
                        w, &result, &error,&num);  
  gsl_integration_cquad_workspace_free (w);
  return result;

}


//An attempt to implement Jacobi matrix calculation

double dTddE(double const& E0, double const& G1, double const& G2, double const& V, double const& E){
  double SumGammaSq = pow(G1 + G2, 2);
  double E0split = (E - (E0 + ((G1 - G2) / (G1 + G2)) * V/2));
  double derivative = -(4* G1*G2)/pow(pow(E0split, 2) + SumGammaSq,2)*2*E0split;
  return derivative;
}

double ddEdE0(double const& E0, double const& G1, double const& G2, double const& V, double const& E){
  return -1.0;
}

double dTdG1(double const& E0, double const& G1, double const& G2, double const& V, double const& E){
  double SumGammaSq = pow((G1 + G2), 2);
  double E0split = pow(E-(E0 + ((G1 - G2)/(G1 + G2))*V/2), 2);
  double derivative = 4*G2/(E0split + SumGammaSq);
  return derivative;
}

double dTdG2(double const& E0, double const& G1, double const& G2, double const& V, double const& E){
  double SumGammaSq = pow((G1 + G2), 2);
  double E0split = pow(E-(E0 + ((G1 - G2)/(G1 + G2))*V/2), 2);
  double derivative = 4*G1/(E0split + SumGammaSq);
  return derivative;
}

double dTdG(double const& E0, double const& G1, double const& G2, double const& V, double const& E){
  double SumGamma = G1 + G2;
  double E0split = pow(E-(E0 + ((G1 - G2)/SumGamma)*V/2), 2);
  double derivative = -4*G1*G2/pow(E0split + pow(SumGamma, 2), 2) * 2*SumGamma;
  return derivative;
}

double dGdG1(double const& E0, double const& G1, double const& G2, double const& V, double const& E){
  return 1.0;
}

double dGdG2(double const& E0, double const& G1, double const& G2, double const& V, double const& E){
  return 1.0;
}

double ddEdG1(double const& E0, double const& G1, double const& G2, double const& V, double const& E){
  double SumGamma = G1 + G2;
  return -G2*V/pow(SumGamma,2);
}
double ddEdG2(double const& E0, double const& G1, double const& G2, double const& V, double const& E){
  double SumGamma = G1 + G2;
  return G1*V/pow(SumGamma, 2);
}
double IntegrateGrad(double (*IntegrandGrad)(double, void*), double const E0,double const G1, double const G2, double const V)
{
  gsl_integration_cquad_workspace * w
    = gsl_integration_cquad_workspace_alloc (650);//was 650

  double result, error;
  size_t num;
  struct int_p alpha = {E0,G1,G2,V};

  gsl_function F;
  F.function = IntegrandGrad;
  F.params = &alpha;

  gsl_integration_cquad (&F, -2, 2, 0, 1e-9,
                        w, &result, &error,&num);  
  gsl_integration_cquad_workspace_free (w);
  return result;

}
std::array <double, 3> current_gradient(double const& E0, double const& G1, double const& G2, double const& V){
  auto full_dTdE0 = [](double x, void* params){struct int_p * parameters = (struct int_p*)params;
                                              double E0 = (parameters->E0);
                                              double G1 = (parameters->G1);
                                              double G2 = (parameters->G2);
                                              double V = (parameters->V);
                                              return (dTddE(E0, G1, G2, V, x) * ddEdE0(E0, G1, G2, V, x))*(Fermi(x, V/2)-Fermi(x, -V/2));};
  auto full_dTdG1 = [](double x, void* params){struct int_p * parameters = (struct int_p*)params;
                                              double E0 = (parameters->E0);
                                              double G1 = (parameters->G1);
                                              double G2 = (parameters->G2);
                                              double V = (parameters->V);
                                              return (dTddE(E0, G1, G2, V, x) * ddEdG1(E0, G1, G2, V, x) +
                                                      dTdG1(E0, G1, G2, V, x) +
                                                      dTdG(E0, G1, G2, V, x)*dGdG1(E0, G1, G2, V, x))*(Fermi(x, V/2)-Fermi(x, -V/2));};
  auto full_dTdG2 = [](double x, void* params){struct int_p * parameters = (struct int_p*)params;
                                              double E0 = (parameters->E0);
                                              double G1 = (parameters->G1);
                                              double G2 = (parameters->G2);
                                              double V = (parameters->V);
                                              return (dTddE(E0, G1, G2, V, x) * ddEdG2(E0, G1, G2, V, x) +
                                                      dTdG2(E0, G1, G2, V, x) +
                                                      dTdG(E0, G1, G2, V, x)*dGdG2(E0, G1, G2, V, x))*(Fermi(x, V/2)-Fermi(x, -V/2));};
  auto full_dIdE0 = IntegrateGrad(full_dTdE0, E0, G1, G2, V) * float(7.74809174e-5);
  auto full_dIdG1 = IntegrateGrad(full_dTdG1, E0, G1, G2, V) * float(7.74809174e-5);
  auto full_dIdG2 = IntegrateGrad(full_dTdG2, E0, G1, G2, V) * float(7.74809174e-5);
  return std::array<double, 3> {full_dIdE0, full_dIdG1, full_dIdG2};
}


int current_df (const gsl_vector * x, void *data, gsl_matrix * J){
  size_t n = ((struct data *)data)->n;
  double *t = ((struct data *)data)->t;
  double *y = ((struct data*)data)->y;

  double E0 = gsl_vector_get (x, 0);
  double G1 = gsl_vector_get (x, 1);
  double G2 = gsl_vector_get (x, 2);

  size_t i;

  for (i = 0; i < n; i++)
    {
      auto derivatives = current_gradient(E0, G1, G2, t[i]); //t[i] is a bias
      gsl_matrix_set (J, i, 0, derivatives[0]);
      gsl_matrix_set (J, i, 1, derivatives[1]);
      gsl_matrix_set (J, i, 2, derivatives[2]);
    }

  return GSL_SUCCESS;
}

//Jacobi calculation ends







int current(const gsl_vector * x, void *data,
        gsl_vector * f)
{
  size_t n = ((struct data *)data)->n;
  double *t = ((struct data *)data)->t;
  double *y = ((struct data *)data)->y;

  double E0 = gsl_vector_get (x, 0);
  double G1 = gsl_vector_get (x, 1);
  double G2 = gsl_vector_get (x, 2);
  size_t i;

  for (i = 0; i < n; i++)
      {
        double Yi= Integrate(E0,G1,G2,t[i])*7.74809174e-5;
        gsl_vector_set (f, i, Yi - y[i]);
      }

    return GSL_SUCCESS;

}

std::array<double,800> calc_current(double const E0,double const  G1,double const  G2, double x[800])
{
std::array<double,800> res;
  for(int i =0; i<800;i++)  {
      res.at(i)=Integrate(E0,G1,G2,x[i])*7.74809174e-5;
  }
return res;
}



void
callback(const size_t iter, void *params,
         const gsl_multifit_nlinear_workspace *w)
{
  gsl_vector *f = gsl_multifit_nlinear_residual(w);
  gsl_vector *x = gsl_multifit_nlinear_position(w);
  double rcond;

  /* compute reciprocal condition number of J(x) */
  gsl_multifit_nlinear_rcond(&rcond, w);
/*
  fprintf(stderr, "iter %2zu: A = %.4f, lambda = %.4f, b = %.4f, cond(J) = %8.4f, |f(x)| = %.4f\n",
          iter,
          gsl_vector_get(x, 0),
          gsl_vector_get(x, 1),
          gsl_vector_get(x, 2),
          1.0 / rcond,
          gsl_blas_dnrm2(f));
 */         
}



double r_weight(double x,double xi,double d)
{
  double w=0;
  
  w=(1-    abs( (x-xi)/d ) * abs( (x-xi)/d ) * abs( (x-xi)/d )       );
  w=w*w*w;
  return w;
}

double median(std::array<double,800> residuals)
{
  sort(residuals.begin(),residuals.end());
  int n=residuals.size();
  if (n % 2 == 0) 
      return (double)(residuals[n/2-1] + residuals[n/2])/2.0;
       
      
  return (double)residuals[(n-1)/2]; 

}

std::array<double,800> 
smooth(std::array<double,800> x,std::array<double,800> y,int number_of_points,int half_window_size)
{
  std::array<double,800> new_y;
  std::array<double,800> regression_weight; /// how to dynamically allocate?
  std::array<double,800> robust_weights;
  std::array<double,800> weight_product;
  robust_weights.fill(1);
  double MAD=0;

/*
  ofstream myfile1, myfile2;
  myfile1.open ("values_X.txt");
  myfile2.open("values_Y_old.txt");
  for(int i=0;i<number_of_points;i++)
  {
  myfile1<<fixed << setw( 11 ) << setprecision( 15 ) << x.at(i)<< "\n";
  myfile2<<fixed << setw( 11 ) << setprecision( 15 ) << y.at(i)<< "\n";
  }
  myfile1.close();
  myfile2.close();
*/
  double c0, c1, cov00, cov01, cov11, chisq;
      int k;
      int n = 2*half_window_size+1;
  /*
  for (int i=0;i<number_of_points;i++)
  {
    cout<<y.at(i)<<endl;
  }
*/


  for(int iteration =0;iteration<no_of_robust_iterations;iteration++)//5
  {
  //regular part

  //generate regression weights
    for(int i=0;i<2*half_window_size+1;i++)
    {   
      regression_weight.at(i)=r_weight(x.at(i),x.at(half_window_size),x.at(0)-x.at(half_window_size));
      //printf("Weight is %.10f\n",regression_weight.at(i));
    }
/*
     for (int i=0;i<number_of_points;i++)
  {
    cout<<"tezina"<<regression_weight.at(i)<<endl;
  } 
*/

  for(int i=half_window_size;i<number_of_points-half_window_size;i++)
  {
      for(int z=0;z<2*half_window_size+1;z++)
      {
        weight_product.at(z)=regression_weight.at(z)*robust_weights.at(i+z-half_window_size);
       // cout<<"producat"<<weight_product.at(z)<<endl;
      }
      gsl_fit_wlinear (&x.at(i-half_window_size), 1, weight_product.begin(), 1, &y.at(i-half_window_size), 1, n, &c0, &c1, &cov00, &cov01, &cov11, &chisq);  
      if(std::isnan(c0)||std::isnan(c1)){cout<<"Error in smoothing procedure"<<endl;}    
      new_y.at(i)=c0+c1*x.at(i);
  }



  for(int i=0;i<half_window_size;i++)
  {
    for(int z=0;z<half_window_size;z++)
      {
        weight_product.at(z)=regression_weight.at(z+half_window_size-i)*robust_weights.at(z);
      }
    gsl_fit_wlinear (&x.at(0), 1, regression_weight.begin(), 1, &y.at(0), 1, n, &c0, &c1, &cov00, &cov01, &cov11, &chisq);
    new_y.at(i)=c0+c1*x.at(i);


  }




  for(int i=number_of_points-half_window_size;i<number_of_points;i++)
  {
    for(int z=0;z<half_window_size;z++)
      {
        weight_product.at(z)=regression_weight.at(z+half_window_size+   -i+(number_of_points-half_window_size))*robust_weights.at(z+number_of_points-half_window_size);
      }
    gsl_fit_wlinear (&x.at(number_of_points-2*half_window_size-1), 1, regression_weight.begin(), 1, &y.at(number_of_points-2*half_window_size-1), 1, n, &c0, &c1, &cov00, &cov01, &cov11, &chisq);
    if(std::isnan(c0)||std::isnan(c1)){cout<<"Error in smoothing procedure"<<endl;}
    new_y.at(i)=c0+c1*x.at(i);





  }
  /*
  ofstream report;
  auto temp_file="/home/filipk/Desktop/"+std::to_string(iteration)+".txt";
  report.open (temp_file);
  //int curve_no=0;
  //report<<"Curve No."<< <<endl; instead use header...
  /*
  std::vector<double> I_fit_n(2000);
  I_fit_n.resize(number_of_points);
  std::copy(y+int((n-number_of_points)/2),y+n-int((n-number_of_points)/2),I_fit_n.begin());
  */
  /*
  for(int i=0;i<800;i++)
  {
  report<<fixed << setw( 11 )<<internal<< setprecision( 15 ) << weight_product.at(i)<<"       "<< setprecision( 15 ) << new_y.at(i)<<"       "<<x.at(i)<<"       "<<y.at(i)<<endl;
  
  }
  report.close();
  */




  if (iteration < no_of_robust_iterations)
  {

    std::array<double,800> abs_residuals;
    for(int i=0;i<number_of_points;i++)
    {
      abs_residuals.at(i)=abs(y.at(i)-new_y.at(i));
    }


    MAD=median(abs_residuals);
    //if (MAD==0) {cout<<"Not good"<<endl;}
    for(int i=0;i<number_of_points;i++)
    {
      if (abs_residuals.at(i)<(6*MAD))
        {
          robust_weights.at(i)=(1-(abs_residuals.at(i)/(6*MAD))*(abs_residuals.at(i)/(6*MAD)) )*(1-(abs_residuals.at(i)/(6*MAD))*(abs_residuals.at(i)/(6*MAD)) );
        }
      else
        {
          robust_weights.at(i)=1e-15;//1e-15;    
        }
      }
      //copy(new_y.begin(),new_y.end(),y.begin());
    }
    std::copy(new_y.begin(),new_y.end(),y.begin());
  }



/*
  ofstream myfile3;
  myfile3.open ("values_Y_new.txt");

  for(int i=0;i<number_of_points;i++)
  {
  myfile3<<fixed << setw( 11 ) << setprecision( 15 ) << new_y.at(i)<< "\n";
  }
  myfile3.close();
*/


  return new_y;
}


int fit(std::string input_file, std::string output_file,int ind,int fitting_range)
{
  
  /* This part should be properly rewritten. We are using workaround by defining weights outside range as 0, without cutting arrays.
  This part of the code can work much faster*/
  int number_of_points=(fitting_range*N) /max_V;
  //cout<<"number_of_points"<<number_of_points;
  ifstream file(input_file, ios_base::in | ios_base::binary);
  filtering_streambuf<input> in;
  in.push(gzip_decompressor());
  in.push(file);  
  //std::vector<int> v={1}; !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  std::istream incoming(&in);

  MeasFile file1(incoming);
  //file1.PrintHeader();
  #if 1

  for(int downup=0;downup<2;downup++)
  {


    const gsl_multifit_nlinear_type *T = gsl_multifit_nlinear_trust;
    gsl_multifit_nlinear_workspace *w;
    gsl_multifit_nlinear_fdf fdf;
    gsl_multifit_nlinear_parameters fdf_params =
          gsl_multifit_nlinear_default_parameters();
    const size_t n = N; //default
    //const size_t n = number_of_points; 
    const size_t p = 3;
    const size_t n_n = number_of_points; 


    ///changing defaults
    //fdf_params.h_df=1;
    //fdf_params.scale=gsl_multifit_nlinear_scale_levenberg;
    fdf_params.solver=gsl_multifit_nlinear_solver_svd;
    fdf_params.fdtype=GSL_MULTIFIT_NLINEAR_CTRDIFF;
  // fdf_params.trs=gsl_multifit_nlinear_trs_dogleg;
    fdf_params.trs = gsl_multifit_nlinear_trs_lmaccel;
    //fdf_params.h_df=1e-18;
    
    
    
    gsl_vector *f;
    gsl_matrix *J;
    gsl_matrix *covar = gsl_matrix_alloc (p, p);
    double t[N], y[N], weights[N];

    std::vector<double> t_n(2000),y_n(2000),weights_n(2000),I_fit_n(2000),error1(2000),error2(2000),error3(2000),error4(2000),sm_n(2000);
    t_n.resize(number_of_points,0);
    y_n.resize(number_of_points,0);
    weights_n.resize(number_of_points,0);
    I_fit_n.resize(number_of_points,0);

    struct data d = { n_n, t_n.data(), y_n.data() };
    ////reuse old parameters in a new fit
    double x_init[3] = { 0.5, 5e-3, 5e-3};
    
    /*

      if (downup==0)
      {
      /* x_init = { E0_d[ind], Gamma1_d[ind], Gamma2_d[ind]};  starting values 
      */
      /*
      x_init[0]=E0_d[ind];
      x_init[1]=Gamma1_d[ind];
      x_init[2]=Gamma2_d[ind];
      }
      else
      {
      x_init[0]=E0_u[ind];
      x_init[1]=Gamma1_u[ind];
      x_init[2]=Gamma2_u[ind];
      }
      */
    
    ////
    gsl_vector_view x = gsl_vector_view_array (x_init, p);
    gsl_vector_view wts = gsl_vector_view_array(weights_n.data(), n_n);
    gsl_rng * r;
    double chisq, chisq0;
    int status, info;
    //size_t i;

    const double xtol = 1e-18;
    const double gtol = 1e-18;
    const double ftol = 1e-18;

    gsl_rng_env_setup();
    r = gsl_rng_alloc(gsl_rng_default);

    /* define the function to be minimized */
    fdf.f = current;
    fdf.df =NULL; //current_df;
    //fdf.df = NULL;   /* set to NULL for finite-difference Jacobian */
    fdf.fvv = NULL;     /* not using geodesic acceleration */
    fdf.n = n_n;
    fdf.p = p;
    fdf.params = &d;
  
    
    std::array<double,800> sm;
    double minel;

    double max_el=0;
    std::array<double,800> CurveSmDiff;

    if (wasthere.at(ind)==0)
    {
      if(downup==0)
        {
        std::copy(file1.data_downsweep_V.begin(),file1.data_downsweep_V.end(),tvec_d.at(ind).begin());//begin(I)
        std::copy(file1.data_downsweep_I.begin(),file1.data_downsweep_I.end(),yvec_d.at(ind).begin());
        //std::copy(file1.data_downsweep_V.begin(),file1.data_downsweep_V.end(),t);//begin(I)
        //std::copy(file1.data_downsweep_I.begin(),file1.data_downsweep_I.end(),y);
        
        sm=smooth(file1.data_downsweep_V,file1.data_downsweep_I,800,20);
        std::copy(sm.begin(),sm.end(),(smvec_d.at(ind)).begin());
        for(int i =0; i<n;i++)
        {
          //CurveSmDiff.at(i) = abs(file1.data_downsweep_I.at(i) - sm.at(i));
          CurveSmDiff.at(i) = abs(y[i] - sm.at(i));
        }
        max_el=*max_element(CurveSmDiff.begin(),CurveSmDiff.end());



        std::fill(weights, weights + n, 0); /*initialize weights as 0*/
        for (int i = 0; i < 800; i++) /*defines region of interest */
        {
          weights[i] =(max_el-CurveSmDiff.at(i))/max_el;
        };
        std::copy(weights,weights+800,weivec_d.at(ind).begin());

        }
      else
        {
        std::copy(file1.data_upsweep_V.begin(),file1.data_upsweep_V.end(),tvec_u.at(ind).begin());
        std::copy(file1.data_upsweep_I.begin(),file1.data_upsweep_I.end(),yvec_u.at(ind).begin()); 
        //std::copy(file1.data_upsweep_V.begin(),file1.data_upsweep_V.end(),t);
        //std::copy(file1.data_upsweep_I.begin(),file1.data_upsweep_I.end(),y); 
        
        sm=smooth(file1.data_upsweep_V,file1.data_upsweep_I,800,40);
        std::copy(sm.begin(),sm.end(),smvec_u.at(ind).begin());
        for(int i =0; i<n;i++)
        {
          //CurveSmDiff.at(i) = abs(file1.data_downsweep_I.at(i) - sm.at(i));
          CurveSmDiff.at(i) = abs(y[i] - sm.at(i));
        }
        max_el=*max_element(CurveSmDiff.begin(),CurveSmDiff.end());



        std::fill(weights, weights + n, 0); /*initialize weights as 0*/
        for (int i = 0; i < 800; i++) /*defines region of interest */
        {
          weights[i] =(max_el-CurveSmDiff.at(i))/max_el;
        };
        std::copy(weights,weights+800,weivec_u.at(ind).begin());
        
        }
    }
    

    if(downup==0)
      {
        std::copy(tvec_d.at(ind).begin(),tvec_d.at(ind).end(),t);
        std::copy(yvec_d.at(ind).begin(),yvec_d.at(ind).end(),y);
        std::copy(weivec_d.at(ind).begin(),weivec_d.at(ind).end(),weights);
        std::copy(smvec_d.at(ind).begin(),smvec_d.at(ind).end(),sm.begin());
      }
    else
      {
        std::copy(tvec_u.at(ind).begin(),tvec_u.at(ind).end(),t);
        std::copy(yvec_u.at(ind).begin(),yvec_u.at(ind).end(),y);
        std::copy(weivec_u.at(ind).begin(),weivec_u.at(ind).end(),weights);
        std::copy(smvec_u.at(ind).begin(),smvec_u.at(ind).end(),sm.begin());
      }
    
    
    /*
    Va=abs(V);
    minel=I[std::min_element(begin(Va),end(Va))-begin(Va)];
    cout<<"File "<<output_file<<minel<<endl;
    for(int i=0;i<800;i++)
    {
      I[i]=I[i]-minel;
    }
    std::copy(std::begin(I),std::end(I),t);


    std::copy(std::begin(V),std::end(V),y);
    */

  
    
    std::copy(weights+int((n-number_of_points)/2),weights+n-int((n-number_of_points)/2),weights_n.begin());
    std::copy(t+int((n-number_of_points)/2),t+n-int((n-number_of_points)/2),t_n.begin());
    std::copy(y+int((n-number_of_points)/2),y+n-int((n-number_of_points)/2),y_n.begin());
    std::copy(sm.begin()+int((n-number_of_points)/2),sm.begin()+n-int((n-number_of_points)/2),sm_n.begin());


    /* allocate workspace with default parameters */
    w = gsl_multifit_nlinear_alloc (T, &fdf_params, n_n, p);

    /* initialize solver with starting point and weights */
    gsl_multifit_nlinear_winit (&x.vector, &wts.vector, &fdf, w);

    /* compute initial cost function */
    f = gsl_multifit_nlinear_residual(w);
    gsl_blas_ddot(f, f, &chisq0);

    
    status = gsl_multifit_nlinear_driver(10000, xtol, gtol, ftol, //was 1000 iteration // last value was 200
                                        callback, NULL, &info, w);

    /* compute covariance of best fit parameters */
    J = gsl_multifit_nlinear_jac(w);
    gsl_multifit_nlinear_covar (J, 0, covar);

    /* compute final cost */
    gsl_blas_ddot(f, f, &chisq);

    #define FIT(i) gsl_vector_get(w->x, i)
    #define ERR(i) sqrt(gsl_matrix_get(covar,i,i))

    //fprintf(stderr, "summary from method '%s/%s'\n",
    //        gsl_multifit_nlinear_name(w),
    //        gsl_multifit_nlinear_trs_name(w));
    //fprintf(stderr, "number of iterations: %zu\n",
           //gsl_multifit_nlinear_niter(w);
    //fprintf(stderr, "function evaluations: %zu\n", fdf.nevalf);
    //fprintf(stderr, "Jacobian evaluations: %zu\n", fdf.nevaldf);
    //fprintf(stderr, "reason for stopping: %s\n",
    //        (info == 1) ? "small step size" : "small gradient");
    //fprintf(stderr, "initial |f(x)| = %f\n", sqrt(chisq0));
    //fprintf(stderr, "final   |f(x)| = %f\n", sqrt(chisq));

    {
      double dof = n_n - p;
      //double c = GSL_MAX_DBL(1, sqrt(chisq / dof));
      double c = chisq / dof;

      //fprintf(stderr, "chisq/dof = %g\n", chisq / dof);

      //fprintf (stderr, "E0     = %.25f +/- %.25f\n", FIT(0), ERR(0));//add c*
      //fprintf (stderr, "Gamma1 = %.25f +/- %.25f\n", FIT(1), ERR(1));
      //fprintf (stderr, "Gamma2 = %.25f +/- %.25f\n", FIT(2), ERR(2));
      //fprintf(stderr,"chi/dof  = %.24f\n",(chisq / dof));
    }

    //fprintf (stderr, "status = %s\n", gsl_strerror (status));
    double p0=FIT(0);
    double p1=FIT(1);
    double p2=FIT(2);
    gsl_multifit_nlinear_free (w);
    gsl_matrix_free (covar);
    gsl_rng_free (r);
    //float k= Integrate(1,2,3,4);
    
    #endif
    std::array<double,800> I_fit;
    
    I_fit=calc_current(p0,p1,p2,t);
    std::copy(I_fit.begin()+int((n-number_of_points)/2),I_fit.begin()+n-int((n-number_of_points)/2),I_fit_n.begin());

    double SSE=0,SST=0,GOF=0;
    double sumw = 0;
    double sumyw = 0;
    double mean = 0;
    /* 
    for(int i=0;i<800;i++)
    {
      sumyw=sumyw+y[i]*weights[i];
      sumw=sumw+weights[i];
    }
    mean=sumyw/sumw;
    for(int i=0;i<800;i++)
    {
      SSE=SSE+weights[i]*(I_fit.at(i)-y[i])*(I_fit.at(i)-y[i]);
      SST=SST+weights[i]*(y[i]-mean)*(y[i]-mean);
    }
    GOF=1-((SSE)*(number_of_points-1))/(SST*(number_of_points-3));
    */
    for(int i=0;i<number_of_points;i++)
      {
      sumyw=sumyw+y_n.at(i)*weights_n.at(i);
      sumw=sumw+weights_n.at(i);
      }
    mean=sumyw/sumw;
    for(int i=0;i<number_of_points;i++)
      {
      SSE=SSE+weights_n.at(i)*(I_fit_n.at(i)-y_n.at(i))*(I_fit_n.at(i)-y_n.at(i));
      error1.at(i)=(I_fit_n.at(i)-y_n.at(i));//(I_fit_n.at(i)-y_n.at(i));
      error2.at(i)=SSE;
      SST=SST+weights_n.at(i)*(y_n.at(i)-mean)*(y_n.at(i)-mean);
      error3.at(i)=(y_n.at(i)-mean);
      error4.at(i)=SST;//(y_n.at(i)-mean);
      }
    GOF=1-((SSE)*(number_of_points-1))/(SST*(number_of_points-3));

    ofstream report;
    //take care of UTF8 encoding size=/= to the number of bytes
    double gamma_swap=0;
    //correct this in light of new things which we know .!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //if (p0<0){p0=std::abs(p0);gamma_swap=p1;p1=p2;p2=gamma_swap;}
    //if (p1<0 && p2<0){p1=std::abs(p1);p2=std::abs(p2);}

    p0=std::abs(p0);
    p1=std::abs(p1);
    p2=std::abs(p2);

    auto temp_file=output_file;
    if(downup==0)
      {
        temp_file.insert(output_file.size()-4,"Down"s);
        GOF_array_d[ind]=GOF;
        E0_d[ind]=p0;
        Gamma1_d[ind]=p1;
        Gamma2_d[ind]=p2;
      }
    else
      {
        temp_file.insert(output_file.size()-4,"Up"s);
        GOF_array_u[ind]=GOF;
        E0_u[ind]=p0;
        Gamma1_u[ind]=p1;
        Gamma2_u[ind]=p2;
      }
    
    report.open (temp_file);
    //int curve_no=0;
    //report<<"Curve No."<< <<endl; instead use header...
    /*
    std::vector<double> I_fit_n(2000);
    I_fit_n.resize(number_of_points);
    std::copy(y+int((n-number_of_points)/2),y+n-int((n-number_of_points)/2),I_fit_n.begin());
    */
    //report<< setw( 11 )<< "E0 is "<< p0 <<"Gamma2 is "<<p1<<"Gamma1 is "<<p2<<"GOF is "<<GOF<<"SST is "<<SST<<"SSR is "<<SSE<<endl;
    //report<<"V[i],"<<"I[i],"<<"I_smooth[i],"<<"Weight[i]"<<endl;
    report<<"V[i],"<<"I[i],"<<"I_smooth[i],"<<"I_fit[i],"<<"Weight[i]"<<endl;
    report<<"#"<<fixed << setw( 11 )<< setprecision( 10 )<< "E0 is "<< p0 <<"     "<<"Gamma1 is "<<p1<<"     "<<"Gamma2 is "<<p2<<"     "<<"GOF is "<<GOF<<"     "<<"SST is "<<SST<<"     "<<"SSR is "<<SSE<<"     "<<"No. iter is "<<gsl_multifit_nlinear_niter(w)<<endl;
    for(int i=0;i<number_of_points;i++)
    {
    /*
    report<<fixed << setw( 11 )<<internal<< setprecision( 15 ) << t_n.at(i)<<"       "<< setprecision( 15 ) <<y_n.at(i)<<"       "<<sm_n.at(i)<<"       "<< setw( 11 )<<internal << setprecision( 15 )<< I_fit_n.at(i)
    <<"       "<< setw( 11 )<<scientific<< setprecision( 15 )<<weights_n.at(i)<<"    "<<error1.at(i)<<"     "<<error2.at(i)<<"    "<<error3.at(i)<<"     "<<error4.at(i)<<endl;
    */
      //report<<scientific<< t_n.at(i)<<","<<y_n.at(i)<<","<<sm_n.at(i)<<","<<weights_n.at(i)<<endl;
      report<<scientific<< t_n.at(i)<<","<<y_n.at(i)<<","<<sm_n.at(i)<<","<<I_fit_n.at(i)<<","<<weights_n.at(i)<<endl;
    }
    report.close();


  }
 progress.fetch_add(1, std::memory_order_relaxed);

 //all important things are calculated. In next function entry with the same index file import and weights will not be done again
 wasthere.at(ind)=1;
}




//namespace fs=std::filesystem


/*initializing the static variable to keep track of active workers*/
int WorkerPool::Active_Workers=0;

int
main (void)
{
cout<<endl<<endl<<endl;
WorkerPool wp;
wp.StartPool();
//cout<<wp.Active_Workers<<endl;

std::vector<std::string> m_paths;
//m_paths.push_back("/home/filipk/Desktop/ALEXs measurements/IV-Sweeps_III_2V_range");

//m_paths.push_back("/home/filipk/Desktop/MCBJ_IV_DATA/Salen/Salen_Mn - Sample 0001/2018-05-04_09.18.09 - [MuMBox]_fluid(Toluol)_molecules(Salen_Mn)");
//m_paths.push_back("/home/filipk/Desktop/MCBJ_IV_DATA/Salen/Salen_Mn - Sample 0001/2018-05-03_09.13.41 - [MuMBox]_fluid(Toluol)_molecules(Salen_Mn)");
//m_paths.push_back("/home/filipk/Desktop/MCBJ_IV_DATA/Salen/Salen_Mn - Sample 0001/2018-04-30_22.06.39 - [MuMBox]_fluid(Toluol)_molecules(Salen_Mn)");
//m_paths.push_back("/home/filipk/Desktop/MCBJ_IV_DATA/Salen/Salen_Mn - Sample 0001/2018-04-28_10.06.43 - [MuMBox]_fluid(Toluol)_molecules(Salen_Mn)");
//m_paths.push_back("/home/filipk/Desktop/MCBJ_IV_DATA/Salen/Salen_Fe - Sample 0001/2018-05-16_08.59.23 - [MuMBox]_Toluol_Salen_FE");
//m_paths.push_back("/home/filipk/Desktop/MCBJ_IV_DATA/Salen/Salen_Fe - Sample 0001/2018-05-15_14.47.35 - [MuMBox]_Toluol_Salen_FE");
//m_paths.push_back("/home/filipk/Desktop/MCBJ_IV_DATA/Salen/Salen_Fe - Sample 0001/2018-05-14_15.34.52 - [MuMBox]_Toluol_Salen_FE");
//m_paths.push_back("/home/filipk/Desktop/MCBJ_IV_DATA/Salen/Salen_Fe - Sample 0001/2018-05-09_15.38.06 - [MuMBox]_Toluol_Salen_FE");
//m_paths.push_back("/home/filipk/Desktop/MCBJ_IV_DATA/Salen/Salen_Co - Sample 0001/2018-06-07_14.07.00 - [MuMBox]_fluid(Toluol)_molecules(Salen_CO)");
//m_paths.push_back("/home/filipk/Desktop/MCBJ_IV_DATA/Salen/Salen_Co - Sample 0001/2018-06-04_14.52.01 - [MuMBox]_fluid(Toluol)_molecules(Salen_CO)");
//m_paths.push_back("/home/filipk/Desktop/MCBJ_IV_DATA/Salen/Salen_Co - Sample 0001/2018-05-31_17.32.47 - [MuMBox]_fluid(Toluol)_molecules(Salen_CO)");

//m_paths.push_back("/home/filipk/Desktop/MCBJ_IV_DATA/Salen/corulene/2018-06-21_12.07.02 - [MuMBox]_fluid(Toluol)_molecules(Corrunulaene)");
//m_paths.push_back("/home/filipk/Desktop/ALEXs measurements/IV-Sweeps_II_2V_range");
//m_paths.push_back("/home/filipk/Desktop/C60_test_fit");
//m_paths.push_back("/home/filipk/Desktop/C60_windsize20");
//m_paths.push_back("/home/filipk/Desktop/sept_25_2019_Alex_low temp_IVs/2019-07-29_09.37.55 - range 500mV");
//m_paths.push_back("/home/filipk/Desktop/sept_25_2019_Alex_low temp_IVs/2019-07-29_10.51.13 - range 500mV");
//m_paths.push_back("/home/filipk/Desktop/sept_25_2019_Alex_low temp_IVs/2019-06-28_15.07.37_20K-290K_2V_AUC60AU");
//m_paths.push_back("/home/filipk/Desktop/C60 IVs 2V sweeps June 2019");
//m_paths.push_back("/home/filipk/Desktop/Fullerene_lowtemp_all/other_model2/2019-10-03_19.44.36 - IV_500mV_bias");
//m_paths.push_back("/home/filipk/Desktop/Fullerene_lowtemp_all/other_model/2019-10-06_22.53.43 - Lowtemp_ 1000mV_bias");
//m_paths.push_back("/home/filipk/Desktop/Fullerene_lowtemp_all/2019-10-07_18.46.33 - C60_low temp_1500mV_Bias");
//////m_paths.push_back("/home/filipk/Desktop/Fullerene_lowtemp_all/other_model/2019-10-07_18.46.33 - C60_low temp_1500mV_Bias_Trans2");
///home/filipk/Desktop/sept_25_2019_Alex_low temp IVs/2019-06-28_15.07.37_20K-290K_2V_AUC60AU
//m_paths.push_back("/home/filipk/Desktop/Fullerene_lowtemp_all/2019-10-05_17.05.27 - 200mV_bias");
//m_paths.push_back("/home/filipk/Desktop/Old salen data/2015-12-06_22.21.03 - [MuMBox]_fluid(oldToluol)_molecules(salen)_MNIV");
//m_paths.push_back("/home/filipk/Desktop/Old salen data/2015-12-19_11.22.43 - [MuMBox]_fluid(Toluol)_salen_co_IV_Curves");
//m_paths.push_back("/home/filipk/Desktop/Old salen data/2018-05-22_18.56.05 - [MuMBox]_fluid(Toluol)_molecules(Salen_CO)");
//m_paths.push_back("/home/filipk/Desktop/Old salen data/2018-05-29_09.43.25 - [MuMBox]_fluid(Toluol)_molecules(Salen_CO)");
//m_paths.push_back("/home/filipk/Desktop/Old salen data/2018-06-18_15.23.33 - [MuMBox]_Liquid(Toluol)_corrunolene");
//m_paths.push_back("/home/filipk/Desktop/Old salen data/2018-06-21_12.07.02 - [MuMBox]_fluid(Toluol)_molecules(Corrunulaene)");
//m_paths.push_back("/home/filipk/Desktop/MCBJ_Coranulene/2019-11-22_21.34.41 - [MuMBox]_fluid(MESITYLENEl)_molecules(CORANULENE)");
//m_paths.push_back("/home/filipk/Desktop/2020-03-05_16.47.13 - unknownEnv");
//m_paths.push_back("/home/filipk/Desktop/2020-03-09_16.46.38 - unknownEnv");
//m_paths.push_back("/home/filipk/Desktop/2020-03-13_01.19.45 - unknownEnv");
//m_paths.push_back("/home/filipk/Desktop/Old salen data/cobalt/Updated IV3");

//m_paths.push_back("/home/filipk/Desktop/Refitinh/Salen_Co - Sample 0001/2018-05-31_17.32.47 - [MuMBox]_fluid(Toluol)_molecules(Salen_CO)");
//m_paths.push_back("/home/filipk/Desktop/Refitinh/Salen_Co - Sample 0001/2018-06-04_14.52.01 - [MuMBox]_fluid(Toluol)_molecules(Salen_CO)");
//m_paths.push_back("/home/filipk/Desktop/Refitinh/Salen_Co - Sample 0001/2018-06-07_14.07.00 - [MuMBox]_fluid(Toluol)_molecules(Salen_CO)");

//m_paths.push_back("/home/filipk/Desktop/Refitinh/Salen_Fe - Sample 0001/2018-05-09_15.38.06 - [MuMBox]_Toluol_Salen_FE");
//m_paths.push_back("/home/filipk/Desktop/Refitinh/Salen_Fe - Sample 0001/2018-05-14_15.34.52 - [MuMBox]_Toluol_Salen_FE");
//m_paths.push_back("/home/filipk/Desktop/Refitinh/Salen_Fe - Sample 0001/2018-05-15_14.47.35 - [MuMBox]_Toluol_Salen_FE");
//m_paths.push_back("/home/filipk/Desktop/Refitinh/Salen_Fe - Sample 0001/2018-05-16_08.59.23 - [MuMBox]_Toluol_Salen_FE");

m_paths.push_back("/home/filipk/Desktop/Refitinh/Salen_Mn - Sample 0001/2018-04-28_10.06.43 - [MuMBox]_fluid(Toluol)_molecules(Salen_Mn)");
m_paths.push_back("/home/filipk/Desktop/Refitinh/Salen_Mn - Sample 0001/2018-04-30_22.06.39 - [MuMBox]_fluid(Toluol)_molecules(Salen_Mn)");
m_paths.push_back("/home/filipk/Desktop/Refitinh/Salen_Mn - Sample 0001/2018-05-03_09.13.41 - [MuMBox]_fluid(Toluol)_molecules(Salen_Mn)");
m_paths.push_back("/home/filipk/Desktop/Refitinh/Salen_Mn - Sample 0001/2018-05-04_09.18.09 - [MuMBox]_fluid(Toluol)_molecules(Salen_Mn)");


cout<<"Number of files to be evaluated  "<<m_paths.size()<<endl;
//for(auto t_path : m_paths){cout<<t_path<<endl;}
for(std::string t_path : m_paths)
{

std::fill(wasthere.begin(),wasthere.end(),0); //why did i comment it out?




//std::fill(&E0_d[0], &E0_d[0]+1000, 0);
//std::fill(&Gamma1_d[0], &Gamma1_d[0]+10000, 5e-3);
//std::fill(&Gamma2_d[0], &Gamma2_d[0]+10000, 5e-3);
//std::fill(&E0_u[0], &E0_u[0]+10000, 0.5);
//std::fill(&Gamma1_u[0], &Gamma1_u[0]+10000, 5e-3);
//std::fill(&Gamma2_u[0], &Gamma2_u[0]+10000, 5e-3);
  //std::fill(wasthere.begin(),wasthere.end(),0);
    
  //cout<<"Input the path to the folder with measurements data";
  //cin<<;
  //gsl_ieee_env_setup ();
  auto start = std::chrono::high_resolution_clock::now();
  auto current_time = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed =start-current_time;
  std::array<std::string,10000> paths;
  std::array<int,10000> indexes;
  std::array<std::string,10000> sorted_paths;
  //sorted_paths.fill("eeeee"); !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
  
  double temp_GOF=0;
  int init_file_number=0;

  std::string temp;
  std::string extr;
  extr.reserve(100);
  std::string lp;
  int counter=0;
  std::string path_folder = t_path; // ="/home/filipk/Desktop/C60 IVs 2V sweeps June 2019"; 
  std::string path = path_folder+"/IV-Sweeps";
  std::vector<double> perf_time;
  
  /*
  //get the list of file paths
  //1508
  //1645
  std::vector<std::string> many_paths(100);
  string line;
  ifstream path_list ("/home/filipk/Desktop/Final_versions_fitting/Fit_program_new_version/ML/many_paths.txt");
  if (path_list.is_open())
  {
    while ( getline (path_list,line) )
    {
      cout << line << '\n';
    }
    path_list.close();
  }
  */
  //for(int vrang=1;vrang<9;vrang++)
  //int vrang =6;
  //{
  /*
  for (const auto & entry : std::experimental::filesystem::directory_iterator(path))
  {
  //std::cout << entry.path() << std::endl;
  paths[counter]=entry.path();
  counter=counter+1;
  }
  for(int i=0;i<counter;i++){
  ifstream file(paths.at(i), ios_base::in | ios_base::binary);
  filtering_streambuf<input> in;
  in.push(gzip_decompressor());
  in.push(file);  
  //std::vector<int> v={1}; !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  std::istream incoming(&in);

  MeasFile file1(incoming);
  //file1.PrintHeader();
  cout<< file1.data_upsweep_V.at(799)<<endl;
  }
  fit("/home/filipk/Desktop/Fit_program/test data/Salen_mn/IV-Sweeps/000315.txt.gz","Output_data/ttt.txt",315);
 */
 counter=0;

#if 1
  for (const auto & entry : boost::filesystem::directory_iterator(path))
  //for (boost::filesystem::directory_iterator itr(path_ss); itr!=boost::filesystem::directory_iterator(); ++itr)
  {
  //std::cout << entry.path() << std::endl;
  //paths[counter]=entry.path();
  //cout<<entry.path().extension().string()<<endl;
  if (entry.path().extension().string()==".gz")
  {
  paths[counter]=entry.path().string();
  //cout<<entry.path().string()<<endl;;
  counter=counter+1;
  }
  else
  {
    cout<<"Non .gz files present "<<endl;
    cout<<entry.path().string()<<endl;
  } 
  }
  for(int i=0;i<counter;i++)
  {
  //std::cout << paths[i] << std::endl;
  temp=paths[i];
  extr.assign(temp);
  extr.assign(extr,temp.length()-13,6);
  indexes.at(i)=stoi(extr);
  //std::copy((paths.at(i)).begin(),(paths.at(i)).end(),&temp);
  //std::cout << extr <<endl;
  }
  init_file_number=*std::min_element(indexes.begin(),indexes.end());
  for(int i=0;i<counter;i++)
  {
   sorted_paths.at(indexes.at(i)-1-init_file_number)=paths.at(i); 
  }


  /*
  for(int i=0;i<counter;i++)
  {
    cout<<sorted_paths.at(i)<<endl;

  }

  cout<< sorted_paths.size()<<endl;
  */


 
  //*std::min_element(indexes.begin(),indexes.end())
  //std::sort(paths.begin(),paths.end(),sorted_paths.begin());
  //for(int i =0;i<counter;i++){cout<<sorted_paths.at(i)<<endl;}
  //removed standard threading scheme
  //std::thread t[counter];
  //adding jobes for Matlab pool

  //check the voltage max
  ifstream file(sorted_paths.at(0), ios_base::in | ios_base::binary);
  filtering_streambuf<input> in;
  in.push(gzip_decompressor());
  in.push(file);  
  //std::vector<int> v={1}; !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  std::istream incoming(&in);
  MeasFile file_check(incoming);
  //std::copy(file1.data_downsweep_V.begin(),file1.data_downsweep_V.end(),tvec_d.at(ind).begin())
  //max_V = std::ceil(*std::max_element(file_check.data_downsweep_V.begin(),file_check.data_downsweep_V.end())*10)*100;
  //cout<< "MAX_V"<<max_V<<endl;
  //max_V=1000;

  max_V=800;
  for(int vrang=8;vrang<int(max_V/100)+1;vrang++)
  {
  //int vrang=20;
  //WorkerPool wp;
 // wp.StartPool(62);
 ////////////////////////////////////////////add logical switch for fitting specific number of curves
 //counter=50;
  for(int i=0;i<counter;i++) //1883
  {
    temp=paths[i];
    extr.assign(temp);
    extr.assign(extr,temp.length()-13,6);
    //cout<<sorted_paths[i]<<endl;
    #if 0
    t[i]=std::thread(fit,sorted_paths.at(i),"Output_data/"+to_string(i)+".txt",i);
    #endif
    //std::string fitdata_folder=""
    std::string command ="mkdir -p ";
    command= command +"\""  + path_folder + "/Output_data_" + to_string(vrang*100)+"mV/"+"\"" ;
    if (std::system(command.data()) == -1) {
               int errsv = errno;
               printf("somecall() failed\n");
               cout<<errsv<<endl;
           }
    std::function<void()> k  = std::bind(fit,sorted_paths.at(i),path_folder+"/Output_data_"+to_string(vrang*100)+"mV/"+to_string(i+1)+".txt",i,vrang*100);
    wp.Add_Job(k);
    //TURN THIS ONN FOR NICE OUTPUT
    //cout<<sorted_paths[i]<<endl;
  }




  //nanosleep((const struct timespec[]){{400, 0}}, NULL);
  std::array<std::string,5> stars={"","*","**","***","****"};
  int looptime=0;
  cout<<endl;
  #if 1
  //execlp("tput","tput", "civis", (char *)NULL);
  //std::system("tput civis");
/*
  try 
    {
      std::system("tput civis"); 
    } 
  catch(const std::system_error& e) 
    {
        std::cout << "Caught system_error with code " << e.code() 
                  << " meaning " << e.what() << '\n';
    }
*/
  if (std::system("tput civis") == -1) {
               int errsv = errno;
               printf("somecall() failed\n");
               cout<<errsv<<endl;
           }



  int prog_temp=0;
  auto offset_time =std::chrono::high_resolution_clock::now();
  while(progress.load(std::memory_order_relaxed)!=counter)
  {
    //cout<<"\t\t\t\t\t\t\t\t";
    current_time = std::chrono::high_resolution_clock::now();
    prog_temp=progress.load(std::memory_order_relaxed);
    elapsed=(current_time-offset_time);
    looptime++;
    cout<<"\r"<<"Progress:   "<<std::right<<setw( 4 )<<stars[looptime%5] <<std::left<< setw( 5 )<<prog_temp<<"out of  "<<counter<<""<<setw( 4 )
    <<stars[looptime%5]<<"Elapsed time:"<<int(elapsed.count())<<"  Time left:"<<((elapsed.count()/prog_temp)*(counter-prog_temp))<< "\t\t";
    std::cout.flush();
    //usleep(1000000);
    nanosleep((const struct timespec[]){{1, 0}}, NULL);
    //sleep(1);
  }
  //execlp("tput","tput", "cnorm", (char *)NULL);
  //std::system("tput cnorm");
  if (std::system("tput cnorm") == -1) {
              int errsv = errno;
              printf("somecall() failed\n");
              cout<<errsv<<endl;
          }
  //wp.StopPool();
#endif
/*
  for(int i=0;i<counter;i++)
  {
    t[i].join();
  }
 */ 

  ofstream myfile3;
  myfile3.open (path_folder+"/GOF"+to_string(vrang*100)+".txt");
  /*write file header*/
  myfile3<<"curve_No"<< 
  "," <<"GOF_array_d"<<"," <<"Transmission_d" << "," <<"E0_d"<<"," << "Gamma1_d"<<"," <<"Gamma2_d"<<
  "," <<"GOF_array_u"<<"," <<"Transmission_u"<< "," <<"E0_u"<<"," <<"Gamma1_u"<<"," <<"Gamma2_u"<<
  "\n";

  for(int i=0;i<counter;i++)
  {
  /*
  myfile3<<fixed << setw( 6 ) << setprecision( 5 ) <<i+1<<"\t"<< 
  setw( 15 ) << setprecision( 5 )<<GOF_array_d[i]<<setw( 15 ) << setprecision( 10 )<<Transmission(0,E0_d[i],Gamma1_d[i],Gamma2_d[i],0) << setw( 15 ) << setprecision( 10 )<<E0_d[i]<<setw( 15 ) << setprecision( 10 )<<Gamma1_d[i]<<setw( 15 ) << setprecision( 10 )<<Gamma2_d[i]<<
  setw( 15 ) << setprecision( 5 )<<GOF_array_u[i]<<setw( 15 ) << setprecision( 10 )<<Transmission(0,E0_u[i],Gamma1_u[i],Gamma2_u[i],0)<< setw( 15 ) << setprecision( 10 )<<E0_u[i]<<setw( 15 ) << setprecision( 10 )<<Gamma1_u[i]<<setw( 15 ) << setprecision( 10 )<<Gamma2_u[i]<<
  "\n";
  */
  myfile3<<fixed << setprecision( 5 ) <<i+1<< 
  "," << setprecision( 5 )<<GOF_array_d[i]<<"," << setprecision( 10 )<<Transmission(0,E0_d[i],Gamma1_d[i],Gamma2_d[i],0) << "," << setprecision( 10 )<<E0_d[i]<<"," << setprecision( 10 )<<Gamma1_d[i]<<"," << setprecision( 10 )<<Gamma2_d[i]<<
  "," << setprecision( 5 )<<GOF_array_u[i]<<"," << setprecision( 10 )<<Transmission(0,E0_u[i],Gamma1_u[i],Gamma2_u[i],0)<< "," << setprecision( 10 )<<E0_u[i]<<"," << setprecision( 10 )<<Gamma1_u[i]<<"," << setprecision( 10 )<<Gamma2_u[i]<<
  "\n";
  }
  myfile3.close();
  progress.store(0,std::memory_order_relaxed);
  //}







  printf("DONE");
  auto finish = std::chrono::high_resolution_clock::now();
  elapsed = finish - start;
  std::cout << "Elapsed time: " << elapsed.count() << " s\n";
  perf_time.push_back(elapsed.count());
#endif 
  }
  std::reverse(perf_time.begin(),perf_time.end());
  
  ofstream myfile4;
  myfile4.open(path_folder+"/performance.txt");
  for(int i=0;i<(max_V/100);i++)
  {
  myfile4<<perf_time.back()<<endl;
  perf_time.pop_back();
  }
  myfile4.close();
  
}
wp.StopPool();
  return 0;
}