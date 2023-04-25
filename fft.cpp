#include <iostream>
#include <complex>
#include <cmath>
#include <iomanip>



using namespace std;

std::complex<double> W(int N,int r)
{
    const double pi = 3.14;
    return  std::complex<double>( cos(  ( ( 2 * pi ) / N ) * r ), sin(  ( ( 2 * pi ) / N ) * r ) ) ;
}




void tab1istab2(complex<double> *t0,complex<double> *t1,int N)
{

    for(int i = 0; i < N; i++)
    {
        t0[i] = t1[i];
    }

}


int bit_revers(int x,int N)
{

    int revers_value = 0;
    int r;

    while(x != 0)
    {
        r = x % 2;
        x>>=1;
        N=N/2;
        revers_value = revers_value + r * N;

    }

    return revers_value;

}


void fft(complex<double> *t0,complex<double> *t1,int N)
{

for(int i = 0;i<N;i++) t1[i] = t0[bit_revers(i,N)];




for(int k = 2; k<=N; k=k*2)
{

 tab1istab2(t0,t1,N);

    for(int i = 0; i < N/k; i++)
    {

        complex<double> w(1 ,0);

         for(int j = 0; j < k/2; j++)
         {

            t1[i*k+j] = t1[i*k+j] + w * t0[i*k+j+k/2];
            t1[i*k+j+k/2] = -w * t1[i*k+j+k/2] + t0[i*k+j];

            w = w/W(N,N/k);

         }

    }



}







}


void reverse_fft(complex<double> *t0,complex<double> *t1,int N)
{
complex<double> pd(2,0);


for(int k = N; k > 1; k=k/2)
{



     tab1istab2(t0,t1,N);

    for(int i = 0; i < N/k; i++)
    {

        complex<double> w(1,0);


        int j = 0;

        t1[i*k+j] = (t1[i*k+j] + w * t0[i*k+j+k/2])/pd;
        t1[i*k+j+k/2] = (-w * t1[i*k+j+k/2] + t0[i*k+j])/pd;

         for(int j = 1; j < k/2; j++)
         {
            t1[i*k+j] = (t1[i*k+j] + t0[i*k+j+k/2])/pd;
            t1[i*k+j+k/2] = (t1[i*k+j+k/2] + t0[i*k+j])/pd;

         }

    }


}




tab1istab2(t0,t1,N);

for(int i = 0;i<N;i++) t1[i] = t0[bit_revers(i,N)];



}




int main()
{
cout << fixed << setprecision(4);
  int N = 16;

   complex<double> n[N] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
   complex<double> n_next[N];





   fft(n,n_next,N);
  reverse_fft(n,n_next,N);



     for(int i = 0;i<N;i++)
   {

       cout<<i<<" ["<<n_next[i]<<"] "<<endl;

   }




    return 0;
}

